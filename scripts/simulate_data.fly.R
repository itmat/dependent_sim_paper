library(dependentsimr)
library(tidyr)
library(dplyr)
library(tibble)
library(DESeq2)
library(rlang)

method <- snakemake@wildcards$method

read_counts <- read.delim("data/GSE81142.counts.txt.gz", row.names=1)
metadata <- read.delim("processed/GSE81142_sample_metadata.txt")

# Produce simulated data
N_SAMPLES <- 100
male_read_counts <- read_counts[,(metadata$sex == "Male") & (metadata$concentration == 0)]
actual_library_sizes <- male_read_counts |> apply(2, sum)
if (method == "pca") {
    set.seed(1)
    # Simulate with the PCA method
    rs <- get_random_structure(list(data=male_read_counts), method="pca", rank=2, types="DESeq2")
} else if (method == "corpcor") {
    set.seed(2)
    # Simulate with the corpcor method
    rs <- get_random_structure(list(data=male_read_counts), method="corpcor", types="DESeq2")
} else if (method == "wishart") {
    set.seed(3)
    # Simulate with the spiked Wishart method
    rs <- get_random_structure(list(data=male_read_counts), rank=11, method="spiked Wishart", types="DESeq2")
} else if (method == "indep") {
    set.seed(4)
    # Simulate without any dependence
    dep_rs <- get_random_structure(list(data=male_read_counts), method="pca", rank=2, type="DESeq2")
    rs <- remove_dependence(dep_rs)
}

library_sizes <- sample(actual_library_sizes / mean(actual_library_sizes), size=N_SAMPLES, replace=TRUE)
generate <- function(rs, name) {
  # Generate and save data+metadata (true values) from a specific random structure
  sim <- draw_from_multivariate_corr(rs, n_samples=N_SAMPLES, size_factors=library_sizes)$data
  colnames(sim) <- paste("sample", 1:ncol(sim), sep='')

  sim %>%
    data.frame %>%
    rownames_to_column("ENSEMBL_ID") %>%
    write.table(file=paste0("simulated_data/", name, ".txt"), sep="\t", row.names=FALSE)

  # metadata
  truth <- data.frame(
    gene_id = rownames(sim),
    mean = rs$marginals$data$q,
    disperion = rs$marginals$data$dispersion
  )
  truth[is.na(truth$mean), 2] <- 0 # NA mean is 0
  write.table(truth, file=paste0("simulated_data/", name, ".true_values.txt"), sep="\t", row.names=FALSE)

  return(sim)
}

# Generate the control fly data
male_sim <- generate(rs, name=paste0("Fly.WholeBody.Male.",method,".Control"))

# Randomly select some genes to be DE
FRAC_DE <- 0.05
LFC_MIN <- 0.2
LFC_MAX <- 2.0
N_DE <- floor(FRAC_DE * nrow(male_sim))
de_genes <- sample(nrow(male_sim), N_DE)
de_lfc <- runif(N_DE, LFC_MIN, LFC_MAX) * sample(c(-1, 1), size=N_DE, replace=TRUE)
## Generate Cases
if (method == "indep") {
    male_de_rs <- dep_rs
    male_de_rs$marginals$data$q[de_genes] <- 2^(de_lfc) * male_de_rs$marginals$data$q[de_genes]
    
    male_indep_de_rs <- remove_dependence(male_de_rs)
    male_indep_de_sim <- generate(male_indep_de_rs, name="Fly.WholeBody.Male.indep.Case")
} else {
    male_de_rs <- rs
    male_de_rs$marginals$data$q[de_genes] <- 2^(de_lfc) * male_de_rs$marginals$data$q[de_genes]
    
    male_de_sim <- generate(male_de_rs, name=paste0("Fly.WholeBody.Male.",method,".Case"))
}
