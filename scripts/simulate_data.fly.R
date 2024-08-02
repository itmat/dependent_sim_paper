library(dependentsimr)
library(tidyr)
library(dplyr)
library(tibble)
library(DESeq2)

read_counts <- read.delim("data/GSE81142.counts.txt.gz", row.names=1)
metadata <- read.delim("processed/GSE81142_sample_metadata.txt")

selected_samples <- (metadata$hour != 0)
selected_read_counts <- read_counts[,selected_samples]
selected_metadata <- metadata[selected_samples,]

# Produce simulated data
N_SAMPLES <- 100
male_read_counts <- read_counts[,(metadata$sex == "Male") & (metadata$concentration == 0)]
actual_library_sizes <- male_read_counts |> apply(2, sum)
male_rs <- get_random_structure(list(data=male_read_counts), rank=2, types="DESeq2")

set.seed(0)
library_sizes <- sample(actual_library_sizes / mean(actual_library_sizes), size=N_SAMPLES, replace=TRUE)
generate <- function(rs, seed, name) {
  # Generate and save data+metadata (true values) from a specific random structure
  set.seed(seed)
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
male_sim <- generate(male_rs, seed=1, name="Fly.WholeBody.Male.k=2.Control")
male_indep_sim <- generate(remove_dependence(male_rs), seed=2, name="Fly.WholeBody.Male.k=0.Control")

# Randomly select some genes to be DE
set.seed(0)
FRAC_DE <- 0.05
LFC_MIN <- 0.2
LFC_MAX <- 2.0
N_DE <- floor(FRAC_DE * nrow(male_sim))
de_genes <- sample(nrow(male_sim), N_DE)
de_lfc <- runif(N_DE, LFC_MIN, LFC_MAX) * sample(c(-1, 1), size=N_DE, replace=TRUE)
male_de_rs <- male_rs
male_de_rs$marginals$data$q[de_genes] <- 2^(de_lfc) * male_de_rs$marginals$data$q[de_genes]

male_de_sim <- generate(male_de_rs, seed=3, name="Fly.WholeBody.Male.k=2.Case")
male_indep_de_sim <- generate(remove_dependence(male_de_rs), seed=4, name="Fly.WholeBody.Male.k=0.Case")
