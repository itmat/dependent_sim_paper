library(dependentsimr)
library(dplyr)
library(tibble)
library(DESeq2)

metadata <- read.delim("processed/GSE151923_sample_metadata.txt", sep="\t")
read_counts <- read.delim("data/GSE151923_metaReadCount_ensembl.txt.gz", sep="\t")

row.names(read_counts) <- read_counts[,1]
read_counts <- read_counts[,6:ncol(read_counts)]

# Drop NAs - the data file is missing part of the very last row, giving NAs
read_counts <- read_counts[!apply(is.na(read_counts), 1, any),]

# Drop a sample that has very low read counts
read_depth <- apply(read_counts, 2, sum)
READ_DEPTH_CUTOFF = 10e6
read_counts <- read_counts[,read_depth > READ_DEPTH_CUTOFF]

# Separate by sex
print(colnames(read_counts))
print(metadata)
sex <- vapply(colnames(read_counts),
    function(id) {
        return(metadata[metadata["sample_id"] == id, "sex"])
    }, "")

male_read_counts <- read_counts[,sex == 'male']
female_read_counts <- read_counts[,sex == 'female']

#### Compute M/F DE genes and log-fold-changes
###dds <- DESeqDataSetFromMatrix(as.matrix(read_counts), data.frame(sex = sex), ~sex)
###dds <- DESeq(dds)
###lfc <- results(dds, c("sex", "female", "male"))$log2FoldChange

# Produce simulated data
male_rs <- get_random_structure(list(data=male_read_counts), method="pca", rank=2, types="DESeq2")

actual_library_sizes <- male_read_counts |> apply(2, sum)

N_SAMPLES <- 100
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

## Generate Controls
# With dependence
male_sim <- generate(male_rs, seed=1, name="Mouse.Cortex.Male.k=2.Control")

# Without dependence (k=0)
male_indep_rs <- remove_dependence(male_rs)
male_indep_sim <- generate(male_indep_rs, seed=2, name="Mouse.Cortex.Male.k=0.Control")

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

## Generate Cases
male_de_sim <- generate(male_de_rs, seed=3, name="Mouse.Cortex.Male.k=2.Case")

male_indep_de_rs <- remove_dependence(male_de_rs)
male_indep_de_sim <- generate(male_indep_de_rs, seed=4, name="Mouse.Cortex.Male.k=0.Case")
