library(dependentsimr)
library(dplyr)
library(tibble)
library(DESeq2)
library(tidyverse)

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
sex <- vapply(colnames(read_counts),
    function(id) {
        return(metadata[metadata["sample_id"] == id, "sex"])
    }, "")

male_read_counts <- read_counts[,sex == 'male']
female_read_counts <- read_counts[,sex == 'female']

# use ZT0 of raw time series data to get dependence structure
male_read_counts <- read.csv("processed/Liver_ZT0-counts.csv") |> 
    column_to_rownames(var="raw_data.EnsemblID")

# Produce simulated data
male_rs <- get_random_structure(list(data=male_read_counts), method="pca", rank=2, types="DESeq2")

actual_library_sizes <- male_read_counts |> apply(2, sum)

N_SAMPLES <- 120
set.seed(0)

generate <- function(rs, seed) {
  # Generate and save data+metadata (true values) from a specific random structure
  set.seed(seed)
  library_sizes <- sample(actual_library_sizes / mean(actual_library_sizes), size=N_SAMPLES, replace=TRUE)
  sim <- draw_from_multivariate_corr(rs, n_samples=N_SAMPLES, size_factors=library_sizes)$data
  colnames(sim) <- paste("sample", 1:ncol(sim), sep='')

  sim <- sim |> 
    data.frame() |> 
    rownames_to_column("ENSEMBL_ID")

  return(sim)
}

## Generate Controls
# With dependence
male_sim <- generate(male_rs, seed=1)

# Without dependence (k=0)
male_indep_rs <- remove_dependence(male_rs)
male_indep_sim <- generate(male_indep_rs, seed=2)

# read time point means
time_point_means <- read.csv("processed/mean_per_time_Liver.csv") |> 
    inner_join(tibble(EnsemblID = male_rs$rownames$data), by='EnsemblID')
time_point_means$mean <- time_point_means[, -1] |> 
    mutate_all(as.numeric) |> 
    rowMeans()

# adjust male_read_counts with time
data_per_time <- list()
indep_data_per_time <- list()
for (i in (2 : (ncol(time_point_means)-1))) {
    time_rs <- male_rs
    factor <- time_point_means[, i]/time_point_means[, 'mean']
    time_rs$marginals$data$q <- factor * time_rs$marginals$data$q
    
    # Generate with dependence
    time_sim <- generate(time_rs, seed=2*i)
    time <- colnames(time_point_means)[i]
    colnames(time_sim) <- c('ENSEMBL_ID', paste0(time, '_sample', 1:120))
    
    # Without dependence (k=0)
    indep_time_rs <- remove_dependence(time_rs)
    indep_time_sim <- generate(indep_time_rs, seed=2*i+1)
    colnames(indep_time_sim) <- c('ENSEMBL_ID', paste0(time, '_sample', 1:120))
    
    data_per_time[[i-1]] <- time_sim
    indep_data_per_time[[i-1]] <- indep_time_sim
}

time_data <- data_per_time |> 
    purrr::reduce(full_join, by="ENSEMBL_ID")
write.csv(time_data, "simulated_data/Liver_120_simulated_time_series_k=2.csv", row.names = FALSE)
norm_time_data <- time_data
read_depth <- apply(time_data[, -1], 2, sum)
norm_time_data[, -1] <- as.tibble(as.matrix(time_data[, -1]) %*% diag(1000000/read_depth))
write.csv(norm_time_data, "simulated_data/Liver_120_normalized_simulated_time_series_k=2.csv", row.names = FALSE)

indep_time_data <- indep_data_per_time |> 
    purrr::reduce(full_join, by="ENSEMBL_ID")
write.csv(indep_time_data, "simulated_data/Liver_120_simulated_time_series_k=0.csv", row.names = FALSE)
norm_indep_time_data <- indep_time_data
read_depth <- apply(indep_time_data[, -1], 2, sum)
norm_indep_time_data[, -1] <- as.tibble(as.matrix(indep_time_data[, -1]) %*% diag(1000000/read_depth))
write.csv(norm_indep_time_data, "simulated_data/Liver_120_normalized_simulated_time_series_k=0.csv", row.names = FALSE)



