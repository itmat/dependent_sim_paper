library(dependentsimr)
library(dplyr)
library(tibble)
library(DESeq2)
library(tidyverse)
library(rlang)

method <- snakemake@wildcards$method

# use ZT0 of raw time series data to get dependence structure
male_read_counts <- read.csv("processed/Cortex_ZT0-counts.csv") |> 
    select(!all_of("X331308_CRM55_WT_CTXR_26WKS_M_ZT0_L1.D704")) |>
    column_to_rownames(var="raw_data.EnsemblID")

# Get dependence structure
if (method == "pca") {
    set.seed(1)
    # Simulate with the PCA method
    male_rs <- get_random_structure(list(data=male_read_counts), method="pca", rank=2, types="DESeq2")
} else if (method == "corpcor") {
    set.seed(2)
    # Simulate with the corpcor method
    male_rs <- get_random_structure(list(data=male_read_counts), method="corpcor", types="DESeq2")
} else if (method == "wishart") {
    set.seed(3)
    # Simulate with the spiked Wishart method
    male_rs <- get_random_structure(list(data=male_read_counts), rank=10, method="spiked Wishart", types="DESeq2")
} else if (method == "indep") {
    set.seed(4)
    # Simulate without any dependence
    dep_rs <- get_random_structure(list(data=male_read_counts), method="pca", rank=2, type="DESeq2")
    male_rs <- remove_dependence(dep_rs)
}

actual_library_sizes <- male_read_counts |> apply(2, sum)

N_SAMPLES <- 120

generate <- function(rs, seed) {
  set.seed(seed)
    # Generate and save data+metadata (true values) from a specific random structure
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
male_sim <- generate(male_rs, .Random.seed)

# read time point means
time_point_means <- read.csv("processed/mean_per_time_Cortex.csv") |> 
    inner_join(tibble(EnsemblID = male_rs$rownames$data), by='EnsemblID')
time_point_means$mean <- time_point_means[, -1] |> 
    mutate_all(as.numeric) |> 
    rowMeans()

# adjust male_read_counts with time
data_per_time <- list()
for (i in (2 : (ncol(time_point_means)-1))) {
    time_rs <- male_rs
    factor <- time_point_means[, i]/time_point_means[, 'mean']
    time_rs$marginals$data$q <- factor * time_rs$marginals$data$q
    
    if (method == "indep") {
        time_rs <- dep_rs
        time_rs$marginals$data$q <- factor * time_rs$marginals$data$q
        time_rs <- remove_dependence(time_rs)
    }
    
    # Generate simulated data at time point
    time_sim <- generate(time_rs, seed=.Random.seed+4*i)
    time <- colnames(time_point_means)[i]
    colnames(time_sim) <- c('ENSEMBL_ID', paste0(time, '_sample', 1:120))
    data_per_time[[i-1]] <- time_sim
}

time_data <- data_per_time |> 
    purrr::reduce(full_join, by="ENSEMBL_ID")
write.csv(time_data, paste0("simulated_data/Cortex_120_simulated_time_series_",method,".csv"), row.names = FALSE)
norm_time_data <- time_data
read_depth <- apply(time_data[, -1], 2, sum)
norm_time_data[, -1] <- as_tibble(as.matrix(time_data[, -1]) %*% diag(1000000/read_depth))
write.csv(norm_time_data, paste0("simulated_data/Cortex_120_normalized_simulated_time_series_",method,".csv"), row.names = FALSE)



