library(mvnfast)
library(tidyverse)
library(rlang)

method <- "mvnfast"
N_GENES <- snakemake@wildcards$n_genes
N_SAMPLES <- 100

# use ZT0 of raw time series data to get dependence structure
male_read_counts <- read.csv("processed/Cortex_ZT0-counts.csv") |> 
    column_to_rownames(var="raw_data.EnsemblID")

selected_rows <- sample.int(nrow(male_read_counts), N_GENES)
input <- male_read_counts[selected_rows,]

## Generate Controls
# With dependence
time <- system.time({
    # Note: add a small amount for positive definiteness
    Sigma <- cov(input |> t()) + diag(N_GENES)*0.01
    data <- rmvn(
        n = N_SAMPLES,
        mu = rowMeans(input),
        sigma = Sigma,
    )
})

message("Elapsed time:")
message(time[['elapsed']])
tibble(
    method = method,
    n_genes = N_GENES,
    time = time[['elapsed']],
) |> write_tsv(paste0("processed/benchmark/", method, "/", N_GENES, ".txt"))
