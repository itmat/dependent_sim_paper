library(dependentsimr)
library(tidyverse)
library(rlang)

method <- snakemake@wildcards$method
N_GENES <- snakemake@wildcards$n_genes
N_SAMPLES <- 100

# use ZT0 of raw time series data to get dependence structure
male_read_counts <- read.csv("processed/Cortex_ZT0-counts.csv") |> 
    select(!all_of("X331308_CRM55_WT_CTXR_26WKS_M_ZT0_L1.D704")) |>
    column_to_rownames(var="raw_data.EnsemblID")

selected_rows <- sample.int(nrow(male_read_counts), N_GENES)
input <- male_read_counts[selected_rows,]
## Generate Controls
# With dependence
time <- system.time({
    actual_library_sizes <- input |> apply(2, sum)
    # Get dependence structure
    if (method == "pca") {
        set.seed(1)
        # Simulate with the PCA method
        rs <- get_random_structure(list(data=input), method="pca", rank=2, types="DESeq2")
    } else if (method == "corpcor") {
        set.seed(2)
        # Simulate with the corpcor method
        rs <- get_random_structure(list(data=input), method="corpcor", types="DESeq2")
    } else if (method == "wishart") {
        set.seed(3)
        # Simulate with the spiked Wishart method
        rs <- get_random_structure(list(data=input), rank=10, method="spiked Wishart", types="DESeq2")
    } else if (method == "indep") {
        set.seed(4)
        # Simulate without any dependence
        dep_rs <- get_random_structure(list(data=input), method="pca", rank=2, type="DESeq2")
        rs <- remove_dependence(dep_rs)
    }

    # Generate simulated data
    library_sizes <- sample(actual_library_sizes / mean(actual_library_sizes), size=N_SAMPLES, replace=TRUE)
    sim <- draw_from_multivariate_corr(rs, n_samples=N_SAMPLES, size_factors=library_sizes)$data
})

message("Elapsed time:")
message(time[['elapsed']])
tibble(
    method = method,
    n_genes = N_GENES,
    time = time[['elapsed']],
) |> write_tsv(paste0("processed/benchmark/", method, "/", N_GENES, ".txt"))
