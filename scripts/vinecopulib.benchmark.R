library(rvinecopulib)
library(tidyverse)
N_GENES <- strtoi(snakemake@wildcards$n_genes)
N_SAMPLES <- 100
method <- "vinecopula"

# Generate some very basic data
X <- runif(N_GENES * N_SAMPLES) |> matrix(N_GENES, N_SAMPLES)
udv <- svd(X)
true_d <- seq(N_SAMPLES, 1)
U <- udv$u
temp <- rnorm(length(true_d) * N_SAMPLES) |> matrix(true_d, N_SAMPLES)
data <- U %*% diag(true_d) %*% temp
norm_data <- pnorm(data) |> t()
print(norm_data |> dim())

if (N_GENES > 2000) {
    # Don't run slow ones
    tibble(
        method = method,
        n_genes = N_GENES,
        time = NA,
    ) |> write_tsv(paste0("processed/benchmark/", method, "/", N_GENES, ".txt"))
} else {
    set.seed(1)
    time <- system.time({
        # Fit a vinecopula model to it
        vc <- rvinecopulib::vinecop(norm_data, family_set = "gaussian")

        # Simulate data using that vinecopula
        sim_data <- rvinecopulib::rvinecop(vc, n= N_SAMPLES)
    })
    print(sim_data |> dim())
    message("Elapsed time:")
    message(time[['elapsed']])
    tibble(
        method = method,
        n_genes = N_GENES,
        time = time[['elapsed']],
    ) |> write_tsv(paste0("processed/benchmark/", method, "/", N_GENES, ".txt"))

}
