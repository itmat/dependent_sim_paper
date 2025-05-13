library(tidyverse)
library(SPsimSeq)

N_SAMPLES <- 100
N_GENES <- snakemake@wildcards$n_genes
method <- "SPsimSeq"

read_counts <- read.csv("processed/Cortex_ZT0-counts.csv") |> 
    column_to_rownames(var="raw_data.EnsemblID")

# Produce simulated data
not_all_zero <- rowSums(read_counts) > 0

input <- read_counts[not_all_zero,] |> as.matrix()
selected_rows <- sample.int(nrow(input), N_GENES, replace=TRUE)
input <- input[selected_rows,]
message(paste("Size of input", input |> nrow()))
set.seed(1)
time <- system.time({
    results <- SPsimSeq(
        n.sim = 1,
        s.data = input,
        n.genes = input |> nrow(),
        batch.config = 1, # only one batch,
        group.config = 1, # only one group
        pDE = 0, # no DE since no groups
        tot.samples = N_SAMPLES,
        model.zero.prob = FALSE, # since bulk RNA
        genewiseCor = TRUE, # gene-gene dependence, notes that it could be SLOW
        log.CPM.transform = TRUE,
        result.format = "list",
        return.details = TRUE,
        verbose = TRUE,
    )
})
message("Elapsed time:")
message(time[['elapsed']])
tibble(
    method = method,
    n_genes = N_GENES,
    time = time[['elapsed']],
) |> write_tsv(paste0("processed/benchmark/", method, "/", N_GENES, ".txt"))
