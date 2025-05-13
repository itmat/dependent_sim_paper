library(tidyverse)
library(SPsimSeq)

N_SAMPLES <- 100

read_counts <- read.csv("processed/Cortex_ZT0-counts.csv") |> 
    column_to_rownames(var="raw_data.EnsemblID")

# Produce simulated data
not_all_zero <- rowSums(read_counts) > 0

input <- read_counts[not_all_zero,] |> as.matrix()
set.seed(1)
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

saveRDS(results, "processed/SPsimSeq/GSE151565.RDS")

counts <- results$sim.data.list[[1]]$counts
# Add back the original gene names onto the simulated data
counts |>
    as.data.frame() |>
    rownames_to_column("temp_id") |>
    left_join(
        results$sim.data.list[[1]]$rowData |>
            rownames_to_column("temp_id"),
        by = join_by(temp_id == temp_id),
    ) |>
    mutate(gene_id = source.ID) |>
    select(-temp_id, -DE.ind, -source.ID) |>
    write_tsv("simulated_data/SPsimSeq.GSE151565.txt")
