library(tidyverse)
library(SPsimSeq)

read_counts <- read.delim("data/GSE81142.counts.txt.gz", row.names=1)
metadata <- read.delim("processed/GSE81142_sample_metadata.txt")

# Produce simulated data
N_SAMPLES <- 100
male_read_counts <- read_counts[,(metadata$sex == "Male") & (metadata$concentration == 0)]

not_all_zero <- rowSums(male_read_counts) > 0

input <- male_read_counts[not_all_zero,] |> as.matrix()
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
    log.CPM.transform = TRUE, #???
    result.format = "list",
    return.details = TRUE,
    verbose = TRUE,
)

saveRDS(results, "processed/SPsimSeq/GSE81142.rds")

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
    write_tsv("simulated_data/SPsimSeq.GSE81142.txt")
