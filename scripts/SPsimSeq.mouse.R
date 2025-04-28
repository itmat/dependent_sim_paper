library(tidyverse)
library(SPsimSeq)

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
    log.CPM.transform = TRUE, #???
    result.format = "list",
    return.details = TRUE,
    verbose = TRUE,
)

saveRDS(results, "processed/SPsimSeq/GSE151923.RDS")

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
    write_tsv("simulated_data/SPsimSeq.GSE151923.txt")
