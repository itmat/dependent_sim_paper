library(tidyverse)
input = snakemake@input$results
benchmarks = snakemake@input$benchmarks
temp <- list()
for (i in seq_along(input)) {
    temp[[length(temp)+1]] <- bind_cols(
        read_tsv(input[[i]]),
        read_tsv(benchmarks[[i]]),
    )
}
bind_rows(temp) |> write_tsv("processed/benchmark/results.txt")
