library(dependentsimr)
library(tidyverse)
library(DESeq2)
library(rlang)

method <- snakemake@wildcards$method

data <- read.csv("data/plasma_nmr.csv")
metadata <- data |> 
    filter(Label == 0) |> 
    select(Samples, Label)
real_data <- data |> 
    filter(Label == 0) |> 
    select(-Label)

transposed_df <- real_data |> 
    pivot_longer(cols = -Samples, names_to = "Metabolites", values_to = "value") |> 
    pivot_wider(names_from = Samples, values_from = value) |> 
    select(-C60)
real_data <- transposed_df |> column_to_rownames(var="Metabolites") |> 
    log()
write_csv(real_data, "processed/Plasma_metabolomics_real.csv")
real_data <- real_data |> 
    as.matrix()

# Produce simulated data
N_SAMPLES <- 59*8
if (method == "pca") {
    set.seed(1)
    # Simulate with the PCA method
    rs <- get_random_structure(list(data=real_data), method="pca", rank=2, types="normal")
} else if (method == "corpcor") {
    set.seed(2)
    # Simulate with the corpcor method
    rs <- get_random_structure(list(data=real_data), method="corpcor", types="normal")
} else if (method == "wishart") {
    set.seed(3)
    # Simulate with the spiked Wishart method
    rs <- get_random_structure(list(data=real_data), rank=11, method="spiked Wishart", types="normal")
} else if (method == "indep") {
    set.seed(4)
    # Simulate without any dependence
    dep_rs <- get_random_structure(list(data=real_data), method="pca", rank=2, type="normal")
    rs <- remove_dependence(dep_rs)
}

generate <- function(rs, name) {
    # Generate and save data+metadata (true values) from a specific random structure
    sim <- draw_from_multivariate_corr(rs, n_samples=N_SAMPLES)$data
    colnames(sim) <- paste("sample", 1:ncol(sim), sep='')

    sim %>%
        data.frame %>%
        rownames_to_column("Metabolites") %>%
        write.table(file=paste0("simulated_data/", name, ".txt"), sep="\t", row.names=FALSE)

    return(sim)
}

# Generate the control fly data
met_sim <- generate(rs, name=paste0("Plasma_metabolomics.",method))
