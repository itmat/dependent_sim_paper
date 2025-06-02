# Comparison of simulated data to a real dataset

library(dependentsimr)
library(ggplot2)
library(ggpointdensity)
library(tidyr)
library(dplyr)
library(readr)
library(patchwork)
library(rlang)

set.seed(55)

# CONFIGURATION
HIGH_EXPR_CUTOFF <- 1
TYPE_ORDER <- c("real", "indep", "pca", "wishart", "corpcor")
TYPE_COLORS <- c("black", RColorBrewer::brewer.pal(length(TYPE_ORDER)-1, "Dark2"))

# Load the real data ---------------------------------
# raw <- read_tsv("../data/GSE151923_metaReadCount_ensembl.txt.gz") |>
#   select(-GeneName, -Description, -Chromosome, -Strand)
# read_counts <- as.matrix(raw[,2:13])
# rownames(read_counts) <- raw[[1]]
# 
# actual_library_sizes <- apply(read_counts, 2, sum)
# 
# # Run dependentsimr on this data --------------------
# N_SAMPLES <- 96
# library_sizes <- (actual_library_sizes + rep(0, N_SAMPLES)) / mean(actual_library_sizes) # Recycle the existing library sizes
# # Simulate with the PCA method
# rs <- get_random_structure(list(counts=read_counts), method="pca", rank=2, type="DESeq2")
# draws <- draw_from_multivariate_corr(rs, n_samples=N_SAMPLES, size_factors=library_sizes)$counts
# 
# # Simulate with the corpcor method
# rs_corpcor <- get_random_structure(list(counts=read_counts), method="corpcor", type="DESeq2")
# draws_corpcor <- draw_from_multivariate_corr(rs_corpcor, n_samples=N_SAMPLES, size_factors=library_sizes)$counts
# 
# # Simulate with the spiked Wishart method
# rs_wishart <- get_random_structure(list(counts=read_counts), rank=11, method="spiked Wishart", type="DESeq2")
# draws_wishart <- draw_from_multivariate_corr(rs_wishart, n_samples=N_SAMPLES, size_factors=library_sizes)$counts
# 
# # Simulate without any dependence
# rs_indep <- remove_dependence(rs)
# indep_draws <- draw_from_multivariate_corr(rs_indep, n_samples=N_SAMPLES, size_factors=library_sizes)$counts
# 
# 
# # Scale to Counts Per Million ------------------------
# cpm <- function(x) { # counts per million
#   return(t(t(x) / apply(x, 2, sum) * 1e6))
# }
# scaled_read_data <- cpm(read_counts)
# scaled_draws <- cpm(draws)
# scaled_draws_corpcor <- cpm(draws_corpcor)
# scaled_draws_wishart <- cpm(draws_wishart)
# scaled_indep_draws <- cpm(indep_draws)

scaled_read_data <- read_csv("data/Plasma_metabolomics_real.csv")
scaled_draws <- read_delim("simulated_data/Plasma_metabolomics.pca.txt", delim = '\t') |> 
    select(-Metabolites) |> 
    as.matrix()
scaled_draws_corpcor <- read_delim("simulated_data/Plasma_metabolomics.corpcor.txt", delim = '\t') |> 
    select(-Metabolites) |> 
    as.matrix()
scaled_draws_wishart <- read_delim("simulated_data/Plasma_metabolomics.wishart.txt", delim = '\t') |> 
    select(-Metabolites) |> 
    as.matrix()
scaled_indep_draws <- read_delim("simulated_data/Plasma_metabolomics.indep.txt", delim = '\t') |> 
    select(-Metabolites) |> 
    as.matrix()

## Generate the figures

simulated_mean <- apply(scaled_draws[,1:ncol(read_counts)], 1, mean)
simulated_variance <- apply(scaled_draws[,1:ncol(read_counts)], 1, var)
real_mean <- apply(scaled_read_data, 1, mean)
real_variance <- apply(scaled_read_data, 1, var)
marg_dist <- data.frame(
  real_mean = real_mean,
  simulated_mean = simulated_mean,
  real_variance = real_variance,
  simulated_variance = simulated_variance
  )
g1<-ggplot(marg_dist|> filter(real_mean > 0.1))+
  geom_pointdensity(aes(x=log(real_mean+1),y=log(simulated_mean+1))) +
  geom_abline(slope=1, intercept=0) +
  scale_color_viridis_c(option="inferno")
g2<-ggplot(marg_dist|> filter(real_mean > 0.1))+
  geom_pointdensity(aes(x=log(real_variance+1),y=log(simulated_variance+1))) +
  geom_abline(slope=1, intercept=0) +
  scale_color_viridis_c(option="inferno")

# Plot the top principal components --------------------
constant_rows <- apply(read_counts, 1, function(x) all(x - mean(x) == 0))
pca <- prcomp(t(scaled_read_data[!constant_rows,]), rank.=2, scale.=TRUE)
projected_read_data <- predict(pca, t(scaled_read_data[!constant_rows, ]))
projected_draws <- predict(pca, t(scaled_draws[!constant_rows, ]))
projected_corpcor <- predict(pca, t(scaled_draws_corpcor[!constant_rows, ]))
projected_wishart <- predict(pca, t(scaled_draws_wishart[!constant_rows, ]))
projected_indep_draws <- predict(pca, t(scaled_indep_draws[!constant_rows, ]))
both_data <- data.frame(list(
    PC1 = c(projected_read_data[,"PC1"], projected_draws[1:12,"PC1"], projected_corpcor[1:12,"PC1"], projected_wishart[1:12,"PC1"], projected_indep_draws[1:12,"PC1"]),
    PC2 = c(projected_read_data[,"PC2"], projected_draws[1:12,"PC2"], projected_corpcor[1:12,"PC2"], projected_wishart[1:12,"PC2"], projected_indep_draws[1:12,"PC2"]),
    type = c(
        rep("real", nrow(projected_read_data)),
        rep("pca", 12),
        rep("corpcor", 12),
        rep("wishart", 12),
        rep("indep", 12)
    )
)) %>%
  mutate(type = factor(type, levels=TYPE_ORDER, ordered=TRUE))
#both_data$type <- factor(both_data$type, levels=c("real", "dependent sim", "independent sim"))
g3<-ggplot(data = both_data, aes(x=PC1, y=PC2,color=type)) +
   geom_point() +
  scale_color_manual(values = TYPE_COLORS, breaks = TYPE_ORDER, name="type")

# Plot gene-gene correlation ---------------------------
high_expr_rows <- which(apply(read_counts, 1, mean) > HIGH_EXPR_CUTOFF)
distinct_cor <- function(x) {
  # correlations of pairs of variables in x
  # except no redundant pairs [e.g., (a,b) and (b,a), or (a,a)]
  mat <- cor(x)
  mat[lower.tri(mat)]
}
corr_data <- list()
for (i in 1:100) {
  # Compute correlation in batches of 100 random rows with randomly chosen simulated samples
  # then those are aggregated to give the whole sample
  genes <- sample(high_expr_rows, size=30, replace=FALSE)
  selected_samples <-sample(ncol(scaled_draws), ncol(read_counts))
  real_corr <- distinct_cor(scaled_read_data[genes,] |> t())
  sim_corr <- distinct_cor(scaled_draws[genes,selected_samples]  |> t())
  corpcor_corr <- distinct_cor(scaled_draws_corpcor[genes,selected_samples]  |> t())
  wishart_corr <- distinct_cor(scaled_draws_wishart[genes,selected_samples]  |> t())
  indep_corr <- distinct_cor(scaled_indep_draws[genes,selected_samples]  |> t())
  corr_data[[length(corr_data)+1]] <- data.frame(
    real = real_corr,
    pca = sim_corr,
    corpcor = corpcor_corr,
    wishart = wishart_corr,
    indep = indep_corr
  )
}
corr_data <- do.call(rbind, corr_data)
corr_data_long <- corr_data |> pivot_longer(everything(), names_to="type", values_to="corr") %>%
  mutate(type=factor(type, levels=TYPE_ORDER, ordered=TRUE))
g4 <- ggplot(data = corr_data_long) +
  geom_density(aes(x = corr, color=type), linewidth=1) +
  scale_color_manual(values = TYPE_COLORS, breaks = TYPE_ORDER, name="type")


# Compare SVD with and without dependence ----------------
PCA <- function(data) {
  constant_rows <- apply(data, 1, function(x) { all(near(x - mean(x), 0)) })
  return(prcomp(t(data[!constant_rows,]), rank.=2, scale.=TRUE))
}
pca_sdevs <- tibble(
  type = "real",
  PC = 1:42,
  sdev = PCA(scaled_read_data)$sdev,
)
for (i in 0:7) {
  pca_sdevs <- rbind(
    pca_sdevs,
    tibble(
      type = "pca",
      PC = 1:12,
      sdev = PCA(scaled_draws[,(1+12*i):(12*(i+1))])$sdev,
    ),
    tibble(
      type = "corpcor",
      PC = 1:12,
      sdev = PCA(scaled_draws_corpcor[,(1+12*i):(12*(i+1))])$sdev,
    ),
    tibble(
      type = "wishart",
      PC = 1:12,
      sdev = PCA(scaled_draws_wishart[,(1+12*i):(12*(i+1))])$sdev,
    ),
    tibble(
      type = "indep",
      PC = 1:12,
      sdev = PCA(scaled_indep_draws[,(1+12*i):(12*(i+1))])$sdev,
    ),
    deparse.level=1
  )
}
pca_sdevs$type <- factor(pca_sdevs$type, levels=TYPE_ORDER, ordered=TRUE)
g5 <- ggplot() +
  facet_wrap(facets=vars(PC), labeller=label_both, ncol=6) +
  geom_hline(data=pca_sdevs |> filter(PC < 12, type == "real"), aes(yintercept=sdev, color="real"), linewidth=1) +
  geom_boxplot(data=pca_sdevs |> filter(PC < 12, type != "real"), aes(x=type, y=sdev, color=type)) +
  scale_color_manual(values=TYPE_COLORS, breaks=TYPE_ORDER, name="type") +
  theme(axis.text.x = element_text(angle=90))

compare_to_real_plot <- (g1 | g2) / (g3 | g4) / (g5) + plot_annotation(tag_levels="a")
