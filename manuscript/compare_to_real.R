
# Comparison of simulated data to a real dataset

library(dependentsimr)
library(ggplot2)
library(ggpointdensity)
library(tidyr)
library(dplyr)
library(patchwork)

set.seed(55)

# CONFIGURATION
N_rows_gene_corr <- 3000
HIGH_EXPR_CUTOFF <- 300

## Generate the simulated data
# Load the real data ---------------------------------
# from GEO: GSE77221
head(Weger18)
weger18_library_sizes <- apply(Weger18, 2, sum)

# Run dependentsimr on this data --------------------
N_SAMPLES <- 96
library_sizes <- (weger18_library_sizes + rep(0, N_SAMPLES)) / mean(weger18_library_sizes) # Recycle the existing library sizes
# Simulate with the PCA method
rs <- get_random_structure(list(counts=Weger18), method="pca", rank=2, type="DESeq2")
draws <- draw_from_multivariate_corr(rs, n_samples=N_SAMPLES, size_factors=library_sizes)$counts

# Simulate with the corpcor method
rs_corpcor <- get_random_structure(list(counts=Weger18), method="corpcor", type="DESeq2")
draws_corpcor <- draw_from_multivariate_corr(rs_corpcor, n_samples=N_SAMPLES, size_factors=library_sizes)$counts

# Simulate without any dependence
rs_indep <- remove_dependence(rs)
indep_draws <- draw_from_multivariate_corr(rs_indep, n_samples=N_SAMPLES, size_factors=library_sizes)$counts


# Scale to Counts Per Million ------------------------
cpm <- function(x) { # counts per million
  return(t(t(x) / apply(x, 2, sum) * 1e6))
}
scaled_read_data <- cpm(Weger18)
scaled_draws <- cpm(draws)
scaled_draws_corpcor <- cpm(draws_corpcor)
scaled_indep_draws <- cpm(indep_draws)

## Generate the figures

simulated_mean <- apply(scaled_draws[,1:ncol(Weger18)], 1, mean)
simulated_variance <- apply(scaled_draws[,1:ncol(Weger18)], 1, var)
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
  geom_abline(slope=1, intercept=0)
g2<-ggplot(marg_dist|> filter(real_mean > 0.1))+
  geom_pointdensity(aes(x=log(real_variance+1),y=log(simulated_variance+1))) +
  geom_abline(slope=1, intercept=0)
  
# Plot gene-gene correlation ---------------------------
high_expr_rows <- which(apply(Weger18, 1, mean) > HIGH_EXPR_CUTOFF)
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
  genes <- sample(high_expr_rows, size=100, replace=FALSE)
  selected_samples <-sample(ncol(draws), ncol(Weger18))
  real_corr <- distinct_cor(scaled_read_data[genes,] |> t())
  sim_corr <- distinct_cor(scaled_draws[genes,selected_samples]  |> t())
  corpcor_corr <- distinct_cor(scaled_draws_corpcor[genes,selected_samples]  |> t())
  indep_corr <- distinct_cor(scaled_indep_draws[genes,selected_samples]  |> t())
  corr_data[[length(corr_data)+1]] <- data.frame(
    real = real_corr,
    pca = sim_corr,
    corpcor = corpcor_corr,
    indep = indep_corr
  )
}
corr_data <- do.call(rbind, corr_data)
corr_data_long <- corr_data |> pivot_longer(everything(), names_to="type", values_to="sim_corr")
quantiles <- seq(0,1,length.out=1001)
corr_quantiles <- corr_data_long |> group_by(type) |> reframe(p=quantiles, correlation=quantile(sim_corr, probs=quantiles))
real_quantiles <- corr_quantiles |> filter(type == "real")
corr_plot_data <- corr_quantiles |>
  select(p, type, sim_correlation = correlation) |>
  left_join(real_quantiles |> select(p, real_correlation = correlation), by = join_by(p)) |>
  filter(type != "real")
g3<-ggplot(data = corr_plot_data) +
  geom_point(aes(x = real_correlation, y = sim_correlation, color=type)) +
  geom_abline(slope=1, intercept=0) +
  labs(
    x = "Correlation coefficient\n(real data)",
    y = "Correlation coefficient\n(simulated data)")

# Plot the top principal components --------------------
constant_rows <- apply(Weger18, 1, function(x) all(x - mean(x) == 0))
pca <- prcomp(t(scaled_read_data[!constant_rows,]), rank.=2, scale.=TRUE)
projected_read_data <- predict(pca, t(scaled_read_data[!constant_rows, ]))
projected_draws <- predict(pca, t(scaled_draws[!constant_rows, ]))
projected_corpcor <- predict(pca, t(scaled_draws_corpcor[!constant_rows, ]))
projected_indep_draws <- predict(pca, t(scaled_indep_draws[!constant_rows, ]))
both_data <- data.frame(list(
  PC1 = c(projected_read_data[,"PC1"], projected_draws[,"PC1"], projected_indep_draws[,"PC1"]),
  PC2 = c(projected_read_data[,"PC2"], projected_draws[,"PC2"], projected_indep_draws[,"PC2"]),
  type = c(
    rep("real", nrow(projected_read_data)),
    rep("dependent sim", nrow(projected_draws)),
    rep("independent sim", nrow(projected_indep_draws))
  )
))
both_data$type <- factor(both_data$type, levels=c("real", "dependent sim", "independent sim"))
g4<-ggplot(data = both_data, aes(x=PC1, y=PC2,color=type)) +
   geom_point() +
  scale_color_manual(values = list(real="blue", `dependent sim`="#00BFC4", `independent sim`="#F8766D"))
# Compare SVD with and without dependence ----------------
PCA <- function(data) {
  constant_rows <- apply(data, 1, function(x) all(x - mean(x) == 0))
  return(prcomp(t(data[!constant_rows,]), rank.=2, scale.=TRUE))
}
pca_sdevs <- tibble(
  type = "real",
  PC = 1:12,
  variance = PCA(scaled_read_data)$sdev,
)
for (i in 0:7) {
  pca_sdevs <- rbind(
    pca_sdevs,
    tibble(
      type = "dependent sim",
      PC = 1:12,
      variance = PCA(scaled_draws[,(1+12*i):(12*(i+1))])$sdev,
    ),
    tibble(
      type = "independent sim",
      PC = 1:12,
      variance = PCA(scaled_indep_draws[,(1+12*i):(12*(i+1))])$sdev,
    ),
    deparse.level=1
  )
}
g5 <- ggplot(data=pca_sdevs |> filter(PC<5)) +
  facet_grid(cols=vars(PC), labeller=label_both) +
  geom_boxplot(aes(x=type, y=variance))

compare_to_real_plot <- (g1 | g2) / (g3 | g4) / (g5) + plot_annotation(tag_levels="a")
