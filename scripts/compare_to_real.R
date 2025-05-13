# Comparison of simulated data to a real dataset

library(dependentsimr)
library(tidyverse)
library(ggpointdensity)
library(patchwork)
library(rlang)

set.seed(55)

#dataset = "GSE81142"
#dataset = "GSE151923"
#dataset = "GSE151565"
dataset <- snakemake@wildcards$dataset

# CONFIGURATION
METHOD_ORDER <- c("real", "indep", "pca", "wishart", "corpcor", "spsimseq")
METHOD_COLORS <- c("black", RColorBrewer::brewer.pal(length(METHOD_ORDER)-1, "Dark2"))
N_SAMPLES <- 96

# Load the real and simulated data ---------------------------------
if (dataset == "GSE151923") {
    HIGH_EXPR_CUTOFF <- 100
    raw <- read_tsv("data/GSE151923_metaReadCount_ensembl.txt.gz") |>
        select(-GeneName, -Description, -Chromosome, -Strand)
    read_counts <- as.matrix(raw[,2:13])
    rownames(read_counts) <- raw[[1]]

    draws_pca <- (read_tsv("simulated_data/Mouse.Cortex.Male.pca.Control.txt") |>
                  column_to_rownames("ENSEMBL_ID"))[,1:N_SAMPLES] |> as.matrix()
    draws_corpcor <- (read_tsv("simulated_data/Mouse.Cortex.Male.corpcor.Control.txt") |>
                  column_to_rownames("ENSEMBL_ID"))[,1:N_SAMPLES] |> as.matrix()
    draws_wishart <- (read_tsv("simulated_data/Mouse.Cortex.Male.wishart.Control.txt") |>
                  column_to_rownames("ENSEMBL_ID"))[,1:N_SAMPLES] |> as.matrix()
    draws_indep <- (read_tsv("simulated_data/Mouse.Cortex.Male.indep.Control.txt") |>
                  column_to_rownames("ENSEMBL_ID"))[,1:N_SAMPLES] |> as.matrix()

    draws_spsimseq <- read_tsv("simulated_data/SPsimSeq.GSE151923.txt") |>
        column_to_rownames("gene_id")
} else if (dataset == "GSE81142") {
    HIGH_EXPR_CUTOFF <- 30
    raw <- read_tsv("data/GSE81142.counts.txt.gz") |> column_to_rownames("GENE_ID")
    metadata <- read_tsv("processed/GSE81142_sample_metadata.txt")
    read_counts <- raw[,(metadata$sex == "Male") & (metadata$concentration == 0)] |> as.matrix()

    draws_pca <- (read_tsv("simulated_data/Fly.WholeBody.Male.pca.Control.txt") |>
                  column_to_rownames("ENSEMBL_ID"))[,1:N_SAMPLES] |> as.matrix()
    draws_wishart <- (read_tsv("simulated_data/Fly.WholeBody.Male.wishart.Control.txt") |>
                  column_to_rownames("ENSEMBL_ID"))[,1:N_SAMPLES] |> as.matrix()
    draws_corpcor <- (read_tsv("simulated_data/Fly.WholeBody.Male.corpcor.Control.txt") |>
                  column_to_rownames("ENSEMBL_ID"))[,1:N_SAMPLES] |> as.matrix()
    draws_indep <- (read_tsv("simulated_data/Fly.WholeBody.Male.indep.Control.txt") |>
                  column_to_rownames("ENSEMBL_ID"))[,1:N_SAMPLES] |> as.matrix()

    draws_spsimseq <- read_tsv("simulated_data/SPsimSeq.GSE81142.txt") |>
        column_to_rownames("gene_id")
} else if (dataset == "GSE151565") {
    HIGH_EXPR_CUTOFF <- 30
    raw <- read.csv("processed/Cortex_ZT0-counts.csv") |> 
        column_to_rownames(var="raw_data.EnsemblID")
    read_counts <- as.matrix(raw)

    draws_pca <- (read_csv("simulated_data/Cortex_120_simulated_time_series_pca.csv") |>
                  select("ENSEMBL_ID", starts_with("ZT0_")) |>
                  column_to_rownames("ENSEMBL_ID"))[,1:N_SAMPLES] |> as.matrix()
    draws_corpcor <- (read_csv("simulated_data/Cortex_120_simulated_time_series_corpcor.csv") |>
                  select("ENSEMBL_ID", starts_with("ZT0_")) |>
                  column_to_rownames("ENSEMBL_ID"))[,1:N_SAMPLES] |> as.matrix()
    draws_wishart <- (read_csv("simulated_data/Cortex_120_simulated_time_series_wishart.csv") |>
                  select("ENSEMBL_ID", starts_with("ZT0_")) |>
                  column_to_rownames("ENSEMBL_ID"))[,1:N_SAMPLES] |> as.matrix()
    draws_indep <- (read_csv("simulated_data/Cortex_120_simulated_time_series_indep.csv") |>
                  select("ENSEMBL_ID", starts_with("ZT0_")) |>
                  column_to_rownames("ENSEMBL_ID"))[,1:N_SAMPLES] |> as.matrix()

    draws_spsimseq <- read_tsv("simulated_data/SPsimSeq.GSE151565.txt") |>
        column_to_rownames("gene_id")
}

N_ORIG_SAMPLES <- ncol(read_counts)

# fill in missing all-zero genes and reorder
missing <- rownames(read_counts)[!(rownames(read_counts) %in% rownames(draws_spsimseq))]
present <- rownames(draws_spsimseq)
zero_genes <- matrix(0, length(missing), ncol(draws_spsimseq))
rownames(zero_genes) <- missing
colnames(zero_genes) <- colnames(draws_spsimseq)
draws_spsimseq <- rbind(
    draws_spsimseq,
    zero_genes
)
draws_spsimseq <- draws_spsimseq[rownames(read_counts),] # put in the right order

# Normalization code --------------------------------
# adapted from edgeR
calcNormFactors.default <- function(object, lib.size=NULL, refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75, ...)
#	Scale normalization of RNA-Seq data, for count matrices
#	Mark Robinson, Gordon Smyth and edgeR team
#	Created 22 October 2009. Last modified 29 Dec 2023.
{
#	Check object
    x <- as.matrix(object)
    if(anyNA(x)) stop("NA counts not permitted")
    nsamples <- ncol(x)

#	Check lib.size
    if(is.null(lib.size)) {
        lib.size <- colSums(x)
    } else {
        if(anyNA(lib.size)) stop("NA lib.sizes not permitted")
        if(length(lib.size) != nsamples) {
            if(length(lib.size) > 1L) warning("calcNormFactors: length(lib.size) doesn't match number of samples",call.=FALSE)
            lib.size <- rep_len(lib.size,nsamples)
        }
    }

#	Remove all zero rows
    allzero <- .rowSums(x>0, nrow(x), nsamples) == 0L
    if(any(allzero)) x <- x[!allzero,,drop=FALSE]

#	Calculate factors
    if( is.null(refColumn) ) {
        f75 <- suppressWarnings(.calcFactorQuantile(data=x, lib.size=lib.size, p=0.75))
        if(median(f75) < 1e-20) {
            refColumn <- which.max(colSums(sqrt(x)))
        } else {
            refColumn <- which.min(abs(f75-mean(f75)))
        }
    }
    f <- rep_len(NA_real_,nsamples)
    for(i in 1:nsamples)
        f[i] <- .calcFactorTMM(obs=x[,i],ref=x[,refColumn], libsize.obs=lib.size[i], libsize.ref=lib.size[refColumn], logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)

#	Factors should multiple to one
    f <- f/exp(mean(log(f)))

#	Output
    names(f) <- colnames(x)
    f
}

.calcFactorTMM <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10)
#	TMM between two libraries
#	Mark Robinson
{
    obs <- as.numeric(obs)
    ref <- as.numeric(ref)

    if( is.null(libsize.obs) ) nO <- sum(obs) else nO <- libsize.obs
    if( is.null(libsize.ref) ) nR <- sum(ref) else nR <- libsize.ref

    logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
    absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
    v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance

#	remove infinite values, cutoff based on A
    fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)

    logR <- logR[fin]
    absE <- absE[fin]
    v <- v[fin]

    if(max(abs(logR)) < 1e-6) return(1)

#	taken from the original mean() function
    n <- length(logR)
    loL <- floor(n * logratioTrim) + 1
    hiL <- n + 1 - loL
    loS <- floor(n * sumTrim) + 1
    hiS <- n + 1 - loS

#	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
#	a fix from leonardo ivan almonacid cardenas, since rank() can return
#	non-integer values when there are a lot of ties
    keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)

    if(doWeighting)
        f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
    else
        f <- mean(logR[keep], na.rm=TRUE)

#	Results will be missing if the two libraries share no features with positive counts
#	In this case, return unity
    if(is.na(f)) f <- 0
    2^f
}

.calcFactorQuantile <- function (data, lib.size, p=0.75)
#	Generalized version of upper-quartile normalization
#	Mark Robinson and Gordon Smyth
#	Created 16 Aug 2010. Last modified 12 Sep 2020.
{
    f <- rep_len(1,ncol(data))
    for (j in seq_len(ncol(data))) f[j] <- quantile(data[,j], probs=p)
    if(min(f)==0) warning("One or more quantiles are zero")
    f / lib.size
}


# Scale to Counts Per Million ------------------------
cpm <- function(x) { # counts per million
    read_depths <- apply(x, 2, sum)

    factors <- calcNormFactors.default(x)
    return(t(t(x) / factors / read_depths * 1e6))
}
scaled_read_data <- cpm(read_counts)
scaled_draws_pca <- cpm(draws_pca)
scaled_draws_corpcor <- cpm(draws_corpcor)
scaled_draws_wishart <- cpm(draws_wishart)
scaled_draws_indep <- cpm(draws_indep)
scaled_draws_spsimseq <- cpm(draws_spsimseq)

## Generate the figures

simulated_mean <- apply(scaled_draws_pca[,1:ncol(read_counts)], 1, mean)
simulated_variance <- apply(scaled_draws_pca[,1:ncol(read_counts)], 1, var)
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
constant_rows <- apply(read_counts, 1, function(x) all(abs(x - mean(x)) < 1e-5))
pca <- prcomp(t(scaled_read_data[!constant_rows,]), rank.=2, scale.=TRUE)
projected_read_data <- predict(pca, t(scaled_read_data[!constant_rows, ]))
projected_pca <- predict(pca, t(scaled_draws_pca[!constant_rows, ]))
projected_corpcor <- predict(pca, t(scaled_draws_corpcor[!constant_rows, ]))
projected_wishart <- predict(pca, t(scaled_draws_wishart[!constant_rows, ]))
projected_indep <- predict(pca, t(scaled_draws_indep[!constant_rows, ]))
projected_spsimseq <- predict(pca, t(scaled_draws_spsimseq[!constant_rows, ]))
both_data <- data.frame(list(
    PC1 = c(
            projected_read_data[1:N_ORIG_SAMPLES,"PC1"],
            projected_pca[1:N_ORIG_SAMPLES,"PC1"],
            projected_corpcor[1:N_ORIG_SAMPLES,"PC1"],
            projected_wishart[1:N_ORIG_SAMPLES,"PC1"],
            projected_indep[1:N_ORIG_SAMPLES,"PC1"],
            projected_spsimseq[1:N_ORIG_SAMPLES, "PC1"]
    ),
    PC2 = c(
            projected_read_data[1:N_ORIG_SAMPLES,"PC2"],
            projected_pca[1:N_ORIG_SAMPLES,"PC2"],
            projected_corpcor[1:N_ORIG_SAMPLES,"PC2"],
            projected_wishart[1:N_ORIG_SAMPLES,"PC2"],
            projected_indep[1:N_ORIG_SAMPLES,"PC2"],
            projected_spsimseq[1:N_ORIG_SAMPLES, "PC2"]
    ),
    method = c(
        rep("real", N_ORIG_SAMPLES),
        rep("pca", N_ORIG_SAMPLES),
        rep("corpcor", N_ORIG_SAMPLES),
        rep("wishart", N_ORIG_SAMPLES),
        rep("indep", N_ORIG_SAMPLES),
        rep("spsimseq", N_ORIG_SAMPLES)
    )
)) %>%
  mutate(method = factor(method, levels=METHOD_ORDER, ordered=TRUE))
#both_data$method <- factor(both_data$method, levels=c("real", "dependent sim", "independent sim"))
g3<-ggplot(data = both_data, aes(x=PC1, y=PC2,color=method, shape=method)) +
   geom_point() +
  scale_color_manual(values = METHOD_COLORS, breaks = METHOD_ORDER, name="method")

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
  genes <- sample(high_expr_rows, size=100, replace=FALSE)
  selected_samples <-sample(ncol(draws_pca), ncol(read_counts))
  real_corr <- distinct_cor(scaled_read_data[genes,] |> t())
  sim_corr <- distinct_cor(scaled_draws_pca[genes,selected_samples]  |> t())
  corpcor_corr <- distinct_cor(scaled_draws_corpcor[genes,selected_samples]  |> t())
  wishart_corr <- distinct_cor(scaled_draws_wishart[genes,selected_samples]  |> t())
  indep_corr <- distinct_cor(scaled_draws_indep[genes,selected_samples]  |> t())
  spsimseq_corr <- distinct_cor(scaled_draws_spsimseq[genes,selected_samples]  |> t())
  corr_data[[length(corr_data)+1]] <- data.frame(
    real = real_corr,
    pca = sim_corr,
    corpcor = corpcor_corr,
    wishart = wishart_corr,
    indep = indep_corr,
    spsimseq = spsimseq_corr
  )
}
corr_data <- do.call(rbind, corr_data)
corr_data_long <- corr_data |> pivot_longer(everything(), names_to="method", values_to="corr") %>%
  mutate(method=factor(method, levels=METHOD_ORDER, ordered=TRUE))
g4 <- ggplot(data = corr_data_long) +
  geom_density(aes(x = corr, color=method), linewidth=1) +
  scale_color_manual(values = METHOD_COLORS, breaks = METHOD_ORDER, name="method")


# Compare SVD with and without dependence ----------------
PCA <- function(data) {
  constant_rows <- apply(data, 1, function(x) { all(near(x - mean(x), 0)) })
  return(prcomp(t(data[!constant_rows,]), rank.=2, scale.=TRUE))
}
pca_sdevs <- tibble(
  method = "real",
  PC = 1:N_ORIG_SAMPLES,
  sdev = PCA(scaled_read_data)$sdev,
)
for (i in 0:7) {
  this_range <- (1+N_ORIG_SAMPLES*i):(N_ORIG_SAMPLES*(i+1))
  pca_sdevs <- rbind(
    pca_sdevs,
    tibble(
      method = "pca",
      PC = 1:N_ORIG_SAMPLES,
      sdev = PCA(scaled_draws_pca[,this_range])$sdev,
    ),
    tibble(
      method = "corpcor",
      PC = 1:N_ORIG_SAMPLES,
      sdev = PCA(scaled_draws_corpcor[,this_range])$sdev,
    ),
    tibble(
      method = "wishart",
      PC = 1:N_ORIG_SAMPLES,
      sdev = PCA(scaled_draws_wishart[,this_range])$sdev,
    ),
    tibble(
      method = "indep",
      PC = 1:N_ORIG_SAMPLES,
      sdev = PCA(scaled_draws_indep[,this_range])$sdev,
    ),
    tibble(
      method = "spsimseq",
      PC = 1:N_ORIG_SAMPLES,
      sdev = PCA(scaled_draws_spsimseq[,this_range])$sdev,
    ),
    deparse.level=1
  )
}
pca_sdevs$method <- factor(pca_sdevs$method, levels=METHOD_ORDER, ordered=TRUE)
g5 <- ggplot() +
  facet_wrap(facets=vars(PC), labeller=label_both, ncol=6) +
  geom_hline(data=pca_sdevs |> filter(PC < N_ORIG_SAMPLES, method == "real"), aes(yintercept=sdev, color="real"), linewidth=1) +
  geom_boxplot(data=pca_sdevs |> filter(PC < N_ORIG_SAMPLES, method != "real"), aes(x=method, y=sdev, color=method)) +
  scale_color_manual(values=METHOD_COLORS, breaks=METHOD_ORDER, name="method") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

compare_to_real_plot <- (g1 | g2) / (g3 | g4) / (g5) + plot_annotation(tag_levels="a")

saveRDS(compare_to_real_plot, paste0("processed/compare_to_real_plot.", dataset, ".RDS"))
ggsave(paste0("results/compare_to_real_plot.", dataset, ".png"), width=10, height=10)
