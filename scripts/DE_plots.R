library(readr)
library(tibble)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)

fdr_df0 <- read_csv("processed/DE/Mouse.Cortex.Male.k=0.fdr.csv")

mouse_cortex_fdr_0 <- ggplot(fdr_df0, aes(cutoff, mean, ymin=lower, ymax=upper)) +
    geom_ribbon(fill="grey") +
    geom_point() +
    geom_abline(slope=1, intercept = 0) +
    labs(title = "k=0") +
    xlab("reported FDR") +
    ylab("true FDR") +
    scale_x_log10() +
    scale_y_log10() +
    coord_cartesian(xlim = c(0.001, 1), ylim= c(0.001,1))

# graph <- ggplot(bind_rows(fdr_lists0, .id="data_frame"),
#                 aes(cutoff, fdr, colour=data_frame)) +
#     geom_point() +
#     scale_x_log10() +
#     scale_y_log10() +
#     labs(title = "k=0") +
#     geom_smooth(color="black", method="lm") +
#     geom_abline(slope=1, intercept = 0)

fdr_df2 <- read_csv("processed/DE/Mouse.Cortex.Male.k=2.fdr.csv")

mouse_cortex_fdr_2 <- ggplot(fdr_df2, aes(cutoff, mean, ymin=lower, ymax=upper)) +
    geom_ribbon(fill="grey") +
    geom_point() +
    geom_abline(slope=1, intercept = 0) +
    labs(title = "k=0") +
    xlab("reported FDR") +
    ylab("true FDR") +
    scale_x_log10() +
    scale_y_log10() +
    coord_cartesian(xlim = c(0.001, 1), ylim= c(0.001,1))


fdr_df0_fly <- read_csv("processed/DE/Fly.WholeBody.Male.k=0.fdr.csv")

fly_fdr_0 <- ggplot(fdr_df0_fly, aes(cutoff, mean, ymin=lower, ymax=upper)) +
    geom_ribbon(fill="grey") +
    geom_point() +
    geom_abline(slope=1, intercept = 0) +
    labs(title = "k=0") +
    xlab("reported FDR") +
    ylab("true FDR") +
    scale_x_log10() +
    scale_y_log10() +
    coord_cartesian(xlim = c(0.001, 1), ylim= c(0.001,1))


fdr_df2_fly <- read_csv("processed/DE/Fly.WholeBody.Male.k=2.fdr.csv")

fly_fdr_2 <- ggplot(fdr_df2_fly, aes(cutoff, mean, ymin=lower, ymax=upper)) +
    geom_ribbon(fill="grey") +
    geom_point() +
    geom_abline(slope=1, intercept = 0) +
    labs(title = "k=0") +
    xlab("reported FDR") +
    ylab("true FDR") +
    scale_x_log10() +
    scale_y_log10() +
    coord_cartesian(xlim = c(0.001, 1), ylim= c(0.001,1))

