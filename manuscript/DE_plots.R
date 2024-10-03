library(readr)
library(tibble)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

mouse_plots = list()
for (method in c("indep", "pca", "wishart", "corpcor")) {
    fdr_df <- read_csv(paste0("../processed/DE/Mouse.Cortex.Male.",method,".fdr.csv"))
    
    mouse_plots[[method]] <- ggplot(fdr_df, aes(cutoff, mean, ymin=lower, ymax=upper)) +
        geom_ribbon(fill="grey") +
        geom_point() +
        geom_abline(slope=1, intercept = 0) +
        ggtitle(method) +
        xlab("reported FDR") +
        ylab(if_else(method=="indep","true FDR","")) +
        scale_x_log10() +
        scale_y_log10() +
        coord_cartesian(xlim = c(0.001, 1), ylim= c(0.001,1))
}
    
# graph <- ggplot(bind_rows(fdr_lists0, .id="data_frame"),
#                 aes(cutoff, fdr, colour=data_frame)) +
#     geom_point() +
#     scale_x_log10() +
#     scale_y_log10() +
#     labs(title = "k=0") +
#     geom_smooth(color="black", method="lm") +
#     geom_abline(slope=1, intercept = 0)

fly_plots = list()
for (method in c("indep", "pca", "wishart", "corpcor")) {
    fdr_df_fly <- read_csv(paste0("../processed/DE/Fly.WholeBody.Male.",method,".fdr.csv"))
    
    fly_plots[[method]] <- ggplot(fdr_df_fly, aes(cutoff, mean, ymin=lower, ymax=upper)) +
        geom_ribbon(fill="grey") +
        geom_point() +
        geom_abline(slope=1, intercept = 0) +
        ggtitle(method) +
        xlab("reported FDR") +
        ylab(if_else(method=="indep","true FDR","")) +
        scale_x_log10() +
        scale_y_log10() +
        coord_cartesian(xlim = c(0.001, 1), ylim= c(0.001,1))
}


DESeq2_fdr_plot <- ((fly_plots$indep | fly_plots$pca | fly_plots$wishart | fly_plots$corpcor) / 
    (mouse_plots$indep | mouse_plots$pca | mouse_plots$wishart | mouse_plots$corpcor)) + plot_annotation(tag_levels="a")
DESeq2_fdr_plot
