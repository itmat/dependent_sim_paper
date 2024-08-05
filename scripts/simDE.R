library(readr)
library(tibble)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)

simDE_analysis <- function(filename, samples) {

    # read the simulated case and control data for the first 5 samples
    case_data <-
        readr::read_tsv(paste0("simulated_data/", filename, ".Case.txt"),
                      col_select = c(1, samples))
    colnames(case_data) <- c("ENSEMBL_ID", "case_S1", "case_S2", "case_S3", "case_S4", "case_S5")
    
    control_data <-
        readr::read_tsv(paste0("simulated_data/", filename, ".Control.txt"),
                        col_select = c(1, samples))
    colnames(control_data) <- c("ENSEMBL_ID", "ctrl_S1", "ctrl_S2", "ctrl_S3", "ctrl_S4", "ctrl_S5")
    
    
    # combine case and control data and specify sample info
    sim_data <- 
        merge(case_data, control_data, by="ENSEMBL_ID")
    
    sample_info <- 
        tibble(sample_id = c(colnames(case_data)[2:6], colnames(control_data)[2:6]),
               group = c(rep("case", 5), rep("control", 5)))
    
    # sim_data as a matrix
    sim_matrix <-
        sim_data |>
        select(-ENSEMBL_ID) |>
        as.matrix()
    
    rownames(sim_matrix) <- sim_data$ENSEMBL_ID
    
    # cluster plot
    clustered_data <-
        t(sim_matrix) |>
        dist() |>
        hclust()
    #plot(clustered_data, main = paste0("Clustering Raw RNA-Seq Quant Data ", filename))
    
    sample_metadata <- 
        sample_info |> column_to_rownames(var = "sample_id")
    
    # create DESeq dataset from the simulated data comparing case and control
    deseq2_data <-
        DESeq2::DESeqDataSetFromMatrix(countData = sim_matrix,
                                       colData = sample_metadata |> 
                                           mutate(group = factor(group, c("case", "control"))),
                                       design = ~group)
    
    # filter genes with high number of zero counts
    keep_genes <- rowSums(counts(deseq2_data) >= 10) >= 3
    
    filtered_deseq2_data <- deseq2_data[keep_genes,]
    
    filtered_deseq2_data <- DESeq(filtered_deseq2_data)
    
    # run DESeq
    deseq2_results <-
        results(filtered_deseq2_data, contrast = c("group", "case", "control"),
                tidy = TRUE, independentFiltering = FALSE) |>
        as_tibble() |>
        rename(ENSEMBL_ID = row, qvalue = padj) |>
        mutate(Comparison = "case_vs_control") |> 
        arrange(pvalue, qvalue)
        
    # create list of differentially expressed genes
    DEgenes <- 
        deseq2_results |>
        select(ENSEMBL_ID, Comparison, log2FoldChange, pvalue, qvalue) |>
        pivot_wider(names_from = Comparison, 
                    values_from = c(log2FoldChange, qvalue),
                    names_sep = ".") |> 
        filter(qvalue.case_vs_control < 0.1)
    
    deseq2_results |> 
        summary(alpha = 0.1)
    
    filtered_results <- 
        deseq2_results |> 
        filter(qvalue < 0.1)
    
    # read true values for the case and control datasets
    case_true_vals <-
        readr::read_tsv(paste0("simulated_data/", filename, ".Case.true_values.txt"),
                        col_select = c("gene_id", "mean"))
    colnames(case_true_vals) <- c("ENSEMBL_ID", "case_mean")
    
    control_true_vals <-
        readr::read_tsv(paste0("simulated_data/", filename, ".Control.true_values.txt"),
                        col_select = c("gene_id", "mean"))
    colnames(control_true_vals) <- c("ENSEMBL_ID", "control_mean")
    
    # calculate log2 fold change for the true values
    true_vals <- 
        merge(case_true_vals, control_true_vals, by="ENSEMBL_ID")
    
    true_vals$true_Log2FoldChange <- log2(true_vals$case_mean/true_vals$control_mean)
    
    # scatter plot comparing log2 fold changes of simulated data and true values
    data_scatter <- 
        merge(deseq2_results, true_vals, by="ENSEMBL_ID") |> 
        select(ENSEMBL_ID, baseMean, pvalue, qvalue, log2FoldChange, true_Log2FoldChange)
    
    #plot(data_scatter$log2FoldChange, data_scatter$true_Log2FoldChange,
    #     main="log2FC_sim_vs_true", xlab="sim", ylab="true")
    
    # sort the genes by their pvalues 
    sorted_data <- 
        data_scatter[order(data_scatter$qvalue),]
    
    # calculate the cumulative sum of truely DE genes
    sorted_data$true_DE <- (sorted_data$true_Log2FoldChange != 0)
    
    sorted_data$cumsum <- cumsum(sorted_data$true_DE)
    
    # for each qvalue cutoff, calculate FDR
    fdr_list <- data.frame(cutoff=10^((-75:-1)*0.04), fdr=rep(1, 75))
    
    for (q in (-75:-1)) {
        filtered_genes <- filter(sorted_data, sorted_data$qvalue < 10^(q*0.04))
        fdr_list$fdr[q+76] <- 1 - max(filtered_genes$cumsum)/nrow(filtered_genes)
    }
    
    fdr_list

}


# run DE analysis on the mouse dataset for k=0 and get fdr values

fdr_lists0 = list()
for (s in (1:20)) {
    fdr_list <- simDE_analysis("Mouse.Cortex.Male.k=0", ((5*s-3):(5*s+1)))
    fdr_lists0[[s]] <- fdr_list
}

# plot mean and confidence interval
fdr_df0 <- bind_rows(fdr_lists0, .id="data_frame") |> 
    group_by(cutoff) |> 
    summarise(mean = mean(fdr), lower = t.test(fdr)$conf.int[1], upper = t.test(fdr)$conf.int[2])

write_csv(fdr_df0, "processed/DE/Mouse.Cortex.Male.k=0.fdr.csv")


# plot fdr graphs for k=2
fdr_lists2 = list()
for (s in (1:20)) {
    #print(((5*s-3):(5*s+1)))
    fdr_list <- simDE_analysis("Mouse.Cortex.Male.k=2", ((5*s-3):(5*s+1)))
    fdr_lists2[[s]] <- fdr_list
}

# plot mean and confidence interval
fdr_df2 <- bind_rows(fdr_lists2, .id="data_frame") |> 
    group_by(cutoff) |> 
    summarise(mean = mean(fdr), lower = t.test(fdr)$conf.int[1], upper = t.test(fdr)$conf.int[2])

write_csv(fdr_df2, "processed/DE/Mouse.Cortex.Male.k=2.fdr.csv")


# run DE analysis on the fly dataset for k=0 and get fdr values

fdr_lists0_fly = list()
for (s in (1:20)) {
    fdr_list <- simDE_analysis("Fly.WholeBody.Male.k=0", ((5*s-3):(5*s+1)))[["fdr_list"]]
    fdr_lists0_fly[[s]] <- fdr_list
}


# plot mean and confidence interval
fdr_df0_fly <- bind_rows(fdr_lists0_fly, .id="data_frame") |> 
    group_by(cutoff) |> 
    summarise(mean = mean(fdr), lower = t.test(fdr)$conf.int[1], upper = t.test(fdr)$conf.int[2])

write_csv(fdr_df0_fly, "processed/DE/Fly.WholeBody.Male.k=0.fdr.csv")


# plot fdr graphs for k=2
fdr_lists2_fly = list()
for (s in (1:20)) {
    fdr_list <- simDE_analysis("Fly.WholeBody.Male.k=2", ((5*s-3):(5*s+1)))[["fdr_list"]]
    fdr_lists2_fly[[s]] <- fdr_list
}

# plot mean and confidence interval
fdr_df2_fly <- bind_rows(fdr_lists2_fly, .id="data_frame") |> 
    group_by(cutoff) |> 
    summarise(mean = mean(fdr), lower = t.test(fdr)$conf.int[1], upper = t.test(fdr)$conf.int[2])

write_csv(fdr_df2_fly, "processed/DE/Fly.WholeBody.Male.k=2.fdr.csv")


# make_qqplot <- function(filename, k) {
#     # distribution of p values
#     pval_lists0 = list()
#     for (s in (1:20)) {
#         pvalue_df <- simDE_analysis(filename, ((5*s-3):(5*s+1)))[["data_scatter"]] |> 
#             filter(true_Log2FoldChange == 0) |> 
#             select(ENSEMBL_ID, pvalue)
#         pval_lists0[[s]] <- pvalue_df
#     }
#     
#     graph <- ggplot(bind_rows(pval_lists0, .id="data_frame"),
#                     aes(sample=pvalue, colour=data_frame)) +
#         stat_qq(distribution=stats::qunif) + 
#         geom_abline(slope=1,intercept=0) + 
#         scale_x_log10() + scale_y_log10() +
#         labs(title = paste0("k=",k))
#     
#     return(graph)
# }
# 
# 
# mouse_qq_0 <- make_qqplot("Mouse.Cortex.Male.k=0", k = 0)
# 
# mouse_qq_2 <- make_qqplot("Mouse.Cortex.Male.k=2", k = 2)


