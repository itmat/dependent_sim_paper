library(readr)
library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

# formula for the circular correlation
get_cir_corr <- function(X, Y) {
    sum1 <- 0
    sum2 <- 0
    sum3 <- 0
    
    if (length(X) != length(Y)) {
        return("X and Y have different lengths")
    } else {
        n <- length(X)
    }
    
    for (i in (1:(n-1)) ) {
        for (j in ((i+1):n)) {
            sum1 <- sum1 + sin(X[i] - X[j])*sin(Y[i] - Y[j])
            sum2 <- sum2 + sin(X[i] - X[j])^2
            sum3 <- sum3 + sin(Y[i] - Y[j])^2
        }
    }
    
    rho <- sum1/(sqrt(sum2)*sqrt(sum3))
    
    return(rho)
        
}

# compute circular correlation between CYCLOP-estimated phases and true_counts_counts_counts phases for all 40 datasets
corrs <- tibble(indep = rep(0, 20), pca = rep(0, 20), wishart = rep(0, 20), corpcor = rep(0, 20))
for (method in c("indep", "pca", "wishart", "corpcor")) {
    for (j in (0:19)) {
        phaselist <- read_csv(paste0("../processed/cyclops/",method,"/batch=", j, "/cyclops_estimated_phaselist.csv"))
        for (i in (1:48)) {
            phaselist[i, 'ID'] <- strsplit(as.character(phaselist[i, 'ID']), "_")[[1]][1] |> 
                substring(3)
        }
        phaselist$new_ID <- as.numeric(phaselist$ID)/12*pi
        
        corrs[[method]][j+1] <- get_cir_corr(phaselist$new_ID, phaselist$phase)
    }
}

# violin plot comparing circular correlation for simulated datasets with and without dependence
corrs_new <- pivot_longer(corrs, cols = 1:4, names_to = "method")
circular_correlation_violin <- ggplot(corrs_new, aes(method, abs(value))) +
    geom_violin() +
    ylab("absolute circular correlation")
circular_correlation_violin

# rainbow plots visualizing true_counts_counts_counts and estimated phases
phaseplots <- list()
for (method in c("indep", "pca", "wishart", "corpcor")) {
    phaselist <- read_csv(paste0("../processed/cyclops/",method,"/batch=", 8, "/cyclops_estimated_phaselist.csv"))
    for (i in (1:48)) {
        phaselist[i, 'ID'] <- strsplit(as.character(phaselist[i, 'ID']), "_")[[1]][1] |> 
            substring(3)
    }
    phaselist$new_ID <- as.numeric(phaselist$ID)/12*pi
    phaselist$true_counts_counts_counts_time <- as.numeric(phaselist$ID)
    
    phaselist$sin_estimated_phase <- sin(phaselist$phase)
    phaselist$cos_estimated_phase <- cos(phaselist$phase)
    phaseplots[[method]] <- ggplot(phaselist,aes(cos_estimated_phase,sin_estimated_phase))+
        scale_color_gradientn(colors = rainbow(8), limits=c(0, 24), name="true time") +
        geom_point(aes(colour=true_counts_counts_counts_time)) +
        xlim(-1, 1) +
        ylim(-1, 1) +
        ggtitle(method)
}

# plot the number of eigengenes used by CYCLOPS
eigen <- read_delim("num_eigengene.txt", col_names=c("num", "file")) |> 
    mutate(method = str_extract(file, "./processed/cyclops/([a-z]+)", group=1))
sim_eigen <- eigen |> filter(method != "real")
real_eigens <- (eigen |> filter(method == "real"))$num
eigen_violin <- ggplot(sim_eigen, aes(method, num)) +
    geom_violin() +
    ylab("number of eigengenes") +
    geom_hline(yintercept = 13, linetype = "dashed", color = "red")

cyclops_plot <- circular_correlation_violin / 
    eigen_violin /
    (phaseplots$indep | phaseplots$pca | phaseplots$wishart | phaseplots$corpcor) + 
    plot_annotation(tag_levels="a")# + plot_layout(guides = 'collect')
cyclops_plot
