library(readr)
library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)

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

# compute circular correlation between CYCLOP-estimated phases and true phases for all 40 datasets
corrs <- tibble(k0 = rep(0, 20), k2 = rep(0, 20))
for (j in (0:19)) {
    #k=0
    phaselist <- read_csv(paste0("processed/cyclops/k=0/batch=", j, "/cyclops_estimated_phaselist.csv"))
    for (i in (1:48)) {
        phaselist[i, 'ID'] <- strsplit(as.character(phaselist[i, 'ID']), "_")[[1]][1] |> 
            substring(3)
    }
    phaselist$new_ID <- as.numeric(phaselist$ID)/12*pi
    
    corrs$k0[j+1] <- get_cir_corr(phaselist$new_ID, phaselist$phase)
    
    #k=2
    phaselist <- read_csv(paste0("processed/cyclops/k=2/batch=", j, "/cyclops_estimated_phaselist.csv"))
    for (i in (1:48)) {
        phaselist[i, 'ID'] <- strsplit(as.character(phaselist[i, 'ID']), "_")[[1]][1] |> 
            substring(3)
    }
    phaselist$new_ID <- as.numeric(phaselist$ID)/12*pi
    
    corrs$k2[j+1] <- get_cir_corr(phaselist$new_ID, phaselist$phase)
}

# violin plot comparing circular correlation for simulated datasets with and without dependence
corrs_new <- pivot_longer(corrs, cols = 1:2, names_to = "k")
circular_correlation_violin <- ggplot(corrs_new, aes(k, abs(value))) +
    geom_violin() +
    ylab("absolute circular correlation")
circular_correlation_violin

# rainbow plots visualizing true and estimated phases
phaselist <- read_csv(paste0("processed/cyclops/k=0/batch=", 5, "/cyclops_estimated_phaselist.csv"))
for (i in (1:48)) {
    phaselist[i, 'ID'] <- strsplit(as.character(phaselist[i, 'ID']), "_")[[1]][1] |> 
        substring(3)
}
phaselist$new_ID <- as.numeric(phaselist$ID)/12*pi
phaselist$true_time <- as.numeric(phaselist$ID)

phaselist$sin_estimated_phase <- sin(phaselist$phase)
phaselist$cos_estimated_phase <- cos(phaselist$phase)
phaseplot_0 <- ggplot(phaselist,aes(cos_estimated_phase,sin_estimated_phase))+
    scale_color_gradientn(colors = rainbow(8), limits=c(0, 24)) +
    geom_point(aes(colour=true_time)) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    ggtitle("without dependence")
phaseplot_0
#ggsave(paste0("plots/k=0_batch=", j, "_scatter.png"),
#       scale = 1, width = 5, height = 4, units = "in")


phaselist <- read_csv(paste0("processed/cyclops/k=2/batch=", 6, "/cyclops_estimated_phaselist.csv"))
for (i in (1:48)) {
    phaselist[i, 'ID'] <- strsplit(as.character(phaselist[i, 'ID']), "_")[[1]][1] |> 
        substring(3)
}
phaselist$new_ID <- as.numeric(phaselist$ID)/12*pi
phaselist$true_time <- as.numeric(phaselist$ID)

phaselist$sin_estimated_phase <- sin(phaselist$phase)
phaselist$cos_estimated_phase <- cos(phaselist$phase)
phaseplot_2 <- ggplot(phaselist,aes(cos_estimated_phase,sin_estimated_phase))+
    scale_color_gradientn(colors = rainbow(8), limits=c(0, 24)) +
    geom_point(aes(colour=true_time)) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    ggtitle("with dependence")
phaseplot_2
#ggsave(paste0("plots/k=2_batch=", j, "_scatter.png"),
#       scale = 1, width = 5, height = 4, units = "in")
