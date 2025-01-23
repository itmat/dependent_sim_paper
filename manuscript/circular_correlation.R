library(readr)
library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

METHOD_ORDER <- c("indep", "pca", "wishart", "corpcor")
METHOD_COLORS <- RColorBrewer::brewer.pal(4, "Dark2")
METHOD_ORDER2 <- c("real", METHOD_ORDER) # Includes "real"
METHOD_COLORS2 <- c("black", METHOD_COLORS)

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
        for (i in (1:32)) {
            phaselist[i, 'ID'] <- strsplit(as.character(phaselist[i, 'ID']), "_")[[1]][1] |>
                substring(3)
        }
        phaselist$new_ID <- as.numeric(phaselist$ID)/12*pi

        corrs[[method]][j+1] <- get_cir_corr(phaselist$new_ID, phaselist$phase)
    }
}

# plot comparing circular correlation for simulated datasets with and without dependence
corrs_new <- pivot_longer(corrs, cols = 1:4, names_to = "method")
corrs_new$method <- factor(corrs_new$method, levels = METHOD_ORDER, ordered = TRUE)
circular_correlation_plot <- ggplot(corrs_new, aes(method, abs(value), color=method)) +
    geom_point(position=position_jitter(width=0.25, height=0)) +
    ylim(0,1) +
    ylab("absolute circular correlation") +
    scale_color_manual(values = METHOD_COLORS, breaks = METHOD_ORDER, name="method", guide="none")
circular_correlation_plot

# rainbow plots visualizing true_counts_counts_counts and estimated phases
phaseplots <- list()
for (method in c("indep", "pca", "wishart", "corpcor")) {
    phaselist <- read_csv(paste0("../processed/cyclops/",method,"/batch=", 5, "/cyclops_estimated_phaselist.csv"))
    for (i in (1:32)) {
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
        ggtitle(method) +
        coord_fixed() +
        labs(x = "cos(phase)", y = "sin(phase)")
}

# plot the number of eigengenes used by CYCLOPS
eigengenes <- list()
for (method in c("indep", "pca", "wishart", "corpcor")) {
    m <- list()
    for (batch in 0:19) {
        n <- read_csv(paste0("../processed/cyclops/",method,"/batch=",batch,"/cyclops_EigengeneExp.csv"), show_col_types = FALSE) |>
            nrow()
        m[[length(m)+1]] <- n
    }
    eigengenes[[length(eigengenes)+1]] <- tibble(
        method = method,
        num = unlist(m),
    )
}
sim_eigen <- bind_rows(eigengenes)
sim_eigen$method <- factor(sim_eigen$method, levels = METHOD_ORDER, ordered = TRUE)
real_eigens <- tibble(
    method = "real",
    num = read_csv("../processed/cyclops/real_data/cyclops_EigengeneExp.csv", show_col_types=FALSE) |> nrow(),
)

eigen_plot <- ggplot() +
    geom_hline(data = real_eigens, aes(yintercept = num, color=method), linetype = "dashed") +
    geom_point(data = sim_eigen, aes(method, num, color=method), position=position_jitter(width=0.25, height=0)) +
    ylab("number of eigengenes") +
    scale_color_manual(values = METHOD_COLORS2, breaks = METHOD_ORDER2, name="method")

cyclops_plot <- (circular_correlation_plot + eigen_plot) /
    (phaseplots$indep | phaseplots$pca | phaseplots$wishart | phaseplots$corpcor) +
    plot_annotation(tag_levels="a") +
    plot_layout(guides = 'collect', axis_titles = 'collect')
cyclops_plot
