library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(tidyverse)

raw_data <- read_csv("data/GSE151565_Cortex-counts.csv.gz")

normalized_data <- raw_data

totals <- summarise_all(raw_data, ~if(is.numeric(.)) sum(.) else "Total")

for (i in (2:ncol(raw_data))) {
    column <- colnames(raw_data)[i]
    total <- as.numeric(30000000/totals[i])
    normalized_data[[column]] <- 
        as.numeric(as.character(normalized_data[[column]]))*total
}

# add row specifying time points
normalized_data <- rbind(rep("", ncol(normalized_data)), normalized_data)
normalized_data[1, 1] <- "time_point"
for (i in (2 : ncol(normalized_data))) {
    normalized_data[1, i] <- unlist(strsplit(colnames(normalized_data)[i], split = "_", fixed = TRUE))[7]
}

time_points <- as.vector(normalized_data[1, ])

raw_data <- rbind(time_points, raw_data)

# calculate normalized data sorted by time and mean per time point
data_per_time = list()
mean_per_time = list()
for (i in (0:7)) {
    data <- cbind(normalized_data$EnsemblID, normalized_data[, time_points==paste0('ZT',i*3)])
    data <- data[-1, ]
    data_per_time[[i+1]] <- data
    data$mean <- data[, -1] |> 
        mutate_all(as.numeric) |> 
        rowMeans()
    mean_per_time[[i+1]] <- data |> 
        select(`normalized_data$EnsemblID`, mean)
}

mean_df <- mean_per_time |> 
    purrr::reduce(full_join, by="normalized_data$EnsemblID")
colnames(mean_df) <- c("EnsemblID", 'ZT0', 'ZT3', 'ZT6', 'ZT9',
                              'ZT12', 'ZT15', 'ZT18', 'ZT21')
write_csv(mean_df, "processed/mean_per_time_Cortex.csv")

processed_data <- data_per_time |> 
    reduce(full_join, by="normalized_data$EnsemblID")
write_csv(processed_data, "processed/Cortex_normalized_time_series_data.csv")

data0 <- cbind(raw_data$EnsemblID, raw_data[, time_points==paste0('ZT',0)])
data0 <- data0[-1, ]
write_csv(data0, "processed/Cortex_ZT0-counts.csv")
