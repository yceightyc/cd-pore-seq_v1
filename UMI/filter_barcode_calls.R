# Load libraries
library(dplyr)
library(tidyr)

# Import argument
args <- commandArgs(trailingOnly = TRUE)
dir <-args[1]
id <- args[2]
output_path <- paste0(dir, "/", id, "_filtered")

# Read base infor at each position in the barcode region
data_9_path <- paste0(dir,"/", id, "_9_merged")
data_9 <- as.data.frame(read.table(data_9_path))
data_10_path <- paste0(dir,"/", id, "_10_merged")
data_10 <- as.data.frame(read.table(data_10_path))
data_11_path <- paste0(dir,"/", id, "_11_merged")
data_11 <- as.data.frame(read.table(data_11_path))
data_12_path <- paste0(dir,"/", id, "_12_merged")
data_12 <- as.data.frame(read.table(data_12_path))
data_13_path <- paste0(dir,"/", id, "_13_merged")
data_13 <- as.data.frame(read.table(data_13_path))

# Merge all position bases
data_merge <- merge (data_9, data_10, by = "V2") %>%
    merge (data_11, by = "V2") %>%
    merge (data_12, by = "V2") %>%
    merge (data_13, by = "V2")

colnames(data_merge) <- c("read_name", "V2", "V3", "V4", "V5", "V6")
data_merge <- data_merge %>%
    unite("context", V2:V6, sep = "")

data_filtered <- data_merge %>%
    filter(context %in% c(".....", ",,,,,"))

write.table(data_filtered, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
