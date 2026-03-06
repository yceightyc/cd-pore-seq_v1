# Load libraries
library(dplyr)
library(tidyr)

# Import argument
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
dir_BarcodeMatchingScore <- args[2] # this dir could be any directory of baseATpositions
output_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_all.tsv")

# Read base infor at each position in the barcode region
channel_rnum_path <- paste0("./", dir_BarcodeMatchingScore, "/", "channel_read_number_all")
channel_rnum_data <- as.data.frame(read.table(channel_rnum_path))
colnames(channel_rnum_data) <- c("rn", "ch", "readnum")
data_68_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_68_merged")
data_68 <- as.data.frame(read.table(data_68_path))
data_69_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_69_merged")
data_69 <- as.data.frame(read.table(data_69_path))
data_70_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_70_merged")
data_70 <- as.data.frame(read.table(data_70_path))
data_71_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_71_merged")
data_71 <- as.data.frame(read.table(data_71_path))
data_72_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_72_merged")
data_72 <- as.data.frame(read.table(data_72_path))
data_89_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_89_merged")
data_89 <- as.data.frame(read.table(data_89_path))
data_90_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_90_merged")
data_90 <- as.data.frame(read.table(data_90_path))
data_91_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_91_merged")
data_91 <- as.data.frame(read.table(data_91_path))
data_92_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_92_merged")
data_92 <- as.data.frame(read.table(data_92_path))
data_93_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_93_merged")
data_93 <- as.data.frame(read.table(data_93_path))
data_94_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_94_merged")
data_94 <- as.data.frame(read.table(data_94_path))

# Merge UMI bases
UMI_merge <- merge (data_68, data_69, by = "V2") %>%
    merge (data_70, by = "V2") %>%
    merge (data_71, by = "V2") %>%
    merge (data_72, by = "V2") %>%
    merge (data_89, by = "V2") %>%
    merge (data_90, by = "V2") %>%
    merge (data_91, by = "V2") %>%
    merge (data_92, by = "V2") %>%
    merge (data_93, by = "V2") %>%
    merge (data_94, by = "V2") 

colnames(UMI_merge) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")
UMI_merge <- UMI_merge %>% 
    unite("UMI", V2:V12, sep = "")

colnames (UMI_merge) <- c("rn", "UMI")

#Merge channel and read number data
final_merged_withChannel <- merge (UMI_merge, channel_rnum_data, by = "rn")

colnames(final_merged_withChannel) <- c("readname", "UMI", "channel", "readnum")

# Save output in the respective directory
write.table(final_merged_withChannel, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
