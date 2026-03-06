# Load libraries
library(dplyr)
library(tidyr)
library(purrr)

# Import argument
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
dir_BarcodeMatchingScore <- args[2] # this dir could be any directory of baseATpositions

# Read base infor at each position in the barcode region
channel_rnum_path <- paste0("./", "channel_read_number_all")
channel_rnum_data <- as.data.frame(read.table(channel_rnum_path))
colnames(channel_rnum_data) <- c("rn", "ch", "readnum")
data_110_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated_rev/", id, "_110_merged")
data_110 <- as.data.frame(read.table(data_110_path))
data_111_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated_rev/", id, "_111_merged")
data_111 <- as.data.frame(read.table(data_111_path))
data_112_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated_rev/", id, "_112_merged")
data_112 <- as.data.frame(read.table(data_112_path))
data_113_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated_rev/", id, "_113_merged")
data_113 <- as.data.frame(read.table(data_113_path))
data_114_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated_rev/", id, "_114_merged")
data_114 <- as.data.frame(read.table(data_114_path))


# Merge UMI bases
UMI_merge <- merge (data_110, data_111, by = "V2") %>%
        merge (data_112, by = "V2") %>%
        merge (data_113, by = "V2") %>%
        merge (data_114, by = "V2")

colnames(UMI_merge) <- c("V1", "V2", "V3", "V4", "V5", "V6")
UMI_merge <- UMI_merge %>%
        unite("UMI", V2:V6, sep = "")

colnames(UMI_merge) <-c("rn", "UMI")

#Merge channel and read number data
final_merged_withChannel <- merge (UMI_merge, channel_rnum_data, by = "rn")

colnames(final_merged_withChannel) <- c("readname", "UMI", "channel", "readnum")

# Save output in the respective directory
output_path <- paste0(dir_BarcodeMatchingScore,  id, "_UMI_bases_all_rev.tsv")
write.table(final_merged_withChannel, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
