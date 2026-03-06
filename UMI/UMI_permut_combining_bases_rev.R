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
data_27_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_27_merged")
data_27 <- as.data.frame(read.table(data_27_path))
data_28_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_28_merged")
data_28 <- as.data.frame(read.table(data_28_path))
data_29_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_29_merged")
data_29 <- as.data.frame(read.table(data_29_path))
data_30_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_30_merged")
data_30 <- as.data.frame(read.table(data_30_path))
data_31_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_31_merged")
data_31 <- as.data.frame(read.table(data_31_path))
data_32_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_32_merged")
data_32 <- as.data.frame(read.table(data_32_path))
data_33_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_33_merged")
data_33 <- as.data.frame(read.table(data_33_path))
data_37_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_37_merged")
data_37 <- as.data.frame(read.table(data_37_path))
data_38_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_38_merged")
data_38 <- as.data.frame(read.table(data_38_path))
data_39_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_39_merged")
data_39 <- as.data.frame(read.table(data_39_path))
data_40_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_40_merged")
data_40 <- as.data.frame(read.table(data_40_path))
data_41_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_41_merged")
data_41 <- as.data.frame(read.table(data_41_path))
data_42_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_42_merged")
data_42 <- as.data.frame(read.table(data_42_path))
data_43_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_43_merged")
data_43 <- as.data.frame(read.table(data_43_path))
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

# Merge Permutation bases
permut_up_merge <- merge (data_27, data_28, by = "V2") %>%
    merge (data_29, by = "V2") %>%
    merge (data_30, by = "V2") %>%
    merge (data_31, by = "V2") %>%
    merge (data_32, by = "V2") %>%
    merge (data_33, by = "V2") 

colnames(permut_up_merge) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")
permut_up_merge <- permut_up_merge %>%
    unite("permut_up", V2:V8, sep = "")

permut_down_merge <- merge (data_37, data_38, by = "V2") %>%
    merge (data_39, by = "V2") %>%
    merge (data_40, by = "V2") %>%
    merge (data_41, by = "V2") %>%
    merge (data_42, by = "V2") %>%
    merge (data_43, by = "V2") 

colnames(permut_down_merge) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")
permut_down_merge <- permut_down_merge %>%
    unite("permut_down", V2:V8, sep = "")

# Merge UMI and permutations
final_merged <- merge (UMI_merge, permut_up_merge, by = "V1") %>%
    merge (permut_down_merge, by = "V1")

colnames (final_merged) <- c("rn", "UMI", "permut_up", "permut_down")

#Merge channel and read number data
final_merged_withChannel <- merge (final_merged, channel_rnum_data, by = "rn")

colnames(final_merged_withChannel) <- c("readname", "UMI", "permut_up", "permut_down", "channel", "readnum")

# Save output in the respective directory
write.table(final_merged_withChannel, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
