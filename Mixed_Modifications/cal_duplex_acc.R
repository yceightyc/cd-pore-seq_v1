# Load libraries
library(dplyr)
library(tidyr)

# Import argument
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
dir_BarcodeMatchingScore <- args[2] # this dir could be any directory of baseATpositions

# umiINreverse
data_132_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_132_merged")
data_132 <- as.data.frame(read.table(data_132_path))
data_133_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_133_merged")
data_133 <- as.data.frame(read.table(data_133_path))
data_134_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_134_merged")
data_134 <- as.data.frame(read.table(data_134_path))

# umiINforward
data_65_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_65_merged")
data_65 <- as.data.frame(read.table(data_65_path))
data_66_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_66_merged")
data_66 <- as.data.frame(read.table(data_66_path))
data_67_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_67_merged")
data_67 <- as.data.frame(read.table(data_67_path))


# Merge UMI bases in the reverse strand
UMI_merge <- merge (data_132, data_133, by = "V2") %>%
    merge (data_134, by = "V2") 
colnames(UMI_merge) <- c("V1", "V2", "V3", "V4")
UMI_merge <- UMI_merge %>% 
    unite("UMI", V2:V4, sep = "")

# Merge UMI bases in the forward strand
UMI_merge_for <- merge (data_65, data_66, by = "V2") %>%
    merge (data_67, by = "V2")
colnames(UMI_merge_for) <- c("V1", "V2", "V3", "V4")
UMI_merge_for <- UMI_merge_for %>%
    unite("UMI", V2:V4, sep = "") %>%
    mutate(UMI = toupper(UMI))

reverse_complement <- function(dna_seq) {
  # Define the complement bases
  complement_map <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C", "*" = "*") 
  
  # Split the sequence into individual bases
  bases <- strsplit(dna_seq, "")[[1]]
  bases <- toupper(bases)
  # Get the complement
  complement_bases <- complement_map[bases]
  
  # Reverse the complement
  reverse_complement_seq <- paste(rev(complement_bases), collapse = "")
  
  return(reverse_complement_seq)
}

# Merge UMIs and permutations
final_merged <- merge (UMI_merge, UMI_merge_for, by = "V1")
colnames (final_merged) <- c("rn", "UMI_rev", "UMI_for")

final_merged_for_included <- final_merged %>%
	mutate(UMI_for_cali = sapply(UMI_rev,reverse_complement), 
	       pos1 = if_else(substr(UMI_for, 1, 1) == substr(UMI_for_cali, 1, 1), 1, 0),
	       pos2 = if_else(substr(UMI_for, 2, 2) == substr(UMI_for_cali, 2, 2), 1, 0),
	       pos3 = if_else(substr(UMI_for, 3, 3) == substr(UMI_for_cali, 3, 3), 1, 0),
	       UMI_all = if_else(UMI_for == UMI_for_cali, 1, 0))

# Save output in the respective directory
output_path <- paste0("UMI_permut_bases_out/",  id, "_raw_UMI_matching_all.tsv")
write.table(final_merged_for_included, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
