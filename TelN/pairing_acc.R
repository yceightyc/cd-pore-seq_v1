# Load libraries
library(dplyr)
library(tidyr)

# Import argument
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
dir_BarcodeMatchingScore <- args[2] # this dir could be any directory of baseATpositions

# Forward barcodes
data_45_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_45_merged")
Base45 <- as.data.frame(read.table(data_45_path))
colnames(Base45) <- c("base45", "rname_for")
Base45 <- Base45 %>%
    filter(base45!= "*")

data_46_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_46_merged")
Base46 <- as.data.frame(read.table(data_46_path))
colnames(Base46) <- c("base46", "rname_for")
Base46 <- Base46 %>%
    filter(base46!= "*")

data_47_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_47_merged")
Base47 <- as.data.frame(read.table(data_47_path))
colnames(Base47) <- c("base47", "rname_for")
Base47 <- Base47 %>%
    filter(base47!= "*")

data_48_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_48_merged")
Base48 <- as.data.frame(read.table(data_48_path))
colnames(Base48) <- c("base48", "rname_for")
Base48 <- Base48 %>%
    filter(base48!= "*")

data_49_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_49_merged")
Base49 <- as.data.frame(read.table(data_49_path))
colnames(Base49) <- c("base49", "rname_for")
Base49 <- Base49 %>%
    filter(base49!= "*")

# Reverse barcodes
data_222_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_222_merged")
Base222 <- as.data.frame(read.table(data_222_path))
colnames(Base222) <- c("base222", "rname_for")
Base222 <- Base222 %>%
    filter(base222!= "*")

data_223_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_223_merged")
Base223 <- as.data.frame(read.table(data_223_path))
colnames(Base223) <- c("base223", "rname_for")
Base223 <- Base223 %>%
    filter(base223!= "*")

data_224_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_224_merged")
Base224 <- as.data.frame(read.table(data_224_path))
colnames(Base224) <- c("base224", "rname_for")
Base224 <- Base224 %>%
    filter(base224!= "*")

data_225_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_225_merged")
Base225 <- as.data.frame(read.table(data_225_path))
colnames(Base225) <- c("base225", "rname_for")
Base225 <- Base225 %>%
    filter(base225!= "*")

data_226_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_226_merged")
Base226 <- as.data.frame(read.table(data_226_path))
colnames(Base226) <- c("base226", "rname_for")
Base226 <- Base226 %>%
    filter(base226!= "*")

# umiINreverse
data_168_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_168_merged")
data_168 <- as.data.frame(read.table(data_168_path))
data_169_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_169_merged")
data_169 <- as.data.frame(read.table(data_169_path))
data_170_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_170_merged")
data_170 <- as.data.frame(read.table(data_170_path))

# umiINforward
data_101_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_101_merged")
data_101 <- as.data.frame(read.table(data_101_path))
data_102_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_102_merged")
data_102 <- as.data.frame(read.table(data_102_path))
data_103_path <- paste0(dir_BarcodeMatchingScore, "/merged_mutated/", id, "_103_merged")
data_103 <- as.data.frame(read.table(data_103_path))


# Merge UMI bases in the reverse strand
UMI_merge <- merge (data_168, data_169, by = "V2") %>%
    merge (data_170, by = "V2") 
colnames(UMI_merge) <- c("V1", "V2", "V3", "V4")
UMI_merge <- UMI_merge %>% 
    unite("UMI", V2:V4, sep = "")

# Merge UMI bases in the forward strand
UMI_merge_for <- merge (data_101, data_102, by = "V2") %>%
    merge (data_103, by = "V2")
colnames(UMI_merge_for) <- c("V1", "V2", "V3", "V4")
UMI_merge_for <- UMI_merge_for %>%
    unite("UMI", V2:V4, sep = "")

# Reverse complementary transformation
reverse_complement <- function(dna_seq) {
  # Define the complement bases
  complement_map <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C", "*" = "*") 
  
  # Split the sequence into individual bases
  bases <- strsplit(dna_seq, "")[[1]]
  
  # Get the complement
  complement_bases <- complement_map[bases]
  
  # Reverse the complement
  reverse_complement_seq <- paste(rev(complement_bases), collapse = "")
  
  return(reverse_complement_seq)
}

# Merge UMIs 
final_merged <- merge (UMI_merge, UMI_merge_for, by = "V1")
colnames (final_merged) <- c("rn", "UMI_rev", "UMI_for")

final_merged_for_included <- final_merged %>%
	mutate(UMI_for_cali = sapply(UMI_rev,reverse_complement), 
	       pos1 = if_else(substr(UMI_for, 1, 1) == substr(UMI_for_cali, 1, 1), 1, 0),
	       pos2 = if_else(substr(UMI_for, 2, 2) == substr(UMI_for_cali, 2, 2), 1, 0),
	       pos3 = if_else(substr(UMI_for, 3, 3) == substr(UMI_for_cali, 3, 3), 1, 0),
	       UMI_all = if_else(UMI_for == UMI_for_cali, 1, 0))

# Filter reads with all accurate basecalls in the forward 5nt barcodes 
barcode_joined <- merge(Base45, Base46, by = "rname_for") %>%
    merge(Base47, by = "rname_for") %>%
    merge(Base48, by = "rname_for") %>%
    merge(Base49, by = "rname_for") %>%
    filter(base45==".") %>%
    filter(base46==".") %>%
    filter(base47==".") %>%
    filter(base48==".") %>%
    filter(base49==".") 

# Merge reverse barcode basecalls
barcode_rev_joined <- merge(Base222, Base223, by = "rname_for") %>%
    merge(Base224, by = "rname_for") %>%
    merge(Base225, by = "rname_for") %>%
    merge(Base226, by = "rname_for") %>%
    mutate(revbar = paste(base222, base223, base224, base225, base226, sep = ""))
revbar_filtered <- merge(barcode_joined, barcode_rev_joined, by = "rname_for") 

# Save UMI matching and pairing accuracy output
output_path_umi <- paste0("./UMI_permut_bases_out/",  id, "_raw_UMI_matching_all.tsv")
write.table(final_merged_for_included, file = output_path_umi, sep = "\t", row.names = FALSE, quote = FALSE)

output_path_acc <- paste0("./data/TelN_pairing_acc.tsv")
write.table(revbar_filtered, file = output_path_acc, sep = "\t", row.names = FALSE, quote = FALSE)