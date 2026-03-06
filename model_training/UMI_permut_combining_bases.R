# Load libraries
library(dplyr)
library(tidyr)

# Import argument
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
dir_BarcodeMatchingScore <- args[2] # this dir could be any directory of baseATpositions

# Permut_upINreverse
data_192_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_192_merged")
data_192 <- as.data.frame(read.table(data_192_path))
data_193_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_193_merged")
data_193 <- as.data.frame(read.table(data_193_path))
data_194_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_194_merged")
data_194 <- as.data.frame(read.table(data_194_path))
data_195_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_195_merged")
data_195 <- as.data.frame(read.table(data_195_path))
data_196_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_196_merged")
data_196 <- as.data.frame(read.table(data_196_path))
data_197_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_197_merged")
data_197 <- as.data.frame(read.table(data_197_path))
data_198_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_198_merged")
data_198 <- as.data.frame(read.table(data_198_path))

# Permut_downINreverse
data_202_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_202_merged")
data_202 <- as.data.frame(read.table(data_202_path))
data_203_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_203_merged")
data_203 <- as.data.frame(read.table(data_203_path))
data_204_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_204_merged")
data_204 <- as.data.frame(read.table(data_204_path))
data_205_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_205_merged")
data_205 <- as.data.frame(read.table(data_205_path))
data_206_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_206_merged")
data_206 <- as.data.frame(read.table(data_206_path))
data_207_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_207_merged")
data_207 <- as.data.frame(read.table(data_207_path))
data_208_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_208_merged")
data_208 <- as.data.frame(read.table(data_208_path))

# umiINreverse
data_168_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_168_merged")
data_168 <- as.data.frame(read.table(data_168_path))
data_169_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_169_merged")
data_169 <- as.data.frame(read.table(data_169_path))
data_170_path <- paste0(dir_BarcodeMatchingScore, "/", id, "_170_merged")
data_170 <- as.data.frame(read.table(data_170_path))

# Merge UMI bases in the reverse strand
UMI_merge <- merge (data_168, data_169, by = "V2") %>%
    merge (data_170, by = "V2") 
colnames(UMI_merge) <- c("V1", "V2", "V3", "V4")
UMI_merge <- UMI_merge %>% 
    unite("UMI", V2:V4, sep = "")

# Merge randomers in the reverse strand
permut_up_merge <- merge (data_192, data_193, by = "V2") %>%
    merge (data_194, by = "V2") %>%
    merge (data_195, by = "V2") %>%
    merge (data_196, by = "V2") %>%
    merge (data_197, by = "V2") %>%
    merge (data_198, by = "V2") 

colnames(permut_up_merge) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")
permut_up_merge <- permut_up_merge %>%
    unite("permut_up", V2:V8, sep = "")

permut_down_merge <- merge (data_202, data_203, by = "V2") %>%
    merge (data_204, by = "V2") %>%
    merge (data_205, by = "V2") %>%
    merge (data_206, by = "V2") %>%
    merge (data_207, by = "V2") %>%
    merge (data_208, by = "V2") 

colnames(permut_down_merge) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")
permut_down_merge <- permut_down_merge %>%
    unite("permut_down", V2:V8, sep = "")

# Compute randomers and UMIs in the forward strand:
## Define a function to compute reverse complement
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

# Merge UMIs and permutations
final_merged <- merge (UMI_merge, permut_up_merge, by = "V1") %>%
    merge (permut_down_merge, by = "V1")
colnames (final_merged) <- c("rn", "UMI_rev", "permut_up_rev", "permut_down_rev")

final_merged_for_included <- final_merged %>%
	mutate(permut_up_for = sapply(permut_down_rev,reverse_complement),
	permut_down_for = sapply(permut_up_rev,reverse_complement),
	UMI_for = sapply(UMI_rev,reverse_complement))

# Save output in the respective directory
output_path <- paste0(dir_BarcodeMatchingScore, "/",  id, "_UMI_permut_bases_all.tsv")
write.table(final_merged_for_included, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
