library(dplyr)
library(tidyr)
library(stringdist)

# Import argument
args <- commandArgs(trailingOnly = TRUE)
for_file <- args[1]
rev_file <- args[2]
output_path <- args[3]

# Read files
for_path <- paste0("./", for_file)
for_reads_all <- as.data.frame(read.table(for_path))
colnames(for_reads_all) <- c("contig_for", "rname_for", "UMI", "channel", "read_num_for")

rev_path <- paste0("./", rev_file)
rev_reads <- as.data.frame(read.table(rev_path))
colnames(rev_reads) <- c("contig_rev","rname_rev","UMI","permut_up","permut_down", "channel","read_num_rev")

# Convert UMI column case in reverse reads
rev_reads$UMI <- toupper(rev_reads$UMI)

# Custom function to count how many letters match between two UMIs
count_matching_letters <- function(umi1, umi2) {
  sum(strsplit(umi1, "")[[1]] == strsplit(umi2, "")[[1]])
}

all_stats <- list()

# Run 10 permutations:
for (i in 1:10) {
 cat("Running iteration", i, "\n")
  
  # Sample 1000 reads
  for_reads <- for_reads_all[sample(nrow(for_reads_all), 1000), ]

# Perform the inner join first, without filtering
merged_data <- for_reads %>%
  inner_join(rev_reads, by = "channel", suffix = c("_for", "_rev"))

# Add a new column for the similarity score 
merged_data <- merged_data %>%
  rowwise() %>%
  mutate(UMI_simi = count_matching_letters(UMI_for, UMI_rev)) %>%
  filter(UMI_simi >=5) %>%
  ungroup()

# Calculate accuracy under different UMI similarity thresholds
stats <- merged_data %>%
  group_by(UMI_simi) %>%
  summarise(
    total = n(),
    contig_match = sum(contig_for == contig_rev),
    match_rate = contig_match / total
  ) %>%
  arrange(desc(UMI_simi)) %>%
  mutate(iteration = i)

all_stats[[i]] <- stats
}
final_stats <- bind_rows(all_stats)

# Save file
write.table(final_stats, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)

