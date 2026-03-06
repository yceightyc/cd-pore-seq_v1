library(dplyr)
library(tidyr)

# Import argument
args <- commandArgs(trailingOnly = TRUE)
for_file <- args[1]
rev_file <- args[2]
out_path <- args[3]

# Read files
for_path <- paste0("./", for_file)
for_reads <- as.data.frame(read.table(for_path))
colnames(for_reads) <- c("rname_for", "UMI", "channel", "read_num_for")

rev_path <- paste0("./", rev_file)
rev_reads <- as.data.frame(read.table(rev_path))
colnames(rev_reads) <- c("rname_rev","UMI","channel","read_num_rev")

# Convert UMI column case in reverse reads
rev_reads$UMI <- toupper(rev_reads$UMI)


# Inner join
merged_data <- merge (for_reads, rev_reads, by = c("UMI", "channel"))

# Save file
write.table(merged_data, file = out_path, sep = "\t", row.names = FALSE, quote = FALSE)

