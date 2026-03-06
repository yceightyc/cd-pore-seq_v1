library(dplyr)
library(tidyr)

# Import argument
args <- commandArgs(trailingOnly = TRUE)
rn_file <- args[1]
umi_file <- args[2]
output_path <- args[3]

# Read files
rn_path <- paste0("./", rn_file)
rn <- as.data.frame(read.table(rn_path))
colnames(rn) <- c("rname", "barcode")
rn <-rn[, "rname"]

umi_path <- paste0("./", umi_file)
umi <- as.data.frame(read.table(umi_path))
colnames(umi) <- c("contig","rname","UMI", "channel","read_num")

# filter by rn
umi_filtered <- umi %>%
	filter(rname %in% rn)

# Save file
write.table(umi_filtered, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)

