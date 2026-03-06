library(tidyr)
library(gtools)
library(dplyr)
library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

inputdir <- args[1]
inputcor <- args[2]
contig_des <- args[3]

contig_info_ribo <- read.table('contig_info_ribo.tsv', sep = "\t", header = FALSE)
colnames(contig_info_ribo) <- c("rname", "contig")
contig_info_ctl <- read.table('contig_info_ctl.tsv', sep = "\t", header = FALSE)
colnames(contig_info_ctl) <- c("rname", "contig")
contig_info <- rbind(contig_info_ribo, contig_info_ctl)

contig_design <- read.table(contig_des, sep = "\t", header = FALSE)
colnames(contig_design) <- c("contig", "l1", "l2", "l3", "tri1", "tri2", "tri3")
contig_tri <- left_join(contig_info, contig_design, by = "contig")

bases <- 63:79
base_list <- lapply(bases, function(pos) {
  df <- read.table(paste0("./", inputdir, "/", pos, "_for_baseofallreads"))
  colnames(df) <- c(paste0("base", pos), "rname")
  df
})

data_for_joined <- Reduce(function(x, y) merge(x, y, by = "rname"), base_list)

bases_rev <- 192:208
base_list_rev <- lapply(bases_rev, function(pos) {
  df <- read.table(paste0("./", inputdir, "/", pos, "_rev_baseofallreads"))
  colnames(df) <- c(paste0("base", pos), "rname")
  df
})

data_rev_joined <- Reduce(function(x, y) merge(x, y, by = "rname"), base_list_rev)

complement_base <- function(base) {
  complement_map <- c(A = "T", T = "A", G = "C", C = "G", "*" = "*")
  return(complement_map[base])
}
cols_to_change <- setdiff(names(data_rev_joined), "rname")
data_rev_joined[cols_to_change] <- lapply(
  data_rev_joined[cols_to_change],
  function(col) complement_base(col)
)

data_correct <- merge(data_for_joined, data_rev_joined, by = "rname") %>%
  left_join(contig_tri, by = "rname") %>%
  mutate(base27correct = ifelse(base63 == base208, 1, 0)) %>%
  mutate(base28correct = ifelse(base64 == base207, 1, 0)) %>%
  mutate(base29correct = ifelse(base65 == base206, 1, 0)) %>%
  mutate(base30correct = ifelse(base66 == base205, 1, 0)) %>%
  mutate(base31correct = ifelse(base67 == base204, 1, 0)) %>%
  mutate(base32correct = ifelse(base68 == base203, 1, 0)) %>%
  mutate(base33correct = ifelse(base69 == base202, 1, 0)) %>%
  mutate(base34correct = ifelse(base70 == tri1, 1, 0)) %>%
  mutate(base35correct = ifelse(base71 == tri2, 1, 0)) %>%
  mutate(base36correct = ifelse(base72 == tri3, 1, 0)) %>%
  mutate(base37correct = ifelse(base73 == base198, 1, 0)) %>%
  mutate(base38correct = ifelse(base74 == base197, 1, 0)) %>%
  mutate(base39correct = ifelse(base75 == base196, 1, 0)) %>%
  mutate(base40correct = ifelse(base76 == base195, 1, 0)) %>%
  mutate(base41correct = ifelse(base77 == base194, 1, 0)) %>%
  mutate(base42correct = ifelse(base78 == base193, 1, 0)) %>%
  mutate(base43correct = ifelse(base79 == base192, 1, 0)) %>%
  mutate_all(~ ifelse(is.na(.),0, .)) %>%
  select(contig,base27correct, base28correct, base29correct, base30correct, base31correct, base32correct, base33correct, base34correct, base35correct, base36correct, base37correct, base38correct, base39correct, base40correct, base41correct, base42correct, base43correct)

data_acc <- data_correct %>%
	group_by(contig) %>%
	summarise(across(everything(), ~ sum(as.numeric(.), na.rm = TRUE) / n(), .names = "acc_{col}"))

colnames(data_acc) <- c("contig","R-8","R-7","R-6","R-5","R-4","R-3","R-2","R-1","R","R+1","R+2","R+3","R+4","R+5","R+6","R+7","R+8")

theme_set(
  theme_classic(base_size = 18, base_family = "Times New Roman") +
    theme(legend.position = "top",
  axis.title = element_text(size = 20, color = "black"),
  axis.text = element_text(size = 18, color = "black"),
  legend.text = element_text(size = 18, color = "black"),
  legend.title = element_text(size = 18, color = "black"),
  axis.ticks.y = element_blank(),
  plot.title = element_text(size = 18, face = "bold", color = "black")
    )
)

data_long <- pivot_longer(data_acc, cols = starts_with("R"), names_to = "Position", values_to = "Accuracy")
data_long$Position <- factor(data_long$Position, levels = unique(data_long$Position))

cali_data <- read.table(inputcor, sep = "\t", header = FALSE)
colnames(cali_data) = c("contig", "Position", "ref", "depth", "basecall", "~", "rname")
cali_data$contig <- sub("_.*", "", cali_data$contig)

pos_labels <- setNames(
  c(paste0("R", -8:-1), "R", paste0("R+", 1:8)),
  63:79
)

stats_cali <- cali_data[cali_data$Position %in% 63:79, ] %>%
    filter(ref != "N") %>%
    group_by(contig,Position) %>%
    mutate(Position = pos_labels[as.character(Position)]) %>%
    summarise(correct = sum(grepl("\\.", basecall)), total = n(), Accuracy = correct / total)

pos_order <- c(paste0("R-", 8:1), "R", paste0("R+", 1:8))
stats_cali$Position <- factor(stats_cali$Position, levels = pos_order)

inputdir_ctl <- args[4]
inputcor_ctl <- args[5]
outpath_plot_raw <- args[6]
outpath_plot_corrected <- args[7]

base_list_ctl <- lapply(bases, function(pos) {
  df <- read.table(paste0("./", inputdir_ctl, "/", pos, "_for_baseofallreads"))
  colnames(df) <- c(paste0("base", pos), "rname")
  df
})

data_for_joined_ctl <- Reduce(function(x, y) merge(x, y, by = "rname"), base_list_ctl)

base_list_rev_ctl <- lapply(bases_rev, function(pos) {
  df <- read.table(paste0("./", inputdir_ctl, "/", pos, "_rev_baseofallreads"))
  colnames(df) <- c(paste0("base", pos), "rname")
  df
})

data_rev_joined_ctl <- Reduce(function(x, y) merge(x, y, by = "rname"), base_list_rev_ctl)

complement_base <- function(base) {
  complement_map <- c(A = "T", T = "A", G = "C", C = "G", "*" = "*")
  return(complement_map[base])
}
cols_to_change <- setdiff(names(data_rev_joined_ctl), "rname")
data_rev_joined_ctl[cols_to_change] <- lapply(
  data_rev_joined_ctl[cols_to_change],
  function(col) complement_base(col)
)

data_correct_ctl <- merge(data_for_joined_ctl, data_rev_joined_ctl, by = "rname") %>%
  left_join(contig_tri, by = "rname") %>%
  mutate(base27correct = ifelse(base63 == base208, 1, 0)) %>%
  mutate(base28correct = ifelse(base64 == base207, 1, 0)) %>%
  mutate(base29correct = ifelse(base65 == base206, 1, 0)) %>%
  mutate(base30correct = ifelse(base66 == base205, 1, 0)) %>%
  mutate(base31correct = ifelse(base67 == base204, 1, 0)) %>%
  mutate(base32correct = ifelse(base68 == base203, 1, 0)) %>%
  mutate(base33correct = ifelse(base69 == base202, 1, 0)) %>%
  mutate(base34correct = ifelse(base70 == tri1, 1, 0)) %>%
  mutate(base35correct = ifelse(base71 == tri2, 1, 0)) %>%
  mutate(base36correct = ifelse(base72 == tri3, 1, 0)) %>%
  mutate(base37correct = ifelse(base73 == base198, 1, 0)) %>%
  mutate(base38correct = ifelse(base74 == base197, 1, 0)) %>%
  mutate(base39correct = ifelse(base75 == base196, 1, 0)) %>%
  mutate(base40correct = ifelse(base76 == base195, 1, 0)) %>%
  mutate(base41correct = ifelse(base77 == base194, 1, 0)) %>%
  mutate(base42correct = ifelse(base78 == base193, 1, 0)) %>%
  mutate(base43correct = ifelse(base79 == base192, 1, 0)) %>%
  mutate_all(~ ifelse(is.na(.),0, .)) %>%
  select(contig,base27correct, base28correct, base29correct, base30correct, base31correct, base32correct, base33correct, base34correct, base35correct, base36correct, base37correct, base38correct, base39correct, base40correct, base41correct, base42correct, base43correct)

data_acc_ctl <- data_correct_ctl %>%
	group_by(contig) %>%
	summarise(across(everything(), ~ sum(as.numeric(.), na.rm = TRUE) / n(), .names = "acc_{col}"))

colnames(data_acc_ctl) <- c("contig","R-8","R-7","R-6","R-5","R-4","R-3","R-2","R-1","R","R+1","R+2","R+3","R+4","R+5","R+6","R+7","R+8")

theme_set(
  theme_classic(base_size = 18, base_family = "Times New Roman") +
    theme(legend.position = "top",
  axis.title = element_text(size = 20, color = "black"),
  axis.text = element_text(size = 18, color = "black"),
  legend.text = element_text(size = 18, color = "black"),
  legend.title = element_text(size = 18, color = "black"),
  axis.ticks.y = element_blank(),
  plot.title = element_text(size = 18, face = "bold", color = "black")
    )
)

data_long_ctl <- pivot_longer(data_acc_ctl, cols = starts_with("R"), names_to = "Position", values_to = "Accuracy")
data_long_ctl$Position <- factor(data_long_ctl$Position, levels = unique(data_long_ctl$Position))

cali_data_ctl <- read.table(inputcor_ctl, sep = "\t", header = FALSE)
colnames(cali_data_ctl) = c("contig", "Position", "ref", "depth", "basecall", "~", "rname")
cali_data_ctl$contig <- sub("_.*", "", cali_data_ctl$contig)

pos_labels <- setNames(
  c(paste0("R", -8:-1), "R", paste0("R+", 1:8)),
  63:79
)

stats_cali_ctl <- cali_data_ctl[cali_data_ctl$Position %in% 63:79, ] %>%
    filter(ref != "N") %>%
    group_by(contig,Position) %>%
    mutate(Position = pos_labels[as.character(Position)]) %>%
    summarise(correct = sum(grepl("\\.", basecall)), total = n(), Accuracy = correct / total)

pos_order <- c(paste0("R-", 8:1), "R", paste0("R+", 1:8))
stats_cali_ctl$Position <- factor(stats_cali_ctl$Position, levels = pos_order)

acc_plot_raw <- ggplot() +
    geom_boxplot(data = data_long_ctl, aes (x = Position, y = Accuracy), color = "#799FA2", alpha = 0.5) +
    stat_summary(data = data_long_ctl, aes (x = Position, y = Accuracy), fun.y = mean, geom = "point", shape = 18, size = 5, color = "#C0B6A4") +
    geom_boxplot(data = data_long, aes (x = Position, y = Accuracy), color = "#799FA2", alpha = 0.5) +
    stat_summary(data = data_long, aes (x = Position, y = Accuracy), fun.y = mean, geom = "point", shape = 18, size = 5, color = "#E75C63") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_cartesian(ylim = c(0, 1))

ggsave(outpath_plot_raw, width = 6, height = 5, plot = acc_plot_raw, device = cairo_pdf)

acc_plot_cor <- ggplot() +
    geom_boxplot(data = stats_cali_ctl, aes (x = Position, y = Accuracy), color = "#799FA2", alpha = 0.5) +
    stat_summary(data = stats_cali_ctl, aes (x = Position, y = Accuracy), fun.y = mean, geom = "point", shape = 18, size = 5, color = "#C0B6A4") +
    geom_boxplot(data = stats_cali, aes (x = Position, y = Accuracy), color = "#799FA2", alpha = 0.5) +
    stat_summary(data = stats_cali, aes (x = Position, y = Accuracy), fun.y = mean, geom = "point", shape = 18, size = 5, color = "#E75C63") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_cartesian(ylim = c(0, 1))

ggsave(outpath_plot_corrected, width = 6, height = 5, plot = acc_plot_cor, device = cairo_pdf)