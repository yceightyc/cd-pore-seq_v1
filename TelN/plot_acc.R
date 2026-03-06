library(dplyr)
library(tidyr)
library(ggplot2)

# Import argument
args <- commandArgs(trailingOnly = TRUE)
teln_path <- args[1]
R10_path <- args[2]
HD_path <- args[3]
mods_path <- args[4]
output_path <- args[5]

teln <- read.table(teln_path, sep = "\t", header = FALSE)
colnames(teln) <- c("contig", "rname_for", "rev_umi","for_umi","for_umi_cali","pos1","pos2","pos3","matching")
mods <- read.table(mods_path, sep = "\t", header = FALSE)
colnames(mods) <- c("contig", "rname_for", "rev_umi","for_umi","for_umi_cali","pos1","pos2","pos3","matching")

r10 <- read.table(R10_path, sep = "\t", header = TRUE)
hd <- read.table(HD_path, sep = "\t", header = TRUE)

teln_stats <- teln %>%
    group_by(contig) %>%
    summarise(matched = sum (matching == 1), total = n(), Matching_Rate = matched / total) %>%
    mutate(Library = "TelN-R10")

mods_stats <- mods %>%
    group_by(contig) %>%
    summarise(matched = sum (matching == 1), total = n(), Matching_Rate = matched / total) %>%
    mutate(Library = "TelN-R10-MM")

hd_stats <- hd %>%
    group_by(contig_rev) %>%
    summarise(matched = sum(contig_for == contig_rev), total = n(), Matching_Rate = matched / total) %>%
    mutate(Library = "UMI-HD")
colnames(hd_stats) = c("contig","matched","total","Matching_Rate", "Library")

r10_stats <- r10 %>%
    group_by(contig_rev) %>%
    summarise(matched = sum(contig_for == contig_rev), total = n(), Matching_Rate = matched / total) %>%
    mutate(Library = "UMI-R10")
colnames(r10_stats) = c("contig","matched","total","Matching_Rate", "Library")

stats_all <- rbind(r10_stats, hd_stats, teln_stats, mods_stats)

theme_set(
  theme_classic(base_size = 18, base_family = "Arial") +
    theme(legend.position = "top",
  axis.title = element_text(size = 20, color = "black"),
  axis.text = element_text(size = 18, color = "black"),
  legend.text = element_text(size = 18, color = "black"),
  legend.title = element_text(size = 18, color = "black"),
  plot.title = element_text(size = 18, face = "bold", color = "black")
    )
)
y_breaks <- seq(0, 1, by = 0.05)

summary_stats <- stats_all %>%
  group_by(Library) %>%
  mutate(
    mean_MR = mean(Matching_Rate, na.rm = TRUE),
    stdv_MR = sd(Matching_Rate, na.rm = TRUE))

acc_plot <- ggplot(data = stats_all, aes(x = Library, y = Matching_Rate)) +
  geom_point(color = "#799FA2") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 5, color = "#E75C63") +
  geom_errorbar(
    data = summary_stats,
    aes(x = Library, ymin = mean_MR - stdv_MR, ymax = mean_MR + stdv_MR),
    width = 0.2,
    color = "#E75C63"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0.5, 1)) +
  xlab("") +
  ylab("Matching Accuracy") +
  scale_y_continuous(breaks = y_breaks)

ggsave(output_path, plot = acc_plot, width = 5, height = 5, dpi = 300, device = cairo_pdf)
