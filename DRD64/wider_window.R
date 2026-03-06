library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(ggridges)
library(patchwork)
library(gtools)
library(knitr)
library(kableExtra)
library(ggrastr)
library(viridis)

theme_set(
  theme_classic(base_size = 18) +
    theme(legend.position = "top",
  axis.title = element_text(size = 20, color = "black"),
  axis.text = element_text(size = 18, color = "black"),
  legend.text = element_text(size = 18, color = "black"),
  legend.title = element_text(size = 18, color = "black"),
  axis.ticks.y = element_blank(),
  plot.title = element_text(size = 18, face = "bold", color = "black")
    )
)

data <- fread("./data/eventalign_output.tsv", sep = "\t", header = TRUE)
data_sum <- data %>%
    group_by(read_name,position) %>%
    summarise(ampl_mean = mean(event_level_mean, na.rm = TRUE), # I keep the amplitude data of kmer regarded as "NNNNNNNNN" to keep all distortions
    ampl_stdv = mean(event_stdv, na.rm = TRUE),
    dwell_time = sum(event_length),
    reference_kmer = first(reference_kmer[model_kmer != "NNNNNNNNN"], default = NA_character_),
    model_kmer = first(model_kmer[model_kmer != "NNNNNNNNN"], default = NA_character_),
    model_mean = mean(model_mean[model_kmer != "NNNNNNNNN"], na.rm = TRUE),
    model_stdv = mean(model_stdv[model_kmer != "NNNNNNNNN"], na.rm = TRUE),
    contig = first(contig) ,
    .groups = "drop"
)

ribo_data <- data_sum[data_sum$position %in% 77:81, c("read_name", "position", "ampl_mean", "ampl_stdv", "dwell_time", "contig")] %>%
    mutate(Type = "Ribo")
contigs <- mixedsort(unique(ribo_data$contig))

all_plots_time_largerWD <- list()
for (x in contigs) {
  data_sum_contig <- data_sum %>%
    filter(contig == x, position %in% 48:102)

  plot <- ggplot() +
    rasterise(geom_boxplot(data = data_sum_contig, aes(x = position, y = dwell_time, group = position), outlier.shape = NA), dpi = 300) +
    coord_cartesian(xlim = c(50, 100), ylim = c(0, 0.015)) +
    geom_vline(xintercept = c(64,67),  linetype = "dashed", color = "#DB7176") +
    geom_vline(xintercept = c(83,86),  linetype = "dashed", color = "#698a9c") +
    geom_segment(
    aes(
      x = 82, xend = 74,
      y = 0, yend = 0
    ),
    inherit.aes = FALSE,
    color = "#DB7176",
    size = 2
  ) +
    geom_segment(
    aes(
      x = 101, xend = 93,
      y = 0, yend = 0
    ),
    inherit.aes = FALSE,
    color = "#698a9c",
    size = 2
  ) + labs(x = "Kmer_Position", y = "Dwell_Time") +
    ggtitle(x)
  all_plots_time_largerWD[[x]] <- plot
}

num_plots <- length(all_plots_time_largerWD)
ncols <- 8
plot_width <- 6  
plot_height <- 6
nrows <- ceiling(num_plots / ncols)

total_width <- plot_width * ncols
total_height <- plot_height * nrows

combined_plot <- wrap_plots(all_plots_time_largerWD, ncol = ncols)
ggsave("./data/all_contigs_dynamics_largerwindow.pdf", plot = combined_plot, width = total_width, height = total_height, device = "pdf")