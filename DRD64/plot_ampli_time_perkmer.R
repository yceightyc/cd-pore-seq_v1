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
library(effectsize)

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
ribo_data$position <- recode(ribo_data$position, "77" = "4th", "78" = "5th", "79" = "6th", "80" = "7th", "81" = "8th")

ctl_data <- data_sum[data_sum$position %in% 96:100, c("read_name", "position", "ampl_mean", "ampl_stdv", "dwell_time", "contig")] %>%
    mutate(Type = "Ctl")
ctl_data$position <- recode(ctl_data$position, "96" = "4th", "97" = "5th", "98" = "6th", "99" = "7th", "100" = "8th")

model_data <- data_sum[data_sum$position %in% 96:100, c("read_name", "position", "model_mean", "model_stdv", "dwell_time", "contig")] %>%
    mutate(Type = "model")
model_data$position <- recode(model_data$position, "96" = "4th", "97" = "5th", "98" = "6th", "99" = "7th", "100" = "8th")
colnames(model_data) = c("read_name", "position", "ampl_mean", "ampl_stdv", "dwell_time", "contig", "Type")

merged_data <- rbind(ribo_data, ctl_data) 

stats_merged <- merged_data %>%
    group_by(position, contig) %>%
    summarise(
    ks_p_value = suppressWarnings(ks.test(
      ampl_mean[Type == "Ctl"],
      ampl_mean[Type == "Ribo"]
    )$p.value),
    U_p_value = suppressWarnings(wilcox.test(
      ampl_mean[Type == "Ctl"],
      ampl_mean[Type == "Ribo"]
    )$p.value),
    ks_p_value_time = suppressWarnings(ks.test(
      dwell_time[Type == "Ctl"],
      dwell_time[Type == "Ribo"]
    )$p.value),
    U_p_value_time = suppressWarnings(wilcox.test(
      dwell_time[Type == "Ctl"],
      dwell_time[Type == "Ribo"]
    )$p.value),    
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_ks = p.adjust(ks_p_value, method = "BH", n = 320),
    significance_ks = case_when(
      p_adj_ks <= 0.001 ~ "***",
      p_adj_ks <= 0.01  ~ "**",
      p_adj_ks <= 0.05  ~ "*",
      TRUE ~ as.character(round(p_adj_ks, 2))
    ),
    p_adj_U = p.adjust(U_p_value, method = "BH", n = 320),
     significance_U = case_when(
      p_adj_U <= 0.001 ~ "***",
      p_adj_U <= 0.01  ~ "**",
      p_adj_U <= 0.05  ~ "*",
      TRUE ~ as.character(round(p_adj_U, 2))
    ),
    p_adj_ks_time = p.adjust(ks_p_value_time, method = "BH", n = 320),
    significance_ks_time = case_when(
      p_adj_ks_time <= 0.001 ~ "***",
      p_adj_ks_time <= 0.01  ~ "**",
      p_adj_ks_time <= 0.05  ~ "*",
      TRUE ~ as.character(round(p_adj_ks_time, 2))
    ),
    p_adj_U_time = p.adjust(U_p_value_time, method = "BH", n = 320),
     significance_U_time = case_when(
      p_adj_U_time <= 0.001 ~ "***",
      p_adj_U_time <= 0.01  ~ "**",
      p_adj_U_time <= 0.05  ~ "*",
      TRUE ~ as.character(round(p_adj_U_time, 2))
    ),
    effect_U <- effectsize::rank_biserial(
     ampl_mean ~ Type,
     data = data.frame(
      ampl_mean = ampl_mean,
      Type = Type
     )
    )$r_rank_biserial ,

    effect_U_time <- effectsize::rank_biserial(
     dwell_time ~ Type,
     data = data.frame(
      dwell_time = dwell_time,
      Type = Type
     )
    )$r_rank_biserial
  )

contigs <- mixedsort(unique(ribo_data$contig))
all_plots <- list()

for (x in contigs) {
plot_data <- merged_data[merged_data$contig == x, c("ampl_mean", "position", "Type")]
colnames (plot_data) = c("Mean_Amplitude", "Kmer_Position", "Type")

model_plot <- model_data[model_data$contig == x, c("ampl_mean", "position", "Type")]  
colnames(model_plot) = c("Mean_Amplitude", "Kmer_Position", "Type")

get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

model_plot <- model_plot %>% group_by(Kmer_Position) %>% summarise(Mean_Amplitudes = get_mode(Mean_Amplitude), type = first(Type)) %>%
 arrange(Kmer_Position) %>%
 mutate(next_position = lead(Kmer_Position))
model_plot[is.na(model_plot$next_position), ]$next_position <- "Inf"

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

custom_colors <- c("Ctl" = "#698a9c", "Ribo" = "#E75C63") 

first_tick_x <- min(pretty(plot_data$Mean_Amplitude))

plot <- ggplot(plot_data, aes(x = Mean_Amplitude, y = Kmer_Position, fill = Type)) +
    geom_density_ridges(
    position = position_points_jitter(height = 0),
    alpha = 0.8,
    scale = 1.2,
    color = NA
  ) +
  geom_segment(
    data = model_plot,
    aes(
      x = Mean_Amplitudes, xend = Mean_Amplitudes,
      y = Kmer_Position, yend = next_position
    ),
    inherit.aes = FALSE,
    color = "#7B7975",
    linetype = "dashed",
    size = 1
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_y_discrete(labels = function(x) ifelse(x == "Inf", "", x)) + 
  coord_cartesian(ylim = c(1.2, length(unique(plot_data$Kmer_Position)) + 0.3)) +
  geom_text(
    data = stats_merged[stats_merged$contig == x, ],
    aes(x = first_tick_x , y = position, label = significance_ks),
    inherit.aes = FALSE,
    angle = 90,          
    hjust = 0.5,         
    vjust = 0.4,         
    size = 6,
    family = "Times New Roman"
  ) + ggtitle(x)
  all_plots[[x]] <- plot
}

num_plots <- length(all_plots)
ncols <- 8
plot_width <- 6 
plot_height <- 6
nrows <- ceiling(num_plots / ncols)

total_width <- plot_width * ncols
total_height <- plot_height * nrows

combined_plot <- wrap_plots(all_plots, ncol = ncols)
ggsave("./data/all_contigs_ampl_density.pdf", plot = combined_plot, width = total_width, height = total_height, device = cairo_pdf)


all_plots_box <- list()
for (x in contigs) {
  plot_data <- merged_data[merged_data$contig == x, c("ampl_mean", "position", "Type")]
  colnames (plot_data) = c("Mean_Amplitude", "Kmer_Position", "Type")
  
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
custom_colors <- c("Ctl" = "#698a9c", "Ribo" = "#DB7176") 

plot <- ggplot(data = plot_data, aes(x = Kmer_Position, y = Mean_Amplitude)) +
    geom_boxplot(color = "#878580", aes(fill = Type), alpha = 1) +
    scale_fill_manual(values = custom_colors) +
    geom_text(
    data = stats_merged[stats_merged$contig == x, ],
    aes(x = position , y = 20, label = significance_U),
    inherit.aes = FALSE,
    hjust = 0.5,         
    vjust = 0.5,         
    size = 8,
    family = "Times New Roman"
  ) + ggtitle(x)
  all_plots_box[[x]] <- plot
    }
                   
num_plots_box <- length(all_plots_box)
ncols <- 8
plot_width <- 6  
plot_height <- 6
nrows_box <- ceiling(num_plots_box / ncols)

total_width <- plot_width * ncols
total_height <- plot_height * nrows_box

combined_plot_box <- wrap_plots(all_plots_box, ncol = ncols)
ggsave("./data/all_contigs_ampl_boxplot.pdf", plot = combined_plot_box, width = total_width, height = total_height, device = cairo_pdf)

all_plots_time <- list()
for (x in contigs) {
  time_data <- merged_data[merged_data$contig == x, c("dwell_time", "position", "Type")]
  colnames (time_data) = c("Dwell_Time", "Kmer_Position", "Type")
  
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
custom_colors <- c("Ctl" = "#698a9c", "Ribo" = "#E75C63") 

plot <- ggplot(time_data, aes(x = Dwell_Time, y = Kmer_Position, fill = Type)) +
  geom_density_ridges(
      position = position_points_jitter(height = 0),
      alpha = 0.8,
      scale = 1.2,
      color = NA
    ) +
  coord_cartesian(xlim = c(0, 0.03)) +
  
  scale_fill_manual(values = custom_colors) +
  geom_text(
    data = stats_merged[stats_merged$contig == x, ],
    aes(x = 0.03 , y = position, label = significance_ks_time),
    inherit.aes = FALSE,
    angle = 90,          
    hjust = 0.5,         
    vjust = 0,         
    size = 6,
    family = "Arial") + ggtitle(x)
  all_plots_time[[x]] <- plot
    }
                   
combined_plot_time <- wrap_plots(all_plots_time, ncol = 4)

ggsave(
  "./data/all_contigs_time_density.pdf",
  plot = combined_plot_time,
  width = 8.27,
  height = 11.69,
  device = cairo_pdf
)


all_plots_box_time <- list()
for (x in contigs) {
  time_data <- merged_data[merged_data$contig == x, c("dwell_time", "position", "Type")]
  colnames (time_data) = c("Dwell_Time", "Kmer_Position", "Type")
  
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
custom_colors <- c("Ctl" = "#698a9c", "Ribo" = "#DB7176") 

plot <- ggplot(data = time_data, aes(x = Kmer_Position, y = Dwell_Time)) +
    geom_boxplot(color = "#878580", aes(fill = Type), alpha = 1) +
  scale_fill_manual(values = custom_colors) +
  geom_text(
    data = stats_merged[stats_merged$contig == x, ],
    aes(x = position , y = 0, label = significance_U),
    inherit.aes = FALSE,
    hjust = 0.5,         
    vjust = 1,         
    size = 8,
    family = "Times New Roman"
  ) +
  coord_cartesian(ylim = c(0, 0.03)) + ggtitle(x)
  all_plots_box_time[[x]] <- plot
    }
                 
num_plots_box_time <- length(all_plots_box_time)
ncols <- 8
plot_width <- 6  
plot_height <- 6
nrows_box_time <- ceiling(num_plots_box_time / ncols)

total_width <- plot_width * ncols
total_height_time <- plot_height * nrows_box_time

combined_plot_box_time <- wrap_plots(all_plots_box_time, ncol = ncols)
ggsave("./data/all_contigs_time_boxplot.pdf", plot = combined_plot_box_time, width = total_width, height = total_height_time, device = cairo_pdf)

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

ribo_offset_time <- data_sum[data_sum$position %in% c(64, 67), c("read_name", "position", "ampl_mean", "ampl_stdv", "dwell_time", "contig")] %>%
    mutate(Type = "Ribo")
ribo_offset_time$position <- recode(ribo_offset_time$position, "64" = "R-18", "67" = "R-15")

ctl_offset_time <- data_sum[data_sum$position %in% c(83, 86), c("read_name", "position", "ampl_mean", "ampl_stdv", "dwell_time", "contig")] %>%
    mutate(Type = "Ctl")
ctl_offset_time$position <- recode(ctl_offset_time$position, "83" = "R-18", "86" = "R-15")

merged_offset_time <- rbind(ribo_offset_time, ctl_offset_time)

stats_merged_offset_time <- merged_offset_time %>%
    group_by(position, contig) %>%
    summarise(
    ks_p_value_time = suppressWarnings(ks.test(
      dwell_time[Type == "Ctl"],
      dwell_time[Type == "Ribo"]
    )$p.value),
    U_p_value_time = suppressWarnings(wilcox.test(
      dwell_time[Type == "Ctl"],
      dwell_time[Type == "Ribo"]
    )$p.value),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_ks_time = p.adjust(ks_p_value_time, method = "BH", n = 192),
    significance_ks_time = case_when(
      p_adj_ks_time <= 0.001 ~ "***",
      p_adj_ks_time <= 0.01  ~ "**",
      p_adj_ks_time <= 0.05  ~ "*",
      TRUE ~ as.character(round(p_adj_ks_time, 2))
    ),
    p_adj_U_time = p.adjust(U_p_value_time, method = "BH", n = 192),
     significance_U_time = case_when(
      p_adj_U_time <= 0.001 ~ "***",
      p_adj_U_time <= 0.01  ~ "**",
      p_adj_U_time <= 0.05  ~ "*",
      TRUE ~ as.character(round(p_adj_U_time, 2))
    )
  )

all_plots_offset_time <- list()
for (x in contigs) {
  offset_time_data <- merged_offset_time[merged_offset_time$contig == x, c("dwell_time", "position", "Type")]
  colnames (offset_time_data) = c("Dwell_Time", "Kmer_Position", "Type")

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
custom_colors <- c("Ctl" = "#698a9c", "Ribo" = "#E75C63")

plot <- ggplot(offset_time_data, aes(x = Dwell_Time, y = Kmer_Position, fill = Type)) +
  geom_density_ridges(
    position = position_points_jitter(height = 0),
    alpha = 0.8,
    scale = 1.2,
    color = NA
  ) +
  coord_cartesian(xlim = c(0, 0.03)) +
  
  scale_fill_manual(values = custom_colors) +
  geom_text(
    data = stats_merged_offset_time[stats_merged_offset_time$contig == x, ],
    aes(x = 0.03 , y = position, label = significance_ks_time),
    inherit.aes = FALSE,
    angle = 90,          
    hjust = 0.5,         
    vjust = 0,         
    size = 6,
    family = "Times New Roman") + ggtitle(x)
  all_plots_offset_time[[x]] <- plot
    }

num_plots_time <- length(all_plots_offset_time)
ncols <- 8
plot_width <- 6  # inches per plot
plot_height <- 6
nrows_time <- ceiling(num_plots_time / ncols)

total_width <- plot_width * ncols
total_height_time <- plot_height * nrows_time

combined_plot_time <- wrap_plots(all_plots_offset_time, ncol = ncols)
ggsave("./data/all_contigs_offset_time_density.pdf", plot = combined_plot_time, width = total_width, height = total_height_time, device = cairo_pdf)


all_plots_box_time_offset <- list()
for (x in contigs) {
  offset_time_data <- merged_offset_time[merged_offset_time$contig == x, c("dwell_time", "position", "Type")]
  colnames (offset_time_data) = c("Dwell_Time", "Kmer_Position", "Type")

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
custom_colors <- c("Ctl" = "#698a9c", "Ribo" = "#DB7176")

plot <- ggplot(data = offset_time_data, aes(x = Kmer_Position, y = Dwell_Time)) +
    geom_boxplot(color = "#878580", aes(fill = Type), alpha = 1) +
  scale_fill_manual(values = custom_colors) +
  geom_text(
    data = stats_merged_offset_time[stats_merged_offset_time$contig == x, ],
    aes(x = position , y = 0, label = significance_U_time),
    inherit.aes = FALSE,
    hjust = 0.5,         
    vjust = 1,         
    size = 8,
    family = "Times New Roman"
  ) +
  coord_cartesian(ylim = c(0, 0.03)) + ggtitle(x)
  all_plots_box_time_offset[[x]] <- plot
    }

num_plots_box_time <- length(all_plots_box_time_offset)
ncols <- 8
plot_width <- 6  
plot_height <- 6
nrows_box_time <- ceiling(num_plots_box_time / ncols)

total_width <- plot_width * ncols
total_height_time <- plot_height * nrows_box_time

combined_plot_box_time <- wrap_plots(all_plots_box_time_offset, ncol = ncols)
ggsave("./data/all_contigs_offset_time_boxplot.pdf", plot = combined_plot_box_time, width = total_width, height = total_height_time, device = cairo_pdf)


rn_ribo <- merged_data[merged_data$Type == "Ribo", "read_name"] %>% mutate(read_name = paste(read_name,"_Ribo", sep = ""))
rn_ctl <- merged_data[merged_data$Type == "Ctl", "read_name"] %>% mutate(read_name = paste(read_name,"_Ctl", sep = ""))
rn <- rbind(rn_ribo, rn_ctl)

model_input <- rbind(merged_data, merged_offset_time) %>%
    left_join(model_data[, c("read_name","position","ampl_mean")], by = c("read_name", "position"))
colnames(model_input) = c("read_name", "position", "ampl_mean", "ampl_stdv", "dwell_time", "contig", "Type", "ampl_model")

model_input <- model_input %>%
    mutate(read_name = paste(read_name,"_",Type, sep = ""), model_deviation = ampl_mean - ampl_model) %>%
    inner_join(rn, by = "read_name") %>% unique()

model_input_final <- model_input[, c("read_name", "position","ampl_mean","dwell_time","contig","Type","model_deviation")] %>%
    pivot_wider(names_from = c(position), values_from = c(ampl_mean, dwell_time, model_deviation))
reorder <- mixedsort(colnames(model_input_final))
model_input_final <- model_input_final %>% select(all_of(reorder))
write.table(model_input_final, "./data/model_input_final.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

contig_trinu <- fread("./input/contig_info.tsv", sep = "\t", header = TRUE)
colnames(contig_trinu) = c("contig","r-1", "r", "r+1", "trinu")

data_withRinfo <- merge(merged_data, contig_trinu, by = "contig")

data_ribo <- data_withRinfo[data_withRinfo$Type == "Ribo", ] %>%
    group_by(contig, r, trinu) %>%
    summarise(ampl_r = mean(ampl_mean), dwell_time_r = mean(dwell_time))
data_ctl <- data_withRinfo[data_withRinfo$Type == "Ctl", ] %>%
    group_by(contig, r, trinu) %>%
    summarise(ampl_c = mean(ampl_mean), dwell_time_c = mean(dwell_time))
rvsc <- merge(data_ribo, data_ctl, by = c("contig","r","trinu")) %>%
    mutate(change_a = abs(ampl_r - ampl_c), change_t =(abs(dwell_time_r - dwell_time_c) / dwell_time_c))
colnames(rvsc) = c("contig","rNMP","trinu","ampl_r","dwell_time_r","ampl_c","dwell_time_c","change_a","change_t")

theme_set(
  theme_classic(base_size = 18, base_family = "serif") +
    theme(legend.position = "top",
  axis.title = element_text(size = 20, color = "black"),
  axis.text = element_text(size = 18, color = "black"),
  legend.text = element_text(size = 18, color = "black"),
  legend.title = element_text(size = 18, color = "black"),
  axis.ticks.y = element_blank(),
  plot.title = element_text(size = 18, face = "bold", color = "black")
    )
)
custom_colors <- c("rA" = "#698A9C", "rC" = "#A7C7E7", "rG" = "#DB7176", "rU" = "#C0B6A4") 

change_plot <- ggplot(rvsc, aes(x = change_a, y = change_t, fill = rNMP)) +
  geom_point(shape = 21, size = 4, color = "white") + 
  geom_text(aes(label = trinu, color = rNMP), vjust = -1, , show.legend = FALSE) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "Change in Amplitude (pA)",
    y = "% Change in Dwell Time"
  )
ggsave("./data/Change_Time_Inten.pdf", plot = change_plot, width = 6, height = 6, device = "pdf")

stats_withRinfo <- merge(stats_merged, contig_trinu, by = "contig")
sig_data <- stats_withRinfo %>%
    select(position, trinu,  p_adj_U,significance_U, p_adj_U_time, significance_U_time)
colnames(sig_data) = c("Kmer","DRD","p_adj_U","significance_U","p_adj_U_time",  "significance_U_time")

stats_fmt <- sig_data %>%
  mutate(
    Amplitude =  p_adj_U,
    Dwell_time =  p_adj_U_time
  ) %>%
  select(DRD, Kmer, Amplitude, Dwell_time)

stats_wide <- stats_fmt %>%
  pivot_wider(
    id_cols = c(DRD),
    names_from = Kmer,
    values_from = c(Amplitude, Dwell_time),
    names_glue = "{Kmer}_{.value}"
  )

stats_counts <- sig_data %>%
  group_by(DRD) %>%
  summarise(Num_significant_Kmer = sum(p_adj_U <= 0.05 | p_adj_U_time <= 0.05),
            .groups = "drop")
stats_wide <- left_join(stats_wide, stats_counts, by = c("DRD"))

id_cols <- c("DRD")
other_cols <- setdiff(names(stats_wide), id_cols)

stats_wide <- stats_wide[, c(id_cols, mixedsort(other_cols))]

pval_cols <- setdiff(names(stats_wide), c("DRD", "Num_significant_Kmer"))
sig_counts <- sapply(stats_wide[pval_cols], function(x) sum(x<=0.05))
summary_row <- c("Num_significant_DRD (out of 64)", as.character(sig_counts), "/")      
final_sig_tab <- rbind(stats_wide, summary_row)

tbl <- final_sig_tab %>%
  kbl(align = "c", booktabs = TRUE,
      col.names = c("DRD",
                  rep(c("Amplitude", "Dwell time"), 5),
                  "Num. Significant Kmer (out of 5)"), escape = FALSE) %>%
  add_header_above(c(" " = 1,
                     "Kmer Position 4" = 2,
                     "Kmer Position 5" = 2,
                     "Kmer Position 6" = 2,
                     "Kmer Position 7" = 2,
                     "Kmer Position 8" = 2,
                     " " = 1)) %>%
  kable_styling(full_width = FALSE, position = "center", font_size = 12)

tbl_html <- as.character(tbl)
html_doc <- paste0(
  "<!doctype html><html><head><meta charset='UTF-8'>",
  "<style>
    body, table { font-family: 'Times New Roman', Times, serif; }
    table { border-collapse: collapse; margin: 24px auto; font-size: 16px}
    th, td { padding: 6px 10px; }
  </style></head><body>",
  tbl_html,
  "</body></html>"
)

effect_data <- stats_withRinfo %>%
  select(
    position,
    trinu,
    effect_U,
    effect_U_time
  )

colnames(effect_data) <- c(
  "Kmer",
  "DRD",
  "effect_U",
  "effect_U_time"
)


effect_fmt <- effect_data %>%
  mutate(
    Amplitude = effect_U,
    Dwell_time = effect_U_time
  ) %>%
  select(
    DRD,
    Kmer,
    Amplitude,
    Dwell_time
  )


effect_wide <- effect_fmt %>%
  pivot_wider(
    id_cols = c(DRD),
    names_from = Kmer,
    values_from = c(
      Amplitude,
      Dwell_time
    ),
    names_glue = "{Kmer}_{.value}"
  )


# Reorder effect-size columns
id_cols_effect <- c("DRD")

other_cols_effect <- setdiff(
  names(effect_wide),
  id_cols_effect
)

effect_wide <- effect_wide[
  ,
  c(
    id_cols_effect,
    mixedsort(other_cols_effect)
  )
]

effect_tbl <- effect_wide %>%
  kbl(
    align = "c",
    booktabs = TRUE,
    col.names = c(
      "DRD",
      rep(
        c(
          "Amplitude",
          "Dwell time"
        ),
        5
      )
    ),
    escape = FALSE
  ) %>%
  add_header_above(
    c(
      " " = 1,
      "Kmer Position 4" = 2,
      "Kmer Position 5" = 2,
      "Kmer Position 6" = 2,
      "Kmer Position 7" = 2,
      "Kmer Position 8" = 2
    )
  ) %>%
  kable_styling(
    full_width = FALSE,
    position = "center",
    font_size = 12
  )


effect_tbl_html <- as.character(
  effect_tbl
)

effect_html_doc <- paste0(
  "<!doctype html><html><head><meta charset='UTF-8'>",
  "<style>
    body, table {
      font-family: Arial, sans-serif;
    }
    table {
      border-collapse: collapse;
      margin: 24px auto;
      font-size: 16px;
    }
    th, td {
      padding: 6px 10px;
    }
  </style></head><body>",
  effect_tbl_html,
  "</body></html>"
)

outfile_effect <- "./data/effect_size_table.html"
writeLines(effect_html_doc, outfile_effect, useBytes = TRUE)

heat_data <- sig_data %>%
  select(DRD, Kmer, p_adj_U, p_adj_U_time) %>%
  pivot_longer(cols = c(p_adj_U, p_adj_U_time),
               names_to = "Type",
               values_to = "p_value") %>%
  mutate(
    Type = recode(Type,
                  p_adj_U = "Amplitude",
                  p_adj_U_time = "Dwell time"),
    significance = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**",
      p_value <= 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Plot heatmap
heat_data$p_log = -log10(heat_data$p_value)
heat_data$p_log[is.infinite(heat_data$p_log)] <- 306
p_value_plot <- ggplot(heat_data, aes(x = DRD, y = Kmer, fill = p_log)) +
  geom_tile(color = "white") +
  facet_wrap(~Type, ncol = 2) +
  scale_fill_gradientn(
  name = "-log10(p)",
  colours = c("white", "#698a9c", "#db7176"),
  values = scales::rescale(c(1.3, 100, 306)),
  limits = c(1.3, 300),
  oob = scales::squish
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Heatmap of Adjusted P-values by DRD and Kmer Position",
    x = "DRD",
    y = "Kmer Position"
  )

ggsave("./data/p_Value.pdf", plot = p_value_plot, device = "pdf", width = 8, height = 3)

outfile <- "./data/signal_stats_table.html"
writeLines(html_doc, outfile, useBytes = TRUE)

theme_set(
  theme_minimal(base_size = 15, base_family = "serif") +
    theme(legend.position = "right", panel.grid = element_blank(),
  axis.title = element_text(size = 18, color = "black"),
  axis.text = element_text(size = 15, color = "black"),
  legend.text = element_text(size = 15, color = "black"),
  legend.title = element_text(size = 15, color = "black"),
  plot.margin = margin(5, 5, 5, 5),
  axis.ticks.y = element_blank(),
  plot.title = element_text(size = 15, face = "bold", color = "black")
    )
)

stats_fmt <- sig_data %>%
  mutate(
    Amplitude = ifelse(p_adj_U <= 0.05,
                       paste0(round(p_adj_U, 3), " (", significance_U, ")"),
                       round(p_adj_U, 3)),
    Dwell_time = ifelse(p_adj_U_time <= 0.05,
                        paste0(round(p_adj_U_time, 3), " (", significance_U_time, ")"),
                        round(p_adj_U_time, 3))
  ) %>%
  mutate(
    Amplitude = ifelse(Amplitude < 0.001, "<0.001 (***)", Amplitude),
    Dwell_time = ifelse(Dwell_time < 0.001, "<0.001 (***)", Dwell_time)
  ) %>%
  select(DRD, Kmer, Amplitude, Dwell_time)

stats_wide <- stats_fmt %>%
  pivot_wider(
    id_cols = c(DRD),
    names_from = Kmer,
    values_from = c(Amplitude, Dwell_time),
    names_glue = "{Kmer}_{.value}"
  )

stats_counts <- sig_data %>%
  group_by(DRD) %>%
  summarise(Num_significant_Kmer = sum(p_adj_U <= 0.05 | p_adj_U_time <= 0.05),
            .groups = "drop")
stats_wide <- left_join(stats_wide, stats_counts, by = c("DRD"))

id_cols <- c("DRD")
other_cols <- setdiff(names(stats_wide), id_cols)

stats_wide <- stats_wide[, c(id_cols, mixedsort(other_cols))]

pval_cols <- setdiff(names(stats_wide), c("DRD", "Num_significant_Kmer"))
sig_counts <- sapply(stats_wide[pval_cols], function(x) sum(grepl("\\*", x)))
summary_row <- c("Num_significant_DRD (out of 64)", as.character(sig_counts), "/")      
final_sig_tab <- rbind(stats_wide, summary_row)

final_sig_long <- final_sig_tab %>%
  select(-Num_significant_Kmer) %>%
  pivot_longer(
    cols = -DRD,
    names_to = "Index_Type",
    values_to = "pval_text"
  ) %>%
  mutate(Significant = grepl("\\*", pval_text)) %>%
  separate(Index_Type, into = c("Index", "Type"), sep = "_")

final_sig_summary <- final_sig_long %>%
  group_by(Index, Type) %>%
  summarise(Signif = sum(Significant), .groups = "drop")

final_sig_combined <- final_sig_long %>%
  group_by(DRD, Index) %>%
  summarise(Any_Signif = any(Significant), .groups = "drop") %>%
  group_by(Index) %>%
  summarise(Signif = sum(Any_Signif), .groups = "drop") %>%
  mutate(Type = "Combined")

final_sig_all <- final_sig_long %>%
  mutate(Index = "All") %>%
  group_by(DRD, Index, Type) %>%
  summarise(Any_Signif = any(Significant), .groups = "drop") %>%
  group_by(Index, Type) %>%
  summarise(Signif = sum(Any_Signif), .groups = "drop")

final_sig_all_combined <- final_sig_long %>%
  mutate(Index = "All") %>%
  group_by(DRD, Index) %>%
  summarise(Any_Signif = any(Significant), .groups = "drop") %>%
  group_by(Index) %>%
  summarise(Signif = sum(Any_Signif), .groups = "drop") %>%
  mutate(Type = "Combined")

final_sig_plot <- bind_rows(final_sig_summary, final_sig_combined, final_sig_all, final_sig_all_combined)
final_sig_plot$Type = gsub("Dwell_time", "Dwell_Time", final_sig_plot$Type)

final_sig_plot$Type <- factor(final_sig_plot$Type, levels = c("Combined", "Dwell_Time","Amplitude"))
final_sig_plot$Index <- factor(final_sig_plot$Index, levels = c("4th", "5th", "6th", "7th", "8th", "All"))

p <- ggplot(final_sig_plot, aes(x = Index, y = Type, fill = Signif)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Signif), color = "black", size = 5) +
  scale_fill_gradientn(
  colours = c("white", "#DB7176"),
  limits = c(40, 64),
  name = "Number\nof DRD"
  ) +
  labs(
    x = "Kmer Position",
    y = "",
    fill = "Significant count",
    title = "Number of Significantly altered DRD"
  )
ggsave("./data/Num_significance_summary.pdf", plot = p, width = 6, height = 4, device = "pdf")

offset_time_withRinfo <- merge(merged_offset_time, contig_trinu, by = "contig")

offset_data_ribo <- offset_time_withRinfo[offset_time_withRinfo$Type == "Ribo", ] %>%
    group_by(contig, r, trinu, position) %>%
    summarise(off_dwell_time_r = mean(dwell_time))
offset_data_ctl <- offset_time_withRinfo[offset_time_withRinfo$Type == "Ctl", ] %>%
    group_by(contig, r, trinu, position) %>%
    summarise(off_dwell_time_c = mean(dwell_time))
offset_rvsc <- merge(offset_data_ribo, offset_data_ctl, by = c("contig","r","trinu", "position")) %>%
    mutate(change_ot =(abs(off_dwell_time_r - off_dwell_time_c) / off_dwell_time_c))
colnames(offset_rvsc) = c("contig","rNMP","trinu","Offset_Position","off_dwell_time_r","off_dwell_time_c","change_ot")

offset_rvsc_plot <- offset_rvsc %>%
  select(rNMP, Offset_Position, change_ot, trinu) %>%
  pivot_wider(
    names_from = Offset_Position,
    values_from = change_ot,
    names_prefix = "Offset_Position_"
  )

theme_set(
  theme_classic(base_size = 18, base_family = "serif") +
    theme(legend.position = "top",
  axis.title = element_text(size = 20, color = "black"),
  axis.text = element_text(size = 18, color = "black"),
  legend.text = element_text(size = 18, color = "black"),
  legend.title = element_text(size = 18, color = "black"),
  axis.ticks.y = element_blank(),
  plot.title = element_text(size = 18, face = "bold", color = "black")
    )
)
custom_colors <- c("rA" = "#698A9C", "rC" = "#A7C7E7", "rG" = "#DB7176", "rU" = "#C0B6A4") 


change_plot <- ggplot(offset_rvsc_plot, aes(x = `Offset_Position_R-15`, y = `Offset_Position_R-18`, fill = rNMP)) +
  geom_point(shape = 21, size = 4, color = "white") + 
  geom_text(aes(label = trinu, color = rNMP), vjust = -1, , show.legend = FALSE) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "% Change in Dwell Time at R-15",
    y = "% Change in Dwell Time at R-18"
  )
ggsave("./data/Change_Offset_Time.pdf", plot = change_plot, width = 6, height = 6, device = "pdf")
