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

data <- fread("./data/eventalign_output.tsv", sep = "\t", header = TRUE)
data_sum <- data %>%
    group_by(read_name,position) %>%
    summarise(ampl_mean = mean(event_level_mean, na.rm = TRUE), # I keep the amplitude data of kmer regarded as "NNNNNNNNN" to keep all distortions
    ampl_stdv = mean(event_stdv, na.rm = TRUE),
    dwell_time = sum(event_length),
    reference_kmer = {
      valid <- !is.na(model_kmer) & model_kmer != "NNNNNNNNN"
      if (any(valid)) first(reference_kmer[valid]) else NA_character_
    },

    model_kmer = {
      valid <- !is.na(model_kmer) & model_kmer != "NNNNNNNNN"
      if (any(valid)) first(model_kmer[valid]) else NA_character_
    },
    model_mean = mean(model_mean[model_kmer != "NNNNNNNNN"], na.rm = TRUE),
    model_stdv = mean(model_stdv[model_kmer != "NNNNNNNNN"], na.rm = TRUE),
    contig = first(contig) ,
    .groups = "drop"
)

# Find the smallest number of contig
min_n <- data_sum %>%
    count(contig) %>%
    summarise(min_n = min(n)) %>%
    pull(min_n)

# Randomly downsample every contig to the same number
set.seed(66)  # makes the random sampling reproducible

data_sum <- data_sum %>%
    group_by(contig) %>%
    slice_sample(n = min_n) %>%
    ungroup()

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
    n_ctl = sum(Type == "Ctl" & !is.na(ampl_mean)),
    n_ribo = sum(Type == "Ribo" & !is.na(ampl_mean)),
    n_ctl_time  = sum(Type == "Ctl" & !is.na(dwell_time)),
    n_ribo_time = sum(Type == "Ribo" & !is.na(dwell_time)),

    ks_p_value = suppressWarnings(ks.test(
      ampl_mean[Type == "Ctl"],
      ampl_mean[Type == "Ribo"]
    )$p.value),
    U_p_value = suppressWarnings(wilcox.test(
      ampl_mean[Type == "Ctl"],
      ampl_mean[Type == "Ribo"],
      exact = FALSE,
      correct = FALSE
    )$p.value),
    ks_p_value_time = suppressWarnings(ks.test(
      dwell_time[Type == "Ctl"],
      dwell_time[Type == "Ribo"]
    )$p.value),
    U_p_value_time = suppressWarnings(wilcox.test(
      dwell_time[Type == "Ctl"],
      dwell_time[Type == "Ribo"],
      exact = FALSE,
      correct = FALSE
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
    z_U = qnorm(U_p_value / 2, lower.tail = FALSE),
    effect_U = abs(z_U) / sqrt(n_ctl + n_ribo),

    z_U_time = qnorm(U_p_value_time / 2, lower.tail = FALSE),
    effect_U_time = abs(z_U_time) / sqrt(n_ctl_time + n_ribo_time)
  )

contig_trinu <- fread("./input/contig_info.tsv", sep = "\t", header = TRUE)
colnames(contig_trinu) = c("contig","r-1", "r", "r+1", "trinu")

stats_withRinfo <- merge(stats_merged, contig_trinu, by = "contig")
sig_data <- stats_withRinfo %>%
    select(position, trinu, p_adj_U, significance_U, p_adj_U_time, significance_U_time, effect_U, effect_U_time)
colnames(sig_data) = c("Kmer","DRD","p_adj_U","significance_U","p_adj_U_time",  "significance_U_time", "effect_U", "effect_U_time")

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
heat_data$p_log[is.infinite(heat_data$p_log)] <-  max(heat_data$p_log)

p_value_plot <- ggplot(heat_data, aes(x = DRD, y = Kmer, fill = p_log)) +
  geom_tile(color = "white") +
  facet_wrap(~Type, ncol = 2) +
  scale_fill_gradientn(
  name = "-log10(p)",
  colours = c("white", "#698a9c", "#db7176"),
  values = scales::rescale(c(0, max(heat_data$p_log)/2, max(heat_data$p_log))),
  limits = c(0, max(heat_data$p_log)),
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

effect_data <- sig_data %>%
  select(DRD, Kmer, effect_U, effect_U_time) %>%
  pivot_longer(
    cols = c(effect_U, effect_U_time),
    names_to = "Type",
    values_to = "effect_size"
  ) %>%
  mutate(
    Type = recode(
      Type,
      effect_U = "Amplitude",
      effect_U_time = "Dwell time"
    )
  )

effect_plot <- ggplot(
  effect_data,
  aes(x = DRD, y = Kmer, fill = effect_size)
) +
  geom_tile(color = "white") +
  facet_wrap(~Type, ncol = 2) +
  
  scale_fill_gradientn(
    name = "Effect size |r|",
    colours = c("white", "#698a9c", "#db7176"),
    limits = c(0, max(effect_data$effect_size, na.rm = TRUE)),
    oob = scales::squish
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.background = element_rect(
      fill = "grey90",
      color = NA
    ),
    strip.text = element_text(face = "bold")
  ) +
  
  labs(
    title = "Wilcoxon Effect Size",
        x = "DRD",
    y = "Kmer Position"
  )

ggsave(
  "./data/effect_Size.pdf",
  plot = effect_plot,
  device = "pdf",
  width = 8,
  height = 3
)
