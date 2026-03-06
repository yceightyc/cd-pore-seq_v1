library(tidyr)
library(dplyr)
library(lsa)
library(ggplot2)
library(data.table)

eventalign_data_test <- as.data.frame(fread("input/eventalign_ctl/real_sig.tsv", header = TRUE, sep = "\t"))
eventalign_data_sim <- as.data.frame(fread("input/eventalign_ctl/sim_sig.tsv", header = TRUE, sep = "\t"))

eventalign_data_test_final <- eventalign_data_test[eventalign_data_test$position %in% 20:233, ] %>%
    group_by(position,read_name) %>%
    summarise(event_level_mean_mean = mean (event_level_mean), dwell_time = sum (event_length), contig = first(contig), position = first(position), read_name = first(read_name))

eventalign_data_sim_final <- eventalign_data_sim[eventalign_data_sim$position %in% 20:233, ] %>%
    group_by(position,read_name) %>%
    summarise(event_level_mean_mean = mean (event_level_mean), dwell_time = sum (event_length), contig = first(contig), position = first(position), read_name = first(read_name))

both_data <- merge (eventalign_data_test_final, eventalign_data_sim_final, by = c("position","read_name"))

sampled_data <- both_data %>%
  filter(read_name %in% sample(unique(read_name), 50))
sampled_data <- na.omit(sampled_data)

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

pearson_corr <- cor(sampled_data$event_level_mean_mean.x, sampled_data$event_level_mean_mean.y, method = "pearson")

cosine_sim <- sum(sampled_data$event_level_mean_mean.x * sampled_data$event_level_mean_mean.y) / (sqrt(sum(sampled_data$event_level_mean_mean.x^2)) * sqrt(sum(sampled_data$event_level_mean_mean.y^2)))

p <- ggplot(sampled_data, aes(x = event_level_mean_mean.x, y = event_level_mean_mean.y)) +
  geom_point(color = "#698A9C", alpha = 0.1) +
  geom_smooth(method = "lm", color = "#DB7176") +
  labs(x = "Mean Amplitude (pA)", y = "Simulated Mean Amplitude (pA)") +
  coord_cartesian(xlim = c(40,180), ylim = c(40,180)) +
  annotate("text",
           x = 50, y = 170,
           label = paste0("Pearson r = ", round(pearson_corr, 3),
                          "\nCosine Sim. = ", round(cosine_sim, 3)),
           hjust = 0, vjust = 1, size = 8, color = "black", family = "Times New Roman")

ggsave("input/eventalign_ctl/Amplitude_simvsreal_50reads.pdf", width = 6, height = 6, dpi = 300, plot = p, device = cairo_pdf)

sampled_data <- na.omit(sampled_data)

spearman_rho <- cor(sampled_data$dwell_time.x, sampled_data$dwell_time.y, method = "spearman")

cosine_sim <- sum(sampled_data$dwell_time.x * sampled_data$dwell_time.y) /
  (sqrt(sum(sampled_data$dwell_time.x^2)) * sqrt(sum(sampled_data$dwell_time.y^2)))

p2 <- ggplot(sampled_data, aes(x = dwell_time.x, y = dwell_time.y)) +
  geom_point(color = "#698A9C", alpha = 0.1) +
  geom_smooth(method = "lm", color = "#DB7176") +
  labs(x = "Dwell Time (s)", y = "Simulated Dwell Time (s)") +
  annotate("text",
           x = 0.005, y = 0.025,
           label = paste0("Spearman ρ = ", round(spearman_rho, 3),
                          "\nCosine Sim. = ", round(cosine_sim, 3)),
           hjust = 0, vjust = 1, size = 8, color = "black", family = "Times New Roman")

ggsave("input/eventalign_ctl/Dwelltime_simvsreal_50reads.pdf", width = 6, height = 6, dpi = 300, plot = p2, device = cairo_pdf)

y_breaks <- seq(0, 1, by = 0.1)

both_data_sum <- both_data %>%
    group_by("read_name") %>%
    mutate(Amplitude_Cos_Sim = sum(event_level_mean_mean.x * event_level_mean_mean.y) / (sqrt(sum(event_level_mean_mean.x^2)) * sqrt(sum(event_level_mean_mean.y^2))),
          Amplitude_r = cor(event_level_mean_mean.x, event_level_mean_mean.y, method = "pearson"),
          Dwell_Time_Cos_Sim = sum(dwell_time.x * dwell_time.y) / (sqrt(sum(dwell_time.x^2)) * sqrt(sum(dwell_time.y^2))),
          Dwell_time_ρ = cor(dwell_time.x, dwell_time.y, method = "spearman"))

df_long <- both_data_sum %>%
  pivot_longer(
    cols = c(Amplitude_Cos_Sim, Amplitude_r, Dwell_Time_Cos_Sim, Dwell_time_ρ),
    names_to = "metric",
    values_to = "value"
  )

df_long$metric <- factor(df_long$metric,
                         levels = c("Amplitude_Cos_Sim", "Amplitude_r", "Dwell_Time_Cos_Sim", "Dwell_time_ρ"))

p3 <- ggplot(df_long, aes(x = metric, y = value)) +
  stat_summary(fun = mean, geom = "bar", fill = "#698A9C", alpha = 0.8, width = 0.6) +
  stat_summary(fun.data = function(x) {
      mean_x <- mean(x)
      sd_x <- sd(x)
      data.frame(y = mean_x, ymin = mean_x - sd_x, ymax = mean_x + sd_x)
    },
    geom = "errorbar",
    width = 0.2,
    color = "black"
  ) +
  scale_y_continuous(breaks = y_breaks) +
  labs(title = "Similarity & Correlation \nacross All Reads (Mean ± SD)", y = "Value", x = NULL) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("input/eventalign_ctl/simtest_allreads.pdf", width = 6, height = 6, dpi = 300, plot = p3, device = cairo_pdf)



