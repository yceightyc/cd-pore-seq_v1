library(ggplot2)
library(tidyr)
library(dplyr)
library(ggridges)
library(knitr)
library(kableExtra)
library(data.table)
library(patchwork)
library(gtools)

event_data <- as.data.frame(read.table("./eventalign_output.tsv", header = TRUE, sep = "\t"))
summarised <- event_data %>%
    filter(model_kmer == reference_kmer) %>%
    group_by(contig,position, read_name) %>%
    summarise(model_kmer = first(model_kmer),
                model_mean = first(model_mean),
                model_stdv = first(model_stdv),
                mean = mean(event_level_mean),
                sum = sum(event_length))

contigs <- mixedsort(unique(summarised$contig))
pdf("signal_dyna/all_contigs_onepage_A4.pdf", width = 8.27, height = 11.69)

# Number of plots to fit
num_plots <- length(contigs)

# Compute number of rows/columns automatically (close to square)
ncol <- ceiling(sqrt(num_plots))
nrow <- ceiling(num_plots / ncol)

# Set plotting layout and margins
par(bg = "white",
    mfrow = c(nrow, ncol),
    mar = c(3, 3, 2, 1),   # smaller margins
    oma = c(2, 2, 2, 2),   # outer margins
    cex.main = 0.8,        # smaller text
    cex.lab = 0.8,
    cex.axis = 0.7)

for (x in contigs) {
  data <- summarised[(summarised$contig == x),] %>%
    group_by(position) %>%
    mutate(mean_mean =  mean(mean),
           q1 = quantile(mean, 0.25),
           q3 = quantile(mean, 0.75),
           sem = sd(mean) / sqrt(length(mean)),
           t_value = qt(0.95, df = length(mean) - 1),
           margin_of_error = t_value * sem,
           lower_ci = mean_mean - margin_of_error,
           upper_ci = mean_mean + margin_of_error)

  data_plot <- distinct(data[c("position","q1", "q3", "model_kmer","model_mean",
                               "mean_mean", "sem", "lower_ci", "upper_ci")], .keep_all = TRUE)
  data_plot_ribo <- data_plot[data_plot$position %in% 77:82, ] %>% arrange(position)
  data_plot_ctl  <- data_plot[data_plot$position %in% 96:101, ] %>% arrange(position)

  actual_value_control <- data_plot_ctl$mean_mean
  upper_bound_control  <- data_plot_ctl$q3
  lower_bound_control  <- data_plot_ctl$q1

  actual_value_ribo <- data_plot_ribo$mean_mean
  upper_bound_ribo  <- data_plot_ribo$q3
  lower_bound_ribo  <- data_plot_ribo$q1

  # Plot
  plot(0, 0,
       ylim = c(50, 180),
       xlim = c(0, 5),
       xlab = "",
       ylab = "",
       main = x,
       xaxt = "n", yaxt = "n",     # suppress both x/y axes
       cex.main = 0.8)

  # Optionally draw tick marks only (no labels)
  axis(side = 1, at = seq(from = 0.5, to = 5.5, by = 1), labels = FALSE, tck = -0.01)
  axis(side = 2, labels = FALSE, tck = -0.01)

  # Draw the lines
  points(c(0:5), actual_value_control, col = "#698a9c", lwd = 2, pch = "-", type = "s")
  points(c(0:5), actual_value_ribo,    col = "#DB7176", lwd = 2, pch = "-", type = "s")

  # CI segments
  for (i in 1:6) {
    segments(x0 = i - 0.5, y0 = upper_bound_control[i],
             x1 = i - 0.5, y1 = lower_bound_control[i],
             col = "black", lwd = 1)
  }
  for (i in 1:6) {
    segments(x0 = i - 0.5, y0 = upper_bound_ribo[i],
             x1 = i - 0.5, y1 = lower_bound_ribo[i],
             col = "black", lwd = 1)
  }

  legend("bottomright", legend = c("Ctl", "Ribo"),
         col = c("#698a9c","#DB7176"),
         fill = c("#698a9c","#DB7176"),
         cex = 0.5, bty = "n")
}

# Add global x/y labels
mtext("Ribonucleotide Containing 9-mer", side = 1, outer = TRUE, line = 1.5, cex = 1)
mtext("Signal Amplitude", side = 2, outer = TRUE, line = 1.5, cex = 1)

dev.off()