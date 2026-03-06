library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(gtools)


# Import argument
args <- commandArgs(trailingOnly = TRUE)
degen <- args[1]
output_path <- args[2]
degen_file <- paste0("./", degen)

data <- read.table(degen_file, header = TRUE, sep = "\t")
theme_set(
	    theme_classic(base_size = 18, base_family = "Times New Roman") +
		        theme(legend.position = "top",
			        axis.title = element_text(size = 20, color = "black"),
				  axis.text = element_text(size = 18, color = "black"),
				  legend.text = element_text(size = 18, color = "black"),
				    legend.title = element_text(size = 18, color = "black"),
				    plot.title = element_text(size = 18, face = "bold", color = "black")
				        )
	    )

summary_stats_ci <- data %>%
	  group_by(UMI_simi) %>%
	    summarise(
		          mean_rate = mean(match_rate),
			      sd_rate = sd(match_rate),
			      n = n(),
			          se = sd_rate / sqrt(n),
			          lower = mean_rate - 1.96 * se,
				      upper = mean_rate + 1.96 * se
				    ) %>%
  arrange(desc(UMI_simi))

x_breaks <- seq(min(summary_stats_ci$UMI_simi), max(summary_stats_ci$UMI_simi), by = 1)

# Finer y-axis breaks (e.g., 0% to 100% every 5%)
y_breaks <- seq(0, 1, by = 0.05)

# Plot using ggplot2
p <- ggplot(summary_stats_ci, aes(x = UMI_simi, y = mean_rate)) +
	  geom_line(size = 1, color = '#698A9C') +
	    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
	      geom_point(size = 3, color = '#DB7176') +
	        labs(
		         x = "UMI Matching Bases",
			     y = "Accuracy (Mean ± 95% CI)"
			   ) +
  scale_x_reverse(breaks = sort(x_breaks, decreasing = TRUE)) +
    scale_y_continuous(breaks = y_breaks)

ggsave(output_path, p,device = cairo_pdf, width = 5, height = 5, dpi = 300)

