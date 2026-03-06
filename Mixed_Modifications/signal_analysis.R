library(ggplot2)
library(tidyr)
library(dplyr)
library(ggridges)
library(knitr)
library(kableExtra)
library(IRdisplay)

event_data <- as.data.frame(read.table("./eventalign_all.tsv", header = TRUE, sep = "\t"))

summarised <- event_data %>%
    filter(model_kmer == reference_kmer) %>%
    group_by(contig,position, read_name) %>%
    summarise(model_kmer = first(model_kmer),
                model_mean = first(model_mean),
                model_stdv = first(model_stdv),
                mean = mean(event_level_mean),
                sum = sum(event_length))

control_data_dA <- summarised[summarised$contig == "dA",] %>%
  group_by(position) %>%
  mutate(mean_mean =  mean(mean),
        sem = sd(mean) / sqrt(length(mean)),
        t_value = qt(0.975, df = length(mean) - 1),
        margin_of_error = t_value * sem,
        lower_ci = mean_mean - margin_of_error,
        upper_ci = mean_mean + margin_of_error)
control_data_dA_plot <- distinct(control_data_dA[c("position","model_kmer","model_mean","mean_mean", "sem", "lower_ci", "upper_ci")], .keep_all = TRUE)
control_data_dA_plot_mod <- control_data_dA_plot[control_data_dA_plot$position %in% 25:36, ] %>% 
  arrange(position)

control_data_dC <- summarised[summarised$contig == "dC",] %>%
  group_by(position) %>%
  mutate(mean_mean =  mean(mean),
        sem = sd(mean) / sqrt(length(mean)),
        t_value = qt(0.975, df = length(mean) - 1),
        margin_of_error = t_value * sem,
        lower_ci = mean_mean - margin_of_error,
        upper_ci = mean_mean + margin_of_error)
control_data_dC_plot <- distinct(control_data_dC[c("position","model_kmer","model_mean","mean_mean", "sem", "lower_ci", "upper_ci")], .keep_all = TRUE)
control_data_dC_plot_mod <- control_data_dC_plot[control_data_dC_plot$position %in% 25:36, ] %>% 
  arrange(position)

control_data_dG <- summarised[summarised$contig == "dG",] %>%
  group_by(position) %>%
  mutate(mean_mean =  mean(mean),
        sem = sd(mean) / sqrt(length(mean)),
        t_value = qt(0.975, df = length(mean) - 1),
        margin_of_error = t_value * sem,
        lower_ci = mean_mean - margin_of_error,
        upper_ci = mean_mean + margin_of_error)
control_data_dG_plot <- distinct(control_data_dG[c("position","model_kmer","model_mean","mean_mean", "sem", "lower_ci", "upper_ci")], .keep_all = TRUE)
control_data_dG_plot_mod <- control_data_dG_plot[control_data_dG_plot$position %in% 25:36, ] %>% 
  arrange(position)

control_data_dT <- summarised[summarised$contig == "dT",] %>%
  group_by(position) %>%
  mutate(mean_mean =  mean(mean),
        sem = sd(mean) / sqrt(length(mean)),
        t_value = qt(0.975, df = length(mean) - 1),
        margin_of_error = t_value * sem,
        lower_ci = mean_mean - margin_of_error,
        upper_ci = mean_mean + margin_of_error)
control_data_dT_plot <- distinct(control_data_dT[c("position","model_kmer","model_mean","mean_mean", "sem", "lower_ci", "upper_ci")], .keep_all = TRUE)
control_data_dT_plot_mod <- control_data_dT_plot[control_data_dT_plot$position %in% 25:36, ] %>% 
  arrange(position)

mod_data_N6mA <- summarised[summarised$contig == "N6mA",] %>%
  group_by(position) %>%
  mutate(mean_mean =  mean(mean),
        sem = sd(mean) / sqrt(length(mean)),
        t_value = qt(0.975, df = length(mean) - 1),
        margin_of_error = t_value * sem,
        lower_ci = mean_mean - margin_of_error,
        upper_ci = mean_mean + margin_of_error)
mod_data_N6mA_plot <- distinct(mod_data_N6mA[c("position","model_kmer","model_mean","mean_mean", "sem", "lower_ci", "upper_ci")], .keep_all = TRUE)
mod_data_N6mA_plot_mod <- mod_data_N6mA_plot[mod_data_N6mA_plot$position %in% 25:36, ] %>% 
  arrange(position)

mod_data_5hmC <- summarised[summarised$contig == "5hmC",] %>%
  group_by(position) %>%
  mutate(mean_mean =  mean(mean),
        sem = sd(mean) / sqrt(length(mean)),
        t_value = qt(0.975, df = length(mean) - 1),
        margin_of_error = t_value * sem,
        lower_ci = mean_mean - margin_of_error,
        upper_ci = mean_mean + margin_of_error)
mod_data_5hmC_plot <- distinct(mod_data_5hmC[c("position","model_kmer","model_mean","mean_mean", "sem", "lower_ci", "upper_ci")], .keep_all = TRUE)
mod_data_5hmC_plot_mod <- mod_data_5hmC_plot[mod_data_5hmC_plot$position %in% 25:36, ] %>% 
  arrange(position)

mod_data_8oxoG <- summarised[summarised$contig == "8oxoG",] %>%
  group_by(position) %>%
  mutate(mean_mean =  mean(mean),
        sem = sd(mean) / sqrt(length(mean)),
        t_value = qt(0.975, df = length(mean) - 1),
        margin_of_error = t_value * sem,
        lower_ci = mean_mean - margin_of_error,
        upper_ci = mean_mean + margin_of_error)
mod_data_8oxoG_plot <- distinct(mod_data_8oxoG[c("position","model_kmer","model_mean","mean_mean", "sem", "lower_ci", "upper_ci")], .keep_all = TRUE)
mod_data_8oxoG_plot_mod <- mod_data_8oxoG_plot[mod_data_8oxoG_plot$position %in% 25:36, ] %>% 
  arrange(position)

mod_data_BrdU <- summarised[summarised$contig == "BrdU",] %>%
  group_by(position) %>%
  mutate(mean_mean =  mean(mean),
        sem = sd(mean) / sqrt(length(mean)),
        t_value = qt(0.975, df = length(mean) - 1),
        margin_of_error = t_value * sem,
        lower_ci = mean_mean - margin_of_error,
        upper_ci = mean_mean + margin_of_error)
mod_data_BrdU_plot <- distinct(mod_data_BrdU[c("position","model_kmer","model_mean","mean_mean", "sem", "lower_ci", "upper_ci")], .keep_all = TRUE)
mod_data_BrdU_plot_mod <- mod_data_BrdU_plot[mod_data_BrdU_plot$position %in% 25:36, ] %>% 
  arrange(position)

mod_data_dU <- summarised[summarised$contig == "dU",] %>%
  group_by(position) %>%
  mutate(mean_mean =  mean(mean),
        sem = sd(mean) / sqrt(length(mean)),
        t_value = qt(0.975, df = length(mean) - 1),
        margin_of_error = t_value * sem,
        lower_ci = mean_mean - margin_of_error,
        upper_ci = mean_mean + margin_of_error)
mod_data_dU_plot <- distinct(mod_data_dU[c("position","model_kmer","model_mean","mean_mean", "sem", "lower_ci", "upper_ci")], .keep_all = TRUE)
mod_data_dU_plot_mod <- mod_data_dU_plot[mod_data_dU_plot$position %in% 25:36, ] %>% 
  arrange(position)

mod_data_dI <- summarised[summarised$contig == "dI",] %>%
  group_by(position) %>%
  mutate(mean_mean =  mean(mean),
        sem = sd(mean) / sqrt(length(mean)),
        t_value = qt(0.975, df = length(mean) - 1),
        margin_of_error = t_value * sem,
        lower_ci = mean_mean - margin_of_error,
        upper_ci = mean_mean + margin_of_error)
mod_data_dI_plot <- distinct(mod_data_dI[c("position","model_kmer","model_mean","mean_mean", "sem", "lower_ci", "upper_ci")], .keep_all = TRUE)
mod_data_dI_plot_mod <- mod_data_dI_plot[mod_data_dI_plot$position %in% 25:36, ] %>% 
  arrange(position)

mod_data_abasic <- summarised[summarised$contig == "didribose",] %>%
  group_by(position) %>%
  mutate(mean_mean =  mean(mean),
        sem = sd(mean) / sqrt(length(mean)),
        t_value = qt(0.975, df = length(mean) - 1),
        margin_of_error = t_value * sem,
        lower_ci = mean_mean - margin_of_error,
        upper_ci = mean_mean + margin_of_error,
        contig = "abasic")
mod_data_abasic_plot <- distinct(mod_data_abasic[c("position","model_kmer","model_mean","mean_mean", "sem", "lower_ci", "upper_ci")], .keep_all = TRUE)
mod_data_abasic_plot_mod <- mod_data_abasic_plot[mod_data_abasic_plot$position %in% 25:36, ] %>% 
  arrange(position)

compare_contigs <- function(mod_name, control_name, n_total = 9) {
  rbind(get(paste0("mod_data_", mod_name)), get(paste0("control_data_",control_name))) %>%
    filter(contig %in% c(mod_name, control_name),
           position %in% 26:34) %>%
    group_by(position) %>%
    summarise(
      ks_p_value = suppressWarnings(
        ks.test(mean[contig == control_name], mean[contig == mod_name])$p.value
      ),
      U_p_value = suppressWarnings(
        wilcox.test(mean[contig == control_name], mean[contig == mod_name])$p.value
      ),
      ks_p_value_time = suppressWarnings(
        ks.test(sum[contig == control_name], sum[contig == mod_name])$p.value
      ),
      U_p_value_time = suppressWarnings(
        wilcox.test(sum[contig == control_name], sum[contig == mod_name])$p.value
      ),
      .groups = "drop"
    ) %>%
    mutate(
      p_adj_ks = p.adjust(ks_p_value, method = "BH", n = n_total),
      significance_ks = case_when(
        p_adj_ks <= 0.001 ~ "***",
        p_adj_ks <= 0.01  ~ "**",
        p_adj_ks <= 0.05  ~ "*",
        TRUE ~ as.character(round(p_adj_ks, 2))
      ),
      p_adj_U = p.adjust(U_p_value, method = "BH", n = n_total),
      significance_U = case_when(
        p_adj_U <= 0.001 ~ "***",
        p_adj_U <= 0.01  ~ "**",
        p_adj_U <= 0.05  ~ "*",
        TRUE ~ as.character(round(p_adj_U, 2))
      ),
      p_adj_ks_time = p.adjust(ks_p_value_time, method = "BH", n = n_total),
      significance_ks_time = case_when(
        p_adj_ks_time <= 0.001 ~ "***",
        p_adj_ks_time <= 0.01  ~ "**",
        p_adj_ks_time <= 0.05  ~ "*",
        TRUE ~ as.character(round(p_adj_ks_time, 2))
      ),
      p_adj_U_time = p.adjust(U_p_value_time, method = "BH", n = n_total),
      significance_U_time = case_when(
        p_adj_U_time <= 0.001 ~ "***",
        p_adj_U_time <= 0.01  ~ "**",
        p_adj_U_time <= 0.05  ~ "*",
        TRUE ~ as.character(round(p_adj_U_time, 2))
      ),
      mod = mod_name,
      control = control_name
    ) %>%
    arrange(position)
}

stats_merged_N6mA <- compare_contigs(mod_name = "N6mA", control_name = "dA")
stats_merged_5hmC <- compare_contigs(mod_name = "5hmC", control_name = "dC")
stats_merged_8oxoG <- compare_contigs(mod_name = "8oxoG", control_name = "dG")
stats_merged_BrdU <- compare_contigs(mod_name = "BrdU", control_name = "dT")
stats_merged_dU <- compare_contigs(mod_name = "dU", control_name = "dT")
stats_merged_dI <- compare_contigs(mod_name = "dI", control_name = "dG")
stats_merged_abasic <- compare_contigs(mod_name = "abasic", control_name = "dT")

actual_value_control = control_data_dA_plot_mod$mean_mean
upper_bound_control =  control_data_dA_plot_mod$upper_ci
lower_bound_control =  control_data_dA_plot_mod$lower_ci

actual_value_mod = mod_data_N6mA_plot_mod$mean_mean
upper_bound_mod =  mod_data_N6mA_plot_mod$upper_ci
lower_bound_mod =  mod_data_N6mA_plot_mod$lower_ci

pdf("./signal_analysis/N6mA_signal.pdf", width = 8, height = 8)
par(bg = "white")
par(cex.main = 1.3, font.main = 8, col.main = "black")   
par(cex.lab = 1.2, col.lab = "black")                 
par(cex.axis = 1, col.axis = "black")                  

plot(0,0, 
     ylim = c(50,150),
     xlim = c(0,11),
     xlab = "Modified Nucleotide Containing 9-mer",
     ylab = "Signal Amplitude",
     cex.axis = 1,
     main = "Signal over Modified/Control Nucleotides",
     xaxt='n')
axis(labels = mod_data_N6mA_plot_mod$model_kmer, at = seq(from = 0.5, to = 11.5, by = 1),side = 1)

points(c(0:11), actual_value_control, col="#698a9c", lwd =3 , pch = "-", type = "s" )
points(c(0:11), actual_value_mod, col="#DB7176", lwd =3 , pch = "-", type = "s" )

for (i in 1:12) {
  segments(x0 = i-0.5, y0 = upper_bound_control[i], 
           x1 = i-0.5, y1 = lower_bound_control[i], col = "black", lwd = 2)
}

for (i in 1:12) {
  segments(x0 = i-0.5, y0 = upper_bound_mod[i], 
           x1 = i-0.5, y1 = lower_bound_mod[i], col = "black", lwd = 2)
}

ypos <- 52  

for (i in 1:9) {
  xpos <- i+1 - 0.5 
  
  text(x = xpos, y = ypos, 
       labels = stats_merged_N6mA$significance_U[i],
       col = "black", cex = 1.2, pos = 3)  
}
legend_labels <- c("dA", "N6mA")
legend("topleft", legend = as.expression(lapply(legend_labels, function(x) bquote(.(x)))),
       col = c("#698a9c","#DB7176"), 
       fill = c("#698a9c","#DB7176"),
       cex = 1.2, bty = "n", text.col = "black", 
       text.font = 5)

dev.off()

actual_value_control = control_data_dC_plot_mod$mean_mean
upper_bound_control =  control_data_dC_plot_mod$upper_ci
lower_bound_control =  control_data_dC_plot_mod$lower_ci

actual_value_mod = mod_data_5hmC_plot_mod$mean_mean
upper_bound_mod =  mod_data_5hmC_plot_mod$upper_ci
lower_bound_mod =  mod_data_5hmC_plot_mod$lower_ci

pdf("./signal_analysis/5hmC_signal.pdf", width = 8, height = 8)
par(bg = "white")
par(cex.main = 1.3, font.main = 8, col.main = "black")   
par(cex.lab = 1.2, col.lab = "black")             
par(cex.axis = 1, col.axis = "black")                

plot(0,0, 
     ylim = c(50,150),
     xlim = c(0,11),
     xlab = "Modified Nucleotide Containing 9-mer",
     ylab = "Signal Amplitude",
     cex.axis = 1,
     main = "Signal over Modified/Control Nucleotides",
     xaxt='n')
axis(labels = mod_data_5hmC_plot_mod$model_kmer, at = seq(from = 0.5, to = 11.5, by = 1),side = 1)

points(c(0:11), actual_value_control, col="#698a9c", lwd =3 , pch = "-", type = "s" )
points(c(0:11), actual_value_mod, col="#DB7176", lwd =3 , pch = "-", type = "s" )

for (i in 1:12) {
  segments(x0 = i-0.5, y0 = upper_bound_control[i], 
           x1 = i-0.5, y1 = lower_bound_control[i], col = "black", lwd = 2)
}

for (i in 1:12) {
  segments(x0 = i-0.5, y0 = upper_bound_mod[i], 
           x1 = i-0.5, y1 = lower_bound_mod[i], col = "black", lwd = 2)
}

ypos <- 52  

for (i in 1:9) {
  xpos <- i+1 - 0.5 
  
  text(x = xpos, y = ypos, 
       labels = stats_merged_5hmC$significance_U[i],
       col = "black", cex = 1.2, pos = 3)  
}
legend_labels <- c("dC", "5hmC")
legend("topleft", legend = as.expression(lapply(legend_labels, function(x) bquote(.(x)))),
       col = c("#698a9c","#DB7176"), 
       fill = c("#698a9c","#DB7176"),
       cex = 1.2, bty = "n", text.col = "black", 
       text.font = 5)

dev.off()

actual_value_control = control_data_dG_plot_mod$mean_mean
upper_bound_control =  control_data_dG_plot_mod$upper_ci
lower_bound_control =  control_data_dG_plot_mod$lower_ci

actual_value_mod = mod_data_8oxoG_plot_mod$mean_mean
upper_bound_mod =  mod_data_8oxoG_plot_mod$upper_ci
lower_bound_mod =  mod_data_8oxoG_plot_mod$lower_ci

pdf("./signal_analysis/8oxoG_signal.pdf", width = 8, height = 8)
par(bg = "white")
par(cex.main = 1.3, font.main = 8, col.main = "black")   
par(cex.lab = 1.2, col.lab = "black")                   
par(cex.axis = 1, col.axis = "black")                 

plot(0,0, 
     ylim = c(50,150),
     xlim = c(0,11),
     xlab = "Modified Nucleotide Containing 9-mer",
     ylab = "Signal Amplitude",
     cex.axis = 1,
     main = "Signal over Modified/Control Nucleotides",
     xaxt='n')
axis(labels = mod_data_8oxoG_plot_mod$model_kmer, at = seq(from = 0.5, to = 11.5, by = 1),side = 1)

points(c(0:11), actual_value_control, col="#698a9c", lwd =3 , pch = "-", type = "s" )
points(c(0:11), actual_value_mod, col="#DB7176", lwd =3 , pch = "-", type = "s" )

for (i in 1:12) {
  segments(x0 = i-0.5, y0 = upper_bound_control[i], 
           x1 = i-0.5, y1 = lower_bound_control[i], col = "black", lwd = 2)
}

for (i in 1:12) {
  segments(x0 = i-0.5, y0 = upper_bound_mod[i], 
           x1 = i-0.5, y1 = lower_bound_mod[i], col = "black", lwd = 2)
}

ypos <- 52  

for (i in 1:9) {
  xpos <- i+1 - 0.5 
  
  text(x = xpos, y = ypos, 
       labels = stats_merged_8oxoG$significance_U[i],
       col = "black", cex = 1.2, pos = 3)  # pos=3 places it *above* ypos
}
legend_labels <- c("dG", "8oxoG")
legend("topleft", legend = as.expression(lapply(legend_labels, function(x) bquote(.(x)))),
       col = c("#698a9c","#DB7176"), 
       fill = c("#698a9c","#DB7176"),
       cex = 1.2, bty = "n", text.col = "black", 
       text.font = 5)

dev.off()

actual_value_control = control_data_dT_plot_mod$mean_mean
upper_bound_control =  control_data_dT_plot_mod$upper_ci
lower_bound_control =  control_data_dT_plot_mod$lower_ci

actual_value_mod = mod_data_BrdU_plot_mod$mean_mean
upper_bound_mod =  mod_data_BrdU_plot_mod$upper_ci
lower_bound_mod =  mod_data_BrdU_plot_mod$lower_ci

pdf("./signal_analysis/BrdU_signal.pdf", width = 8, height = 8)
par(bg = "white")
par(cex.main = 1.3, font.main = 8, col.main = "black")  
par(cex.lab = 1.2, col.lab = "black")                 
par(cex.axis = 1, col.axis = "black")               

plot(0,0, 
     ylim = c(50,150),
     xlim = c(0,11),
     xlab = "Modified Nucleotide Containing 9-mer",
     ylab = "Signal Amplitude",
     cex.axis = 1,
     main = "Signal over Modified/Control Nucleotides",
     xaxt='n')
axis(labels = mod_data_BrdU_plot_mod$model_kmer, at = seq(from = 0.5, to = 11.5, by = 1),side = 1)

points(c(0:11), actual_value_control, col="#698a9c", lwd =3 , pch = "-", type = "s" )
points(c(0:11), actual_value_mod, col="#DB7176", lwd =3 , pch = "-", type = "s" )

for (i in 1:12) {
  segments(x0 = i-0.5, y0 = upper_bound_control[i], 
           x1 = i-0.5, y1 = lower_bound_control[i], col = "black", lwd = 2)
}

for (i in 1:12) {
  segments(x0 = i-0.5, y0 = upper_bound_mod[i], 
           x1 = i-0.5, y1 = lower_bound_mod[i], col = "black", lwd = 2)
}

ypos <- 52  

for (i in 1:9) {
  xpos <- i+1 - 0.5 
  
  text(x = xpos, y = ypos, 
       labels = stats_merged_BrdU$significance_U[i],
       col = "black", cex = 1.2, pos = 3)  # pos=3 places it *above* ypos
}
legend_labels <- c("dT", "BrdU")
legend("topleft", legend = as.expression(lapply(legend_labels, function(x) bquote(.(x)))),
       col = c("#698a9c","#DB7176"), 
       fill = c("#698a9c","#DB7176"),
       cex = 1.2, bty = "n", text.col = "black", 
       text.font = 5)

dev.off()

actual_value_control = control_data_dT_plot_mod$mean_mean
upper_bound_control =  control_data_dT_plot_mod$upper_ci
lower_bound_control =  control_data_dT_plot_mod$lower_ci

actual_value_mod = mod_data_dU_plot_mod$mean_mean
upper_bound_mod =  mod_data_dU_plot_mod$upper_ci
lower_bound_mod =  mod_data_dU_plot_mod$lower_ci

pdf("./signal_analysis/dU_signal.pdf", width = 8, height = 8)
par(bg = "white")
par(cex.main = 1.3, font.main = 8, col.main = "black")  
par(cex.lab = 1.2, col.lab = "black")                  
par(cex.axis = 1, col.axis = "black")                

plot(0,0, 
     ylim = c(50,150),
     xlim = c(0,11),
     xlab = "Modified Nucleotide Containing 9-mer",
     ylab = "Signal Amplitude",
     cex.axis = 1,
     main = "Signal over Modified/Control Nucleotides",
     xaxt='n')
axis(labels = mod_data_dU_plot_mod$model_kmer, at = seq(from = 0.5, to = 11.5, by = 1),side = 1)

points(c(0:11), actual_value_control, col="#698a9c", lwd =3 , pch = "-", type = "s" )
points(c(0:11), actual_value_mod, col="#DB7176", lwd =3 , pch = "-", type = "s" )

for (i in 1:12) {
  segments(x0 = i-0.5, y0 = upper_bound_control[i], 
           x1 = i-0.5, y1 = lower_bound_control[i], col = "black", lwd = 2)
}

for (i in 1:12) {
  segments(x0 = i-0.5, y0 = upper_bound_mod[i], 
           x1 = i-0.5, y1 = lower_bound_mod[i], col = "black", lwd = 2)
}

ypos <- 52  

for (i in 1:9) {
  xpos <- i+1 - 0.5 
  
  text(x = xpos, y = ypos, 
       labels = stats_merged_dU$significance_U[i],
       col = "black", cex = 1.2, pos = 3)  # pos=3 places it *above* ypos
}
legend_labels <- c("dT", "dU")
legend("topleft", legend = as.expression(lapply(legend_labels, function(x) bquote(.(x)))),
       col = c("#698a9c","#DB7176"), 
       fill = c("#698a9c","#DB7176"),
       cex = 1.2, bty = "n", text.col = "black", 
       text.font = 5)

dev.off()

actual_value_control = control_data_dT_plot_mod$mean_mean
upper_bound_control =  control_data_dT_plot_mod$upper_ci
lower_bound_control =  control_data_dT_plot_mod$lower_ci

actual_value_mod = mod_data_abasic_plot_mod$mean_mean
upper_bound_mod =  mod_data_abasic_plot_mod$upper_ci
lower_bound_mod =  mod_data_abasic_plot_mod$lower_ci

pdf("./signal_analysis/abasic_signal.pdf", width = 8, height = 8)
par(bg = "white")
par(cex.main = 1.3, font.main = 8, col.main = "black")   
par(cex.lab = 1.2, col.lab = "black")                   
par(cex.axis = 1, col.axis = "black")                 

plot(0,0, 
     ylim = c(50,150),
     xlim = c(0,11),
     xlab = "Modified Nucleotide Containing 9-mer",
     ylab = "Signal Amplitude",
     cex.axis = 1,
     main = "Signal over Modified/Control Nucleotides",
     xaxt='n')
axis(labels = mod_data_abasic_plot_mod$model_kmer, at = seq(from = 0.5, to = 11.5, by = 1),side = 1)

points(c(0:11), actual_value_control, col="#698a9c", lwd =3 , pch = "-", type = "s" )
points(c(0:11), actual_value_mod, col="#DB7176", lwd =3 , pch = "-", type = "s" )

for (i in 1:12) {
  segments(x0 = i-0.5, y0 = upper_bound_control[i], 
           x1 = i-0.5, y1 = lower_bound_control[i], col = "black", lwd = 2)
}

for (i in 1:12) {
  segments(x0 = i-0.5, y0 = upper_bound_mod[i], 
           x1 = i-0.5, y1 = lower_bound_mod[i], col = "black", lwd = 2)
}

ypos <- 52  

for (i in 1:9) {
  xpos <- i+1 - 0.5 
  
  text(x = xpos, y = ypos, 
       labels = stats_merged_dU$significance_U[i],
       col = "black", cex = 1.2, pos = 3)  
}
legend_labels <- c("dT", "abasic")
legend("topleft", legend = as.expression(lapply(legend_labels, function(x) bquote(.(x)))),
       col = c("#698a9c","#DB7176"), 
       fill = c("#698a9c","#DB7176"),
       cex = 1.2, bty = "n", text.col = "black", 
       text.font = 5)

dev.off()

actual_value_control = control_data_dG_plot_mod$mean_mean
upper_bound_control =  control_data_dG_plot_mod$upper_ci
lower_bound_control =  control_data_dG_plot_mod$lower_ci

actual_value_mod = mod_data_dI_plot_mod$mean_mean
upper_bound_mod =  mod_data_dI_plot_mod$upper_ci
lower_bound_mod =  mod_data_dI_plot_mod$lower_ci

pdf("./signal_analysis/dI_signal.pdf", width = 8, height = 8)
par(bg = "white")
par(cex.main = 1.3, font.main = 8, col.main = "black")   
par(cex.lab = 1.2, col.lab = "black")                    
par(cex.axis = 1, col.axis = "black")                

plot(0,0, 
     ylim = c(50,150),
     xlim = c(0,11),
     xlab = "Modified Nucleotide Containing 9-mer",
     ylab = "Signal Amplitude",
     cex.axis = 1,
     main = "Signal over Modified/Control Nucleotides",
     xaxt='n')
axis(labels = mod_data_dI_plot_mod$model_kmer, at = seq(from = 0.5, to = 11.5, by = 1),side = 1)

points(c(0:11), actual_value_control, col="#698a9c", lwd =3 , pch = "-", type = "s" )
points(c(0:11), actual_value_mod, col="#DB7176", lwd =3 , pch = "-", type = "s" )

for (i in 1:12) {
  segments(x0 = i-0.5, y0 = upper_bound_control[i], 
           x1 = i-0.5, y1 = lower_bound_control[i], col = "black", lwd = 2)
}

for (i in 1:12) {
  segments(x0 = i-0.5, y0 = upper_bound_mod[i], 
           x1 = i-0.5, y1 = lower_bound_mod[i], col = "black", lwd = 2)
}

ypos <- 52  

for (i in 1:9) {
  xpos <- i+1 - 0.5 
  
  text(x = xpos, y = ypos, 
       labels = stats_merged_dI$significance_U[i],
       col = "black", cex = 1.2, pos = 3) 
}
legend_labels <- c("dG", "dI")
legend("topleft", legend = as.expression(lapply(legend_labels, function(x) bquote(.(x)))),
       col = c("#698a9c","#DB7176"), 
       fill = c("#698a9c","#DB7176"),
       cex = 1.2, bty = "n", text.col = "black", 
       text.font = 5)

dev.off()

dA_time_data_plot <- control_data_dA %>%
    filter(position %in% 26:34) %>%
    select(c("contig","position","sum")) %>%
    mutate(Kmer_Position = position - 25, Dwell_Time = sum)

N6mA_time_data_plot <- mod_data_N6mA %>%
    filter(position %in% 26:34) %>%
    select(c("contig","position","sum")) %>%
    mutate(Kmer_Position = position - 25, Dwell_Time = sum)

N6mA_dwell_time_plot <- rbind(dA_time_data_plot, N6mA_time_data_plot) %>%
    arrange(Kmer_Position) %>%
    mutate(Kmer_Position = factor(Kmer_Position), Type = contig)

summarise_dwell_time <- function(mod_name, control_name) {
  rbind(get(paste0("mod_data_", mod_name)), get(paste0("control_data_",control_name))) %>%
    filter(position %in% 26:34) %>%
    select(c("contig","position","sum")) %>%
    mutate(Kmer_Position = position - 25, Dwell_Time = sum) %>%
    arrange(Kmer_Position) %>%
    mutate(Kmer_Position = factor(Kmer_Position), Type = contig)
}

dwell_time_plot_N6mA <- summarise_dwell_time(mod_name ="N6mA", control_name = "dA")

dwell_time_plot_5hmC <- summarise_dwell_time(mod_name ="5hmC", control_name = "dC")

dwell_time_plot_8oxoG <- summarise_dwell_time(mod_name ="8oxoG", control_name = "dG")

dwell_time_plot_dI <- summarise_dwell_time(mod_name ="dI", control_name = "dG")

dwell_time_plot_BrdU <- summarise_dwell_time(mod_name ="BrdU", control_name = "dT")

dwell_time_plot_dU <- summarise_dwell_time(mod_name ="dU", control_name = "dT")

dwell_time_plot_abasic <- summarise_dwell_time(mod_name ="abasic", control_name = "dT")

theme_set(
  theme_classic(base_size = 18, base_family = "serif") +
    theme(legend.position = "top",
  axis.title = element_text(size = 20, color = "black"),
  axis.text = element_text(size = 18, color = "black"),
  legend.text = element_text(size = 18, color = "black"),
  legend.title = element_text(size = 18, color = "black"),
  axis.ticks.y = element_blank(),
  plot.title = element_text(size = 18, face = "bold", color = "black")
    ))

custom_colors <- c("dA" = "#698a9c", "N6mA" = "#E75C63")
plot <- ggplot(dwell_time_plot_N6mA, aes(x = Dwell_Time, y = Kmer_Position, fill = Type)) +
  geom_density_ridges(
    position = position_points_jitter(height = 0),
    alpha = 0.8,
    scale = 1.2,
    color = NA
  ) +
  coord_cartesian(xlim = c(0, 0.03)) +
  scale_fill_manual(values = custom_colors) +
  geom_text(
    data = stats_merged_N6mA,
    aes(x = 0.03 , y = position-25, label = significance_ks_time),
    inherit.aes = FALSE,
    angle = 90,      
    hjust = 0.5,        
    vjust = 0,      
    size = 6,
    family = "serif") 

ggsave("./signal_analysis/N6mA_dwell_time.pdf", plot = plot)

dwell_time_plot_5hmC <- dwell_time_plot_5hmC %>%
  mutate(Type = factor(Type, levels = c("dC", "5hmC")))
  
custom_colors <- c("dC" = "#698a9c", "5hmC" = "#E75C63")
plot <- ggplot(dwell_time_plot_5hmC, aes(x = Dwell_Time, y = Kmer_Position, fill = Type)) +
  geom_density_ridges(
    position = position_points_jitter(height = 0),
    alpha = 0.8,
    scale = 1.2,
    color = NA
  ) +
  coord_cartesian(xlim = c(0, 0.03)) +
  scale_fill_manual(values = custom_colors) +
  geom_text(
    data = stats_merged_5hmC,
    aes(x = 0.03 , y = position-25, label = significance_ks_time),
    inherit.aes = FALSE,
    angle = 90,         
    hjust = 0.5,         
    vjust = 0,        
    size = 6,
    family = "serif") 

ggsave("./signal_analysis/5hmC_dwell_time.pdf", plot = plot)

dwell_time_plot_8oxoG <- dwell_time_plot_8oxoG %>%
  mutate(Type = factor(Type, levels = c("dG", "8oxoG")))
  
custom_colors <- c("dG" = "#698a9c", "8oxoG" = "#E75C63")
plot <- ggplot(dwell_time_plot_8oxoG, aes(x = Dwell_Time, y = Kmer_Position, fill = Type)) +
  geom_density_ridges(
    position = position_points_jitter(height = 0),
    alpha = 0.8,
    scale = 1.2,
    color = NA
  ) +
  coord_cartesian(xlim = c(0, 0.03)) +
  scale_fill_manual(values = custom_colors) +
  geom_text(
    data = stats_merged_8oxoG,
    aes(x = 0.03 , y = position-25, label = significance_ks_time),
    inherit.aes = FALSE,
    angle = 90,        
    hjust = 0.5,        
    vjust = 0,      
    size = 6,
    family = "serif") 

ggsave("./signal_analysis/8oxoG_dwell_time.pdf", plot = plot)

dwell_time_plot_dI <- dwell_time_plot_dI %>%
  mutate(Type = factor(Type, levels = c("dG", "dI")))
  
custom_colors <- c("dG" = "#698a9c", "dI" = "#E75C63")
plot <- ggplot(dwell_time_plot_dI, aes(x = Dwell_Time, y = Kmer_Position, fill = Type)) +
  geom_density_ridges(
    position = position_points_jitter(height = 0),
    alpha = 0.8,
    scale = 1.2,
    color = NA
  ) +
  coord_cartesian(xlim = c(0, 0.03)) +
  scale_fill_manual(values = custom_colors) +
  geom_text(
    data = stats_merged_dI,
    aes(x = 0.03 , y = position-25, label = significance_ks_time),
    inherit.aes = FALSE,
    angle = 90,         
    hjust = 0.5,     
    vjust = 0,         
    size = 6,
    family = "serif") 

ggsave("./signal_analysis/dI_dwell_time.pdf", plot = plot)

dwell_time_plot_BrdU <- dwell_time_plot_BrdU %>%
  mutate(Type = factor(Type, levels = c("dT", "BrdU")))
  
custom_colors <- c("dT" = "#698a9c", "BrdU" = "#E75C63")
plot <- ggplot(dwell_time_plot_BrdU, aes(x = Dwell_Time, y = Kmer_Position, fill = Type)) +
  geom_density_ridges(
    position = position_points_jitter(height = 0),
    alpha = 0.8,
    scale = 1.2,
    color = NA
  ) +
  coord_cartesian(xlim = c(0, 0.03)) +
  scale_fill_manual(values = custom_colors) +
  geom_text(
    data = stats_merged_BrdU,
    aes(x = 0.03 , y = position-25, label = significance_ks_time),
    inherit.aes = FALSE,
    angle = 90,        
    hjust = 0.5,       
    vjust = 0,        
    size = 6,
    family = "serif") 
ggsave("./signal_analysis/BrdU_dwell_time.pdf", plot = plot)

dwell_time_plot_dU <- dwell_time_plot_dU %>%
  mutate(Type = factor(Type, levels = c("dT", "dU")))
  
custom_colors <- c("dT" = "#698a9c", "dU" = "#E75C63")
plot <- ggplot(dwell_time_plot_dU, aes(x = Dwell_Time, y = Kmer_Position, fill = Type)) +
  geom_density_ridges(
    position = position_points_jitter(height = 0),
    alpha = 0.8,
    scale = 1.2,
    color = NA
  ) +
  coord_cartesian(xlim = c(0, 0.03)) +
  scale_fill_manual(values = custom_colors) +
  geom_text(
    data = stats_merged_dU,
    aes(x = 0.03 , y = position-25, label = significance_ks_time),
    inherit.aes = FALSE,
    angle = 90,        
    hjust = 0.5,       
    vjust = 0,         
    size = 6,
    family = "serif") 

ggsave("./signal_analysis/dU_dwell_time.pdf", plot = plot)

dwell_time_plot_abasic <- dwell_time_plot_abasic %>%
  mutate(Type = factor(Type, levels = c("dT", "abasic")))
  
custom_colors <- c("dT" = "#698a9c", "abasic" = "#E75C63")
plot <- ggplot(dwell_time_plot_abasic, aes(x = Dwell_Time, y = Kmer_Position, fill = Type)) +
  geom_density_ridges(
    position = position_points_jitter(height = 0),
    alpha = 0.8,
    scale = 1.2,
    color = NA
  ) +
  coord_cartesian(xlim = c(0, 0.03)) +
  scale_fill_manual(values = custom_colors) +
  geom_text(
    data = stats_merged_abasic,
    aes(x = 0.03 , y = position-25, label = significance_ks_time),
    inherit.aes = FALSE,
    angle = 90,     
    hjust = 0.5,        
    vjust = 0,     
    size = 6,
    family = "serif") 

ggsave("./signal_analysis/abasic_dwell_time.pdf", plot = plot)

stats_all <- rbind(stats_merged_N6mA, stats_merged_5hmC, stats_merged_8oxoG, stats_merged_dI, stats_merged_BrdU, stats_merged_dU, stats_merged_abasic)

stats_fmt <- stats_all %>%
  mutate(
    Kmer = position - 25,
    Amplitude = p_adj_U ,
    Dwell_time = p_adj_U_time 
  ) %>%
  select(mod, control, Kmer, Amplitude, Dwell_time)

stats_wide <- stats_fmt %>%
  pivot_wider(
    id_cols = c(mod, control),
    names_from = Kmer,
    values_from = c(Amplitude, Dwell_time),
    names_glue = "{Kmer}_{.value}"
  )

stats_counts <- stats_all %>%
  mutate(Kmer = position - 25) %>%
  group_by(mod, control) %>%
  summarise(Num_significant = sum(p_adj_U <= 0.05 | p_adj_U_time <= 0.05),
            .groups = "drop")

stats_wide <- left_join(stats_wide, stats_counts, by = c("mod","control"))

tbl <- stats_wide %>%
  kbl(align = "c", booktabs = TRUE,
      col.names = c("Modification", "Control",
                  rep(c("Amplitude", "Dwell time"), 9),
                  "Num. Significant (≤0.05)"), escape = FALSE) %>%
  add_header_above(c(" " = 2,
                     "Kmer Position 1" = 2,
                     "Kmer Position 2" = 2,
                     "Kmer Position 3" = 2,
                     "Kmer Position 4" = 2,
                     "Kmer Position 5" = 2,
                     "Kmer Position 6" = 2,
                     "Kmer Position 7" = 2,
                     "Kmer Position 8" = 2,
                     "Kmer Position 9" = 2,
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
outfile <- "./signal_analysis/signal_stats_table.html"
writeLines(html_doc, outfile, useBytes = TRUE)

summarised_new <- summarised %>%
    mutate(model_deviation = mean - model_mean) %>%
    filter(position %in% 26:34)

colnames(summarised_new) = c("contig", "position", "read_name", "model_kmer", "model_mean", "model_stdv", "ampl_mean", "dwell_time", "model_deviation")
summarised_new$position <- recode(summarised_new$position, "26" = "1st", "27" = "2nd", "28" = "3rd", "29" = "4th", "30" = "5th", "31" = "6th", "32" = "7th", "33" = "8th", "34" = "9th")
pivoted_all <- summarised_new[, c("contig", "position","read_name","ampl_mean","dwell_time","model_deviation")] %>%
    arrange(read_name) %>%
    pivot_wider(names_from = c(position), values_from = c(ampl_mean, dwell_time, model_deviation))
write.table(pivoted_all, "./signal_analysis/training_input.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
