library(dplyr)
library(tidyr)
library(data.table)
library(gtools)

args <- commandArgs(trailingOnly = TRUE)

input_batch_real <- args[1]
input_batch_sim <- args[2]
output_path <- args[3]

teln_ribo_real <- fread(
  paste0("./", input_batch_real),
  sep = "\t",
  header = TRUE
)

teln_ribo_sim <- fread(
  paste0("./", input_batch_sim),
  sep = "\t",
  header = TRUE
)


teln_ribo_real_mut <- teln_ribo_real %>%
    filter(model_kmer != "NNNNNNNNN", position %in% 20:233) %>%
    mutate(model_deviation = event_level_mean - model_mean,
    contig_short = sub("_.*", "", contig)
  ) %>%
  as.data.table()
#teln_ribo_real_mut[, model_deviation_filtered := ifelse(model_kmer != "NNNNNNNNN", model_deviation, NA)]

teln_ribo_real_final <- teln_ribo_real_mut[, .(
  event_level_mean_mean = mean(event_level_mean),
  model_deviation_mean = mean(model_deviation),
  dwell_time = sum(event_length),
  contig = first(contig_short)
), by = .(position, read_name)]

teln_ribo_sim_mut <- teln_ribo_sim %>%
    filter(model_kmer != "NNNNNNNNN", position %in% 20:233) %>%
    mutate(model_deviation = event_level_mean - model_mean,
    contig_short = sub("_.*", "", contig)
  ) %>%
  as.data.table()
#teln_ribo_real_mut[, model_deviation_filtered := ifelse(model_kmer != "NNNNNNNNN", model_deviation, NA)]

teln_ribo_sim_final <- teln_ribo_sim_mut[, .(
  event_level_mean_mean = mean(event_level_mean),
  model_deviation_mean = mean(model_deviation),
  dwell_time = sum(event_length),
  contig = first(contig_short)
), by = .(position, read_name)]

amplitude_ribo_real <- teln_ribo_real_final %>%
    filter (position %in% 62:70) %>%
    select (position, read_name, event_level_mean_mean, model_deviation_mean, contig)

dwell_time_offset_ribo_real <- teln_ribo_real_final %>%
    filter (position %in% c(52, 55)) %>%
    select (position, read_name,dwell_time, contig)

amplitude_wide_ribo_real <- amplitude_ribo_real %>%
    pivot_wider(names_from = position, values_from = c(event_level_mean_mean, model_deviation_mean)) %>%
    mutate(ribo = 1)

colnames(amplitude_wide_ribo_real) = c("read_name","contig","level_mean_26","level_mean_27","level_mean_28","level_mean_29","level_mean_30","level_mean_31","level_mean_32","level_mean_33","level_mean_34","model_deviation_mean_26","model_deviation_mean_27","model_deviation_mean_28","model_deviation_mean_29","model_deviation_mean_30","model_deviation_mean_31","model_deviation_mean_32","model_deviation_mean_33","model_deviation_mean_34","ribo")

time_wide_ribo_real <- dwell_time_offset_ribo_real %>%
    pivot_wider(names_from = position, values_from = c(dwell_time))

colnames(time_wide_ribo_real) = c("read_name","contig","dwell_time_-18","dwell_time_-15")

all_ribo_real <- left_join(amplitude_wide_ribo_real, time_wide_ribo_real, by = c("read_name", "contig"))

amplitude_riboup1_real <- teln_ribo_real_final %>%
    filter (position %in% 61:69) %>%
    select (position, read_name, event_level_mean_mean, model_deviation_mean, contig)

dwell_time_offset_riboup1_real <- teln_ribo_real_final %>%
    filter (position %in% c(51, 54)) %>%
    select (position, read_name,dwell_time, contig)

amplitude_wide_riboup1_real <- amplitude_riboup1_real %>%
    pivot_wider(names_from = position, values_from = c(event_level_mean_mean, model_deviation_mean)) %>%
    mutate(ribo = 0) %>%
    select(c(read_name, contig), mixedsort(setdiff(names(.), c("read_name", "contig"))), ribo)

colnames(amplitude_wide_riboup1_real) = c("read_name","contig","level_mean_26","level_mean_27","level_mean_28","level_mean_29","level_mean_30","level_mean_31","level_mean_32","level_mean_33","level_mean_34","model_deviation_mean_26","model_deviation_mean_27","model_deviation_mean_28","model_deviation_mean_29","model_deviation_mean_30","model_deviation_mean_31","model_deviation_mean_32","model_deviation_mean_33","model_deviation_mean_34","ribo")

time_wide_riboup1_real <- dwell_time_offset_riboup1_real %>%
    pivot_wider(names_from = position, values_from = c(dwell_time)) %>%
    select(c(read_name, contig), mixedsort(setdiff(names(.), c("read_name", "contig"))))

colnames(time_wide_riboup1_real) = c("read_name","contig","dwell_time_-18","dwell_time_-15")

all_riboup1_real <- left_join(amplitude_wide_riboup1_real, time_wide_riboup1_real, by = c("read_name", "contig"))
all_riboup1_real$read_name <- paste0(all_riboup1_real$read_name, "_upa")

amplitude_ribodown1_real <- teln_ribo_real_final %>%
    filter (position %in% 63:71) %>%
    select (position, read_name, event_level_mean_mean, model_deviation_mean, contig)

dwell_time_offset_ribodown1_real <- teln_ribo_real_final %>%
    filter (position %in% c(53, 56)) %>%
    select (position, read_name,dwell_time, contig)

amplitude_wide_ribodown1_real <- amplitude_ribodown1_real %>%
    pivot_wider(names_from = position, values_from = c(event_level_mean_mean, model_deviation_mean)) %>%
    mutate(ribo = 0) %>%
    select(c(read_name, contig), mixedsort(setdiff(names(.), c("read_name", "contig"))), ribo)


colnames(amplitude_wide_ribodown1_real) = c("read_name","contig","level_mean_26","level_mean_27","level_mean_28","level_mean_29","level_mean_30","level_mean_31","level_mean_32","level_mean_33","level_mean_34","model_deviation_mean_26","model_deviation_mean_27","model_deviation_mean_28","model_deviation_mean_29","model_deviation_mean_30","model_deviation_mean_31","model_deviation_mean_32","model_deviation_mean_33","model_deviation_mean_34","ribo")

time_wide_ribodown1_real <- dwell_time_offset_ribodown1_real %>%
    pivot_wider(names_from = position, values_from = c(dwell_time)) %>%
    select(c(read_name, contig), mixedsort(setdiff(names(.), c("read_name", "contig"))))


colnames(time_wide_ribodown1_real) = c("read_name","contig","dwell_time_-18","dwell_time_-15")


all_ribodown1_real <- left_join(amplitude_wide_ribodown1_real, time_wide_ribodown1_real, by = c("read_name", "contig"))
all_ribodown1_real$read_name <- paste0(all_ribodown1_real$read_name, "_downa")

amplitude_ribodown2_real <- teln_ribo_real_final %>%
    filter (position %in% 64:72) %>%
    select (position, read_name, event_level_mean_mean, model_deviation_mean, contig)

dwell_time_offset_ribodown2_real <- teln_ribo_real_final %>%
    filter (position %in% c(54, 57)) %>%
    select (position, read_name,dwell_time, contig)

amplitude_wide_ribodown2_real <- amplitude_ribodown2_real %>%
    pivot_wider(names_from = position, values_from = c(event_level_mean_mean, model_deviation_mean)) %>%
    mutate(ribo = 0) %>%
    select(c(read_name, contig), mixedsort(setdiff(names(.), c("read_name", "contig"))), ribo)


colnames(amplitude_wide_ribodown2_real) = c("read_name","contig","level_mean_26","level_mean_27","level_mean_28","level_mean_29","level_mean_30","level_mean_31","level_mean_32","level_mean_33","level_mean_34","model_deviation_mean_26","model_deviation_mean_27","model_deviation_mean_28","model_deviation_mean_29","model_deviation_mean_30","model_deviation_mean_31","model_deviation_mean_32","model_deviation_mean_33","model_deviation_mean_34","ribo")

time_wide_ribodown2_real <- dwell_time_offset_ribodown2_real %>%
    pivot_wider(names_from = position, values_from = c(dwell_time)) %>%
    select(c(read_name, contig), mixedsort(setdiff(names(.), c("read_name", "contig"))))


colnames(time_wide_ribodown2_real) = c("read_name","contig","dwell_time_-18","dwell_time_-15")


all_ribodown2_real <- left_join(amplitude_wide_ribodown2_real, time_wide_ribodown2_real, by = c("read_name", "contig"))
all_ribodown2_real$read_name <- paste0(all_ribodown2_real$read_name, "_downb")

amplitude_ribodown3_real <- teln_ribo_real_final %>%
    filter (position %in% 65:73) %>%
    select (position, read_name, event_level_mean_mean, model_deviation_mean, contig)

dwell_time_offset_ribodown3_real <- teln_ribo_real_final %>%
    filter (position %in% c(55, 58)) %>%
    select (position, read_name,dwell_time, contig)

amplitude_wide_ribodown3_real <- amplitude_ribodown3_real %>%
    pivot_wider(names_from = position, values_from = c(event_level_mean_mean, model_deviation_mean)) %>%
    mutate(ribo = 0) %>%
    select(c(read_name, contig), mixedsort(setdiff(names(.), c("read_name", "contig"))), ribo)


colnames(amplitude_wide_ribodown3_real) = c("read_name","contig","level_mean_26","level_mean_27","level_mean_28","level_mean_29","level_mean_30","level_mean_31","level_mean_32","level_mean_33","level_mean_34","model_deviation_mean_26","model_deviation_mean_27","model_deviation_mean_28","model_deviation_mean_29","model_deviation_mean_30","model_deviation_mean_31","model_deviation_mean_32","model_deviation_mean_33","model_deviation_mean_34","ribo")

time_wide_ribodown3_real <- dwell_time_offset_ribodown3_real %>%
    pivot_wider(names_from = position, values_from = c(dwell_time)) %>%
    select(c(read_name, contig), mixedsort(setdiff(names(.), c("read_name", "contig"))))


colnames(time_wide_ribodown3_real) = c("read_name","contig","dwell_time_-18","dwell_time_-15")
all_ribodown3_real <- left_join(amplitude_wide_ribodown3_real, time_wide_ribodown3_real, by = c("read_name", "contig"))
all_ribodown3_real$read_name <- paste0(all_ribodown3_real$read_name, "_downc")
surrounding_deoxy <- rbind(all_riboup1_real,all_ribodown1_real,all_ribodown2_real, all_ribodown3_real)

amplitude_ribo_sim <- teln_ribo_sim_final %>%
    filter (position %in% 62:70) %>%
    select (position, read_name, event_level_mean_mean, model_deviation_mean, contig)

dwell_time_offset_ribo_sim <- teln_ribo_sim_final %>%
    filter (position %in% c(181, 184)) %>%
    select (position, read_name,dwell_time, contig)

amplitude_wide_ribo_sim <- amplitude_ribo_sim %>%
    pivot_wider(names_from = position, values_from = c(event_level_mean_mean, model_deviation_mean)) %>%
    mutate(ribo = 0) %>%
    select(c(read_name, contig), mixedsort(setdiff(names(.), c("read_name", "contig"))), ribo)

colnames(amplitude_wide_ribo_sim) = c("read_name","contig","level_mean_26","level_mean_27","level_mean_28","level_mean_29","level_mean_30","level_mean_31","level_mean_32","level_mean_33","level_mean_34","model_deviation_mean_26","model_deviation_mean_27","model_deviation_mean_28","model_deviation_mean_29","model_deviation_mean_30","model_deviation_mean_31","model_deviation_mean_32","model_deviation_mean_33","model_deviation_mean_34","ribo")

time_wide_ribo_sim <- dwell_time_offset_ribo_sim %>%
    pivot_wider(names_from = position, values_from = c(dwell_time)) %>%
    select(c(read_name, contig), mixedsort(setdiff(names(.), c("read_name", "contig"))))
colnames(time_wide_ribo_sim) = c("read_name","contig","dwell_time_-18","dwell_time_-15")

all_ribo_sim <- left_join(amplitude_wide_ribo_sim, time_wide_ribo_sim, by = c("read_name", "contig"))

new_RF_input_nosimsurrounding <- rbind(all_ribo_sim, all_ribo_real, surrounding_deoxy)
write.table(new_RF_input_nosimsurrounding, file = file.path(output_path), row.names = FALSE, quote = FALSE, sep = "\t")
