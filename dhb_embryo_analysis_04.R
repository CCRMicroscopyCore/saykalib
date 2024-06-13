# DHB embryo analysis for Bechara Saykali 
# Set 03
# Use output from dhb_embryo_analysis.py 
# 03/21/23 

# Author: Andy D. Tran, CCR Microscopy Core, LCBG, CCR, NCI 

# Libraries and themes---------------------------------------------------------

library(tidyverse)
library(ggsci)

theme <-  theme(panel.grid.major=element_blank(), 
                panel.grid.minor=element_blank(),
                panel.background=element_blank(),
                axis.line=element_line(color="black", linewidth=1), 
                axis.ticks=element_line(color="black", linewidth=1),
                text=element_text(size=18), 
                plot.title=element_text(size=18), 
                axis.text=element_text(color="black"),
                plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

# Define functions-------------------------------------------------------------

nuc_open <- function(img_name){
  img_path <- file.path(set_path, img_name) 
  csv_list <- list.files(img_path, pattern = '_nuc_stats') 
  for(csv in csv_list){
    csv_path <- file.path(img_path, csv) 
    name_tmp <- str_replace(csv, '_nuc_stats_', '-')
    tmp <- read_csv(csv_path) %>% 
      mutate(image = str_split(name_tmp, '-')[[1]][1]) %>% 
      mutate(z_plane = str_replace(str_split(name_tmp, '-')[[1]][2], '.csv', '')) %>%
      mutate(z_plane = as.integer(z_plane)) %>% 
      mutate(id = as.integer(id), area = as.numeric(area),
             ruby_int = as.numeric(ruby_int), clover_int = as.numeric(clover_int),
             ch405_int = as.numeric(ch405_int),
             ch647_int = as.numeric(ch647_int)) %>% 
      select(image, id, z_plane, area, ruby_int, clover_int, ch405_int, ch647_int) 
    df_tmp <- full_join(df_tmp, tmp)
  }
  return(df_tmp)
}

cyto_open <- function(img_name){
  img_path <- file.path(set_path, img_name)
  csv_list <- list.files(img_path, pattern = '_cyto_stats') 
  for(csv in csv_list){
    csv_path <- file.path(img_path, csv) 
    name_tmp <- str_replace(csv, '_cyto_stats_', '-')
    tmp <- read_csv(csv_path) %>% 
      mutate(image = str_split(name_tmp, '-')[[1]][1]) %>% 
      mutate(z_plane = str_replace(str_split(name_tmp, '-')[[1]][2], '.csv', '')) %>%
      mutate(z_plane = as.integer(z_plane)) %>% 
      rename(cyto_area = area) %>% 
      rename(cyto_clover_int = clover_int) %>% 
      rename(cyto_ruby_int = ruby_int) %>% 
      rename(cyto_ch405_int = ch405_int) %>% 
      rename(cyto_ch647_int = ch647_int) %>% 
      mutate(id = as.integer(id), cyto_area = as.numeric(cyto_area),
             cyto_clover_int = as.numeric(cyto_clover_int),
             cyto_ruby_int = as.numeric(cyto_ruby_int),
             cyto_ch405_int = as.numeric(cyto_ch405_int),
             cyto_ch647_int = as.numeric(cyto_ch647_int)) %>% 
      select(image, id, z_plane, cyto_area, cyto_ruby_int, cyto_clover_int,
             cyto_ch405_int, cyto_ch647_int)
    df_tmp <- full_join(df_tmp, tmp)
  }
  return(df_tmp)
}

embryo_open <- function(img_name){
  img_path <- file.path(set_path, img_name)
  csv_list <- list.files(img_path, pattern = '_embryo') 
  for(csv in csv_list){
    csv_path <- file.path(img_path, csv) 
    name_tmp <- str_replace(csv, '_embryo_', '-')
    tmp <- read_csv(csv_path) %>% 
      mutate(image = str_split(name_tmp, '-')[[1]][1]) %>% 
      mutate(z_plane = str_replace(str_split(name_tmp, '-')[[1]][2], '.csv', '')) %>%
      mutate(z_plane = as.integer(z_plane)) %>% 
      rename(id = roi_id) %>% 
      mutate(id = as.integer(id), embryo_id = as.integer(embryo_id)) %>% 
      select(image, id, embryo_id, z_plane)
    df_tmp <- full_join(df_tmp, tmp)
  }
  return(df_tmp)
}

# Define paths-----------------------------------------------------------------

base_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Bechara Saykali/output/dhb_embryo_04'

output_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Bechara Saykali/results/dhb_embryo_04'

# Upload data------------------------------------------------------------------

set_list <- list.files(base_path)

nuc <- tibble(image = character(), set = character(), id = integer(), 
              z_plane = integer(), area = numeric(), ruby_int = numeric(),
              clover_int = numeric(), ch405_int = numeric(), 
              ch647_int = numeric())
df_tmp <- tibble(image = character(), id = integer(), z_plane = integer(),
                 area = numeric(), ruby_int = numeric(), clover_int = numeric(),
                 ch405_int = numeric(), ch647_int = numeric())

for(set_name in set_list){
  set_path <- file.path(base_path, set_name)
  img_list <- list.files(set_path) 
  for(img in img_list){
    nuc_tmp <- nuc_open(img) %>% 
      mutate(set = set_name)
    nuc <- full_join(nuc, nuc_tmp)
  }
}

nuc <- nuc %>% 
  arrange(set, image, z_plane)

any(is.na(nuc))
unique(nuc$set)

cyto <- tibble(image = character(), set = character(), id = integer(),
               z_plane = integer(), cyto_area = numeric(),
               cyto_ruby_int = numeric(), cyto_clover_int = numeric(),
               cyto_ch405_int = numeric(), cyto_ch647_int = numeric())
df_tmp <- tibble(image = character(), id = integer(), z_plane = integer(),
                 cyto_area = numeric(), cyto_ruby_int = numeric(),
                 cyto_clover_int = numeric(), cyto_ch405_int = numeric(),
                 cyto_ch647_int = numeric())

for(set_name in set_list){
  set_path <- file.path(base_path, set_name)
  img_list <- list.files(set_path)
  for(img in img_list){
    cyto_tmp <- cyto_open(img) %>% 
      mutate(set = set_name)
    cyto <- full_join(cyto, cyto_tmp)
  }
}

any(is.na(cyto))
unique(cyto$set)

embryo <- tibble(image = character(), set = character(), id = integer(), 
                 embryo_id = integer(), z_plane = integer())
df_tmp <- tibble(image = character(), id = integer(), 
                 embryo_id = integer(), z_plane = integer())

for(set_name in set_list){
  set_path <- file.path(base_path, set_name)
  img_list <- list.files(set_path)
  for(img in img_list){
    embryo_tmp <- embryo_open(img) %>% 
      mutate(set = set_name)
    embryo <- full_join(embryo, embryo_tmp)
  }
}

any(is.na(embryo))
unique(embryo$set)

# Data cleaning----------------------------------------------------------------

df_clean <- nuc %>% 
  group_by(set, image, id) %>% 
  summarise(z = n()) %>% 
  mutate(keep = ifelse(z > 1, 1, 0)) %>% 
  select(image, id, keep)

stage_csv_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Bechara Saykali/notes/dhb_embryo_03.csv'
stage_csv <- read_csv(stage_csv_path)

df <- full_join(nuc, cyto) %>% 
  full_join(embryo) %>% 
  full_join(df_clean) %>% 
  #full_join(stage_csv) %>% 
  filter(keep > 0) %>% 
  filter(!is.na(area)) %>% 
  filter(!is.na(cyto_area)) %>% 
  mutate(total_ruby_int = area * ruby_int) %>% 
  mutate(total_cyto_ruby_int = cyto_area * cyto_ruby_int) %>% 
  mutate(total_clover_int = area * clover_int) %>% 
  mutate(total_cyto_clover_int = cyto_area * cyto_clover_int) %>% 
  mutate(total_ch405_int = area * ch405_int) %>% 
  mutate(total_cyto_ch405_int = cyto_area * cyto_ch405_int) %>% 
  mutate(total_ch647_int = area * ch647_int) %>% 
  mutate(total_cyto_ch647_int = cyto_area * cyto_ch647_int) %>% 
  group_by(set, image, id, embryo_id) %>% 
  summarise(area = sum(area), cyto_area = sum(cyto_area),
            total_ruby_int = sum(total_ruby_int),
            total_cyto_ruby_int = sum(total_cyto_ruby_int),
            total_clover_int = sum(total_clover_int),
            total_cyto_clover_int = sum(total_cyto_clover_int),
            total_ch405_int = sum(total_ch405_int),
            total_cyto_ch405_int = sum(total_cyto_ch405_int),
            total_ch647_int = sum(total_ch647_int),
            total_cyto_ch647_int = sum(total_cyto_ch647_int)) %>% 
  mutate(ruby_int = total_ruby_int / area,
         cyto_ruby_int = total_cyto_ruby_int / cyto_area,
         clover_int = total_clover_int / area,
         cyto_clover_int = total_cyto_clover_int / cyto_area,
         ch405_int = total_ch405_int / area,
         cyto_ch405_int = total_cyto_ch405_int / cyto_area,
         ch647_int = total_ch647_int / area,
         cyto_ch647_int = total_cyto_ch647_int / cyto_area) %>% 
  mutate(clover_norm = clover_int / ruby_int) %>% 
  mutate(clover_ratio = clover_int / cyto_clover_int) %>% 
  mutate(ruby_ratio = ruby_int / cyto_ruby_int) %>% 
  mutate(ch405_ratio = ch405_int / cyto_ch405_int) %>% 
  mutate(ch647_ratio = ch647_int / cyto_ch647_int) %>% 
  select(image, set, id, embryo_id, area, cyto_area, ruby_int, clover_int, 
         clover_norm, clover_ratio, ruby_ratio, 
         ch405_ratio, ch647_ratio, ch405_int, ch647_int) %>% 
  dplyr::rowwise() %>% 
  mutate(stain = case_when(
    grepl('CDX', image) ~ 'CDX2_Sox2',
    grepl('Gata', image) ~ 'Gata6_Nanog',
    TRUE ~ 'NA'
  ))

any(is.na(df))
unique(df$set)
unique(df$stage)
unique(df$stain)

# df_test <- df %>% 
#   filter(is.na(embryo_id))

df_count <- df %>% 
  group_by(image, set, embryo_id) %>% 
  summarise(nuc_count = n())

# Determine nuclear staining cutoffs-------------------------------------------
# CDX2 / Sox2
df_cdx <- df %>% 
  filter(stain == 'CDX2_Sox2')

#cdx2_threshold <- quantile(df_cdx$ch405_int, probs = c(0.85))
#sox2_threshold <- quantile(df_cdx$ch647_int, probs = c(0.85))

cdx2_threshold <- c(500)
sox2_threshold <- c(300)

ggplot(df_cdx, aes(x = ch405_ratio, y = ch647_ratio)) +
  geom_point(size = 1, alpha = 0.75, show.legend = F) +
  geom_density2d() +
  geom_vline(xintercept = 1, linetype = 2, col = 'red') +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  theme +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks(base = 10, sides = 'bl') +
  labs(title = 'CDX2 vs Sox2 intensity', x = 'CDX2 intensity [AU]',
       y = 'Sox2 intensity [AU]')

ggplot(df_cdx, aes(x = ch405_int, y = ch647_int)) +
  geom_point(size = 1, alpha = 0.75, show.legend = F) +
  geom_density2d() +
  geom_vline(xintercept = cdx2_threshold, linetype = 2, col = 'red') +
  geom_hline(yintercept = sox2_threshold, linetype = 2, col = 'red') +
  theme +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks(base = 10, sides = 'bl') +
  labs(title = 'CDX2 vs Sox2 intensity', x = 'CDX2 intensity [AU]',
       y = 'Sox2 intensity [AU]')

plot_name <- 'cdx2_sox2_intensity.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 16,
       height = 11, dpi = 150)

ggplot(df_cdx, aes(x = ch405_int)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, 
                 fill = 'white', col = 'black') +
  geom_density(lwd = 2) +
  geom_vline(xintercept = cdx2_threshold, 
             linetype = 2, col = 'red') +
  annotate(geom = 'label', x = cdx2_threshold, y = 0.002, 
           label = signif(cdx2_threshold, digits = 6)) +
  theme +
  labs(title = 'CDX2 intensity', x = 'Intensity [AU]')

plot_name <- 'cdx2_intensity_histogram_500.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)

ggplot(df_cdx, aes(x = ch647_int)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, 
                 fill = 'white', col = 'black') +
  geom_density(lwd = 2) +
  geom_vline(xintercept = sox2_threshold, 
             linetype = 2, col = 'red') +
  annotate(geom = 'label', x = sox2_threshold, y = 0.002, 
           label = signif(sox2_threshold, digits = 6)) +
  theme +
  labs(title = 'Sox2 intensity', x = 'Intensity [AU]')

plot_name <- 'sox2_intensity_histogram.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)

df_cdx <- df_cdx %>% 
  mutate(ch405_class = ifelse(ch405_int >= cdx2_threshold, 
                              'cdx2_pos', 'cdx2_neg')) %>% 
  mutate(ch647_class = ifelse(ch647_int >= sox2_threshold,
                              'sox2_pos', 'sox2_neg')) %>% 
  mutate(marker_class = case_when(
    ch405_class == 'cdx2_neg' & ch647_class == 'sox2_neg' ~ 'dbl_neg',
    ch405_class == 'cdx2_pos' & ch647_class == 'sox2_neg' ~ 'cdx2_pos_sox2_neg',
    ch405_class == 'cdx2_pos' & ch647_class == 'sox2_pos' ~ 'dbl_pos',
    ch405_class == 'cdx2_neg' & ch647_class == 'sox2_pos' ~ 'cdx2_neg_sox2_pos',
    TRUE ~ 'NA'
  ))

unique(df_cdx$ch405_class)  
unique(df_cdx$ch647_class)
unique(df_cdx$marker_class)


# Merge dataframe

df_total <- df_cdx %>% 
  full_join(df_count)

any(is.na(df_total))
unique(df_total$stage)

df_total$stage <- factor(df_total$stage, 
                         levels = c('Morula', 'E3.5', 'E4.0', 'E4.5'))

# Plot Ruby nuclear intensity--------------------------------------------------

tmp_plot <- df_total

ggplot(tmp_plot, aes(x = stage, y = ruby_int, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 60, linetype = 2, col = 'red', lwd = 1) +
  theme +
  scale_y_log10() +
  annotation_logticks(base = 10, sides = 'l') +
  labs(title = 'mRuby nuclear intensity', x = 'Embryo stage',
       y = 'Mean intensity [AU]')

plot_name <- 'ruby_nuclear_intensity_total.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Plot Clover nuclear intensity------------------------------------------------

tmp_plot <- df_total

ggplot(tmp_plot, aes(x = stage, y = clover_int, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, col = 'black', show.legend = F,
               outlier.size = -1) +
  theme +
  scale_y_log10() +
  annotation_logticks(base = 10, sides = 'l') +
  labs(title = 'Clover nuclear intensity', x = 'Embryo stage',
       y = 'Mean intensity [AU]')

plot_name <- 'clover_nuclear_intensity_total.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Plot normalized Clover intensity---------------------------------------------

tmp_plot <- df_total %>% 
  filter(ruby_int > 60)

ggplot(tmp_plot, aes(x = stage, y = clover_norm, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, col = 'black', show.legend = F,
               outlier.size = -1) +
  #geom_hline(yintercept = 1, linetype = 2, col = 'red', lwd = 1) +
  theme +
  #ylim(0, 1.1) +
  labs(title = 'Clover normalized intensity', x = 'Embryo stage', 
       y = 'Normalized Intensity [clover/ruby]')

plot_name <- 'clover_norm_intensity_total.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Plot Clover ratio------------------------------------------------------------

tmp_plot <- df_total %>% 
  filter(ruby_int > 60) 

ggplot(tmp_plot, aes(x = stage, y = clover_ratio, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red', lwd = 1) +
  theme +
  #ylim(0, 1.1) +
  labs(title = 'Clover nuclear/cytoplasimc ratio', x = 'Embryo stage', 
       y = 'Intensity ratio [nuclear/cytoplasmic]')

plot_name <- 'clover_ratio_intensity_total.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Separate stainings-----------------------------------------------------------
# CDX2 / Sox2 

tmp_plot <- df_total %>% 
  filter(ruby_int > 60) %>% 
  filter(stain == 'CDX2_Sox2')

ggplot(tmp_plot, aes(x = stage, y = clover_norm, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  #geom_hline(yintercept = 1, linetype = 2, col = 'red', lwd = 1) +
  theme +
  facet_grid(rows = vars(ch405_class)) +
  labs(title = 'Clover normalized intensity', x = 'Embryo stage', 
       y = 'Normalized Intensity [clover/ruby]')

plot_name <- 'clover_norm_intensity_cdx2.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)  

ggplot(tmp_plot, aes(x = stage, y = clover_ratio, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red', lwd = 1) +
  theme +
  #ylim(0, 1.1) +
  facet_grid(rows = vars(ch405_class)) +
  labs(title = 'Clover nuclear/cytoplasimc ratio', x = 'Embryo stage', 
       y = 'Intensity ratio [nuclear/cytoplasmic]')

plot_name <- 'clover_ratio_intensity_cdx2.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

ggplot(tmp_plot, aes(x = stage, y = clover_norm, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  #geom_hline(yintercept = 1, linetype = 2, col = 'red', lwd = 1) +
  theme +
  facet_grid(rows = vars(ch647_class)) +
  labs(title = 'Clover normalized intensity', x = 'Embryo stage', 
       y = 'Normalized Intensity [clover/ruby]')

plot_name <- 'clover_norm_intensity_sox2.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)   

ggplot(tmp_plot, aes(x = stage, y = clover_ratio, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red', lwd = 1) +
  theme +
  #ylim(0, 1.1) +
  facet_grid(rows = vars(ch647_class)) +
  labs(title = 'Clover nuclear/cytoplasimc ratio', x = 'Embryo stage', 
       y = 'Intensity ratio [nuclear/cytoplasmic]')

plot_name <- 'clover_ratio_intensity_sox2.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Plot selected stages for staining--------------------------------------------
# CDX2 

tmp_plot <- df_total %>% 
  filter(ruby_int > 60) %>% 
  filter(stain == 'CDX2_Sox2') %>% 
  filter(stage == 'BlastE' | stage == 'E3.5late' | stage == 'E4'| stage == 'E4.5')

ggplot(tmp_plot, aes(x = ch405_class, y = clover_norm, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  #geom_hline(yintercept = 1, linetype = 2, col = 'red', lwd = 1) +
  theme +
  facet_grid(cols = vars(stage)) +
  labs(title = 'Clover normalized intensity', x = 'Embryo stage', 
       y = 'Normalized Intensity [clover/ruby]')

plot_name <- 'clover_norm_intensity_cdx2_stage.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)   

ggplot(tmp_plot, aes(x = ch405_class, y = clover_ratio, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col ='red', lwd = 1) +
  theme +
  facet_grid(cols = vars(stage)) +
  labs(title = 'Clover nuclear/cytoplasimc ratio', x = 'Embryo stage', 
       y = 'Intensity ratio [nuclear/cytoplasmic]')

plot_name <- 'clover_ratio_intensity_cdx2_stage.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)   

# Sox2 

tmp_plot <- df_total %>% 
  filter(ruby_int > 60) %>% 
  filter(stain == 'CDX2_Sox2') %>% 
  filter(stage == 'BlastE' | stage == 'E3.5late' | stage == 'E4' | stage == 'E4.5')

ggplot(tmp_plot, aes(x = ch647_class, y = clover_norm, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  #geom_hline(yintercept = 1, linetype = 2, col = 'red', lwd = 1) +
  theme +
  facet_grid(cols = vars(stage)) +
  labs(title = 'Clover normalized intensity', x = 'Embryo stage', 
       y = 'Normalized Intensity [clover/ruby]')

plot_name <- 'clover_norm_intensity_sox2_stage.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)   

ggplot(tmp_plot, aes(x = ch647_class, y = clover_ratio, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col ='red', lwd = 1) +
  theme +
  facet_grid(cols = vars(stage)) +
  labs(title = 'Clover nuclear/cytoplasimc ratio', x = 'Embryo stage', 
       y = 'Intensity ratio [nuclear/cytoplasmic]')

plot_name <- 'clover_ratio_intensity_sox2_stage.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150) 

# Narrow analysis--------------------------------------------------------------
# CDX2 / Sox2------------------------------------------------------------------

df2 <- df_total

unique(df2$image)

df2$marker_class <- factor(df2$marker_class, 
                           levels = c('dbl_neg', 'cdx2_pos_sox2_neg',
                                      'dbl_pos', 'cdx2_neg_sox2_pos'))

unique(df2$marker_class)

# CDX2-------------------------------------------------------------------------

tmp_plot <- df2 %>% 
  filter(ruby_int > 60)

ggplot(tmp_plot, aes(x = ch405_class, y = clover_norm, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  #geom_hline(yintercept = 1, linetype = 2, col = 'red', lwd = 1) +
  theme +
  facet_grid(cols = vars(stage)) +
  labs(title = 'Clover normalized intensity', x = 'Embryo stage', 
       y = 'Normalized Intensity [clover/ruby]')

plot_name <- 'narrow_clover_norm_intensity_cdx2_stage_500.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)   

ggplot(tmp_plot, aes(x = ch405_class, y = clover_ratio, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col ='red', lwd = 1) +
  theme +
  facet_grid(cols = vars(stage)) +
  labs(title = 'Clover nuclear/cytoplasmic ratio', x = 'Embryo stage', 
       y = 'Intensity ratio [nuclear/cytoplasmic]')

plot_name <- 'narrow_clover_ratio_intensity_cdx2_stage_500.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150) 

# Sox2------------------------------------------------------------------------- 

tmp_plot <- df2 %>% 
  filter(ruby_int > 60)

ggplot(tmp_plot, aes(x = ch647_class, y = clover_norm, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  #geom_hline(yintercept = 1, linetype = 2, col = 'red', lwd = 1) +
  theme +
  facet_grid(cols = vars(stage)) +
  labs(title = 'Clover normalized intensity', x = 'Embryo stage', 
       y = 'Normalized Intensity [clover/ruby]')

plot_name <- 'narrow_clover_norm_intensity_sox2_stage.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)   

ggplot(tmp_plot, aes(x = ch647_class, y = clover_ratio, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col ='red', lwd = 1) +
  theme +
  facet_grid(cols = vars(stage)) +
  labs(title = 'Clover nuclear/cytoplasimc ratio', x = 'Embryo stage', 
       y = 'Intensity ratio [nuclear/cytoplasmic]')

plot_name <- 'narrow_clover_ratio_intensity_sox2_stage.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150) 

# CDX2 marker class------------------------------------------------------------

ggplot(tmp_plot, aes(x = marker_class, y = clover_norm, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  #geom_hline(yintercept = 1, linetype = 2, col = 'red', lwd = 1) +
  theme +
  theme(axis.text = element_text(angle = 90)) +
  facet_grid(cols = vars(stage)) +
  labs(title = 'Clover normalized intensity', x = 'Embryo stage', 
       y = 'Normalized Intensity [clover/ruby]')

plot_name <- 'narrow_clover_norm_intensity_cdx2_marker_500.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)   

ggplot(tmp_plot, aes(x = marker_class, y = clover_ratio, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col ='red', lwd = 1) +
  theme +
  theme(axis.text = element_text(angle = 90)) +
  facet_grid(cols = vars(stage)) +
  labs(title = 'Clover nuclear/cytoplasmic ratio', x = 'Embryo stage', 
       y = 'Intensity ratio [nuclear/cytoplasmic]')

plot_name <- 'narrow_clover_ratio_intensity_cdx2_marker_500.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)   

# CDX2 marker class, 2 class---------------------------------------------------

tmp_plot2 <- tmp_plot %>% 
  filter(grepl('cdx2', marker_class))

ggplot(tmp_plot2, aes(x = marker_class, y = clover_norm, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  #geom_hline(yintercept = 1, linetype = 2, col = 'red', lwd = 1) +
  theme +
  theme(axis.text = element_text(angle = 90)) +
  facet_grid(cols = vars(stage)) +
  labs(title = 'Clover normalized intensity', x = 'Embryo stage', 
       y = 'Normalized Intensity [clover/ruby]')

plot_name <- 'narrow_clover_norm_intensity_cdx2_marker2_500.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)   

ggplot(tmp_plot2, aes(x = marker_class, y = clover_ratio, fill = stage)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col ='red', lwd = 1) +
  theme +
  theme(axis.text = element_text(angle = 90)) +
  facet_grid(cols = vars(stage)) +
  labs(title = 'Clover nuclear/cytoplasmic ratio', x = 'Embryo stage', 
       y = 'Intensity ratio [nuclear/cytoplasmic]')

plot_name <- 'narrow_clover_ratio_intensity_cdx2_marker2_500.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)   

# Stratify by cell count-------------------------------------------------------
# Ruby nuclear intensity-------------------------------------------------------

tmp_plot <- df_total %>% 
  filter(embryo_id > 0)

ggplot(tmp_plot, aes(x = nuc_count, y = ruby_int)) +
  geom_point(size = 2, alpha = 0.25) +
  theme +
  labs(title = 'mRuby nuclear intensity vs embryo nuclear count',
       y = 'Mean intensity [AU]', x = 'Embryo nuclear count')

plot_name <- 'embryo_count_ruby_intensity.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)  

# Clover nuclear intensity-----------------------------------------------------

tmp_plot <- df_total %>% 
  filter(embryo_id > 0)

ggplot(tmp_plot, aes(x = nuc_count, y = clover_int)) +
  geom_point(size = 2, alpha = 0.25) +
  theme +
  labs(title = 'Clover nuclear intensity vs embryo nuclear count',
       y = 'Mean intensity [AU]', x = 'Embryo nuclear count')

plot_name <- 'embryo_count_clover_intensity.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)  

# Clover ratio-----------------------------------------------------------------

tmp_plot <- df_total %>% 
  filter(embryo_id > 0)

ggplot(tmp_plot, aes(x = nuc_count, y = clover_ratio)) +
  geom_point(size = 2, alpha = 0.25) +
  theme +
  labs(title = 'Clover nuclear cytoplasmic ratio vs embryo nuclear count',
       y = 'Intensity ratio [nucleus/cytoplasmic]', 
       x = 'Embryo nuclear count')

plot_name <- 'embryo_count_clover_ratio.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)  

# Clover ratio CDX2 separated--------------------------------------------------

tmp_plot <- df_total %>% 
  filter(embryo_id > 0)

ggplot(tmp_plot, aes(x = nuc_count, y = clover_ratio)) +
  geom_point(size = 2, alpha = 0.25) +
  facet_grid(rows = vars(ch405_class)) +
  theme +
  labs(title = 'Clover nuclear cytoplasmic ratio vs embryo nuclear count',
       y = 'Intensity ratio [nucleus/cytoplasmic]', 
       x = 'Embryo nuclear count')

plot_name <- 'embryo_count_clover_ratio_cdx.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)  

# Clover ratio Sox2 separated--------------------------------------------------

tmp_plot <- df_total %>% 
  filter(embryo_id > 0)

ggplot(tmp_plot, aes(x = nuc_count, y = clover_ratio)) +
  geom_point(size = 2, alpha = 0.25) +
  facet_grid(rows = vars(ch647_class)) +
  theme +
  labs(title = 'Clover nuclear cytoplasmic ratio vs embryo nuclear count',
       y = 'Intensity ratio [nucleus/cytoplasmic]', 
       x = 'Embryo nuclear count')

plot_name <- 'embryo_count_clover_ratio_sox2.png'
ggsave(plot_name, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)  

# Separate for each individual embryo------------------------------------------
# Clover ratio-----------------------------------------------------------------

tmp_plot <- df_total %>% 
  filter(embryo_id > 0) %>% 
  mutate(unique_embryo_id = str_c(image, as.character(embryo_id), sep = '_'))

tmp_plot_summary <- tmp_plot %>% 
  group_by(unique_embryo_id) %>% 
  summarise

embryo_list <- tmp_plot_summary$unique_embryo_id

output_path2 <- file.path(output_path, 'individual_embryo')

for(embryo_name in embryo_list){
  print(embryo_name)
  tmp_plot2 <- tmp_plot %>% 
    filter(unique_embryo_id == embryo_name)
  
  ggplot(tmp_plot2, aes(x = ch405_class, y = clover_ratio)) +
    geom_jitter(width = 0.25, size = 1, alpha = 0.75, show.legend = F) +
    geom_violin(width = 0.5, alpha = 0.1, col = 'black', show.legend = F) +
    geom_boxplot(width = 0.1, alpha = 0.1, col = 'black', show.legend = F,
                 outlier.size = -1) +
    geom_hline(yintercept = 1, linetype = 2, col = 'red', linewidth = 1) +
    theme +
    labs(title = embryo_name,
         y = 'Intensity ratio [nucleus/cytoplasmic]', 
         x = '')
  
  tmp_plotname <- str_c(embryo_name, '_cdx2.png', sep = '')
  ggsave(tmp_plotname, plot = last_plot(), path = output_path2, width = 8,
         height = 11, dpi = 150)
  
  ggplot(tmp_plot2, aes(x = ch647_class, y = clover_ratio)) +
    geom_jitter(width = 0.25, size = 1, alpha = 0.75, show.legend = F) +
    geom_violin(width = 0.5, alpha = 0.1, col = 'black', show.legend = F) +
    geom_boxplot(width = 0.1, alpha = 0.1, col = 'black', show.legend = F,
                 outlier.size = -1) +
    geom_hline(yintercept = 1, linetype = 2, col = 'red', linewidth = 1) +
    theme +
    labs(title = embryo_name,
         y = 'Intensity ratio [nucleus/cytoplasmic]', 
         x = '')
  
  tmp_plotname <- str_c(embryo_name, '_sox2.png', sep = '')
  ggsave(tmp_plotname, plot = last_plot(), path = output_path2, width = 8,
         height = 11, dpi = 150)
}

# Export spreadsheet-----------------------------------------------------------

csv_name <- 'dhb_embryo_04_output.csv'
csv_path <- file.path(output_path, csv_name)

write_csv(df_total, csv_path)
