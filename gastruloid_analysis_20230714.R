# Gastruloid analysis for Bechara Saykali
# 07/14/2023 set
# 07/26/23

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

csv_open <- function(csv_pattern){
  csv_list <- list.files(base_path, pattern = csv_pattern)
  for(csv in csv_list){
    csv_path <- file.path(base_path, csv)
    name_tmp <- str_replace(csv, csv_pattern, '')
    tmp <- read_csv(csv_path) %>% 
      mutate(image = name_tmp) %>% 
      mutate(roi_id = as.integer(roi_id), area = as.integer(area), 
             gast_id = as.integer(gast_id), dist_mean = as.numeric(dist_mean),
             dist_max = as.integer(dist_max), centroid_y = as.numeric(centroid_y),
             centroid_x = as.numeric(centroid_x),
             intensity_405 = as.numeric(intensity_405),
             intensity_488 = as.numeric(intensity_488),
             intensity_568 = as.numeric(intensity_568),
             intensity_647 = as.numeric(intensity_647)) %>% 
      select(image, roi_id, gast_id, area, dist_mean, dist_max, 
             centroid_y, centroid_x,
             intensity_405, intensity_488, intensity_568, intensity_647)
    df_tmp <- full_join(df_tmp, tmp)
  }
  return(df_tmp)
}

# Define paths-----------------------------------------------------------------

base_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Bechara Saykali/images/gastruloids/2023/07_14'
output_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Bechara Saykali/results/gastruloids/2023/07_14'

# Upload data------------------------------------------------------------------

df_tmp <- tibble(image = character(), roi_id = integer(), gast_id = integer(),
                 area = integer(), dist_mean = numeric(), dist_max = integer(),
                 centroid_y = numeric(), centroid_x = numeric(),
                 intensity_405 = numeric(), intensity_488 = numeric(),
                 intensity_568 = numeric(), intensity_647 = numeric())

gast_df_tmp <- csv_open('_gastruloid.csv')
gast_df <- gast_df_tmp %>% 
  rename(gast_area = area) %>%
  rename(gast_dist_max = dist_max) %>% 
  rename(gast_centroid_y = centroid_y, gast_centroid_x = centroid_x) %>% 
  select(image, gast_id, gast_area, gast_dist_max, 
         gast_centroid_y, gast_centroid_x)

nuc_df_tmp <- csv_open('_nuclei.csv')
nuc_df <- nuc_df_tmp %>% 
  rename(nuc_area = area) %>% 
  rename(nuc_intensity_405 = intensity_405) %>% 
  rename(nuc_intensity_488 = intensity_488) %>% 
  rename(nuc_intensity_568 = intensity_568) %>% 
  rename(nuc_intensity_647 = intensity_647) %>% 
  select(image, roi_id, gast_id, nuc_area, dist_mean, centroid_y, centroid_x,
         nuc_intensity_405, nuc_intensity_488, 
         nuc_intensity_568, nuc_intensity_647)

cyto_df_tmp <- csv_open('_cyto.csv')
cyto_df <- cyto_df_tmp %>% 
  rename(cyto_area = area) %>% 
  rename(cyto_intensity_405 = intensity_405) %>% 
  rename(cyto_intensity_488 = intensity_488) %>% 
  rename(cyto_intensity_568 = intensity_568) %>% 
  rename(cyto_intensity_647 = intensity_647) %>% 
  select(image, roi_id, cyto_area,
         cyto_intensity_405, cyto_intensity_488,
         cyto_intensity_568, cyto_intensity_647)

# Data cleaning----------------------------------------------------------------

px_per_um <- 0.324567986740054

df <- full_join(nuc_df, cyto_df) %>% 
  full_join(gast_df) %>% 
  filter(gast_id > 0) %>% 
  mutate(dist_from_center_px = gast_dist_max - dist_mean) %>% 
  mutate(dist_from_center_um = dist_from_center_px * px_per_um) %>% 
  mutate(corr_centroid_y = centroid_y - gast_centroid_y,
         corr_centroid_x = centroid_x - gast_centroid_x) %>% 
  mutate(dhb_ratio = nuc_intensity_488 / cyto_intensity_488) %>% 
  mutate(ring = case_when(
    dist_from_center_um < 100 ~ 'A',
    dist_from_center_um >= 100 & dist_from_center_um < 200 ~ 'B',
    dist_from_center_um >= 200 & dist_from_center_um < 300 ~ 'C',
    dist_from_center_um >= 300 & dist_from_center_um < 400 ~ 'D',
    dist_from_center_um >= 400 ~ 'E'
  )) %>% 
  mutate(radial_distance = dist_from_center_px / gast_dist_max) %>% 
  mutate(experiment = ifelse(grepl('plv', image), 'Expt2', 'Expt1')) %>% 
  mutate(bmp = ifelse(grepl('noBMP4', image), 'noBMP4', 'BMP4'))

any(is.na(df))
unique(df$image)
max(df$dist_from_center_um)
unique(df$ring)

df$ring <- factor(df$ring, levels = c('A', 'B', 'C', 'D', 'E'))

unique(df$bmp)
df$bmp <- factor(df$bmp, levels = c('noBMP4', 'BMP4'))

unique(df$experiment)

# Experiment 1-----------------------------------------------------------------
# CDX2 / Ncad------------------------------------------------------------------
# Split dataframe--------------------------------------------------------------

df_exp1 <- df %>% 
  filter(experiment == 'Expt1') %>% 
  filter(grepl('CDX2', image))

# Plot CDX2 and Ncad intensity and cutoff--------------------------------------

tmp_plot <- df_exp1 %>% 
  select(roi_id, nuc_intensity_405, nuc_intensity_647)

ncad_threshold = 115

ggplot(tmp_plot, aes(x = nuc_intensity_647)) +
  geom_histogram(aes(y = after_stat(density)), col = 'black', fill = 'white') +
  geom_density() +
  geom_vline(xintercept = ncad_threshold, linetype = 2, col = 'red') +
  annotate(geom = 'label', x = ncad_threshold, y = 20, label = ncad_threshold) +
  scale_x_log10() +
  annotation_logticks(base = 10, sides = 'b') +
  theme +
  labs(title = 'Ncad nuclear intensity', x = 'Intensity [AU]')

plotname <- 'expt1_ncad_intensity.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)

cdx_threshold = 160

ggplot(tmp_plot, aes(x = nuc_intensity_405)) +
  geom_histogram(aes(y = after_stat(density)), col = 'black', fill = 'white') +
  geom_density() +
  geom_vline(xintercept = cdx_threshold, linetype = 2, col = 'red') +
  annotate(geom = 'label', x = cdx_threshold, y = 10, label = cdx_threshold) +
  scale_x_log10() +
  annotation_logticks(base = 10, sides = 'b') +
  theme +
  labs(title = 'CDX2 nuclear intensity', x = 'Intensity [AU]')

plotname <- 'expt1_cdx2_intensity.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)

ggplot(tmp_plot, aes(x = nuc_intensity_647, y = nuc_intensity_405)) +
  geom_point(size = 1, alpha = 0.08) +
  geom_density2d() +
  geom_vline(xintercept = ncad_threshold, linetype = 2, col = 'red') +
  annotate(geom = 'label', x = ncad_threshold, y = 250, label = ncad_threshold) +
  geom_hline(yintercept = cdx_threshold, linetype = 2, col = 'red') +
  annotate(geom = 'label', y = cdx_threshold, x = 250, label = cdx_threshold) +
  scale_x_log10() +
  scale_y_log10() +
  theme +
  labs(title = 'Ncad vs CDX2 nuclear intensity', x = 'Ncad intensity [AU]',
       y = 'CDX2 intensity [AU]')

plotname <- 'expt1_ncad_cdx2_intensity.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 8,
       height = 11, dpi = 150)

df_exp1_2 <- df_exp1 %>% 
  mutate(germ_layer = ifelse(nuc_intensity_405 > cdx_threshold, 'trophoectoderm', 
                             'non-trophoectoderm'))

df_exp1_2$germ_layer <- factor(df_exp1_2$germ_layer, levels = c('non-trophoectoderm',
                                                    'trophoectoderm'))

df_exp1_2 <- df_exp1_2 %>% 
  mutate(cdx_status = ifelse(nuc_intensity_405 > cdx_threshold, 'CDX2_pos',
                             'CDX2_neg')) %>% 
  mutate(ncad_status = ifelse(nuc_intensity_647 > ncad_threshold, 'Ncad_pos',
                              'Ncad_neg'))

df_exp1_2$cdx_status <- factor(df_exp1_2$cdx_status, levels = c('CDX2_neg',
                                                                'CDX2_pos'))

df_exp1_2$ncad_status <- factor(df_exp1_2$ncad_status, levels = c('Ncad_neg',
                                                                'Ncad_pos'))

# Plot DHB ratio in total cells------------------------------------------------

ggplot(df_exp1_2, aes(x = germ_layer, y = dhb_ratio, col = germ_layer)) +
  geom_jitter(width = 0.35, size = 1, alpha = 0.1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() +
  facet_grid(cols = vars(bmp)) +
  ylim(0, 3) +
  theme +
  labs(title = 'DHB ratio', x = 'Germ layer', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt1_ncad_cdx2_dhb_ratio_total.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Plot DHB ratio divided by ring-----------------------------------------------

ggplot(df_exp1_2, aes(x = ring, y = dhb_ratio, col = ring)) +
  geom_jitter(width = 0.35, size = 1, alpha = 0.1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() + 
  facet_grid(cols = vars(bmp)) +
  ylim(0, 3) +
  theme + 
  labs(title = 'DHB ratio vs concentric ring', 
       x = 'Ring', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt1_ncad_cdx2_dhb_ratio_ring.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Plot DHB ratio vs radial distance--------------------------------------------

ggplot(df_exp1_2, aes(x = radial_distance, y = dhb_ratio, col = ring)) +
  geom_point(size = 1, alpha = 0.1, show.legend = F) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() +
  facet_grid(rows = vars(bmp)) +
  ylim(0, 3) +
  theme +
  labs(title = 'DHB ratio vs radial distance',
       x = 'Radial distance', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt1_ncad_cdx2_dhb_ratio_radial.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Plot DHB ratio by Ncad status------------------------------------------------

ggplot(df_exp1_2, aes(x = ncad_status, y = dhb_ratio, col = ncad_status)) +
  geom_jitter(width = 0.35, size = 1, alpha = 0.1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() +
  facet_grid(cols = vars(bmp)) +
  ylim(0, 3) +
  theme +
  labs(title = 'DHB ratio', x = '', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt1_ncad_cdx2_dhb_ratio_ncad.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Plot DHB ratio by CDX2 status------------------------------------------------

ggplot(df_exp1_2, aes(x = cdx_status, y = dhb_ratio, col = cdx_status)) +
  geom_jitter(width = 0.35, size = 1, alpha = 0.1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() +
  facet_grid(cols = vars(bmp)) +
  ylim(0, 3) +
  theme +
  labs(title = 'DHB ratio', x = '', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt1_ncad_cdx2_dhb_ratio_cdx2.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Pseudo-map DHB ratio---------------------------------------------------------
rng <- seq(from = 0, to = 1.75, by = 0.25)

ggplot(df_exp1_2, aes(x = corr_centroid_x, y = corr_centroid_y, fill = dhb_ratio)) +
  geom_point(size = 1, pch = 21, col = 'white', alpha = 0.25) +
  scale_fill_gradient2(low = 'white', mid = 'lightgray', high = 'black',
                       midpoint = mean(rng),
                       limits = c(min(rng), max(rng)), name = '') +
  ylim(-2000, 2000) +
  xlim(-2000, 2000) +
  facet_grid(cols = vars(bmp)) +
  theme +
  labs(title = 'DHB ratio pseudo map', x = '', y = '')

plotname <- 'expt1_ncad_cdx2_dhb_ratio_pseudomap.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

max(df_exp1_2$nuc_intensity_405)
rng <- seq(from = 25, to = 500, by = 25)

ggplot(df_exp1_2, aes(x = corr_centroid_x, y = corr_centroid_y, 
                fill = nuc_intensity_405)) +
  geom_point(size = 1, pch = 21, col = 'white', alpha = 0.25) +
  scale_fill_gradient2(low = 'white', mid = 'cyan', high = 'royalblue4',
                       midpoint = mean(rng),
                       limits = c(min(rng), max(rng)), name = '') +
  ylim(-2000, 2000) +
  xlim(-2000, 2000) +
  facet_grid(cols = vars(bmp)) +
  theme +
  labs(title = 'CDX2 intensity pseudo map', x = '', y = '')

plotname <- 'expt1_ncad_cdx2_cdx2_intensity_pseudomap.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

max(df_exp1_2$nuc_intensity_647)
rng <- seq(from = 25, to = 150, by = 25)

ggplot(df_exp1_2, aes(x = corr_centroid_x, y = corr_centroid_y, 
                fill = nuc_intensity_647)) +
  geom_point(size = 1, pch = 21, col = 'white', alpha = 0.25) +
  scale_fill_gradient2(low = 'white', mid = 'pink', high = 'darkred',
                       midpoint = mean(rng),
                       limits = c(min(rng), max(rng)), name = '') +
  ylim(-2000, 2000) +
  xlim(-2000, 2000) +
  facet_grid(cols = vars(bmp)) +
  theme +
  labs(title = 'Ncad intensity pseudo map', x = '', y = '')

plotname <- 'expt1_ncad_cdx2_ncad_intensity_pseudomap.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Sox2 / Tfap2c----------------------------------------------------------------
# Split dataframe--------------------------------------------------------------
unique(df$image)
df_exp1 <- df %>% 
  filter(experiment == 'Expt1') %>% 
  filter(grepl('SOX2', image))

# Plot Sox2 and Tfap2c intensity and cutoff--------------------------------------

tmp_plot <- df_exp1 %>% 
  select(roi_id, nuc_intensity_405, nuc_intensity_647)

tfap_threshold = 140

ggplot(tmp_plot, aes(x = nuc_intensity_647)) +
  geom_histogram(aes(y = after_stat(density)), col = 'black', fill = 'white') +
  geom_density() +
  geom_vline(xintercept = tfap_threshold, linetype = 2, col = 'red') +
  annotate(geom = 'label', x = tfap_threshold, y = 20, label = tfap_threshold) +
  scale_x_log10() +
  annotation_logticks(base = 10, sides = 'b') +
  theme +
  labs(title = 'Tfap2c nuclear intensity', x = 'Intensity [AU]')

plotname <- 'expt1_tfap2c_intensity.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)

sox2_threshold = 160

ggplot(tmp_plot, aes(x = nuc_intensity_405)) +
  geom_histogram(aes(y = after_stat(density)), col = 'black', fill = 'white') +
  geom_density() +
  geom_vline(xintercept = sox2_threshold, linetype = 2, col = 'red') +
  annotate(geom = 'label', x = sox2_threshold, y = 10, label = sox2_threshold) +
  scale_x_log10() +
  annotation_logticks(base = 10, sides = 'b') +
  theme +
  labs(title = 'Sox2 nuclear intensity', x = 'Intensity [AU]')

plotname <- 'expt1_sox2_intensity.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)

ggplot(tmp_plot, aes(x = nuc_intensity_647, y = nuc_intensity_405)) +
  geom_point(size = 1, alpha = 0.08) +
  geom_density2d() +
  geom_vline(xintercept = tfap_threshold, linetype = 2, col = 'red') +
  annotate(geom = 'label', x = tfap_threshold, y = 250, label = tfap_threshold) +
  geom_hline(yintercept = sox2_threshold, linetype = 2, col = 'red') +
  annotate(geom = 'label', y = sox2_threshold, x = 250, label = sox2_threshold) +
  scale_x_log10() +
  scale_y_log10() +
  theme +
  labs(title = 'Tfap2c vs Sox2 nuclear intensity', x = 'Tfap2c intensity [AU]',
       y = 'Sox2 intensity [AU]')

plotname <- 'expt1_tfap_sox2_intensity.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 8,
       height = 11, dpi = 150)

df_exp1_2 <- df_exp1 %>% 
  mutate(germ_layer = ifelse(nuc_intensity_405 > sox2_threshold, 'trophoectoderm', 
                             'non-trophoectoderm'))

df_exp1_2$germ_layer <- factor(df_exp1_2$germ_layer, levels = c('non-trophoectoderm',
                                                                'trophoectoderm'))

df_exp1_2 <- df_exp1_2 %>% 
  mutate(sox2_status = ifelse(nuc_intensity_405 > sox2_threshold, 'SOX2_pos',
                             'SOX2_neg')) %>% 
  mutate(tfap_status = ifelse(nuc_intensity_647 > tfap_threshold, 'TFAP_pos',
                              'TFAP_neg'))

df_exp1_2$sox2_status <- factor(df_exp1_2$sox2_status, levels = c('SOX2_neg',
                                                                'SOX2_pos'))

df_exp1_2$tfap_status <- factor(df_exp1_2$tfap_status, levels = c('TFAP_neg',
                                                                  'TFAP_pos'))

# Plot DHB ratio in total cells------------------------------------------------

ggplot(df_exp1_2, aes(x = germ_layer, y = dhb_ratio, col = germ_layer)) +
  geom_jitter(width = 0.35, size = 1, alpha = 0.1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() +
  facet_grid(cols = vars(bmp)) +
  ylim(0, 3) +
  theme +
  labs(title = 'DHB ratio', x = 'Germ layer', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt1_tfap_sox2_dhb_ratio_total.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Plot DHB ratio divided by ring-----------------------------------------------

ggplot(df_exp1_2, aes(x = ring, y = dhb_ratio, col = ring)) +
  geom_jitter(width = 0.35, size = 1, alpha = 0.1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() + 
  facet_grid(cols = vars(bmp)) +
  ylim(0, 3) +
  theme + 
  labs(title = 'DHB ratio vs concentric ring', 
       x = 'Ring', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt1_tfap_sox2_dhb_ratio_ring.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Plot DHB ratio vs radial distance--------------------------------------------

ggplot(df_exp1_2, aes(x = radial_distance, y = dhb_ratio, col = ring)) +
  geom_point(size = 1, alpha = 0.1, show.legend = F) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() +
  facet_grid(rows = vars(bmp)) +
  ylim(0, 3) +
  theme +
  labs(title = 'DHB ratio vs radial distance',
       x = 'Radial distance', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt1_tfap_sox2_dhb_ratio_radial.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Plot DHB ratio by TFAP status------------------------------------------------

ggplot(df_exp1_2, aes(x = tfap_status, y = dhb_ratio, col = tfap_status)) +
  geom_jitter(width = 0.35, size = 1, alpha = 0.1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() +
  facet_grid(cols = vars(bmp)) +
  ylim(0, 3) +
  theme +
  labs(title = 'DHB ratio', x = '', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt1_tfap_sox2_dhb_ratio_tfap.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Plot DHB ratio by SOX2 status------------------------------------------------

ggplot(df_exp1_2, aes(x = sox2_status, y = dhb_ratio, col = sox2_status)) +
  geom_jitter(width = 0.35, size = 1, alpha = 0.1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() +
  facet_grid(cols = vars(bmp)) +
  ylim(0, 3) +
  theme +
  labs(title = 'DHB ratio', x = '', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt1_tfap_sox2_dhb_ratio_sox2.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Pseudo-map DHB ratio---------------------------------------------------------
rng <- seq(from = 0, to = 1.75, by = 0.25)

ggplot(df_exp1_2, aes(x = corr_centroid_x, y = corr_centroid_y, fill = dhb_ratio)) +
  geom_point(size = 1, pch = 21, col = 'white', alpha = 0.25) +
  scale_fill_gradient2(low = 'white', mid = 'lightgray', high = 'black',
                       midpoint = mean(rng),
                       limits = c(min(rng), max(rng)), name = '') +
  ylim(-2000, 2000) +
  xlim(-2000, 2000) +
  facet_grid(cols = vars(bmp)) +
  theme +
  labs(title = 'DHB ratio pseudo map', x = '', y = '')

plotname <- 'expt1_tfap_sox2_dhb_ratio_pseudomap.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

max(df_exp1_2$nuc_intensity_405)
rng <- seq(from = 25, to = 225, by = 25)

ggplot(df_exp1_2, aes(x = corr_centroid_x, y = corr_centroid_y, 
                      fill = nuc_intensity_405)) +
  geom_point(size = 1, pch = 21, col = 'white', alpha = 0.25) +
  scale_fill_gradient2(low = 'white', mid = 'cyan', high = 'royalblue4',
                       midpoint = mean(rng),
                       limits = c(min(rng), max(rng)), name = '') +
  ylim(-2000, 2000) +
  xlim(-2000, 2000) +
  facet_grid(cols = vars(bmp)) +
  theme +
  labs(title = 'Sox2 intensity pseudo map', x = '', y = '')

plotname <- 'expt1_tfap_sox2_sox2_intensity_pseudomap.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

max(df_exp1_2$nuc_intensity_647)
rng <- seq(from = 25, to = 300, by = 25)

ggplot(df_exp1_2, aes(x = corr_centroid_x, y = corr_centroid_y, 
                      fill = nuc_intensity_647)) +
  geom_point(size = 1, pch = 21, col = 'white', alpha = 0.25) +
  scale_fill_gradient2(low = 'white', mid = 'pink', high = 'darkred',
                       midpoint = mean(rng),
                       limits = c(min(rng), max(rng)), name = '') +
  ylim(-2000, 2000) +
  xlim(-2000, 2000) +
  facet_grid(cols = vars(bmp)) +
  theme +
  labs(title = 'Tfap2c intensity pseudo map', x = '', y = '')

plotname <- 'expt1_tfap_sox2_tfap_intensity_pseudomap.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 8, dpi = 150)

# Experiment 2-----------------------------------------------------------------
# Split dataframe--------------------------------------------------------------

df_exp2 <- df %>% 
  filter(experiment == 'Expt2') %>% 
  mutate(cdk = case_when(
    grepl('CDK24C8', image) ~ 'CDK24C8',
    grepl('CDK24C1', image) ~ 'CDK24C1',
    grepl('CDK2', image) ~ 'CDK2',
    TRUE ~ 'control'
  ))

df_exp2$cdk <- factor(df_exp2$cdk, levels = c('control', 'CDK2', 
                                              'CDK24C1', 'CDK24C8'))

# Plot CDX2 and Ncad intensity and cutoff--------------------------------------

tmp_plot <- df_exp2 %>% 
  select(roi_id, nuc_intensity_405, nuc_intensity_647)

ncad_threshold = 115

ggplot(tmp_plot, aes(x = nuc_intensity_647)) +
  geom_histogram(aes(y = after_stat(density)), col = 'black', fill = 'white') +
  geom_density() +
  geom_vline(xintercept = ncad_threshold, linetype = 2, col = 'red') +
  annotate(geom = 'label', x = ncad_threshold, y = 20, label = ncad_threshold) +
  scale_x_log10() +
  annotation_logticks(base = 10, sides = 'b') +
  theme +
  labs(title = 'Ncad nuclear intensity', x = 'Intensity [AU]')

plotname <- 'expt2_ncad_intensity.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)

cdx_threshold = 160

ggplot(tmp_plot, aes(x = nuc_intensity_405)) +
  geom_histogram(aes(y = after_stat(density)), col = 'black', fill = 'white') +
  geom_density() +
  geom_vline(xintercept = cdx_threshold, linetype = 2, col = 'red') +
  annotate(geom = 'label', x = cdx_threshold, y = 10, label = cdx_threshold) +
  scale_x_log10() +
  annotation_logticks(base = 10, sides = 'b') +
  theme +
  labs(title = 'CDX2 nuclear intensity', x = 'Intensity [AU]')

plotname <- 'expt2_cdx2_intensity.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 11,
       height = 8, dpi = 150)

ggplot(tmp_plot, aes(x = nuc_intensity_647, y = nuc_intensity_405)) +
  geom_point(size = 1, alpha = 0.08) +
  geom_density2d() +
  geom_vline(xintercept = ncad_threshold, linetype = 2, col = 'red') +
  annotate(geom = 'label', x = ncad_threshold, y = 250, label = ncad_threshold) +
  geom_hline(yintercept = cdx_threshold, linetype = 2, col = 'red') +
  annotate(geom = 'label', y = cdx_threshold, x = 250, label = cdx_threshold) +
  scale_x_log10() +
  scale_y_log10() +
  theme +
  labs(title = 'Ncad vs CDX2 nuclear intensity', x = 'Ncad intensity [AU]',
       y = 'CDX2 intensity [AU]')

plotname <- 'expt2_ncad_cdx2_intensity.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 8,
       height = 11, dpi = 150)

df_exp2_2 <- df_exp2 %>% 
  mutate(germ_layer = ifelse(nuc_intensity_405 > cdx_threshold, 'trophoectoderm', 
                             'non-trophoectoderm'))

df_exp2_2$germ_layer <- factor(df_exp2_2$germ_layer, levels = c('non-trophoectoderm',
                                                                'trophoectoderm'))

df_exp2_2 <- df_exp2_2 %>% 
  mutate(cdx_status = ifelse(nuc_intensity_405 > cdx_threshold, 'CDX2_pos',
                             'CDX2_neg')) %>% 
  mutate(ncad_status = ifelse(nuc_intensity_647 > ncad_threshold, 'Ncad_pos',
                              'Ncad_neg'))

df_exp2_2$cdx_status <- factor(df_exp2_2$cdx_status, levels = c('CDX2_neg',
                                                                'CDX2_pos'))

df_exp2_2$ncad_status <- factor(df_exp2_2$ncad_status, levels = c('Ncad_neg',
                                                                  'Ncad_pos'))

# Plot DHB ratio in total cells------------------------------------------------

ggplot(df_exp2_2, aes(x = germ_layer, y = dhb_ratio, col = germ_layer)) +
  geom_jitter(width = 0.35, size = 1, alpha = 0.1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() +
  facet_grid(cols = vars(bmp), rows = vars(cdk)) +
  ylim(0, 3) +
  theme +
  labs(title = 'DHB ratio', x = 'Germ layer', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt2_ncad_cdx2_dhb_ratio_total.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 11,
       height = 16, dpi = 150)

# Plot DHB ratio divided by ring-----------------------------------------------

ggplot(df_exp2_2, aes(x = ring, y = dhb_ratio, col = ring)) +
  geom_jitter(width = 0.35, size = 1, alpha = 0.1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() + 
  facet_grid(cols = vars(bmp), rows = vars(cdk)) +
  ylim(0, 3) +
  theme + 
  labs(title = 'DHB ratio vs concentric ring', 
       x = 'Ring', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt2_ncad_cdx2_dhb_ratio_ring.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 11, dpi = 150)

# Plot DHB ratio vs radial distance--------------------------------------------

ggplot(df_exp2_2, aes(x = radial_distance, y = dhb_ratio, col = ring)) +
  geom_point(size = 1, alpha = 0.1, show.legend = F) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() +
  facet_grid(cols = vars(bmp), rows = vars(cdk)) +
  ylim(0, 3) +
  theme +
  labs(title = 'DHB ratio vs radial distance',
       x = 'Radial distance', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt2_ncad_cdx2_dhb_ratio_radial.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 11, dpi = 150)

# Plot DHB ratio by Ncad status------------------------------------------------

ggplot(df_exp2_2, aes(x = ncad_status, y = dhb_ratio, col = ncad_status)) +
  geom_jitter(width = 0.35, size = 1, alpha = 0.1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() +
  facet_grid(cols = vars(bmp), rows = vars(cdk)) +
  ylim(0, 3) +
  theme +
  labs(title = 'DHB ratio', x = '', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt2_ncad_cdx2_dhb_ratio_ncad.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 11, dpi = 150)

# Plot DHB ratio by CDX2 status------------------------------------------------

ggplot(df_exp2_2, aes(x = cdx_status, y = dhb_ratio, col = cdx_status)) +
  geom_jitter(width = 0.35, size = 1, alpha = 0.1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = 'black', show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = 'black', show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, linetype = 2, col = 'red') +
  scale_color_aaas() +
  facet_grid(cols = vars(bmp), rows = vars(cdk)) +
  ylim(0, 3) +
  theme +
  labs(title = 'DHB ratio', x = '', y = 'DHB ratio [nuc/cyto]')

plotname <- 'expt2_ncad_cdx2_dhb_ratio_cdx2.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 11, dpi = 150)

# Pseudo-map DHB ratio---------------------------------------------------------
max(df_exp2_2$dhb_ratio)
rng <- seq(from = 0, to = 2, by = 0.25)

ggplot(df_exp2_2, aes(x = corr_centroid_x, y = corr_centroid_y, fill = dhb_ratio)) +
  geom_point(size = 1, pch = 21, col = 'white', alpha = 0.25) +
  scale_fill_gradient2(low = 'white', mid = 'lightgray', high = 'black',
                       midpoint = mean(rng),
                       limits = c(min(rng), max(rng)), name = '') +
  ylim(-2000, 2000) +
  xlim(-2000, 2000) +
  facet_grid(cols = vars(bmp), rows = vars(cdk)) +
  theme +
  labs(title = 'DHB ratio pseudo map', x = '', y = '')

plotname <- 'expt2_ncad_cdx2_dhb_ratio_pseudomap.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 11,
       height = 16, dpi = 150)

max(df_exp1_2$nuc_intensity_405)
rng <- seq(from = 25, to = 200, by = 25)

ggplot(df_exp2_2, aes(x = corr_centroid_x, y = corr_centroid_y, 
                      fill = nuc_intensity_405)) +
  geom_point(size = 1, pch = 21, col = 'white', alpha = 0.25) +
  scale_fill_gradient2(low = 'white', mid = 'cyan', high = 'royalblue4',
                       midpoint = mean(rng),
                       limits = c(min(rng), max(rng)), name = '') +
  ylim(-2000, 2000) +
  xlim(-2000, 2000) +
  facet_grid(cols = vars(bmp), rows = vars(cdk)) +
  theme +
  labs(title = 'CDX2 intensity pseudo map', x = '', y = '')

plotname <- 'expt2_ncad_cdx2_cdx2_intensity_pseudomap.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 11,
       height = 16, dpi = 150)

max(df_exp2_2$nuc_intensity_647)
rng <- seq(from = 25, to = 175, by = 25)

ggplot(df_exp2_2, aes(x = corr_centroid_x, y = corr_centroid_y, 
                      fill = nuc_intensity_647)) +
  geom_point(size = 1, pch = 21, col = 'white', alpha = 0.25) +
  scale_fill_gradient2(low = 'white', mid = 'pink', high = 'darkred',
                       midpoint = mean(rng),
                       limits = c(min(rng), max(rng)), name = '') +
  ylim(-2000, 2000) +
  xlim(-2000, 2000) +
  facet_grid(cols = vars(bmp), rows = vars(cdk)) +
  theme +
  labs(title = 'Ncad intensity pseudo map', x = '', y = '')

plotname <- 'expt2_ncad_cdx2_ncad_intensity_pseudomap.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 11,
       height = 16, dpi = 150)

# Export spreadsheet-----------------------------------------------------------
# Experiment 1 - CDX2----------------------------------------------------------

df_exp1 <- df %>% 
  filter(experiment == 'Expt1') %>% 
  filter(grepl('CDX2', image))

df_exp1_2 <- df_exp1 %>% 
  mutate(germ_layer = ifelse(nuc_intensity_405 > cdx_threshold, 'trophoectoderm', 
                             'non-trophoectoderm'))

df_exp1_2$germ_layer <- factor(df_exp1_2$germ_layer, levels = c('non-trophoectoderm',
                                                                'trophoectoderm'))

df_exp1_2 <- df_exp1_2 %>% 
  mutate(cdx_status = ifelse(nuc_intensity_405 > cdx_threshold, 'CDX2_pos',
                             'CDX2_neg')) %>% 
  mutate(ncad_status = ifelse(nuc_intensity_647 > ncad_threshold, 'Ncad_pos',
                              'Ncad_neg'))

df_exp1_2$cdx_status <- factor(df_exp1_2$cdx_status, levels = c('CDX2_neg',
                                                                'CDX2_pos'))

df_exp1_2$ncad_status <- factor(df_exp1_2$ncad_status, levels = c('Ncad_neg',
                                                                  'Ncad_pos'))


csvname <- 'expt1_cdx2_ncad_output.csv'
csv_path <- file.path(output_path, csvname)

write_csv(df_exp1_2, csv_path)

# Experiment 1 - SOX2----------------------------------------------------------

df_exp1 <- df %>% 
  filter(experiment == 'Expt1') %>% 
  filter(grepl('SOX2', image))

df_exp1_2 <- df_exp1 %>% 
  mutate(germ_layer = ifelse(nuc_intensity_405 > sox2_threshold, 'trophoectoderm', 
                             'non-trophoectoderm'))

df_exp1_2$germ_layer <- factor(df_exp1_2$germ_layer, levels = c('non-trophoectoderm',
                                                                'trophoectoderm'))

df_exp1_2 <- df_exp1_2 %>% 
  mutate(sox2_status = ifelse(nuc_intensity_405 > sox2_threshold, 'SOX2_pos',
                              'SOX2_neg')) %>% 
  mutate(tfap_status = ifelse(nuc_intensity_647 > tfap_threshold, 'TFAP_pos',
                              'TFAP_neg'))

df_exp1_2$sox2_status <- factor(df_exp1_2$sox2_status, levels = c('SOX2_neg',
                                                                  'SOX2_pos'))

df_exp1_2$tfap_status <- factor(df_exp1_2$tfap_status, levels = c('TFAP_neg',
                                                                  'TFAP_pos'))


csvname <- 'expt1_sox2_tfap_output.csv'
csv_path <- file.path(output_path, csvname)

write_csv(df_exp1_2, csv_path)

# Experiment 2-----------------------------------------------------------------

df_exp2 <- df %>% 
  filter(experiment == 'Expt2') %>% 
  mutate(cdk = case_when(
    grepl('CDK24C8', image) ~ 'CDK24C8',
    grepl('CDK24C1', image) ~ 'CDK24C1',
    grepl('CDK2', image) ~ 'CDK2',
    TRUE ~ 'control'
  ))

df_exp2$cdk <- factor(df_exp2$cdk, levels = c('control', 'CDK2', 
                                              'CDK24C1', 'CDK24C8'))

df_exp2_2 <- df_exp2 %>% 
  mutate(germ_layer = ifelse(nuc_intensity_405 > cdx_threshold, 'trophoectoderm', 
                             'non-trophoectoderm'))

df_exp2_2$germ_layer <- factor(df_exp2_2$germ_layer, levels = c('non-trophoectoderm',
                                                                'trophoectoderm'))

df_exp2_2 <- df_exp2_2 %>% 
  mutate(cdx_status = ifelse(nuc_intensity_405 > cdx_threshold, 'CDX2_pos',
                             'CDX2_neg')) %>% 
  mutate(ncad_status = ifelse(nuc_intensity_647 > ncad_threshold, 'Ncad_pos',
                              'Ncad_neg'))

df_exp2_2$cdx_status <- factor(df_exp2_2$cdx_status, levels = c('CDX2_neg',
                                                                'CDX2_pos'))

df_exp2_2$ncad_status <- factor(df_exp2_2$ncad_status, levels = c('Ncad_neg',
                                                                  'Ncad_pos'))

csvname <- 'expt2_cdx2_ncad_output.csv'
csv_path <- file.path(output_path, csvname)

write_csv(df_exp2_2, csv_path)
