# FGF4 CDX2 embryo analysis
# 05/31/24

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
  for(csv_name in csv_list){
    print(csv_name)
    csv_path <- file.path(base_path, csv_name)
    tmp <- read_csv(csv_path) %>% 
      mutate(image = str_replace(csv_name, csv_pattern, '')) %>% 
      mutate(id = as.integer(id),
             embryo_id = as.integer(embryo_id),
             z_plane = as.integer(z_plane),
             area = as.integer(area),
             ch405_int = as.numeric(ch405_int),
             ch488_int = as.numeric(ch488_int),
             ch561_int = as.numeric(ch561_int)) %>% 
      select(image, id, embryo_id, z_plane, area, 
             ch405_int, ch488_int, ch561_int)
    df_tmp <- full_join(df_tmp, tmp)
  }
  return(df_tmp)
}

# Define paths-----------------------------------------------------------------

base_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Bechara Saykali/images/fgf4_cdx2_embryo/05_17'

output_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Bechara Saykali/results/fgf4_cdx2_embryo/05_17'

# Upload data------------------------------------------------------------------

df_tmp <- tibble(image = character(), id = integer(), embryo_id = integer(), 
                 z_plane = integer(), area = integer(),
                 ch405_int = numeric(), ch488_int = numeric(), 
                 ch561_int = numeric())

df_nuc <- csv_open('_output.csv')
df_nuc <- df_nuc %>% 
  rename(nuc_area = area,
         nuc_ch405_int = ch405_int,
         nuc_ch488_int = ch488_int,
         nuc_ch561_int = ch561_int)

df_cyto <- csv_open('_cyto.csv')
df_cyto <- df_cyto %>% 
  rename(cyto_area = area,
         cyto_ch405_int = ch405_int,
         cyto_ch488_int = ch488_int,
         cyto_ch561_int = ch561_int)

# Data cleaning----------------------------------------------------------------

df_tmp2 <- full_join(df_nuc, df_cyto) %>% 
  filter(embryo_id > 0) %>% 
  filter(!is.na(cyto_area))

any(is.na(df_tmp2))

df_clean <- df_tmp2 %>% 
  group_by(image, id) %>% 
  summarise(z = n()) %>% 
  mutate(keep = ifelse(z > 1, 1, 0)) %>% 
  select(image, id, keep)

df <- full_join(df_tmp2, df_clean) %>% 
  filter(keep > 0) %>% 
  mutate(total_nuc_ch405_int = nuc_area * nuc_ch405_int,
         total_cyto_ch405_int = cyto_area * cyto_ch405_int,
         total_nuc_ch488_int = nuc_area * nuc_ch488_int,
         total_cyto_ch488_int = cyto_area * cyto_ch488_int,
         total_nuc_ch561_int = nuc_area * nuc_ch561_int,
         total_cyto_ch561_int = cyto_area * cyto_ch561_int) %>% 
  group_by(image, id, embryo_id) %>% 
  summarise(nuc_area = sum(nuc_area), cyto_area = sum(cyto_area),
            total_nuc_ch405_int = sum(total_nuc_ch405_int),
            total_cyto_ch405_int = sum(total_cyto_ch405_int),
            total_nuc_ch488_int = sum(total_nuc_ch488_int),
            total_cyto_ch488_int = sum(total_cyto_ch488_int),
            total_nuc_ch561_int = sum(total_nuc_ch561_int),
            total_cyto_ch561_int = sum(total_cyto_ch561_int)) %>% 
  mutate(nuc_ch405_int = total_nuc_ch405_int / nuc_area,
         cyto_ch405_int = total_cyto_ch405_int / cyto_area,
         nuc_ch488_int = total_nuc_ch488_int / nuc_area,
         cyto_ch488_int = total_cyto_ch488_int / cyto_area,
         nuc_ch561_int = total_nuc_ch561_int / nuc_area,
         cyto_ch561_int = total_cyto_ch561_int / cyto_area) %>% 
  mutate(ch405_ratio = cyto_ch405_int / nuc_ch405_int,
         ch488_ratio = cyto_ch488_int / nuc_ch488_int,
         ch561_ratio = cyto_ch561_int / nuc_ch561_int) %>% 
  dplyr::rowwise() %>% 
  mutate(treatment = case_when(
    grepl('nofgf', image) ~ 'No FGF4',
    grepl('fgf', image) ~ 'FGF4'
  )) %>% 
  select(image, treatment, id, embryo_id, nuc_area, cyto_area,
         nuc_ch405_int, cyto_ch405_int, ch405_ratio,
         nuc_ch488_int, cyto_ch488_int, ch488_ratio,
         nuc_ch561_int, cyto_ch561_int, ch561_ratio)

df$treatment <- factor(df$treatment, levels = c('No FGF4', 'FGF4'))
unique(df$treatment)

any(is.na(df))

# Classify CDX2 status---------------------------------------------------------

cdx2_threshold <- c(105)

ggplot(df, aes(x = nuc_ch405_int)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, 
                 fill = 'white', col = 'black') +
  geom_density(linewidth = 2) +
  geom_vline(xintercept = cdx2_threshold, 
             linetype = 2, col = 'red') +
  annotate(geom = 'label', x = cdx2_threshold, y = 0.002, 
           label = signif(cdx2_threshold, digits = 6)) +
  theme +
  labs(title = 'CDX2 intensity', x = 'Intensity [AU]')

plotname <- 'cdx2_intensity_histogram.png'
ggsave(plotname, plot = last_plot(), path = output_path, 
       width = 11, height = 8, dpi = 150)

df <- df %>% 
  mutate(cdx2_status = ifelse(nuc_ch405_int > cdx2_threshold,
                              'cdx2_pos', 'cdx2_neg'))

df$cdx2_status <- factor(df$cdx2_status, levels = c('cdx2_pos', 'cdx2_neg'))

# Plot DHB ratio---------------------------------------------------------------

ggplot(df, aes(x = cdx2_status, y = ch488_ratio, fill = cdx2_status)) +
  geom_jitter(width = 0.25, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.5, show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.5, show.legend = F,
               outlier.size = -1) +
  geom_hline(yintercept = 1, col = 'red', linetype = 2, linewidth = 1) +
  theme +
  scale_y_log10() +
  annotation_logticks(base = 10, sides = 'l') +
  facet_grid(cols = vars(treatment)) +
  labs(title = 'Clover intensity ratio', x = 'CDX2 status', 
       y = 'Clover ratio [cytoplasmic/nuclear]')

plotname <- 'clover_ratio.png'
ggsave(plotname, plot = last_plot(), path = output_path, 
       width = 11, height = 8, dpi = 150)

# Export spreadsheet-----------------------------------------------------------

csvname <- 'fgf4_embryo_output.csv'
csvpath <- file.path(output_path, csvname)

write_csv(df, csvpath)
