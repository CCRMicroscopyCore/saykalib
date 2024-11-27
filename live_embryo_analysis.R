# Live embryo analysis for saykalib
# 11/19/24
# First live dataset includes 9_12 - 9_27

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
                plot.title=element_text(size=24), 
                axis.text=element_text(color="black"),
                plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

# Compile data-----------------------------------------------------------------
### Merge cytoplasmic and nuclear data-----------------------------------------

img_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Bechara Saykali/images/live_embryo/cyto_nuc_output/240913_liDHB_live_day0_mtg_M4'

cyto_list <- list.files(img_path, pattern = 'cyto_output')

df_tmp <- tibble(nuc_id = integer(),
                 area = integer(),
                 intensity = numeric(),
                 time_point = integer(), 
                 z_plane = integer())

iter <- 1

for(cyto_name in cyto_list){
  print(paste(as.character(iter), as.character(length(cyto_list)), sep = '/'))
  cyto_path <- file.path(img_path, cyto_name)
  cyto_tmp <- read_csv(cyto_path) %>% 
    mutate(nuc_id = as.integer(nuc_id),
           area = as.integer(area),
           intensity = as.numeric(intensity),
           time_point = as.integer(time_point),
           z_plane = as.integer(z_plane)) %>% 
    select(nuc_id, area, intensity, time_point, z_plane)
  df_tmp <- full_join(df_tmp, cyto_tmp)
  iter <- iter + 1
}

cyto_df <- df_tmp %>% 
  mutate(total_intensity = area * intensity) %>% 
  group_by(nuc_id, time_point) %>% 
  summarise(total_area = sum(area),
            total_intensity = sum(total_intensity)) %>% 
  mutate(cyto_intensity = total_intensity / total_area) %>%
  select(nuc_id, time_point, cyto_intensity)

nuc_list <- list.files(img_path, pattern = 'nuc_output')

df_tmp <- tibble(nuc_id = integer(),
                 area = integer(),
                 intensity = numeric(),
                 time_point = integer(), 
                 z_plane = integer())

iter <- 1

for(nuc_name in nuc_list){
  print(paste(as.character(iter), as.character(length(nuc_list)), sep = '/'))
  nuc_path <- file.path(img_path, nuc_name)
  nuc_tmp <- read_csv(nuc_path) %>% 
    mutate(nuc_id = as.integer(nuc_id),
           area = as.integer(area),
           intensity = as.numeric(intensity),
           time_point = as.integer(time_point),
           z_plane = as.integer(z_plane)) %>% 
    select(nuc_id, area, intensity, time_point, z_plane)
  df_tmp <- full_join(df_tmp, nuc_tmp)
  iter <- iter + 1
}

nuc_df <- df_tmp %>% 
  mutate(total_intensity = area * intensity) %>% 
  group_by(nuc_id, time_point) %>% 
  summarise(total_area = sum(area),
            total_intensity = sum(total_intensity)) %>% 
  mutate(nuc_intensity = total_intensity / total_area) %>%
  select(nuc_id, time_point, nuc_intensity)

df <- full_join(nuc_df, cyto_df) %>% 
  mutate(cyto_nuc_ratio = cyto_intensity / nuc_intensity) %>% 
  filter(!is.na(cyto_nuc_ratio))

### Download track and classification data-------------------------------------

img_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Bechara Saykali/images/live_embryo/output2/240913_liDHB_live_day0_mtg_M5_Statistics'

csv_list <- list.files(img_path, pattern = 'Mean_Ch=3_Img=1')
csv_name <- csv_list[1]
csv_path <- file.path(img_path, csv_name)

csv_df <- read_csv(csv_path, skip = 3) %>% 
  rename(time_point = contains('Time'),
         nuc_id = contains('Intensity'),
         track_id = contains('TrackID'),
         class = contains('Set')) %>% 
  mutate(time_point = as.integer(time_point),
         nuc_id = as.integer(nuc_id),
         track_id = as.integer(track_id)) %>% 
  mutate(time_point = time_point - 1) %>% 
  replace_na(list(track_id = 0, 
                  class = 'cdx2')) %>% 
  select(time_point, nuc_id, track_id, class)

### Merge all data-------------------------------------------------------------

df_tmp <- full_join(csv_df, df) %>% 
  filter(!is.na(track_id)) %>% 
  filter(track_id > 0)

df_track <- df_tmp %>% 
  group_by(track_id) %>% 
  summarise(initial_time_point = min(time_point),
            last_time_point = max(time_point)) %>% 
  mutate(track_length = last_time_point - initial_time_point) %>% 
  select(track_id, track_length)

df2 <- full_join(df_tmp, df_track) %>% 
  filter(track_length > 20) %>% 
  filter(!is.na(cyto_nuc_ratio))

any(is.na(df2))

df2$class <- factor(df2$class, levels = c('sox2',
                                          'cdx2'))

### Data analysis--------------------------------------------------------------

output_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Bechara Saykali/results/live_embryo/240913_liDHB_live_day0_mtg_M5'

output_df_path <- file.path(output_path, 'dhb_intensity.csv')
write_csv(df2, output_df_path)

### Plot data------------------------------------------------------------------

plot_df <- df2 %>% 
  group_by(time_point, track_id, class) %>% 
  summarise(cyto_nuc_ratio = mean(cyto_nuc_ratio))

ggplot(plot_df, aes(x = time_point, y = cyto_nuc_ratio, col = track_id)) +
  geom_line(aes(group = track_id),
            alpha = 0.75,
            show.legend = F) +
  geom_hline(yintercept = 1, col = 'red', linetype = 2) +
  facet_grid(rows = vars(class)) +
  theme +
  labs(title = 'DHB C/N ratio',
       y = 'Cytoplasmic / nuclear ratio',
       x = 'Time point')

plotname <- 'dhb_ratio.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 11, dpi = 150)

# # Log2 data--------------------------------------------------------------------
# 
# csv_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Bechara Saykali/results/live_embryo/240927_iDHB_day0_mtg_2024_09_27__16_14_33(1)/dhb_intensity.csv'
# output_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Bechara Saykali/results/live_embryo/240927_iDHB_day0_mtg_2024_09_27__16_14_33(1)'

# df <- read_csv(csv_path)
# 
# df$class <- factor(df$class, levels = c('sox2', 'cdx2'))

df_out <- df2 %>% 
  mutate(log2_ratio = log2(cyto_nuc_ratio)) %>% 
  group_by(time_point, track_id, class) %>% 
  summarise(log2_ratio = mean(log2_ratio))

output_df_path <- file.path(output_path, 'dhb_log2_ratio.csv')
write_csv(df_out, output_df_path)

ggplot(df_out, aes(x = time_point, y = log2_ratio, col = track_id)) +
  geom_line(aes(group = track_id),
            alpha = 0.75,
            show.legend = F) +
  geom_hline(yintercept = 0, col = 'red', linetype = 2) +
  facet_grid(rows = vars(class)) +
  theme +
  labs(title = 'DHB Log2 C/N ratio',
       y = 'Log2 Cytoplasmic / nuclear ratio',
       x = 'Time point')

plotname <- 'dhb_log2_ratio.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 11, dpi = 150)

# Plot branching---------------------------------------------------------------

plot_df <- df2 %>% 
  group_by(time_point, track_id, class) %>% 
  summarise(branches = n())

ggplot(plot_df, aes(x = time_point, y = branches, col = track_id)) +
  geom_line(aes(group = track_id),
            alpha = 0.75,
            show.legend = F) +
  facet_grid(rows = vars(class)) +
  theme +
  labs(title = 'Track branch number',
       y = '# of branches',
       x = 'Time point')

plotname <- 'track_branch.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 11, dpi = 150)

# Temp Zeiss files-------------------------------------------------------------

csv_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Bechara Saykali/images/live_embryo/output2/240912_iDHB_live_day0_mtg_2024_09_12__15_39_18(4)_Statistics/merge_0_0_0_Intensity_Mean_Ch=3_Img=1.csv'

output_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Bechara Saykali/results/live_embryo/240912_iDHB_live_day0_mtg_2024_09_12__15_39_18(4)'

df <- read_csv(csv_path, skip = 3) %>% 
  rename(time_point = contains('Time'),
         nuc_id = contains('Intensity'),
         track_id = contains('TrackID'),
         class = contains('Set')) %>% 
  mutate(time_point = as.integer(time_point),
         nuc_id = as.integer(nuc_id),
         track_id = as.integer(track_id)) %>% 
  mutate(time_point = time_point - 1) %>% 
  replace_na(list(track_id = 0, 
                  class = 'cdx2')) %>% 
  select(time_point, nuc_id, track_id, class)

df$class <- factor(df$class, levels = c('sox2', 'cdx2'))

plot_df <- df %>% 
  filter(track_id > 0) %>% 
  group_by(time_point, track_id, class) %>% 
  summarise(branches = n())

ggplot(plot_df, aes(x = time_point, y = branches, col = track_id)) +
  geom_line(aes(group = track_id),
            alpha = 0.75,
            show.legend = F) +
  facet_grid(rows = vars(class)) +
  theme +
  labs(title = 'Track branch number',
       y = '# of branches',
       x = 'Time point')

plotname <- 'track_branch.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16,
       height = 11, dpi = 150)
