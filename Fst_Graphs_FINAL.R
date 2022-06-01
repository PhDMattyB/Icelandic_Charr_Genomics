##############################
## Fst Outlier graphs combined
##
## Matt Brachmann (PhDMattyB)
##
## 2020-05-01
##
##############################

library(patchwork)
library(tidyverse)

AC04q.1_29_split = function(data){
  AC04q.1_29 = data %>% 
    filter(AC_CHR == 'AC04q.1:29') %>% 
    mutate(AC_CHR = as.factor(case_when(
      win_mid > '40000000' ~ 'AC04q.1',
      win_mid < '40000000' ~ 'AC29')))
  
  data = data %>% 
    filter(AC_CHR != 'AC04q.1:29')
  
  data = bind_rows(data, 
                   AC04q.1_29)
}

Mb_Conversion = function(data){
  data %>% 
    group_by(AC_CHR) %>% 
    mutate(win_mid_mb = win_mid/1000000)
}


setwd('~/PhD/SNP Demographic modelling/Outliers_directory/Fst_RDA_Colocalization')

# Both_outlier = read_csv('TSBPL_BOTH_Colocalization_Fst_Outliers_200Kb.csv') %>%
Both_outlier = read_csv('TLGBPL_BOTH_Colocalization_Fst_Outliers_200Kb.csv') %>% 
  filter(FST_n > 3) %>% 
  mutate(AC_CHR = as.factor(case_when(
    CHR == '1' ~ 'AC01',
    CHR == '2' ~ 'AC02',
    CHR == '3' ~ 'AC03',
    CHR == '4' ~ 'AC04p',
    CHR == '5' ~ 'AC04q.1:29',
    CHR == '6' ~ 'AC04q.2',
    CHR == '7' ~ 'AC05',
    CHR == '8' ~ 'AC06',
    CHR == '9' ~ 'AC06',
    CHR == '10' ~ 'AC07',
    CHR == '11' ~ 'AC08',
    CHR == '12' ~ 'AC09',
    CHR == '13' ~ 'AC10',
    CHR == '14' ~ 'AC11',
    CHR == '15' ~ 'AC12',
    CHR == '16' ~ 'AC13',
    CHR == '17' ~ 'AC14',
    CHR == '18' ~ 'AC15',
    CHR == '19' ~ 'AC16',
    CHR == '20' ~ 'AC17',
    CHR == '21' ~ 'AC18',
    CHR == '22' ~ 'AC19',
    CHR == '23' ~ 'AC20',
    CHR == '24' ~ 'AC21',
    CHR == '25' ~ 'AC22',
    CHR == '26' ~ 'AC23',
    CHR == '27' ~ 'AC24',
    CHR == '28' ~ 'AC25',
    CHR == '29' ~ 'AC26',
    CHR == '30' ~ 'AC27',
    CHR == '31' ~ 'AC28',
    CHR == '32' ~ 'AC30',
    CHR == '33' ~ 'AC31',
    CHR == '34' ~ 'AC32',
    CHR == '35' ~ 'AC33',
    CHR == '36' ~ 'AC34',
    CHR == '37' ~ 'AC35',
    CHR == '38' ~ 'AC36',
    CHR == '39' ~ 'AC37',
    CHR == '0' ~ 'Contigs')))

iso_outlier = read_csv('TLGBPL_ISO_Colocalization_Fst_Outliers_200Kb.csv') %>% 
  filter(FST_n > 3) %>%
  mutate(AC_CHR = as.factor(case_when(
    CHR == '1' ~ 'AC01',
    CHR == '2' ~ 'AC02',
    CHR == '3' ~ 'AC03',
    CHR == '4' ~ 'AC04p',
    CHR == '5' ~ 'AC04q.1:29',
    CHR == '6' ~ 'AC04q.2',
    CHR == '7' ~ 'AC05',
    CHR == '8' ~ 'AC06',
    CHR == '9' ~ 'AC06',
    CHR == '10' ~ 'AC07',
    CHR == '11' ~ 'AC08',
    CHR == '12' ~ 'AC09',
    CHR == '13' ~ 'AC10',
    CHR == '14' ~ 'AC11',
    CHR == '15' ~ 'AC12',
    CHR == '16' ~ 'AC13',
    CHR == '17' ~ 'AC14',
    CHR == '18' ~ 'AC15',
    CHR == '19' ~ 'AC16',
    CHR == '20' ~ 'AC17',
    CHR == '21' ~ 'AC18',
    CHR == '22' ~ 'AC19',
    CHR == '23' ~ 'AC20',
    CHR == '24' ~ 'AC21',
    CHR == '25' ~ 'AC22',
    CHR == '26' ~ 'AC23',
    CHR == '27' ~ 'AC24',
    CHR == '28' ~ 'AC25',
    CHR == '29' ~ 'AC26',
    CHR == '30' ~ 'AC27',
    CHR == '31' ~ 'AC28',
    CHR == '32' ~ 'AC30',
    CHR == '33' ~ 'AC31',
    CHR == '34' ~ 'AC32',
    CHR == '35' ~ 'AC33',
    CHR == '36' ~ 'AC34',
    CHR == '37' ~ 'AC35',
    CHR == '38' ~ 'AC36',
    CHR == '39' ~ 'AC37',
    CHR == '0' ~ 'Contigs')))

morpho_outlier = read_csv('TLGBPL_MORPHO_Colocalization_Fst_Outliers_200Kb.csv') %>% 
  filter(FST_n > 3) %>%
  mutate(AC_CHR = as.factor(case_when(
    CHR == '1' ~ 'AC01',
    CHR == '2' ~ 'AC02',
    CHR == '3' ~ 'AC03',
    CHR == '4' ~ 'AC04p',
    CHR == '5' ~ 'AC04q.1:29',
    CHR == '6' ~ 'AC04q.2',
    CHR == '7' ~ 'AC05',
    CHR == '8' ~ 'AC06',
    CHR == '9' ~ 'AC06',
    CHR == '10' ~ 'AC07',
    CHR == '11' ~ 'AC08',
    CHR == '12' ~ 'AC09',
    CHR == '13' ~ 'AC10',
    CHR == '14' ~ 'AC11',
    CHR == '15' ~ 'AC12',
    CHR == '16' ~ 'AC13',
    CHR == '17' ~ 'AC14',
    CHR == '18' ~ 'AC15',
    CHR == '19' ~ 'AC16',
    CHR == '20' ~ 'AC17',
    CHR == '21' ~ 'AC18',
    CHR == '22' ~ 'AC19',
    CHR == '23' ~ 'AC20',
    CHR == '24' ~ 'AC21',
    CHR == '25' ~ 'AC22',
    CHR == '26' ~ 'AC23',
    CHR == '27' ~ 'AC24',
    CHR == '28' ~ 'AC25',
    CHR == '29' ~ 'AC26',
    CHR == '30' ~ 'AC27',
    CHR == '31' ~ 'AC28',
    CHR == '32' ~ 'AC30',
    CHR == '33' ~ 'AC31',
    CHR == '34' ~ 'AC32',
    CHR == '35' ~ 'AC33',
    CHR == '36' ~ 'AC34',
    CHR == '37' ~ 'AC35',
    CHR == '38' ~ 'AC36',
    CHR == '39' ~ 'AC37',
    CHR == '0' ~ 'Contigs')))

normal = read_tsv('TLGBPL_Fst_neutral_200Kb_window.txt') %>%
  filter(FST_n > 3) %>% 
  mutate(AC_CHR = as.factor(case_when(
    CHR == '1' ~ 'AC01',
    CHR == '2' ~ 'AC02',
    CHR == '3' ~ 'AC03',
    CHR == '4' ~ 'AC04p',
    CHR == '5' ~ 'AC04q.1:29',
    CHR == '6' ~ 'AC04q.2',
    CHR == '7' ~ 'AC05',
    CHR == '8' ~ 'AC06',
    CHR == '9' ~ 'AC06',
    CHR == '10' ~ 'AC07',
    CHR == '11' ~ 'AC08',
    CHR == '12' ~ 'AC09',
    CHR == '13' ~ 'AC10',
    CHR == '14' ~ 'AC11',
    CHR == '15' ~ 'AC12',
    CHR == '16' ~ 'AC13',
    CHR == '17' ~ 'AC14',
    CHR == '18' ~ 'AC15',
    CHR == '19' ~ 'AC16',
    CHR == '20' ~ 'AC17',
    CHR == '21' ~ 'AC18',
    CHR == '22' ~ 'AC19',
    CHR == '23' ~ 'AC20',
    CHR == '24' ~ 'AC21',
    CHR == '25' ~ 'AC22',
    CHR == '26' ~ 'AC23',
    CHR == '27' ~ 'AC24',
    CHR == '28' ~ 'AC25',
    CHR == '29' ~ 'AC26',
    CHR == '30' ~ 'AC27',
    CHR == '31' ~ 'AC28',
    CHR == '32' ~ 'AC30',
    CHR == '33' ~ 'AC31',
    CHR == '34' ~ 'AC32',
    CHR == '35' ~ 'AC33',
    CHR == '36' ~ 'AC34',
    CHR == '37' ~ 'AC35',
    CHR == '38' ~ 'AC36',
    CHR == '39' ~ 'AC37',
    CHR == '0' ~ 'Contigs')))

# fst_outlier = read_tsv('SLGPEL_Fst_outliers_200Kb_window.txt') %>%
  fst_outlier = read_tsv('TLGBPL_Fst_outliers_200Kb_window.txt') %>% 
filter(FST_n > 3) %>% 
  mutate(AC_CHR = as.factor(case_when(
    CHR == '1' ~ 'AC01',
    CHR == '2' ~ 'AC02',
    CHR == '3' ~ 'AC03',
    CHR == '4' ~ 'AC04p',
    CHR == '5' ~ 'AC04q.1:29',
    CHR == '6' ~ 'AC04q.2',
    CHR == '7' ~ 'AC05',
    CHR == '8' ~ 'AC06',
    CHR == '9' ~ 'AC06',
    CHR == '10' ~ 'AC07',
    CHR == '11' ~ 'AC08',
    CHR == '12' ~ 'AC09',
    CHR == '13' ~ 'AC10',
    CHR == '14' ~ 'AC11',
    CHR == '15' ~ 'AC12',
    CHR == '16' ~ 'AC13',
    CHR == '17' ~ 'AC14',
    CHR == '18' ~ 'AC15',
    CHR == '19' ~ 'AC16',
    CHR == '20' ~ 'AC17',
    CHR == '21' ~ 'AC18',
    CHR == '22' ~ 'AC19',
    CHR == '23' ~ 'AC20',
    CHR == '24' ~ 'AC21',
    CHR == '25' ~ 'AC22',
    CHR == '26' ~ 'AC23',
    CHR == '27' ~ 'AC24',
    CHR == '28' ~ 'AC25',
    CHR == '29' ~ 'AC26',
    CHR == '30' ~ 'AC27',
    CHR == '31' ~ 'AC28',
    CHR == '32' ~ 'AC30',
    CHR == '33' ~ 'AC31',
    CHR == '34' ~ 'AC32',
    CHR == '35' ~ 'AC33',
    CHR == '36' ~ 'AC34',
    CHR == '37' ~ 'AC35',
    CHR == '38' ~ 'AC36',
    CHR == '39' ~ 'AC37',
    CHR == '0' ~ 'Contigs')))

# Split AC04q.1:29 --------------------------------------------------------

Both_outlier = AC04q.1_29_split(data = Both_outlier)
iso_outlier = AC04q.1_29_split(data = iso_outlier)
morpho_outlier = AC04q.1_29_split(data = morpho_outlier)
normal = AC04q.1_29_split(data = normal)
fst_outlier = AC04q.1_29_split(data = fst_outlier)

# Conversion --------------------------------------------------------------

Both_outlier = Mb_Conversion(data = Both_outlier)
iso_outlier = Mb_Conversion(data = iso_outlier)
morpho_outlier = Mb_Conversion(data = morpho_outlier)
normal = Mb_Conversion(normal)
fst_outlier = Mb_Conversion(data = fst_outlier)

##
# Fst outlier zoom 1-5 ----------------------------------------------------
# test = normal %>% 
#   filter(AC_CHR != 'Contigs', 
#          AC_CHR != 'AC01', 
#          AC_CHR != 'AC02', 
#          AC_CHR != 'AC03', 
#          AC_CHR != 'AC04p', 
#          AC_CHR != 'AC04q.2')
# View(test)
normal_sub = normal %>% 
  filter(AC_CHR %in% c('AC01', 
                       'AC02', 
                       'AC03', 
                       'AC04p', 
                       'AC04q.2', 
                       'AC05'))
Fst_sub = fst_outlier %>% 
  filter(AC_CHR %in% c('AC01', 
                       'AC02', 
                       'AC03', 
                       'AC04p', 
                       'AC04q.2', 
                       'AC05'))
overlap_sub = overlap %>% 
  filter(AC_CHR %in% c('AC01', 
                       'AC02', 
                       'AC03', 
                       'AC04p', 
                       'AC04q.2', 
                       'AC05'))
morpho_sub = morpho_outlier %>% 
  filter(AC_CHR %in% c('AC01', 
                       'AC02', 
                       'AC03', 
                       'AC04p', 
                       'AC04q.2', 
                       'AC05'))
iso_sub = iso_outlier %>% 
  filter(AC_CHR %in% c('AC01', 
                       'AC02', 
                       'AC03', 
                       'AC04p', 
                       'AC04q.2', 
                       'AC05'))
Both_sub = Both_outlier %>% 
  filter(AC_CHR %in% c('AC01', 
                       'AC02', 
                       'AC03', 
                       'AC04p', 
                       'AC04q.2', 
                       'AC05'))

## Make the plot
VBRSIL_Fst_Zoom =
  ggplot(data = normal_sub, 
         aes(x = win_mid_mb,
             y = FST_mean, 
             group = AC_CHR)) +
  geom_point(col = '#8C8C8C', 
             size = 2)+
  geom_point(data = Fst_sub,
             aes(x = win_mid_mb, 
                 y = FST_mean, 
                 group = AC_CHR),
             col = '#5488E3', 
             size = 2)+
  # geom_point(data = overlap_sub,
  #            aes(x = win_mid_mb,
  #                y = FST_mean,
  #                group =AC_CHR),
  #            col = '#3D1673',
  #            size = 2)+
  # geom_point(data = morpho_sub,
  #            aes(x = win_mid_mb,
  #                y = FST_mean,
  #                group = AC_CHR),
  #            col = '#274F73',
  #            size = 2)+
  # geom_point(data = iso_sub,
  #            aes(x = win_mid_mb,
  #                y = FST_mean,
  #                group = AC_CHR),
  #            col = '#58F252',
  #            size = 2)+
  # geom_point(data = Both_sub,
  #            aes(x = win_mid_mb,
  #                y = FST_mean,
  #                group = AC_CHR),
  #            col = '#F28705',
  #            size = 2)+
  facet_grid(~AC_CHR, 
             scales = 'free')+
  labs(x = 'Chromosomal position (Mb)', 
       y = 'Fst',
       title = 'Vatnshlidarvatn')+
  ylim(0.00, 1.00)+
  theme(panel.grid = element_blank(), 
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 12,
                                   angle = 90,
                                   hjust = 1),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        # axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face = 'bold',
                                  size = 14),
        strip.background = element_rect(fill = 'white',
                                        colour = 'black'), 
        plot.title = element_text(size = 18,
                                  face = 'bold'))

## Need to combine the plots for all morph pairs
SLGBPEL_Fst_Zoom
SLGBPEL_Fst_Zoom
TLGBPL_Fst_Zoom
TSBPL_Fst_Zoom
VBRSIL_Fst_Zoom

zoom = SLGBPEL_Fst_Zoom/VBRSIL_Fst_Zoom
zoom2 = SLGBPEL_Fst_Zoom/SLGBPEL_Fst_Zoom/TLGBPL_Fst_Zoom/TSBPL_Fst_Zoom/VBRSIL_Fst_Zoom

ggsave('Manhattan_plot_zoom_10.12.2020.tiff',
       plot = zoom, 
       units = 'cm', 
       dpi = 'retina')

ggsave('Manhattan_plot_zoom_allpops_10.12.2020.tiff',
       plot = zoom2, 
       units = 'cm', 
       dpi = 'retina', 
       height = 25, 
       width = 30)

##
# Fst Outlier Zoom --------------------------------------------------------
## Need to subset the datasets for AC24 and AC25
normal_sub = normal %>% 
  filter(AC_CHR %in% c('AC24', 'AC25'))
Fst_sub = fst_outlier %>% 
  filter(AC_CHR %in% c('AC24', 'AC25'))
overlap_sub = overlap %>% 
  filter(AC_CHR %in% c('AC24', 'AC25'))

## Make the plot
VBRSIL_Fst_Zoom =
  ggplot(data = normal_sub, 
         aes(x = win_mid_mb,
             y = FST_mean, 
             group = AC_CHR)) +
  geom_point(col = '#8C8C8C', 
             size = 2)+
  geom_point(data = Fst_sub,
             aes(x = win_mid_mb, 
                 y = FST_mean, 
                 group = AC_CHR),
             col = '#5488E3', 
             size = 2)+
  geom_point(data = overlap_sub, 
             aes(x = win_mid_mb, 
                 y = FST_mean, 
                 group =AC_CHR), 
             col = '#3D1673',
             size = 2)+
  facet_grid(~AC_CHR, 
             scales = 'free')+
  labs(x = 'Chromosomal position (Mb)', 
       y = 'Fst',
       title = 'Vatnshlidarvatn')+
  ylim(0.00, 1.00)+
  theme(panel.grid = element_blank(), 
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 12,
                                   angle = 90,
                                   hjust = 1),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        # axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face = 'bold',
                                  size = 14),
        strip.background = element_rect(fill = 'white',
                                        colour = 'black'), 
        plot.title = element_text(size = 18,
                                  face = 'bold'))

## Need to combine the plots for all morph pairs
SLGBPEL_Fst_Zoom
SLGBPEL_Fst_Zoom
TLGBPL_Fst_Zoom
TSBPL_Fst_Zoom
VBRSIL_Fst_Zoom

zoom_combo_big = SLGBPEL_Fst_Zoom/SLGBPEL_Fst_Zoom/TLGBPL_Fst_Zoom/TSBPL_Fst_Zoom/VBRSIL_Fst_Zoom

zoom_combo_small = SLGBPEL_Fst_Zoom/VBRSIL_Fst_Zoom

ggsave('Common_Fst_outliers_GV_26.11.2020.tiff',
       plot = zoom_combo_small, 
       units = 'cm', 
       dpi = 'retina')

##
# BIG Graph -------------------------------------------------------------------

## set the theme for all of the plots
theme_set(theme_bw())

overlap = read_tsv('TLGBPL_Fst_overlaping_regions.txt')
## Drop windows with less than 3 SNPs

TLGBPL_plot_fst2 =
  ggplot(data = normal, 
                  aes(x = win_mid_mb,
             y = FST_mean, 
             group = AC_CHR)) +
  geom_point(col = '#8C8C8C', 
             size = 2)+
  geom_point(data = fst_outlier,
             aes(x = win_mid_mb, 
                 y = FST_mean, 
                 group = AC_CHR),
             col = '#5488E3', 
             size = 2)+
  geom_point(data = iso_outlier,
             aes(x = win_mid_mb, 
                 y = FST_mean, 
                 group = AC_CHR),
             col = '#2A8C39', 
             size = 2)+
  geom_point(data = morpho_outlier,
             aes(x = win_mid_mb, 
                 y = FST_mean, 
                 group = AC_CHR),
             col = '#DB0808', 
             size = 2)+
  geom_point(data = Both_outlier,
             aes(x = win_mid_mb, 
                 y = FST_mean, 
                 group = AC_CHR),
             col = '#F2CD13', 
             size = 2)+
  geom_point(data = overlap, 
             aes(x = win_mid_mb, 
                 y = FST_mean, 
                 group =AC_CHR), 
             col = '#3D1673',
             size = 2)+
    facet_grid(~ AC_CHR, 
             scales = 'free')+
  labs(x = 'Chromosomal position (Mb)', 
       y = 'Fst',
       title = 'Thingvallavatn LGB-PL')+
  # ylim(0.00, 1.00)+
  theme(panel.grid = element_blank(), 
        axis.text.x = element_blank(),
        # axis.text.x = element_text(size = 12,
        #                            angle = 90,
        #                            hjust = 1),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16),
        # axis.title.x = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face = 'bold',
                                  size = 14),
        strip.background = element_rect(fill = 'white',
                                        colour = 'black'), 
        plot.title = element_text(size = 18,
                                  face = 'bold'))

GSBPI_plot_fst2
SLGBPEL_plot_fst2
TLGBPL_plot_fst2
TSBPL_plot_fst2


VBRSIL_plot_fst2 =
  ggplot(data = normal, 
         aes(x = win_mid_mb,
             y = FST_mean, 
             group = AC_CHR)) +
  geom_point(col = '#8C8C8C', 
             size = 2)+
  geom_point(data = fst_outlier,
             aes(x = win_mid_mb, 
                 y = FST_mean, 
                 group = AC_CHR),
             col = '#5488E3', 
             size = 2)+
  # geom_point(data = Both_outlier,
  #            aes(x = win_mid_mb, 
  #                y = FST_mean, 
  #                group = AC_CHR),
  #            col = '#EF0E78', 
  #            size = 2)+
  # geom_point(data = iso_outlier,
  #            aes(x = win_mid_mb, 
  #                y = FST_mean, 
  #                group = AC_CHR),
  #            col = '#53C997', 
  #            size = 2)+
  # geom_point(data = morpho_outlier,
  #            aes(x = win_mid_mb, 
  #                y = FST_mean, 
  #                group = AC_CHR),
  #            col = '#F25D27', 
  #            size = 2)+
  geom_point(data = overlap, 
             aes(x = win_mid_mb, 
                 y = FST_mean, 
                 group =AC_CHR), 
             col = '#3D1673', 
             size = 2)+
  facet_grid(~ AC_CHR, 
             scales = 'free')+
  labs(x = 'Chromosomal position (Mb)', 
       y = 'Fst',
       title = 'Vatnshlidarvatn')+
  # ylim(0.00, 1.00)+
  theme(panel.grid = element_blank(), 
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 12,
                                   angle = 90,
                                   hjust = 1),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        # axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face = 'bold',
                                  size = 14),
        strip.background = element_rect(fill = 'white',
                                        colour = 'black'), 
        plot.title = element_text(size = 18,
                                  face = 'bold'))

VBRSIL_plot_fst2
#
# Bar graph ---------------------------------------------------------------

setwd('~/Fst_sliding_window/')
outliers = read_tsv('Polypops_colocalization_data_fig5.txt') %>% 
  filter(value == 'Outlier') %>% 
  group_by(population)

# outliers$population

near_outliers = read_tsv('Polypops_colocalization_data_fig5.txt') %>% 
  filter(value == 'Non-Outlier') %>% 
  group_by(population)
# near_outliers$population

neutral = read_tsv('Polypops_neutral_snps_data.txt') %>% 
  group_by(population)

# setwd('~/Fst_sliding_window/SLGBPI/')
# outliers = read_tsv('SLGBPI_Colocalization_data.txt') %>% 
#   filter(value == 'Outlier')
# near_outliers = read_tsv('SLGBPI_Colocalization_data.txt') %>% 
#   filter(value == 'Non-Outlier')
# 
# neutral = read_tsv('SLGBPI_Neutral_snps.txt')
# 
# neutral = anti_join(neutral, 
#                     near_outliers, 
#                     by = 'SNP')

class1 = summarise(outliers, 
                   num_snps = n(),
                   FST_mean = mean(FST), 
                   FST_sd = sd(FST),
                   FST_se = sd(FST)/sqrt(n()),
                   FST_min = min(FST), 
                   FST_max = max(FST))

class2 = summarise(near_outliers,
                   num_snps = n(),
                   FST_mean = mean(FST),
                   FST_sd = sd(FST),
                   FST_se = sd(FST)/sqrt(n()),
                   FST_min = min(FST), 
                   FST_max = max(FST))

class3 = summarise(neutral, 
                   num_snps = n(),
                   FST_mean = mean(FST),
                   FST_sd = sd(FST),
                   FST_se = sd(FST)/sqrt(n()),
                   FST_min = min(FST), 
                   FST_max = max(FST))

class_data = bind_rows(class1, 
                       class2, 
                       class3) %>% 
  rowid_to_column() %>% 
mutate(class = as.factor(case_when(
  rowid == '1' ~ 'Outlier loci',
  rowid == '2' ~ 'Outlier loci',
  rowid == '3' ~ 'Outlier loci',
  rowid == '4' ~ 'Outlier loci',
  rowid == '5' ~ 'Outlier loci',
  # rowid == '6' ~ 'Outlier loci',
  rowid == '6' ~ 'Near outlier loci',
  rowid == '7' ~ 'Near outlier loci',
  rowid == '8' ~ 'Near outlier loci',
  rowid == '9' ~ 'Near outlier loci',
  rowid == '10' ~ 'Near outlier loci',
  # rowid == '12' ~ 'Near outlier loci',
  rowid == '11' ~ 'Neutral loci',
  rowid == '12' ~ 'Neutral loci',
  rowid == '13' ~ 'Neutral loci',
  rowid == '14' ~ 'Neutral loci',
  rowid == '15' ~ 'Neutral loci',
  # rowid == '18' ~ 'Neutral loci',

))) %>%
  mutate(full_morph = as.factor(case_when(
    population == 'VBRSIL' ~ 'G: Benthic - Pelagic',
    population == 'VBRSIL' ~ 'S: Benthic - Pelagic',
    # population == 'SLGBPL' ~ 'S: Benthic - Pelagic 1',
    # population == 'SLGBPI' ~ 'S: Benthic - Pelagic 2',
    population == 'VBRSIL' ~ 'T: Benthic 1 - Pelagic',
    population == 'VBRSIL' ~ 'T: Benthic 2 - Pelagic',
    population == 'VBRSIL' ~ 'V: Benthic - Pelagic'))) %>%
arrange(population)
  # mutate(class = as.factor(case_when(
  #   rowid == '1' ~ 'Class 1',
  #   rowid == '2' ~ 'Class 2',
  #   rowid == '3' ~ 'Class 3')))

  

colours = c('#FF3721', 
            '#E2FF3B', 
            '#5097FF')

class_data$class = factor(class_data$class, 
                          levels=c("Outlier loci", 
                                   "Near outlier loci", 
                                   "Neutral loci"))
bar_graph = class_data %>% 
  ggplot(aes(x = class, 
             y = FST_mean))+
  geom_bar(stat = 'identity', 
           aes(col = class, 
               fill = class))+
  geom_errorbar(aes(ymin = FST_mean - FST_se, 
                    ymax = FST_mean + FST_se), 
                width = 0.09, 
                position = position_dodge(0.9))+
  ylim(0.0, 1.00)+
  facet_grid(~full_morph, 
             scales = 'free')+
  scale_color_manual(values = colours)+
  scale_fill_manual(values = colours) +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text.x = element_text(size = 12, 
                                   angle = 90,
                                   hjust = 1),
        axis.text.y = element_text(size = 12), 
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none',
        # axis.title.y = element_blank(), 
        strip.text = element_text(face = 'bold',
                                  size = 10),
        strip.background = element_rect(fill = 'white',
                                        colour = 'black'))+
  labs(y = 'Fst') 

bar_graph

ggsave('Fig5_islands_differentiation.tiff', 
       width = 10,
       plot = bar_graph)
# combo and save ----------------------------------------------------------

# combo = VBRSIL_plot_fst2 + VBRSIL_bar + plot_layout(widths = c(3, 1))
# 
# ggsave('VBRSIL_Colocalization_FstOutliers_200kb_plotstyle2.tiff',
#        plot = combo, 
#        width = 30,
#        height = 10,
#        limitsize = FALSE)

# VBRSIL1 = VBRSIL_plot_fst1 
# VBRSIL1 = VBRSIL_plot_fst1 + VBRSIL_bar + plot_layout(widths = c(3, 1))
# VBRSIL1 = VBRSIL_plot_fst1 + VBRSIL_bar + plot_layout(widths = c(3, 1))
# VBRSIL1 = VBRSIL_plot_fst1 + VBRSIL_bar + plot_layout(widths = c(3, 1))
# SLGBPL1 = SLGBPL_plot_fst1 + SLGBPL_bar + plot_layout(widths = c(3, 1))
# SLGBPI1 = SLGBPI_plot_fst1 + SLGBPI_bar + plot_layout(widths = c(3, 1))
# 
# total_style1 = VBRSIL1/SLGBPL1/SLGBPI1/VBRSIL1/VBRSIL1/VBRSIL1
# 
# ggsave('ALLPOPS_Colocalization_FstOutliers_200kb_plotstyle1.tiff',
#       plot = total_style1, 
#       width = 30,
#       height = 15,
#       limitsize = FALSE)

GSBPI2 = GSBPI_plot_fst2
SLGBPEL2 = SLGBPEL_plot_fst2
TLGBPL2 = TLGBPL_plot_fst2
TSBPL2 = TSBPL_plot_fst2
# SLGBPL =  SLGBPL_plot_fst2
# SLGBPI2 = SLGBPI_plot_fst2
VBRSIL2 = VBRSIL_plot_fst2

total_style2 = GSBPI2 / SLGBPEL2 / TLGBPL2 / TSBPL2 / VBRSIL2

ggsave('ALLPOPS_Colocalization_FstOutliers_200kb_13.04.2021.tiff',
       plot = total_style2, 
       width = 33,
       height = 20,
       dpi = 'retina',
       limitsize = FALSE)

##
# Gene_graph --------------------------------------------------------------
theme_set(theme_bw())

outlier %>% 
  group_by(AC_CHR) %>% 
  top_n(5, FST_mean) %>% 
  distinct(FST_mean, 
           AC_CHR) %>% 
  arrange(AC_CHR) %>% 
  View()

outs = outlier %>%
  filter(AC_CHR %in% c('AC01', 
                       'AC02', 
                       'AC04p', 
                       'AC04q.2', 
                       'AC05'))
# filter(AC_CHR %in% c('AC02', 
#                      'AC15'))
# filter(AC_CHR %in% c('AC02',
#                      'AC04p',
#                      'AC05',
#                      'AC07',
#                      'AC14',
#                      'AC15',
#                      'AC25',
#                      'AC32',
#                      'AC33',
#                      'AC36'))
normy = normal %>%
  filter(AC_CHR %in% c('AC01', 
                       'AC02', 
                       'AC04p', 
                       'AC04q.2', 
                       'AC05'))
  # filter(AC_CHR %in% c('AC02', 
  #                      'AC15'))
# filter(AC_CHR %in% c('AC02',
#                      'AC04p',
#                      'AC05',
#                      'AC07',
#                      'AC14',
#                      'AC15',
#                      'AC25',
#                      'AC32',
#                      'AC33',
#                      'AC36'))

Gene_data = read_csv('Polypops_Common_genes.csv') %>% 
  # rename(AC_CHR = X6) %>% 
  # select(-AKA) %>% 
  group_by(AC_CHR) %>% 
  mutate(start_mb = start/1000000,
         end_mb = end/1000000) %>%
  select(salvelinus_id,
         AC_CHR,
         start, 
         end,
         start_mb,
         end_mb, 
         distance) %>%
  rename(xmin = start_mb,
         xmax = end_mb) %>%
  filter(distance == '100') %>% 
  group_by(AC_CHR) %>%
  arrange(AC_CHR) 
  # ungroup()


VBRSIL_Gene_100 = ggplot() +
  geom_point(data = normy,
             aes(x = win_mid_mb,
                 y = FST_mean,
                 group = AC_CHR),
             col = 'grey49')+
  geom_smooth(data = normy,
              aes(x = win_mid_mb,
                  y = FST_mean,
                  group = AC_CHR),
              col = '#1F2440',
              size = 2)+
  geom_point(data = outs,
             aes(x = win_mid_mb,
                 y = FST_mean,
                 group = AC_CHR),
             col = '#FF3721')+
    # geom_segment(data = Gene_data,
    #              aes(x = xmin,
    #                xend = xmax,
    #                y = 0.4,
    #                yend = 0.5),
    #              # alpha = 0.1,
    #              fill = "#B8004B",
    #              col = '#B8004B',
    #              inherit.aes = FALSE)+
  facet_grid(~ AC_CHR, 
             scales = 'free')+
  labs(x = 'Chromosomal position (Mb)', 
       y = 'Fst',
       title = 'G: Benthic - Pelagic')+
  ylim(0.00, 1.00)+
   theme(panel.grid = element_blank(), 
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 12,
                                   angle = 90,
                                   hjust = 1),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        # axis.title.x = element_blank(),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(face = 'bold',
                                  size = 14),
        strip.background = element_rect(fill = 'white',
                                        colour = 'black'), 
        plot.title = element_text(size = 18,
                                  face = 'bold'))

VBRSIL_Gene_100
SLGBPL_Gene_100
SLGBPL_Gene_100
VBRSIL_Gene_100
VBRSIL_Gene_100
VBRSIL_Gene_100


Total_Graph = VBRSIL_Gene_100/SLGBPL_Gene_100/SLGBPI_Gene_100/VBRSIL_Gene_100/VBRSIL_Gene_100/VBRSIL_Gene_100

ggsave('Polypops_Shared_Genes_100kb.tiff',
       plot = Total_Graph, 
       width = 15,
       height = 15,
       limitsize = FALSE)


# Chromosome Zoom ---------------------------------------------------------

theme_set(theme_bw())

outlier %>% 
  group_by(AC_CHR) %>% 
  top_n(5, FST_mean) %>% 
  distinct(FST_mean, 
           AC_CHR) %>% 
  arrange(-FST_mean) %>% 
  ungroup() %>% 
  slice(1:5) 

outs = outlier %>%
  # filter(AC_CHR %in% c('ACO1', 
  #                      'AC06', 
  #                      'AC09',
  #                      'AC20', 
  #                      'AC29'))
  # filter(AC_CHR %in% c('AC04q.2', 
  #                      'AC19', 
  #                      'AC20',
  #                      'AC26', 
  #                      'AC30'))
  # filter(AC_CHR %in% c('AC01', 
  #                      'AC06', 
  #                      'AC09',
  #                      'AC20', 
  #                      'AC29'))
  # filter(AC_CHR %in% c('AC02', 
  #                    'AC04p', 
  #                    'AC13',
  #                    'AC15', 
  #                    'AC20'))
# filter(AC_CHR %in% c('AC15', 
#                      'AC19', 
#                      'AC20',
#                      'AC21', 
#                      'AC35'))
  filter(AC_CHR %in% c('AC06', 
                     'AC10', 
                     'AC24',
                     'AC25', 
                     'AC26'))
normy = normal %>%
  filter(AC_CHR %in% c('AC06', 
                       'AC10', 
                       'AC24',
                       'AC25', 
                       'AC26'))

ggplot() +
  geom_point(data = normy,
             aes(x = win_mid_mb,
                 y = FST_mean,
                 group = AC_CHR),
             col = 'grey49')+
  geom_smooth(data = normy,
              aes(x = win_mid_mb,
                  y = FST_mean,
                  group = AC_CHR),
              col = '#1F2440',
              size = 2)+
  geom_point(data = outs,
             aes(x = win_mid_mb,
                 y = FST_mean,
                 group = AC_CHR),
             col = '#FF3721')+
  facet_grid(~ AC_CHR, 
             scales = 'free')+
  labs(x = 'Chromosomal position (Mb)', 
       y = 'Fst',
       title = 'V: Benthic - Pelagic')+
  ylim(0.00, 1.00)+
  theme(panel.grid = element_blank(), 
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 12,
                                   angle = 90,
                                   hjust = 1),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        # axis.title.x = element_blank(),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(face = 'bold',
                                  size = 14),
        strip.background = element_rect(fill = 'white',
                                        colour = 'black'), 
        plot.title = element_text(size = 18,
                                  face = 'bold'))

ggsave('VBRSIL_Top5_Div_Chrom.tiff',
       plot = last_plot())



# Fst_distribution --------------------------------------------------------

## ploting the distribution of Fst values within each morph pair
## Here we are not using the distribution of sliding windows
## we have that graph and I think this may be more informative

setwd('~/PhD/SNP Demographic modelling/Outliers_directory/Fst_data')

fst_data = read_tsv('VBRSIL_fst.FST')%>% 
  mutate(AC_CHR = as.factor(case_when(
    CHR == '1' ~ 'AC01',
    CHR == '2' ~ 'AC02',
    CHR == '3' ~ 'AC03',
    CHR == '4' ~ 'AC04p',
    CHR == '5' ~ 'AC04q.1:29',
    CHR == '6' ~ 'AC04q.2',
    CHR == '7' ~ 'AC05',
    CHR == '8' ~ 'AC06',
    CHR == '9' ~ 'AC06',
    CHR == '10' ~ 'AC07',
    CHR == '11' ~ 'AC08',
    CHR == '12' ~ 'AC09',
    CHR == '13' ~ 'AC10',
    CHR == '14' ~ 'AC11',
    CHR == '15' ~ 'AC12',
    CHR == '16' ~ 'AC13',
    CHR == '17' ~ 'AC14',
    CHR == '18' ~ 'AC15',
    CHR == '19' ~ 'AC16',
    CHR == '20' ~ 'AC17',
    CHR == '21' ~ 'AC18',
    CHR == '22' ~ 'AC19',
    CHR == '23' ~ 'AC20',
    CHR == '24' ~ 'AC21',
    CHR == '25' ~ 'AC22',
    CHR == '26' ~ 'AC23',
    CHR == '27' ~ 'AC24',
    CHR == '28' ~ 'AC25',
    CHR == '29' ~ 'AC26',
    CHR == '30' ~ 'AC27',
    CHR == '31' ~ 'AC28',
    CHR == '32' ~ 'AC30',
    CHR == '33' ~ 'AC31',
    CHR == '34' ~ 'AC32',
    CHR == '35' ~ 'AC33',
    CHR == '36' ~ 'AC34',
    CHR == '37' ~ 'AC35',
    CHR == '38' ~ 'AC36',
    CHR == '39' ~ 'AC37',
    CHR == '0' ~ 'Contigs')))

iso_outlier = read_tsv('SLGBPEL_ISOoutliers_Fst.txt')
Morpho_outliers = read_tsv('SLGBPEL_MORPHOoutliers_Fst.txt')
Overlap_outliers = read_tsv('SLGBPEL_BOTHoutliers_Fst.txt')


theme_set(theme_bw())

vbrsil = fst_data %>% 
  ggplot(aes(x = FST))+
  geom_density(fill = '#BF1B1B',
               alpha = 0.8)+
  # geom_density(data = iso_outlier, 
  #              aes(x = FST), 
  #              fill = '#58F252', 
  #              alpha = 0.8)+
  # geom_density(data = Morpho_outliers, 
  #               aes(x = FST), 
  #               fill = '#274F73', 
  #               alpha = 0.8)+
  # geom_density(data = Overlap_outliers, 
  #              aes(x = FST), 
  #              fill = '#F28705', 
  #              alpha = 0.8)+
  labs(x = 'Fst value', 
       y = 'Frequency', 
       title = 'Vatnshlidarvatn')+
  # scale_fill_manual(values = '#BF1B1B')+
  theme(panel.grid = element_blank(), 
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12), 
        axis.ticks = element_line(size = 1), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold', 
                                  size = 12), 
        legend.position = 'none')
gsbpi
slgbpel
tlgbpl
tsbpl
vbrsil

fst_outliers = (gsbpi + slgbpel)/(tlgbpl + tsbpl) 

fst_allplots = fst_outliers | vbrsil

ggsave('Fst_distributions_allpops.tiff', 
       plot = fst_allplots, 
       units = 'cm', 
       width = 27,
       height = 13.5,
       dpi = 'retina')



# Genomic parallelism windows ---------------------------------------------

Gal_both_outlier = Both_outlier
Gal_ISO_outlier = iso_outlier
Gal_morpho_outlier = morpho_outlier
Gal_normal = normal
Gal_fst_outlier = fst_outlier


Svin_both_outlier = Both_outlier
Svin_ISO_outlier = iso_outlier
Svin_morpho_outlier = morpho_outlier
Svin_normal = normal
Svin_fst_outlier = fst_outlier

TLGBPL_both_outlier = Both_outlier
TLGBPL_ISO_outlier = iso_outlier
TLGBPL_morpho_outlier = morpho_outlier
TLGBPL_normal = normal
TLGBPL_fst_outlier = fst_outlier

TSBPL_both_outlier = Both_outlier
TSBPL_ISO_outlier = iso_outlier
TSBPL_morpho_outlier = morpho_outlier
TSBPL_normal = normal
TSBPL_fst_outlier = fst_outlier

## 

intersect(Gal_fst_outlier$win_mid, 
          Svin_fst_outlier$win_mid)

g_df = 
  # Gal_both_outlier %>% 
  # Gal_morpho_outlier %>% 
  Gal_ISO_outlier %>%
  # Gal_fst_outlier %>%
  distinct(
    # AC_CHR, 
           # win_mid, 
           FST_mean,
           .keep_all = TRUE) %>% 
  select(CHR, 
         AC_CHR,
         win_mid, 
         win_mid_mb)

s_df = 
  # Svin_both_outlier %>% 
  # Svin_morpho_outlier %>% 
  # Svin_ISO_outlier %>%
  Svin_fst_outlier %>%
  distinct(FST_mean, 
           .keep_all = TRUE) %>% 
  filter(CHR != 0)%>% 
  select(CHR, 
         AC_CHR,
         win_mid, 
         win_mid_mb)

t1_df = 
  # TLGBPL_both_outlier %>% 
  # TLGBPL_morpho_outlier %>% 
  # TLGBPL_ISO_outlier %>%
  TLGBPL_fst_outlier %>%
  distinct(
    # AC_CHR, 
    # win_mid, 
    FST_mean,
    .keep_all = TRUE) %>% 
  select(CHR, 
         AC_CHR,
         win_mid, 
         win_mid_mb)

t2_df = 
  # TSBPL_both_outlier %>% 
  # TSBPL_morpho_outlier %>% 
  # TSBPL_ISO_outlier %>%
  TSBPL_fst_outlier %>%
  distinct(
    # AC_CHR, 
    # win_mid, 
    FST_mean,
    .keep_all = TRUE) %>% 
  select(CHR, 
         AC_CHR,
         win_mid, 
         win_mid_mb)

  intersect(g_df, 
            s_df)

  intersect(g_df, 
            t1_df)
  
  intersect(g_df, 
            t2_df)
  
  intersect(s_df, 
            t1_df)
  
  intersect(s_df, 
            t2_df)
  
  intersect(t1_df, 
            t2_df)
  
test = inner_join(g_df, 
             s_df)
  
test = inner_join(test, 
                  t1_df)

test = inner_join(test, 
                  t2_df)
## 4 windows that are parallel across all gst1t2


# Table 2 window calculations ---------------------------------------------
win_size = 190000
Gal_both_outlier%>% 
  select(AC_CHR, 
         # value, 
         win_start, 
         win_end, 
         win_mid, 
         win_mid_mb,
         FST_n, 
         FST_mean) %>% 
  # View()
  group_by(AC_CHR) %>%
  arrange(AC_CHR) %>%
  mutate(round_win = round(win_mid_mb,
                           digits = 2)) %>%
  mutate(distinct_win = factor(round_win)) %>%
  distinct(distinct_win,
           .keep_all = T) %>%
  ungroup() %>%
  group_by(AC_CHR,
           FST_n,
           FST_mean) %>%
  summarise(num_win = n())

Gal_ISO_outlier %>% 
  select(AC_CHR, 
         value, 
         win_start, 
         win_end, 
         win_mid, 
         win_mid_mb,
         FST_n, 
         FST_mean) %>% 
  # View()
  group_by(AC_CHR) %>%
  arrange(AC_CHR) %>%
  mutate(round_win = round(win_mid_mb,
                              digits = 2)) %>%
  mutate(distinct_win = factor(round_win)) %>%
  distinct(distinct_win,
           .keep_all = T) %>%
  ungroup() %>%
  group_by(AC_CHR,
           FST_n,
           FST_mean) %>%
  summarise(num_win = n())
  
Gal_morpho_outlier%>% 
  select(AC_CHR, 
         value, 
         win_start, 
         win_end, 
         win_mid, 
         win_mid_mb,
         FST_n, 
         FST_mean) %>% 
  # View()
  group_by(AC_CHR) %>%
  arrange(AC_CHR) %>%
  mutate(round_win = round(win_mid_mb,
                           digits = 2)) %>%
  mutate(distinct_win = factor(round_win)) %>%
  distinct(distinct_win,
           .keep_all = T) %>%
  ungroup() %>%
  group_by(AC_CHR,
           FST_n,
           FST_mean) %>%
  summarise(num_win = n())

# Gal_normal
# Gal_fst_outlier




Svin_both_outlier %>% 
  select(AC_CHR, 
         # value, 
         win_start, 
         win_end, 
         win_mid, 
         win_mid_mb,
         FST_n, 
         FST_mean) %>% 
  # View()
  group_by(AC_CHR) %>%
  arrange(AC_CHR) %>%
  mutate(round_win = round(win_mid_mb,
                           digits = 2)) %>%
  mutate(distinct_win = factor(round_win)) %>%
  distinct(distinct_win,
           .keep_all = T) %>%
  ungroup() %>%
  group_by(AC_CHR,
           FST_n,
           FST_mean) %>%
  summarise(num_win = n())
Svin_ISO_outlier%>% 
  select(AC_CHR, 
         # value, 
         win_start, 
         win_end, 
         win_mid, 
         win_mid_mb,
         FST_n, 
         FST_mean) %>% 
  # View()
  group_by(AC_CHR) %>%
  arrange(AC_CHR) %>%
  mutate(round_win = round(win_mid_mb,
                           digits = 2)) %>%
  mutate(distinct_win = factor(round_win)) %>%
  distinct(distinct_win,
           .keep_all = T) %>%
  ungroup() %>%
  group_by(AC_CHR,
           FST_n,
           FST_mean) %>%
  summarise(num_win = n())

Svin_morpho_outlier %>% 
  select(AC_CHR, 
         # value, 
         win_start, 
         win_end, 
         win_mid, 
         win_mid_mb,
         FST_n, 
         FST_mean) %>% 
  # View()
  group_by(AC_CHR) %>%
  arrange(AC_CHR) %>%
  mutate(round_win = round(win_mid_mb,
                           digits = 2)) %>%
  mutate(distinct_win = factor(round_win)) %>%
  distinct(distinct_win,
           .keep_all = T) %>%
  ungroup() %>%
  group_by(AC_CHR,
           FST_n,
           FST_mean) %>%
  summarise(num_win = n())

Svin_normal = normal
Svin_fst_outlier = fst_outlier

TLGBPL_both_outlier %>% 
  select(AC_CHR, 
         # value, 
         win_start, 
         win_end, 
         win_mid, 
         win_mid_mb,
         FST_n, 
         FST_mean) %>% 
  # View()
  group_by(AC_CHR) %>%
  arrange(AC_CHR) %>%
  mutate(round_win = round(win_mid_mb,
                           digits = 2)) %>%
  mutate(distinct_win = factor(round_win)) %>%
  distinct(distinct_win,
           .keep_all = T) %>%
  ungroup() %>%
  group_by(AC_CHR,
           FST_n,
           FST_mean) %>%
  summarise(num_win = n())
TLGBPL_ISO_outlier %>% 
  select(AC_CHR, 
         # value, 
         win_start, 
         win_end, 
         win_mid, 
         win_mid_mb,
         FST_n, 
         FST_mean) %>% 
  # View()
  group_by(AC_CHR) %>%
  arrange(AC_CHR) %>%
  mutate(round_win = round(win_mid_mb,
                           digits = 2)) %>%
  mutate(distinct_win = factor(round_win)) %>%
  distinct(distinct_win,
           .keep_all = T) %>%
  ungroup() %>%
  group_by(AC_CHR,
           FST_n,
           FST_mean) %>%
  summarise(num_win = n())
TLGBPL_morpho_outlier%>% 
  select(AC_CHR, 
         # value, 
         win_start, 
         win_end, 
         win_mid, 
         win_mid_mb,
         FST_n, 
         FST_mean) %>% 
  # View()
  group_by(AC_CHR) %>%
  arrange(AC_CHR) %>%
  mutate(round_win = round(win_mid_mb,
                           digits = 2)) %>%
  mutate(distinct_win = factor(round_win)) %>%
  distinct(distinct_win,
           .keep_all = T) %>%
  ungroup() %>%
  group_by(AC_CHR,
           FST_n,
           FST_mean) %>%
  summarise(num_win = n())
TLGBPL_normal = normal
TLGBPL_fst_outlier = fst_outlier

TSBPL_both_outlier %>% 
  select(AC_CHR, 
         # value, 
         win_start, 
         win_end, 
         win_mid, 
         win_mid_mb,
         FST_n, 
         FST_mean) %>% 
  # View()
  group_by(AC_CHR) %>%
  arrange(AC_CHR) %>%
  mutate(round_win = round(win_mid_mb,
                           digits = 2)) %>%
  mutate(distinct_win = factor(round_win)) %>%
  distinct(distinct_win,
           .keep_all = T) %>%
  ungroup() %>%
  group_by(AC_CHR,
           FST_n,
           FST_mean) %>%
  summarise(num_win = n())

TSBPL_ISO_outlier %>% 
  select(AC_CHR, 
         # value, 
         win_start, 
         win_end, 
         win_mid, 
         win_mid_mb,
         FST_n, 
         FST_mean) %>% 
  # View()
  group_by(AC_CHR) %>%
  arrange(AC_CHR) %>%
  mutate(round_win = round(win_mid_mb,
                           digits = 2)) %>%
  mutate(distinct_win = factor(round_win)) %>%
  distinct(distinct_win,
           .keep_all = T) %>%
  ungroup() %>%
  group_by(AC_CHR,
           FST_n,
           FST_mean) %>%
  summarise(num_win = n())

TSBPL_morpho_outlier %>% 
  select(AC_CHR, 
         # value, 
         win_start, 
         win_end, 
         win_mid, 
         win_mid_mb,
         FST_n, 
         FST_mean) %>% 
  # View()
  group_by(AC_CHR) %>%
  arrange(AC_CHR) %>%
  mutate(round_win = round(win_mid_mb,
                           digits = 2)) %>%
  mutate(distinct_win = factor(round_win)) %>%
  distinct(distinct_win,
           .keep_all = T) %>%
  ungroup() %>%
  group_by(AC_CHR,
           FST_n,
           FST_mean) %>%
  summarise(num_win = n())
TSBPL_normal  %>% 
  select(AC_CHR, 
         # value, 
         win_start, 
         win_end, 
         win_mid, 
         win_mid_mb,
         FST_n, 
         FST_mean) %>% 
  # View()
  group_by(AC_CHR) %>%
  arrange(AC_CHR) %>%
  mutate(round_win = round(win_mid_mb,
                           digits = 2)) %>%
  mutate(distinct_win = factor(round_win)) %>%
  distinct(distinct_win,
           .keep_all = T) %>%
  ungroup() %>%
  group_by(AC_CHR,
           FST_n,
           FST_mean) %>%
  summarise(num_win = n())

TSBPL_fst_outlier  %>% 
  select(AC_CHR, 
         # value, 
         win_start, 
         win_end, 
         win_mid, 
         win_mid_mb,
         FST_n, 
         FST_mean) %>% 
  # View()
  group_by(AC_CHR) %>%
  arrange(AC_CHR) %>%
  mutate(round_win = round(win_mid_mb,
                           digits = 2)) %>%
  mutate(distinct_win = factor(round_win)) %>%
  distinct(distinct_win,
           .keep_all = T) %>%
  ungroup() %>%
  group_by(AC_CHR,
           FST_n,
           FST_mean) %>%
  summarise(num_win = n())


