##############################
## Fst sliding window
##
## Matt Brachmann (PhDMattyB)
##
## 2020-04-01
##
##############################
# dir.create('~/Fst_sliding_window/Galtabol_fst/')

# library(devtools)
library(windowscanr)
library(patchwork)
library(tidyverse)

# setwd('~/Fst_sliding_window')
setwd('~/PhD_Genomics_Chapter3/Fst_Iceland_pops/')
# setwd('~/Fst_sliding_window/Galtabol_fst/')
# setwd('~/Fst_sliding_window/TLGBPL/')
# setwd('~/Fst_sliding_window/TSBPL/')
# setwd('~/Fst_sliding_window/VBRSIL/')
# setwd('~/Fst_sliding_window/SLGBPL/')
# setwd('~/Fst_sliding_window/SLGBPI/')

# data --------------------------------------------------------------------
## This is the data needed for the sliding window analysis. 
## We get this data from PLINK 1.9 --fst --within flags
## Need to make sure we specify the --chr-set to 39
## and make sure the contigs or unplaced regions are 
## set to 0
GSBPI_data = read_tsv('Galtabol_chr_fix.fst') %>% 
  group_by(CHR) %>% 
  filter(CHR != 0)
TLGBPL_data = read_tsv('TLGBPL_fst.fst')%>% 
  group_by(CHR)%>% 
  filter(CHR != 0)
TSBPL_data = read_tsv('TSBPL_fst.fst')%>% 
  group_by(CHR)%>% 
  filter(CHR != 0)
VBRSIL_data = read_tsv('VBRSIL_fst.fst')%>% 
  group_by(CHR)%>% 
  filter(CHR != 0)
# data = read_tsv('SLGBPL_fst.fst')
# data = read_tsv('SLGBPI_fst.fst')
SLGBPEL_data = read_tsv('SLGBPEL_fst.fst')%>% 
  group_by(CHR)%>% 
  filter(CHR != 0)



# FST outliers SNP --------------------------------------------------------
GSBPI_top5 = GSBPI_data[GSBPI_data$FST > quantile(GSBPI_data$FST, 
                                                  prob = 1-5/100),]
SLGBPEL_top5 = SLGBPEL_data[SLGBPEL_data$FST > quantile(SLGBPEL_data$FST, 
                                                  prob = 1-5/100),]
TLGBPL_top5 = TLGBPL_data[TLGBPL_data$FST > quantile(TLGBPL_data$FST, 
                                                  prob = 1-5/100),]
TSBPL_top5 = TSBPL_data[TSBPL_data$FST > quantile(TSBPL_data$FST, 
                                                  prob = 1-5/100),]
VBRSIL_top5 = VBRSIL_data[VBRSIL_data$FST > quantile(VBRSIL_data$FST, 
                                                  prob = 1-5/100),]

intersect(GSBPI_top5$SNP, 
          SLGBPEL_top5$SNP) %>% 
  as_tibble()
intersect(GSBPI_top5$SNP, 
          TLGBPL_top5$SNP) %>% 
  as_tibble()
intersect(GSBPI_top5$SNP, 
          TSBPL_top5$SNP) %>% 
  as_tibble()
intersect(GSBPI_top5$SNP, 
          VBRSIL_top5$SNP) %>% 
  as_tibble()

intersect(SLGBPEL_top5$SNP, 
          TLGBPL_top5$SNP) %>% 
  as_tibble()
intersect(SLGBPEL_top5$SNP, 
          TSBPL_top5$SNP) %>% 
  as_tibble()
intersect(SLGBPEL_top5$SNP, 
          VBRSIL_top5$SNP) %>% 
  as_tibble()

intersect(TLGBPL_top5$SNP, 
          TSBPL_top5$SNP) %>% 
  as_tibble()

intersect(TLGBPL_top5$SNP, 
          VBRSIL_top5$SNP) %>% 
  as_tibble()
intersect(TSBPL_top5$SNP, 
          VBRSIL_top5$SNP) %>% 
  as_tibble()


parallel = inner_join(GSBPI_top5, 
                      SLGBPEL_top5, 
                      by = 'SNP')
parallel = inner_join(parallel, 
                      TLGBPL_top5, 
                      by = 'SNP')
parallel = inner_join(parallel, 
                      TSBPL_top5, 
                      by = 'SNP')
parallel = inner_join(parallel, 
                      VBRSIL_top5, 
                      by = 'SNP')

theme_set(theme_bw())
GSBPI_Outlier_Dist = ggplot(data = GSBPI_top5, 
                            aes(x = FST))+
  geom_histogram(col = 'black',
                 fill = '#1A5173')+
  # geom_dotplot(col = '#1A5173', 
  #                fill = '#1A5173',
  #              method = 'histodot',
  #              binwidth = 1/100)+
  xlim(0, 1.00)+
  labs(x = 'Fst', 
       y = 'Frequncy', 
       title = 'Galtabol')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

SLGBPEL_Outlier_Dist = ggplot(data = SLGBPEL_top5, 
                              aes(x = FST))+
  geom_histogram(col = 'black',
                 fill = '#3CA661')+
  # geom_dotplot(col = '#3CA661', 
  #              fill = '#3CA661',
  #              method = 'histodot',
  #              binwidth = 1/100)+
  xlim(0,1.00)+
  labs(x = 'Fst', 
       y = 'Frequncy', 
       title = 'Svinavatn')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

TLGBPL_Outlier_Dist = ggplot(data = TLGBPL_top5, 
                             aes(x = FST))+
  geom_histogram(col = 'black',
                 fill = '#BF3030')+
  # geom_dotplot(col = '#F2C230',
  #              fill = '#F2C230',
  #              method = 'histodot',
  #              binwidth = 1/100)+
  xlim(0,1.00)+
  labs(x = 'Fst', 
       y = 'Frequncy', 
       title = 'Thingvallavatn 1')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

TSBPL_Outlier_Dist = ggplot(data = TSBPL_top5, 
                            aes(x = FST))+
  geom_histogram(col = 'black',
                 fill = '#F28A80')+
  # geom_dotplot(col = '#F28A80', 
  #              fill = '#F28A80',
  #              method = 'histodot',
  #              binwidth = 1/100)+
  xlim(0,1.00)+
  labs(x = 'Fst', 
       y = 'Frequncy', 
       title = 'Thingvallavatn 2')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

VBRSIL_Outlier_Dist = ggplot(data = VBRSIL_top5, 
                             aes(x = FST))+
  geom_histogram(col = 'black',
                 fill = '#F2C230')+
  # geom_dotplot(col = '#BF3030', 
  #              fill = '#BF3030',
  #              method = 'histodot',
  #              binwidth = 1/200)+
  xlim(0,1.00)+
  labs(x = 'Fst', 
       y = 'Frequncy', 
       title = 'Vatnshlidarvatn')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))


combo = GSBPI_Outlier_Dist/SLGBPEL_Outlier_Dist/TLGBPL_Outlier_Dist/TSBPL_Outlier_Dist/VBRSIL_Outlier_Dist

ggsave('FST_outlier_distributions.tiff', 
       plot = combo, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 25)


# FST outliers 200Kb regions ----------------------------------------------

GSBPI = read_tsv('GSBPI_Fst_200Kb_window.txt') %>% 
  filter(CHR != 0, 
         FST_n >3) %>% 
  Mb_Conversion()
GSBPI_200kb = read_csv('GSBPI_Fst_200kb_outliers_19.04.2021.csv') %>% 
  filter(AC_CHR != 'Contigs')

SLGBPEL = read_tsv('SLGBPEL_Fst_200Kb_window.txt')%>% 
  filter(CHR != 0, 
         FST_n >3)%>% 
  Mb_Conversion()
SLGBPEL_200kb = read_csv('SLGBPEL_Fst_200kb_outliers_19.04.2021.csv') %>% 
  filter(AC_CHR != 'Contigs')

TLGBPL = read_tsv('TLGBPL_Fst_200Kb_window.txt')%>% 
  filter(CHR != 0, 
         FST_n >3)%>% 
  Mb_Conversion()
TLGBPL_200kb = read_csv('TLGBPL_Fst_200kb_outliers_19.04.2021.csv') %>% 
  filter(AC_CHR != 'Contigs')

TSBPL = read_tsv('TSBPL_Fst_200Kb_window.txt') %>% 
  filter(CHR != 0, 
         FST_n > 3)%>% 
  Mb_Conversion()
TSBPL_200kb = read_csv('TSBPL_Fst_200kb_outliers_19.04.2021.csv') %>% 
  filter(AC_CHR != 'Contigs')

VBRSIL = read_tsv('VBRSIL_Fst_200Kb_window.txt') %>% 
  filter(CHR != 0, 
         FST_n > 3)%>% 
  Mb_Conversion()
VBRSIL_200kb = read_csv('VBRSIL_Fst_200kb_outliers_19.04.2021.csv') %>% 
  filter(AC_CHR != 'Contigs')

## Need to pull out the Fst outliers from the Neutral dataset!!

GSBPI = anti_join(GSBPI, 
                  GSBPI_200kb)
SLGBPEL = anti_join(SLGBPEL, 
                    SLGBPEL_200kb)
TLGBPL = anti_join(TLGBPL, 
                   TLGBPL_200kb)
TSBPL = anti_join(TSBPL, 
                  TSBPL_200kb)
VBRSIL = anti_join(VBRSIL, 
                   VBRSIL_200kb)

Fst_neutral = ggplot()+
  geom_point(data = GSBPI, 
             aes(x = win_mid_mb, 
                 y = FST_mean), 
             col = '#1A5173', 
             size = 3) +
  geom_point(data = SLGBPEL, 
             aes(x = win_mid_mb, 
                 y = FST_mean), 
             col = '#3CA661', 
             size = 3) +
  geom_point(data = TLGBPL, 
             aes(x = win_mid_mb, 
                 y = FST_mean), 
             col = '#BF3030', 
             size = 3) +
  geom_point(data = TSBPL, 
             aes(x = win_mid_mb, 
                 y = FST_mean), 
             col = '#F28A80', 
             size = 3) +
  geom_point(data = VBRSIL, 
             aes(x = win_mid_mb, 
                 y = FST_mean), 
             col = '#F2C230', 
             size = 3) +
  facet_grid(~AC_CHR)+
  labs(y = 'Mean Fst', 
       x = 'Distance along chromosome (Mb)', 
       title = 'A)')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 8, 
                                   angle = 45, 
                                   hjust = 1, 
                                   vjust = 1), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(colour = 'black'))

Fst_outliers = ggplot()+
  geom_point(data = GSBPI_200kb, 
             aes(x = win_mid_mb, 
                 y = FST_mean), 
             col = '#1A5173', 
             size = 3) +
  geom_point(data = SLGBPEL_200kb, 
             aes(x = win_mid_mb, 
                 y = FST_mean), 
             col = '#3CA661', 
             size = 3) +
  geom_point(data = TLGBPL_200kb, 
             aes(x = win_mid_mb, 
                 y = FST_mean), 
             col = '#BF3030', 
             size = 3) +
  geom_point(data = TSBPL_200kb, 
             aes(x = win_mid_mb, 
                 y = FST_mean), 
             col = '#F28A80', 
             size = 3) +
  geom_point(data = VBRSIL_200kb, 
             aes(x = win_mid_mb, 
                 y = FST_mean), 
             col = '#F2C230', 
             size = 3) +
  facet_grid(~AC_CHR)+
  labs(y = 'Mean Fst', 
       x = 'Distance along chromosome (Mb)', 
       title = 'B)')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 8, 
                                   angle = 45, 
                                   hjust = 1, 
                                   vjust = 1), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(colour = 'black'))

Fst_combo = Fst_neutral/Fst_outliers

ggsave('~/PhD_Genomics_Chapter3/Fst_Iceland_pops/MeanFst_Neutral_Outliers.tiff', 
       plot = Fst_combo, 
       dpi = 'retina', 
       units = 'cm', 
       width = 45, 
       height = 15)

# Sliding window analysis -------------------------------------------------

## This is performed on the lab computer
## it makes mine get emotional. 
## Use the data above, from PLINK, to perform the analysis

fst_position = winScan(x = data, 
                       groups = 'CHR', 
                       position = 'POS',
                       values = 'FST', 
                       win_size = 200000, 
                       win_step = 199000, 
                       funs = c('mean', 'sd'))
fst_position = fst_position %>%
  as_tibble()
# write avg fst per window ------------------------------------------------
## Write the txt file for each window size. 
## Need to compare the different window sizes to see which one
## is the most appropriate. 
## small window size == greater chance for false positives
## large window size == less chance to find differences

write_tsv(fst_position, 
          'SLGBPEL_Fst_200Kb_window.txt')

# Plot data ---------------------------------------------------------------

# setwd('~/Fst_sliding_window/Galtabol_fst/')
# setwd('~/Fst_sliding_window/SLGBPL/')
# setwd('~/Fst_sliding_window/')
# setwd('~/Fst_sliding_window/TLGBPL/') 
# setwd('~/Fst_sliding_window/TSBPL/')
# setwd('~/Fst_sliding_window/VBRSIL/')

# setwd('~/PhD_Genomics_Chapter3/Fst_Iceland_pops/')

data = read_tsv('VBRSIL_Fst_200Kb_window.txt') %>% 
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



# Conversion --------------------------------------------------------------

Mb_Conversion = function(data){
  data %>% 
    group_by(AC_CHR) %>% 
    mutate(win_mid_mb = win_mid/1000000)
}

data = Mb_Conversion(data)

#

#
# Fst outliers 200Kb regions------------------------------------------------------------

data = data %>% 
  filter(FST_n > 3) %>% 
  na.omit() %>% 
  ungroup() %>% 
  filter(CHR != 0)

top5 = data[data$FST_mean > quantile(data$FST_mean, 
                                     prob = 1-5/100),]

write_csv(top5, 
         'VBRSIL_200Kb_Fst_outlier_23.07.2021.csv')
## need to combine the overlap and fst results for each popn
top5_overlap = read_tsv('Fst_window_overlap_allBPpairs.txt')

overlap_regions = data %>% 
  filter(AC_CHR %in% c('AC24', 'AC25'))

overlap = inner_join(overlap_regions, 
                     top5_overlap, 
                     by = 'win_mid') %>% 
  select(-AC_CHR.y) %>% 
  rename(AC_CHR = AC_CHR.x) %>% 
  group_by(AC_CHR) 

overlap
write_tsv(overlap, 
          'VBRSIL_Fst_overlaping_regions.txt')
#
# Avg fst outlier per chromosome ------------------------------------------
##filter the data
data = data %>% 
  filter(FST_n > 3) %>% 
  na.omit() %>% 
  ungroup()


# Fst outlier comparison --------------------------------------------------

GSBPI_outlier = read_csv('GSBPI_200Kb_Fst_outlier_23.07.2021.csv')
SLGBPEL_outlier = read_csv('SLGBPEL_200Kb_Fst_outlier_23.07.2021.csv')
TLGBPL_outlier = read_csv('TLGBPL_200Kb_Fst_outlier_23.07.2021.csv')
TSBPL_outlier = read_csv('TSBPL_200Kb_Fst_outlier_23.07.2021.csv')
VBRSIL_outlier = read_csv('VBRSIL_200Kb_Fst_outlier_23.07.2021.csv')

inner_join(GSBPI_outlier, 
           TSBPL_outlier, 
           by = 'win_mid')

##
# Fst bar graph neutral ---------------------------------------------------
data = data %>% 
  filter(FST_n > 3) %>% 
  na.omit() %>% 
  ungroup()

T_top5 = data[data$FST_mean > quantile(data$FST_mean, 
                                       prob = 1-5/100),]

write_tsv(T_top5, 
          'VBRSIL_Fst_outliers_200Kb_window.txt')
T_neutral = anti_join(data, 
                      T_top5, 
                      by = c('win_start', 
                             'win_end', 
                             'win_mid', 
                             'FST_n',
                             'FST_mean', 
                             'FST_sd', 
                             'AC_CHR', 
                             'win_mid_mb'))
write_tsv(T_neutral, 
          'VBRSIL_Fst_neutral_200Kb_window.txt')

G_neutral_mean = G_neutral %>% 
  group_by(AC_CHR) %>% 
  summarise(avg_neutral_fst = mean(FST_mean))

label = rep('V: Benthic - pelagic', 
            length(V_neutral_mean$avg_neutral_fst)) %>% 
  as_tibble()
## Bind that back together
V_neutral_data = bind_cols(V_neutral_mean, 
                           label) %>% 
  rename(Morph_compare = value)


G_neutral_data
S_neutral_data
T1_neutral_data
T2_neutral_data
V_neutral_data

neutral_data = bind_rows(G_neutral_data, 
                         S_neutral_data, 
                         T1_neutral_data, 
                         T2_neutral_data, 
                         V_neutral_data)
new_colors = c('#88A825',
               '#35203B',
               '#911146',
               '#CF4A30',
               '#ED8C2B')

Avg_fst_neutral = neutral_data %>% 
  ggplot()+
  geom_bar(aes(x = AC_CHR, 
               y = avg_neutral_fst, 
               # col = Morph_compare, 
               fill = Morph_compare),
           stat = 'identity', 
           position = 'dodge')+
  # scale_color_manual(values = new_colors)+
  scale_fill_manual(values = new_colors) +
  labs(title = 'A)', 
       x = 'Chromosome', 
       y = 'Average Fst')+
  guides(fill = guide_legend(title = 'Morph Pair'))+
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 9,
                                   angle = 90,
                                   hjust = 1), 
        axis.ticks = element_line(size = 2), 
        panel.grid = element_blank())


# Fst bar graph - outliers ------------------------------------------------
G_top5 = read_csv('GSBPI_Fst_200kb_outliers_19.04.2021.csv')
S_top5 = read_csv('SLGBPEL_Fst_200kb_outliers_19.04.2021.csv')
T1_top5 = read_csv('TLGBPL_Fst_200kb_outliers_19.04.2021.csv')
T2_top5 = read_csv('TSBPL_Fst_200kb_outliers_19.04.2021.csv')
V_top5 = read_csv('VBRSIL_Fst_200kb_outliers_19.04.2021.csv')


## find top 5% Fst outliers
G_top5 = data[data$FST_mean > quantile(data$FST_mean, 
                                       prob = 1-5/100),]

write_tsv(G_top5, 
          'GSBPI_Fst_outliers_200Kb_window.txt')
## summarise avg fst per chromosome
V_mean = V_top5 %>% 
  group_by(AC_CHR) %>% 
  summarise(avg_fst = mean(FST_mean))
## need a population label for the data frame
label = rep('V: Benthic - pelagic', 
            length(V_mean$avg_fst)) %>% 
  as_tibble()
## Bind that back together
V_data = bind_cols(V_mean, 
                   label) %>% 
  rename(Morph_compare = value)

G_data
S_data 
T1_data
T2_data
V_data

outlier_data = bind_rows(G_data, 
                         S_data, 
                         T1_data, 
                         T2_data, 
                         V_data)

# colors = c('#467BB3',
#                '#FF1E0C',
#                '#108565',
#                '#D1BA0A')

new_colors = c('#88A825',
               '#35203B',
               '#911146',
               '#CF4A30',
               '#ED8C2B')

Avg_fst_outliers = outlier_data %>% 
  ggplot()+
  geom_bar(aes(x = AC_CHR, 
               y = avg_fst, 
               # col = Morph_compare, 
               fill = Morph_compare),
           stat = 'identity', 
           position = 'dodge')+
  # scale_color_manual(values = new_colors)+
  scale_fill_manual(values = new_colors) +
  labs(title = 'B)', 
       x = 'Chromosome', 
       y = 'Average Fst')+
  guides(fill = guide_legend(title = 'Morph Pair'))+
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 9,
                                   angle = 90,
                                   hjust = 1), 
        axis.ticks = element_line(size = 2), 
        panel.grid = element_blank(), 
        legend.position = 'none')




Avg_fst_neutral/Avg_fst_outliers

ggsave('Fst_window_neutral_outliers_combined.tiff', 
       plot = last_plot())

#
# histogram of outliers ---------------------------------------------------


## histogram of mean Fst values per window
G_hist = ggplot(data = data)+
  geom_histogram(aes(x = FST_mean),
                 col = '#3E423A',
                 fill = '#3E423A')+
  geom_histogram(data = top5, 
                 aes(x = FST_mean),
                 col = '#02808E', 
                 fill = '#02808E') +
  labs(title = 'A)', 
       x = 'Mean Fst', 
       y = 'Number of windows')+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 12,
                                  face = 'bold'),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.y = element_text(size = 12),
        axis.ticks = element_line(size = 2),
        # axis.title.x = element_blank()
  )

V_hist 

#
# Fst Outlier plot --------------------------------------------------------

## plot of Fst including outliers per chromosome
V_fst = ggplot(data = data,
               aes(x = win_mid_mb,
                   y = FST_mean, 
                   group = AC_CHR)) +
  geom_point(col = '#3E423A')+
  # geom_smooth(col = '#1F2440',
  #             size = 2)+
  geom_point(data = top5, 
             aes(x = win_mid_mb, 
                 y = FST_mean, 
                 group = AC_CHR), 
             col = '#02808E')+
  geom_point(data = overlap, 
             aes(x = win_mid_mb, 
                 y = FST_mean, 
                 group = AC_CHR), 
             col = '#D9042B')+
  facet_grid(~ AC_CHR, 
             scales = 'free')+
  labs(title = 'E)', 
       x = 'Chromosomal position (Mb)', 
       y = 'Fst')+
  ylim(0.00, 1.00)+
  theme(panel.grid = element_blank(), 
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 9,
                                   angle = 90,
                                   hjust = 1),
        # axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face = 'bold',
                                  size = 10),
        strip.background = element_rect(fill = 'white',
                                        colour = 'black'), 
        plot.title = element_text(size = 12,
                                  face = 'bold'))

G_fst


# Combine plots -----------------------------------------------------------

library(patchwork)

G_hist/S_hist/T1_hist/T2_hist/V_hist

ggsave('FigS3_Fst_outlier_histograms.tiff',
       plot = last_plot(), 
       height = 9)


G_fst/S_fst/T1_fst/T2_fst/V_fst

ggsave('Fig3_Fst_200kb_window_overlap.tiff',
       plot = last_plot(), 
       width = 25,
       height = 9,
       limitsize = FALSE)


# Fst window outlier overlap ----------------------------------------------

filter_data = data %>% 
  filter(FST_n > 3) %>% 
  na.omit() %>% 
  ungroup()

write_csv(filter_data, 
          'VBRSIL_Fst_Neutral_200Kb_windows.19.04.2021.csv')
# filter_data %>%
#   group_by(AC_CHR) %>% 
#   summarise(n = n()) %>% 
#   View()


V_top5 = filter_data[filter_data$FST_mean > quantile(filter_data$FST_mean, 
                                                     prob = 1-5/100),]

write_csv(V_top5, 
          'VBRSIL_Fst_200kb_outliers_19.04.2021.csv')
# View(G_top5)


## need to figure out the number of distinct outliers 
## THis is due to the high overlap between windows
V_distinct_outliers = V_top5 %>% 
  group_by(AC_CHR) %>% 
  distinct(FST_mean)

G_distinct_outliers
S_distinct_outliers
T1_distinct_outliers
T2_distinct_outliers
V_distinct_outliers

# write_csv(S_distinct_outliers, 
#           'SLGBPEL_Fst_outliers_19.04.2021.csv')

## sort through the top 5% Fst outliers for each population
G_top5 = G_top5 %>% 
  group_by(AC_CHR) 
S_top5 = S_top5 %>% 
  group_by(AC_CHR) 
T1_top5 = T1_top5 %>% 
  group_by(AC_CHR) 
T2_top5 = T2_top5 %>% 
  group_by(AC_CHR) 
V_top5 = V_top5 %>% 
  group_by(AC_CHR) 

dplyr::intersect(G_top5,
                 S_top5,
                 by = 'win_mid')

dplyr::intersect(G_top5,
                 T1_top5,
                 by = 'win_mid')

dplyr::intersect(G_top5,
                 T2_top5,
                 by = 'win_mid')

dplyr::intersect(G_top5,
                 V_top5,
                 by = 'win_mid')
dplyr::intersect(S_top5,
                 T1_top5,
                 by = 'win_mid')

dplyr::intersect(S_top5,
                 T2_top5,
                 by = 'win_mid')
dplyr::intersect(S_top5,
                 V_top5,
                 by = 'win_mid')
dplyr::intersect(T1_top5,
                 T2_top5,
                 by = 'win_mid')
dplyr::intersect(T1_top5,
                 V_top5,
                 by = 'win_mid')
dplyr::intersect(T2_top5,
                 V_top5,
                 by = 'win_mid')


# tests = intersect(tests, 
#                   T1_top5, 
#                   by = 'win_mid')
# 
# tests = intersect(tests, 
#                   T2_top5, 
#                   by = 'win_mid')
# 
# tests = intersect(tests, 
#                   V_top5,
#                   by = 'win_mid')
# 
# write_tsv(tests, 
#           'Fst_window_overlap_allBPpairs.txt')
#
# SNP density plot -------------------------------------------------------------

SNP_density_plot = data %>% 
  group_by(AC_CHR) %>% 
  mutate(window = row_number()) %>% 
  # filter(AC_CHR == 'AC29') %>%
  filter(AC_CHR != 'Contigs') %>%
  ggplot(aes(x = window,
             y = FST_n))+
  geom_step(stat = 'identity', 
            size = 1)+
  scale_y_continuous(breaks = seq(0, 5000, by = 1))+
  # scale_x_continuous(breaks = seq(0, 1000000, by = 1000))+
  geom_hline(yintercept = 3, 
             col = 'red',
             size = 3)+
  # geom_density(stat = 'identity')
  facet_grid(~AC_CHR, 
             scales = 'free')+
  labs(x = 'Window number per chromosome',
       y = 'Number of SNPS',
       title = 'Number of SNPS per window')+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16),
        axis.ticks = element_line(size = 1), 
        axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = 'white',
                                        colour = 'black'),
        plot.title = element_text(size = 18,
                                  face = 'bold'))

# setwd('~/Fst_sliding_window/')
ggsave('SNP_density_per_window_19.04.2021.tiff',
       plot = SNP_density_plot, 
       width = 30,
       height = 10,
       limitsize = FALSE)


# Plot the results --------------------------------------------------------

## set the theme for all of the plots
theme_set(theme_bw())

## Drop windows with less than 3 SNPs

plot_fst = data %>% 
  filter(FST_n > 3) %>% 
  ggplot(aes(x = win_mid_mb,
             y = FST_mean, 
             group = AC_CHR)) +
  geom_point(col = 'grey49')+
  # geom_step(aes(x = win_mid,
  #               y = FST_mean),
  #           direction = 'mid',
  #           size = 2,
  #           col = '#1F2440')+
  geom_smooth(col = '#1F2440',
              # col = '#02E084',
              size = 2)+
  facet_grid(~ AC_CHR, 
             scales = 'free')+
  labs(title = 'SLGBPEL Fst 200kb window', 
       x = 'Chromosomal position (Mb)', 
       y = 'Fst')+
  ylim(0.00, 1.00)+
  theme(panel.grid = element_blank(), 
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 9,
                                   angle = 90,
                                   hjust = 1),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face = 'bold',
                                  size = 10),
        strip.background = element_rect(fill = 'white',
                                        colour = 'black'), 
        plot.title = element_text(size = 12,
                                  face = 'bold'))

plot_fst
# final = plot_10kb/plot_50kb/plot_100kb/plot_150kb/plot_200kb/plot_250kb

ggsave('SLGBPEL_Fst_200kb_window.tiff',
       plot = plot_fst, 
       width = 25,
       height = 10,
       limitsize = FALSE)


