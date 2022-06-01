##############################
## Co-localization of outliers and Fst peaks
##
## Matt Brachmann (PhDMattyB)
##
## 2020-04-28
##
##############################

library(devtools)
library(windowscanr)
library(patchwork)
library(tidyverse)

# setwd('~/Fst_sliding_window')
# setwd('~/Fst_sliding_window/Galtabol_fst/')
# setwd('~/Fst_sliding_window/TLGBPL/')
# setwd('~/Fst_sliding_window/TSBPL/')
# setwd('~/Fst_sliding_window/VBRSIL/')
# setwd('~/Fst_sliding_window/SLGBPL/')
# setwd('~/Fst_sliding_window/SLGBPI/')
# setwd('~/Fst_sliding_window/')

setwd('~/PhD_Genomics_Chapter3/Fst_Iceland_pops/')

# Data manipulation -------------------------------------------------------
data = read_tsv('Galtabol_chr_fix.fst')
data = read_tsv('SLGBPEL_fst.fst')
data = read_tsv('TLGBPL_fst.fst')
data = read_tsv('TSBPL_fst.fst')


data %>% 
  group_by(CHR)

outliers = read_csv('TSBPL_Outlier_data.csv') %>% 
  rename(SNP = `Marker ID`) %>% 
  select(SNP, 
         `#Chromosome`,
         `Genetic distance`,
         `Physical position`, 
         CHROME3, 
         value) %>% 
  arrange(SNP)


# RDA outlier Fst comparison ----------------------------------------------
RDA_outliers = inner_join(outliers, 
                          data, 
                          by = 'SNP')

RDA_outliers %>% 
  ggplot(aes(x = value, 
             y = FST))+
  geom_violin()

RDA_outliers %>% 
  group_by(value) %>% 
  summarise(mean_trait = mean(FST))

RDA_loci_diffs = aov(data = RDA_outliers, 
                     FST ~ value)

summary(RDA_loci_diffs)
TukeyHSD(RDA_loci_diffs, 
         conf.level = 0.99)
plot(TukeyHSD(RDA_loci_diffs))

# RDA outliers vs other loci ----------------------------------------------

RDA_outliers = inner_join(outliers, 
           data, 
           by = 'SNP')
RDA_label = rep('RDA_outlier', 
                nrow(RDA_outliers)) %>% 
  as_tibble()
RDA_outliers = bind_cols(RDA_outliers, 
                         RDA_label)%>% 
  select(CHR, 
         SNP, 
         POS, 
         NMISS, 
         FST, 
         value...11) %>% 
  rename(value = value...11)

Neutral_snps = anti_join(data, 
          outliers, 
          by = 'SNP')
neutral_label = rep('Neutral', 
                    nrow(Neutral_snps)) %>% 
  as_tibble()

Neutral_snps = bind_cols(Neutral_snps, 
                         neutral_label)

ultimate_data = bind_rows(RDA_outliers, 
                          Neutral_snps)

ultimate_data %>% 
  ggplot(aes(x = value, 
             y = FST))+
  geom_violin()

ultimate_data %>% 
  group_by(value) %>% 
  summarise(mean_FST = mean(FST))

aov_test = aov(data = ultimate_data, 
               FST ~ value)
summary(aov_test)

##
# Old shit we don't need anymore ------------------------------------------


# poly_outliers = read_csv('Polypop_Outlier_data.csv')%>% 
#   rename(SNP = `Marker ID`) %>% 
#   select(SNP, 
#          `#Chromosome`, 
#          `Genetic distance`, 
#          `Physical position`, 
#          CHROME3) %>% 
#   arrange(SNP)
# 
# poly_outliers_new = anti_join(poly_outliers,
#           outliers,
#           by = c('SNP',
#                  '#Chromosome',
#                  'Genetic distance',
#                  'Physical position',
#                  'CHROME3'))
# 
# outliers = full_join(outliers, 
#                  poly_outliers) %>% 
#   arrange(SNP) %>% 
#   distinct()

  # View()

# outliers = left_join(outliers, 
#                      poly_outliers, 
#                      by = c('SNP')) %>% 
#   select(SNP, 
#          `#Chromosome.x.x`, 
#          `Genetic distance.x.x`, 
#          `Physical position.x.x`, 
#          CHROME3.x.x)
# 
# inner_join(outliers,
#            poly_outliers,
#            by = c('SNP',
#                   '#Chromosome.x',
#                   'Genetic distance.x',
#                   'Physical position.x',
#                   'CHROME3.x'))



# ## VBRSIL special case
# not_outliers = anti_join(data,
#                          poly_outliers,
#                          by = 'SNP')
# ####

not_outliers = anti_join(data,
                         outliers,
                         by = 'SNP')

grouping = rep('Non-Outlier',
               length(not_outliers$FST)) %>%
  as_tibble()

not_outliers = bind_cols(not_outliers,
                         grouping)

write_tsv(not_outliers,
          'SLGBPEL_Neutral_snps.txt', col_names = TRUE)

# outliers$value
overlap_outliers = outliers %>% 
  filter(value == 'Overlapping locus')

BothTraits_Fst = left_join(overlap_outliers, 
                           data, 
                           by = 'SNP') %>% 
  select(SNP, CHR, POS, FST, CHROME3, value)


iso_trait = outliers %>% 
  filter(value == 'Isotope locus')
Iso_fst = left_join(iso_trait, 
                    data, 
                    by = 'SNP') %>% 
  select(SNP, CHR, POS, FST, CHROME3, value)

morpho_trait = outliers %>% 
  filter(value == 'Morphological locus')

morpho_fst = left_join(morpho_trait, 
                       data, 
                       by = 'SNP') %>% 
  select(SNP, CHR, POS, FST, CHROME3, value)

# Fst_outliers = left_join(outliers,
#                    data,
#                    by = 'SNP') %>%
#   select(CHR,
#          SNP,
#          POS,
#          NMISS,
#          FST)

# ##Vatnshlidarvatn isssues!!
# Fst_outliers = left_join(poly_outliers,
#                          data,
#                          by = 'SNP') %>%
#   select(CHR,
#          SNP,
#          POS,
#          NMISS,
#          FST) %>% 
#   distinct()

# grouping = rep('Both Outliers', 
#                length(BothTraits_Fst$FST)) %>% 
#   as_tibble()
# 
# BothTraits_Fst = bind_cols(BothTraits_Fst, 
#                            grouping)
# grouping = rep('Iso outlier', 
#                length(Iso_fst$FST)) %>% 
#   as_tibble()
# Iso_fst = bind_cols(Iso_fst, 
#                     grouping)
# 
# grouping = rep('Morpho outlier', 
#                length(morpho_fst$FST)) %>% 
#   as_tibble()
# morpho_fst = bind_cols(morpho_fst, 
#                        grouping)
# 
# grouping = rep('Outlier',
#                length(Fst_outliers$FST)) %>% 
#   as_tibble()
# 
# Fst_outliers = bind_cols(Fst_outliers, 
#                          grouping)

## write the three files as csvs
write_tsv(BothTraits_Fst, 
          'SLGBPEL_BOTHoutliers_Fst.txt')

write_tsv(Iso_fst, 
          'SLGBPEL_ISOoutliers_Fst.txt')

write_tsv(morpho_fst, 
          'SLGBPEL_MORPHOoutliers_Fst.txt')

## read in the neutral data
df_nonout = read_tsv('TSBPL_Neutral_snps.txt')

##read in the outlier data
# df_outliers = read_tsv('TSBPL_outliers_Fst.txt')

overlap_outs = read_tsv('TSBPL_BOTHoutliers_Fst.txt') %>% 
  filter(CHR != 0)
iso_outs = read_tsv('TSBPL_ISOoutliers_Fst.txt') %>% 
  filter(CHR != 0)
morpho_outs = read_tsv('TSBPL_MORPHOoutliers_Fst.txt') %>% 
  filter(CHR != 0)

range_cal = function(data){
  data %>% 
    arrange(CHR, 
            POS) %>% 
    group_by(SNP, CHR) %>% 
    mutate(low = POS - 200000,
           high = POS + 200000)
}

# df_outliers = range_cal(df_outliers)
overlap_outs = range_cal(overlap_outs)
iso_outs = range_cal(iso_outs)
morpho_outs = range_cal(morpho_outs)


# df_outliers = range_cal(df_outliers)
# View(df_outliers)
## CAM came in clutch and showed us how to automate this  
  in_range = function(val, min, max){
    if((val >= min) & (val <= max)){
      return(TRUE)
    }
    return(FALSE)
  }

list_of_dfs = list()

df_outliers = overlap_outs
df_outliers = iso_outs
df_outliers = morpho_outs


for(i in 1:nrow(df_outliers)){
  
  #moved the subsetting out of lapply and here so its easier to keep track of
  id = df_outliers[['SNP']][[i]]
  min = df_outliers[['low']][[i]]
  max = df_outliers[['high']][[i]]
  
  #we are working with multiple chromosomes so should subset on that first
  chr = df_outliers[['CHR']][[i]]
  #get all snps for the chr of the outlier snp
  chr_snps = df_nonout[df_nonout$CHR == chr,]
  
  #simplified. here we are iterating over the position column for the chr, not full df
  #looks a little cleaner than my initial description because we make the min and max
  #variables above, not in the lapply
  keep_vec = unlist(lapply(chr_snps[['POS']], function(x){
    in_range(x, min, max)
  }))
  
  #subset the chr non-outliers based off the boolean vector we just made
  new_df = chr_snps[keep_vec,]
  #add a new member to the list, name for each list member is the outlier snp,
  #value for each list member is the dataframe corresponding to the snps in the
  #window
  list_of_dfs[[id]] = new_df
}

# list_of_dfs

#double checks
length(list_of_dfs) == nrow(df_outliers)
names(list_of_dfs) == df_outliers[['SNP']]
#look at it to understand the structure



## Need to make a big dataframe from all of the 163 dataframes
## in list_of_dfs

neutral_range_snps = bind_rows(list_of_dfs, .id = "column_label") %>% 
  select(-column_label)

outs = df_outliers %>% 
  select(-low, -high)

outs = outs %>% 
  dplyr::select(CHR, 
                SNP, 
                POS, 
                FST, 
                value...7)

neutral_range_snps = neutral_range_snps %>% 
  dplyr::select(-NMISS)


# bind_rows(neutral_range_snps, 
#                    outs) %>% 
#   arrange(CHR, POS) %>% 
#   write_tsv('SLGBPEL_MORPHO_Outliers_Colocalization_data.txt')


# bind_rows(neutral_range_snps, 
#           outs) %>% 
#   arrange(CHR, POS) %>% 
#   group_by(CHR)

# Cal peak size -----------------------------------------------------------


peak_cal = function(data){
  data %>% 
    summarise(last_pos = last(POS), 
            first_pos = first(POS), 
            avg_Fst = mean(FST)) %>% 
    mutate(peak_dist = last_pos - first_pos) %>% 
    mutate(Mb_peak_dist = peak_dist/1000000) 
}

peak_size = lapply(list_of_dfs, 
       peak_cal)

peak_df = bind_rows(peak_size, 
          .id = "column_label") %>%
  na.omit() %>% 
  summarise(mean_peak_dist = mean(peak_dist), 
            mean_peak_Mb = mean(Mb_peak_dist),
            sd_peak_Mb = sd(Mb_peak_dist),
            mean_Fst_peak = mean(avg_Fst))
peak_df

##  
# Sliding window ----------------------------------------------------------
# Both = read_tsv('GSBPI_BOTHOutliers_Colocalization_data.txt')
# iso = read_tsv('GSBPI_ISOOutliers_Colocalization_data.txt')
# morpho = read_tsv('GSBPI_MORPHOOutliers_Colocalization_data.txt')

Both = read_tsv('SLGBPEL_BOTH_Outliers_Colocalization_data.txt')
iso = read_tsv('SLGBPEL_ISO_Outliers_Colocalization_data.txt')
morpho = read_tsv('SLGBPEL_MORPHO_Outliers_Colocalization_data.txt')


fst_window = winScan(x = morpho,
                     groups = 'CHR',
                     # groups = 'value',
                     position = 'POS',
                     values = 'FST',
                     win_size = 200000,
                     win_step = 199000,
                     funs = c('mean', 'sd'))

fst_window %>%
  as_tibble() %>%
  filter(FST_n >3) %>%
  write_csv('SLGBPEL_MORPHO_Colocalization_Fst_Outliers_200Kb.csv',
            col_names = T)


# fst_window_flag = winScan(x = Both,
#                      groups = c('CHR', 'value'),
#                      # groups = 'value',
#                      position = 'POS',
#                      values = 'FST',
#                      win_size = 200000,
#                      win_step = 1000,
#                      funs = c('mean', 'sd'))
# 
# fst_window_flag %>% 
#   as_tibble() %>%
#   filter(FST_n >3) %>% 
#   write_csv('TSBPL_BOTH_Colocalization_Fst_Outliers_200Kb.csv',
#             col_names = T)
# 
# read_csv('GSBPI_BOTH_Colocalization_Fst_Outliers_200Kb.csv') %>% 
#   filter(value == 'Overlapping locus')


# Fst peaks and data filter -----------------------------------------------

## Read in the data for the outliers
TLGBPL_iso = read_tsv('TLGBPL_ISO_Outliers_Colocalization_data.txt')
TLGBPL_morpho = read_tsv('TLGBPL_MORPHO_Outliers_Colocalization_data.txt')
TLGBPL_both = read_tsv('TLGBPL_BOTH_Outliers_Colocalization_data.txt')

## make labels for each group of snps
iso_outlier = TLGBPL_iso %>% 
  # filter(value == 'Isotope locus')
  filter(value == 'Non-Outlier')

# label = rep('Outlier', 
#             length(iso_outlier$value)) %>% 
#   as_tibble() %>% 
#   rename(label = value)
label = rep('Near Outlier',
            length(iso_outlier$value)) %>%
  as_tibble() %>%
  rename(label = value)
iso_outlier = bind_cols(iso_outlier, 
                        label)


morpho_outlier = TLGBPL_morpho %>% 
  # filter(value == 'Morphological locus')
  filter(value == 'Non-Outlier')
label = rep('Near Outlier', 
            length(morpho_outlier$value)) %>% 
  as_tibble() %>% 
  rename(label = value)
morpho_outlier = bind_cols(morpho_outlier, 
                        label)

both_outlier = TLGBPL_both %>% 
  # filter(value == 'Overlapping locus')
  filter(value == 'Non-Outlier')
label = rep('Near Outlier', 
            length(both_outlier$value)) %>% 
  as_tibble() %>% 
  rename(label = value)
both_outlier = bind_cols(both_outlier, 
                        label)
## write the csv for each population
bind_rows(iso_outlier, 
          morpho_outlier, 
          both_outlier) %>% 
  write_csv('TLGBPL_Near_Outlier_data_21.04.2021.csv')

## Need to read in data and create a population column
# GSBPI = read_csv('GSBPI_Outlier_Colocalization_data_21.04.2021.csv')
GSBPI = read_csv('GSBPI_Near_Outlier_data_21.04.2021.csv')

population = rep('GSBPI', length(GSBPI$label)) %>% 
  as_tibble() %>% 
  rename(Population = value)
GSBPI = bind_cols(GSBPI, 
                  population)

# SLGBPEL = read_csv('SLGBPEL_Outlier_Colocalization_data_21.04.2021.csv')
SLGBPEL = read_csv('SLGBPEL_Near_Outlier_data_21.04.2021.csv')
population = rep('SLGBPEL', length(SLGBPEL$label)) %>% 
  as_tibble() %>% 
  rename(Population = value)
SLGBPEL = bind_cols(SLGBPEL, 
                  population)

# TLGBPL = read_csv('TLGBPL_Outlier_Colocalization_data_21.04.2021.csv')
TLGBPL = read_csv('TLGBPL_Near_Outlier_data_21.04.2021.csv')
population = rep('TLGBPL', length(TLGBPL$label)) %>% 
  as_tibble() %>% 
  rename(Population = value)
TLGBPL = bind_cols(TLGBPL, 
                  population)

# TSBPL = read_csv('TSBPL_Outlier_Colocalization_data_21.04.2021.csv')
TSBPL = read_csv('TSBPL_Near_Outlier_data_21.04.2021.csv')
population = rep('TSBPL', length(TSBPL$label)) %>% 
  as_tibble() %>% 
  rename(Population = value)
TSBPL = bind_cols(TSBPL, 
                  population)
## write the csv
bind_rows(GSBPI, 
          SLGBPEL, 
          TLGBPL, 
          TSBPL) %>% 
  write_csv('AllPops_Near_Outlier_data.21.04.2021.csv')


## Need to find the neutral snps

GSBPI_neutral = read_tsv('GSBPI_Neutral_snps.txt')
label = rep('Neutral',
            length(GSBPI_neutral$value)) %>% 
  as_tibble() %>% 
  rename(label = value)
population = rep('GSBPI', 
                 length(GSBPI_neutral$value)) %>% 
  as_tibble() %>% 
  rename(Population = value)

GSBPI_neutral = bind_cols(GSBPI_neutral, 
          label, 
          population)


SLGBPEL_neutral = read_tsv('SLGBPEL_Neutral_snps.txt')
label = rep('Neutral',
            length(SLGBPEL_neutral$value)) %>% 
  as_tibble() %>% 
  rename(label = value)
population = rep('SLGBPEL', 
                 length(SLGBPEL_neutral$value)) %>% 
  as_tibble() %>% 
  rename(Population = value)
SLGBPEL_neutral = bind_cols(SLGBPEL_neutral, 
                          label, 
                          population)

TLGBPL_neutral = read_tsv('TLGBPL_Neutral_snps.txt')
label = rep('Neutral',
            length(TLGBPL_neutral$value)) %>% 
  as_tibble() %>% 
  rename(label = value)
population = rep('TLGBPL', 
                 length(TLGBPL_neutral$value)) %>% 
  as_tibble() %>% 
  rename(Population = value)
TLGBPL_neutral = bind_cols(TLGBPL_neutral, 
                            label, 
                            population)

TSBPL_neutral = read_tsv('TSBPL_Neutral_snps.txt')
label = rep('Neutral',
            length(TSBPL_neutral$value)) %>% 
  as_tibble() %>% 
  rename(label = value)
population = rep('TSBPL', 
                 length(TSBPL_neutral$value)) %>% 
  as_tibble() %>% 
  rename(Population = value)
TSBPL_neutral = bind_cols(TSBPL_neutral, 
                            label, 
                            population)

bind_rows(GSBPI_neutral, 
          SLGBPEL_neutral, 
          TLGBPL_neutral, 
          TSBPL_neutral) %>% 
  write_csv('Allpops_neutral_snps_21.4.2021.csv')

# Shared and non-shared regions -------------------------------------------
# read_csv('GSBPI_Fst_Outliers_Combined_200Kb.csv')
# read_csv('GSBPI_Fst_Outliers_Combined_200Kb.csv')

GSBPI = read_csv('GSBPI_Fst_200kb_outliers_19.04.2021.csv') 
  # filter(value == 'Outlier', 
  #        FST_n > 0) %>%
  # filter(row_number() %% 101 == 1) %>% 
  # filter(row_number() %% 2 == 1)

# SLGBPL = read_csv('SLGBPL_Colocalization_Fst_Outliers_200Kb.csv') %>% 
SLGBPEL  
filter(value == 'Outlier', 
         FST_n > 0) %>%
  filter(row_number() %% 101 == 1) %>% 
  filter(row_number() %% 2 == 1)

GSBPI = read_csv('.csv') %>% 
  filter(value == 'Outlier', 
         FST_n > 0) %>%
  filter(row_number() %% 101 == 1) %>% 
  filter(row_number() %% 2 == 1)

GSBPI = read_csv('.csv') %>% 
  filter(value == 'Outlier', 
         FST_n > 0) %>%
  filter(row_number() %% 101 == 1) %>% 
  filter(row_number() %% 2 == 1)

GSBPI = read_csv('.csv') %>% 
  filter(value == 'Outlier', 
         FST_n > 0) %>%
  filter(row_number() %% 101 == 1) %>% 
  filter(row_number() %% 2 == 1)

GSBPI = read_csv('.csv') %>% 
  filter(value == 'Outlier', 
         FST_n > 0) %>%
  filter(row_number() %% 101 == 1) %>% 
  filter(row_number() %% 2 == 1)


# GRAPHS!! ----------------------------------------------------------------

data = read_csv('SLGBPEL_Colocalization_Fst_Outliers_200Kb.csv') %>% 
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

data = AC04q.1_29_split(data = data)

# Conversion --------------------------------------------------------------

Mb_Conversion = function(data){
  data %>% 
    group_by(AC_CHR) %>% 
    mutate(win_mid_mb = win_mid/1000000)
}

data = Mb_Conversion(data)

# FST peak colocalization GRAPH -------------------------------------------------------------------
## set the theme for all of the plots
theme_set(theme_bw())

## Drop windows with less than 3 SNPs

plot_fst = data %>%  
  ggplot(aes(x = win_mid_mb,
             y = FST_mean, 
             group = AC_CHR)) +
  geom_point(col = 'grey49')+
  # geom_step(aes(x = win_mid,
  #               y = FST_mean),
  #           direction = 'mid',
  #           size = 2,
  #           col = '#1F2440')+
  geom_smooth(col = '#FF3721',
              # col = '#02E084',
              size = 2)+
  facet_grid(~ AC_CHR, 
             scales = 'free')+
  labs(title = 'VBRSIL Colocalization Fst outliers 200kb window', 
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

ggsave(plot = last_plot(), 
       'VBRSIL_Outlier_Fst_peaks_colocalization.tiff')
# Colocalization table and bar chart ------------------------------------------------
# setwd('~/Fst_sliding_window')

GSBPI_out = read_tsv('GSBPI_Colocalization_data.txt')
value = rep('GSBPI', length(GSBPI_out$SNP)) %>% 
  as_tibble() %>% 
  rename(population = value)
GSBPI_out = bind_cols(GSBPI_out, 
                      value) %>% 
  write_tsv('GSBPI_colocalization_data_fig5.txt')

SLGBPEL_out = read_tsv('SLGBPEL_Colocalization_data.txt')
value = rep('SLGBPEL', length(SLGBPEL_out$SNP)) %>%
  as_tibble() %>%
  rename(population = value)
SLGBPEL_out = bind_cols(SLGBPEL_out,
                      value) %>%
  write_tsv('SLGBPEL_colocalization_data_fig5.txt')


# SLGBPL_out = read_tsv('SLGBPL_Colocalization_data.txt')
# value = rep('SLGBPL', length(SLGBPL_out$SNP)) %>% 
#   as_tibble() %>% 
#   rename(population = value)
# SLGBPL_out = bind_cols(SLGBPL_out, 
#                       value) %>% 
#   write_tsv('SLGBPL_colocalization_data_fig5.txt')
# 
# SLGBPI_out = read_tsv('SLGBPI_Colocalization_data.txt')
# value = rep('SLGBPI', length(SLGBPI_out$SNP)) %>% 
#   as_tibble() %>% 
#   rename(population = value)
# SLGBPI_out = bind_cols(SLGBPI_out, 
#                       value) %>% 
#   write_tsv('SLGBPI_colocalization_data_fig5.txt')

TLGBPL_out = read_tsv('TLGBPL_Colocalization_data.txt')
value = rep('TLGBPL', length(TLGBPL_out$SNP)) %>% 
  as_tibble() %>% 
  rename(population = value)
TLGBPL_out = bind_cols(TLGBPL_out, 
                      value) %>% 
  write_tsv('TLGBPL_colocalization_data_fig5.txt')

TSBPL_out = read_tsv('TSBPL_Colocalization_data.txt')
value = rep('TSBPL', length(TSBPL_out$SNP)) %>% 
  as_tibble() %>% 
  rename(population = value)
TSBPL_out = bind_cols(TSBPL_out, 
                      value) %>% 
  write_tsv('TSBPL_colocalization_data_fig5.txt')

VBRSIL_out = read_tsv('VBRSIL_Colocalization_data.txt')
value = rep('VBRSIL', length(VBRSIL_out$SNP)) %>% 
  as_tibble() %>% 
  rename(population = value)
VBRSIL_out = bind_cols(VBRSIL_out, 
                      value) %>% 
  write_tsv('VBRSIL_colocalization_data_fig5.txt')

GSBPI = read_tsv('GSBPI_colocalization_data_fig5.txt')
SLGBPEL = read_tsv('SLGBPEL_colocalization_data_fig5.txt')
# SLGBPL = read_tsv('SLGBPL_colocalization_data_fig5.txt')
# SLGBPI = read_tsv('SLGBPI_colocalization_data_fig5.txt')
TLGBPL = read_tsv('TLGBPL_colocalization_data_fig5.txt')
TSBPL = read_tsv('TSBPL_colocalization_data_fig5.txt')
VBRSIL = read_tsv('VBRSIL_colocalization_data_fig5.txt')

data = bind_rows(GSBPI, 
                 SLGBPEL,
                     # SLGBPL, 
                     # SLGBPI, 
                     TLGBPL, 
                     TSBPL, 
                     VBRSIL) %>% 
  write_tsv('Polypops_colocalization_data_fig5.txt')

outliers = read_tsv('Polypops_colocalization_data_fig5.txt') %>% 
  filter(value == 'Outlier') %>% 
  group_by(population)

near_outliers = read_tsv('Polypops_colocalization_data_fig5.txt') %>% 
  filter(value == 'Non-Outlier') %>% 
  group_by(population)

GSBPI_neutral = read_tsv('GSBPI_neutral_SNPS.txt')
value = rep('GSBPI', length(GSBPI_neutral$SNP)) %>% 
  as_tibble() %>% 
  rename(population = value)
GSBPI_neutral = bind_cols(GSBPI_neutral, 
                      value) %>% 
  write_tsv('GSBPI_neutral_data_fig5.txt')

SLGBPEL_neutral = read_tsv('SLGBPEL_Neutral_snps.txt')
value = rep('SLGBPEL', length(SLGBPEL_neutral$SNP)) %>% 
  as_tibble() %>% 
  rename(population = value)
SLGBPEL_neutral = bind_cols(SLGBPEL_neutral, 
                           value) %>% 
  write_tsv('SLGBPEL_neutral_data_fig5.txt')

# SLGBPL_neutral = read_tsv('SLGBPL_Neutral_snps.txt')
# value = rep('SLGBPL', length(SLGBPL_neutral$SNP)) %>% 
#   as_tibble() %>% 
#   rename(population = value)
# SLGBPL_neutral = bind_cols(SLGBPL_neutral, 
#                           value) %>% 
#   write_tsv('SLGBPL_neutral_data_fig5.txt')
# 
# SLGBPI_neutral = read_tsv('SLGBPI_Neutral_snps.txt')
# value = rep('SLGBPI', length(SLGBPI_neutral$SNP)) %>% 
#   as_tibble() %>% 
#   rename(population = value)
# SLGBPI_neutral = bind_cols(SLGBPI_neutral, 
#                            value) %>% 
#   write_tsv('SLGBPI_neutral_data_fig5.txt')
TLGBPL_neutral = read_tsv('TLGBPL_Neutral_snps.txt')
value = rep('TLGBPL', length(TLGBPL_neutral$SNP)) %>% 
  as_tibble() %>% 
  rename(population = value)
TLGBPL_neutral = bind_cols(TLGBPL_neutral, 
                           value) %>% 
  write_tsv('TLGBPL_neutral_data_fig5.txt')
TSBPL_neutral = read_tsv('TSBPL_Neutral_SNPS.txt')
value = rep('TSBPL', length(TSBPL_neutral$SNP)) %>% 
  as_tibble() %>% 
  rename(population = value)
TSBPL_neutral = bind_cols(TSBPL_neutral, 
                           value) %>% 
  write_tsv('TSBPL_neutral_data_fig5.txt')
VBRSIL_neutral = read_tsv('VBRSIL_neutral_SNPS.txt')
value = rep('VBRSIL', length(VBRSIL_neutral$SNP)) %>% 
  as_tibble() %>% 
  rename(population = value)
VBRSIL_neutral = bind_cols(VBRSIL_neutral, 
                          value) %>% 
  write_tsv('VBRSIL_neutral_data_fig5.txt')

GSBPI_neutral = read_tsv('GSBPI_neutral_data_fig5.txt')
SLGBPEL_neutral = read_tsv('SLGBPEL_neutral_data_fig5.txt')
# SLGBPL_neutral = read_tsv('SLGBPL_neutral_data_fig5.txt')
# SLGBPI_neutral = read_tsv('SLGBPI_neutral_data_fig5.txt')
TLGBPL_neutral = read_tsv('TLGBPL_neutral_data_fig5.txt')
TSBPL_neutral = read_tsv('TSBPL_neutral_data_fig5.txt')
VBRSIL_neutral = read_tsv('VBRSIL_neutral_data_fig5.txt')

neutral_data = bind_rows(GSBPI_neutral,
                         SLGBPEL_neutral,
                 # SLGBPL_neutral, 
                 # SLGBPI_neutral, 
                 TLGBPL_neutral, 
                 TSBPL_neutral, 
                 VBRSIL_neutral) 

neutral = anti_join(neutral_data, 
                    near_outliers, 

                    by = c('CHR',
                           'SNP',
                           'POS',
                           'NMISS',
                           'FST',
                           'value')) %>% 
  group_by(population) %>% 
  write_tsv('Polypops_neutral_snps_data.txt')

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
    rowid == '6' ~ 'Outlier loci',
    rowid == '7' ~ 'Near outlier loci',
    rowid == '8' ~ 'Near outlier loci',
    rowid == '9' ~ 'Near outlier loci',
    rowid == '10' ~ 'Near outlier loci',
    rowid == '11' ~ 'Near outlier loci',
    rowid == '12' ~ 'Near outlier loci',
    rowid == '13' ~ 'Neutral loci',
    rowid == '14' ~ 'Neutral loci',
    rowid == '15' ~ 'Neutral loci',
    rowid == '16' ~ 'Neutral loci',
    rowid == '17' ~ 'Neutral loci',
    rowid == '18' ~ 'Neutral loci',
    
  ))) 

colours = c('#FF3721', 
            '#E2FF3B', 
            '#5097FF')

class_data %>%
  group_by(population) %>% 
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
  facet_grid(~population, 
             scales = 'free')+
  scale_color_manual(values = colours)+
  scale_fill_manual(values = colours) +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        axis.title.x = element_blank(), 
        legend.position = 'none')+
  labs(y = 'Mean Fst')

ggsave(plot = last_plot(), 
       'VBRSIL_Fst_loci_types.tiff')
# FST outlier colocalization aov -----------------------------------------------------------
setwd('~/Fst_sliding_window/Galtabol_fst/')
setwd('~/Fst_sliding_window/TLGBPL/')
setwd('~/Fst_sliding_window/TSBPL/')
setwd('~/Fst_sliding_window/VBRSIL/')
setwd('~/Fst_sliding_window/SLGBPL/')
setwd('~/Fst_sliding_window/SLGBPI/')

outliers = read_tsv('VBRSIL_Colocalization_data.txt') %>% 
  filter(value == 'Outlier')

near_outliers = read_tsv('VBRSIL_Colocalization_data.txt') %>% 
  filter(value == 'Non-Outlier')

neutral = read_tsv('VBRSIL_Neutral_snps.txt')

neutral = anti_join(neutral, 
                    near_outliers, 
                    
                    by = c('CHR',
                           'SNP',
                           'POS',
                           'NMISS',
                           'FST',
                           'value'))
near_outliers = near_outliers %>% 
  mutate(value = as.factor(case_when(
    value == 'Non-Outlier' ~ 'Near outlier'
  )))

data = bind_rows(outliers, 
                 near_outliers, 
                 neutral)
fst_aov = aov(FST ~ value, data = data)
summary(fst_aov)

pairwise.t.test(data$FST, 
                data$value,
                p.adjust.method = 'bonferroni')




# Differentiation peak range cal ------------------------------------------

# Galtabol peak sizes -----------------------------------------------------

G_Both = read_tsv('GSBPI_BOTHOutliers_Colocalization_data.txt')
G_iso = read_tsv('GSBPI_ISOOutliers_Colocalization_data.txt')
G_morpho = read_tsv('GSBPI_MORPHOOutliers_Colocalization_data.txt')

View(G_Both)

G_peak1 = G_Both %>% 
  slice(1:9)
G_peak1_size = (G_peak1[9,3]-G_peak1[1,3])/1000000

G_peak2 = G_Both %>% 
  slice(10:13)
G_peak2_size = (G_peak2[4,3]-G_peak1[1,3])/1000000

G_peak3 = G_Both %>% 
  slice(14:18)
G_peak3_size = (G_peak3[5,3]-G_peak3[1,3])/1000000

G_peak4 = G_Both %>% 
  slice(19:20)
G_peak4_size = (G_peak4[2,3]-G_peak4[1,3])/1000000

G_peak5 = G_Both %>% 
  slice(21:26)
G_peak5_size = (G_peak5[6,3]-G_peak5[1,3])/1000000

G_peak6 = G_Both %>% 
  slice(28:29)
G_peak6_size = (G_peak6[2,3]-G_peak6[1,3])/1000000

G_peak7 = G_Both %>% 
  slice(30:40)
G_peak7_size = (G_peak7[11,3]-G_peak7[1,3])/1000000

G_peak8 = G_Both %>% 
  slice(41:44)
G_peak8_size = (G_peak8[4,3]-G_peak8[1,3])/1000000

G_peak9 = G_Both %>% 
  slice(45:49)
G_peak9_size = (G_peak9[5,3]-G_peak9[1,3])/1000000

G_peak10 = G_Both %>% 
  slice(50:53)
G_peak10_size = (G_peak10[4,3]-G_peak10[1,3])/1000000

G_peak11 = G_Both %>% 
  slice(54:56)
G_peak11_size = (G_peak11[3,3]-G_peak11[1,3])/1000000

G_peak12 = G_Both %>% 
  slice(57:59)
G_peak12_size = (G_peak12[3,3]-G_peak12[1,3])/1000000

G_peak13 = G_Both %>% 
  slice(61:64)
G_peak13_size = (G_peak13[4,3]-G_peak13[1,3])/1000000

G_peak14 = G_Both %>% 
  slice(65:68)
G_peak14_size = (G_peak14[4,3]-G_peak14[1,3])/1000000

G_peak15 = G_Both %>% 
  slice(70:72)
G_peak15_size = (G_peak15[3,3]-G_peak15[1,3])/1000000

G_peak16 = G_Both %>% 
  slice(73:75)
G_peak16_size = (G_peak16[3,3]-G_peak16[1,3])/1000000

G_peak17 = G_Both %>% 
  slice(76:80)
G_peak17_size = (G_peak17[5,3]-G_peak17[1,3])/1000000

G_peak18 = G_Both %>% 
  slice(81:84)
G_peak18_size = (G_peak18[4,3]-G_peak18[1,3])/1000000

G_peak19 = G_Both %>% 
  slice(85:92)
G_peak19_size = (G_peak19[8,3]-G_peak19[1,3])/1000000

G_peak20 = G_Both %>% 
  slice(93:95)
G_peak20_size = (G_peak20[3,3]-G_peak20[1,3])/1000000

G_peak21 = G_Both %>% 
  slice(96:98)
G_peak21_size = (G_peak21[3,3]-G_peak21[1,3])/1000000

G_peak22 = G_Both %>% 
  slice(99:103)
G_peak22_size = (G_peak22[5,3]-G_peak22[1,3])/1000000

G_peak23 = G_Both %>% 
  slice(104:105)
G_peak23_size = (G_peak23[2,3]-G_peak23[1,3])/1000000

G_peak24 = G_Both %>% 
  slice(106:120)
G_peak24_size = (G_peak24[15,3]-G_peak24[1,3])/1000000

G_peak25 = G_Both %>% 
  slice(121:125)
G_peak25_size = (G_peak25[5,3]-G_peak25[1,3])/1000000

G_peak26 = G_Both %>% 
  slice(126:128)
G_peak26_size = (G_peak26[3,3]-G_peak26[1,3])/1000000

G_peak27 = G_Both %>% 
  slice(129:134)
G_peak27_size = (G_peak27[6,3]-G_peak27[1,3])/1000000

G_peak28 = G_Both %>% 
  slice(135:138)
G_peak28_size = (G_peak28[4,3]-G_peak28[1,3])/1000000

G_peak29 = G_Both %>% 
  slice(139:144)
G_peak29_size = (G_peak29[6,3]-G_peak29[1,3])/1000000

G_peak30 = G_Both %>% 
  slice(140:146)
G_peak30_size = (G_peak30[7,3]-G_peak30[1,3])/1000000

G_peak31 = G_Both %>% 
  slice(140:146)
G_peak31_size = (G_peak31[7,3]-G_peak31[1,3])/1000000
