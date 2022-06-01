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
setwd('~/Fst_sliding_window/Galtabol_fst/')
setwd('~/Fst_sliding_window/TLGBPL/')
setwd('~/Fst_sliding_window/TSBPL/')
setwd('~/Fst_sliding_window/VBRSIL/')
setwd('~/Fst_sliding_window/SLGBPL/')
setwd('~/Fst_sliding_window/SLGBPI/')
setwd('~/Fst_sliding_window/')

setwd('~/PhD_Genomics_Chapter3/Fst_Iceland_pops/')

# Data manipulation -------------------------------------------------------
data = read_tsv('GSBPI_fst.fst')

data = read_tsv('Galtabol_fst_NA.fst')

outliers = read_csv('GSBPI_Outlier_data2.csv') %>% 
  rename(SNP = `Marker ID`) %>% 
  select(SNP, 
         `#Chromosome`,
         `Genetic distance`,
         `Physical position`, 
         CHROME3, 
         value) %>% 
  arrange(SNP)

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
          'GSBPI_Neutral_snps.txt', col_names = TRUE)

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

grouping = rep('Both Outliers', 
               length(BothTraits_Fst$FST)) %>% 
  as_tibble()

BothTraits_Fst = bind_cols(BothTraits_Fst, 
                           grouping)
grouping = rep('Iso outlier', 
               length(Iso_fst$FST)) %>% 
  as_tibble()
Iso_fst = bind_cols(Iso_fst, 
                    grouping)

grouping = rep('Morpho outlier', 
               length(morpho_fst$FST)) %>% 
  as_tibble()
morpho_fst = bind_cols(morpho_fst, 
                       grouping)
# 
# grouping = rep('Outlier',
#                length(Fst_outliers$FST)) %>% 
#   as_tibble()
# 
# Fst_outliers = bind_cols(Fst_outliers, 
#                          grouping)

## write the three files as csvs
write_tsv(BothTraits_Fst, 
          'GSBPI_BOTHoutliers_Fst.txt')

write_tsv(Iso_fst, 
          'GSBPI_ISOoutliers_Fst.txt')

write_tsv(morpho_fst, 
          'GSBPI_MORPHOoutliers_Fst.txt')

## read in the neutral data
df_nonout = read_tsv('GSBPI_Neutral_snps.txt')

##read in the outlier data
# df_outliers = read_tsv('GSBPI_outliers_Fst.txt')

overlap_outs = read_tsv('GSBPI_BOTHoutliers_Fst.txt') %>% 
  filter(CHR != 0)
iso_outs = read_tsv('GSBPI_ISOoutliers_Fst.txt') %>% 
  filter(CHR != 0)
morpho_outs = read_tsv('GSBPI_MORPHOoutliers_Fst.txt') %>% 
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
list_of_dfs
names(list_of_dfs) #outlier IDS
list_of_dfs[1]	#df with outlier label
list_of_dfs[[1]]#just the df



## Need to make a big dataframe from all of the 163 dataframes
## in list_of_dfs

neutral_range_snps = bind_rows(list_of_dfs, .id = "column_label") %>% 
  select(-column_label)

outs = df_outliers %>% 
  select(-low, -high)

bind_rows(neutral_range_snps, 
                   outs) %>% 
  arrange(CHR, POS) %>% 
  write_tsv('GSBPI_BOTHOutliers_Colocalization_data.txt')

##  
# Sliding window ----------------------------------------------------------
Both = read_tsv('GSBPI_BOTHOutliers_Colocalization_data.txt')
iso = read_tsv('GSBPI_ISOOutliers_Colocalization_data.txt')
morpho = read_tsv('GSBPI_MORPHOOutliers_Colocalization_data.txt')


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
  write_csv('GSBPI_MORPHO_Colocalization_Fst_Outliers_200Kb.csv',
            col_names = T)

read_csv('TSBPL_MORPHO_Colocalization_Fst_Outliers_200Kb.csv') 

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



# Shared and non-shared regions -------------------------------------------
# read_csv('GSBPI_Fst_Outliers_Combined_200Kb.csv')
# read_csv('GSBPI_Fst_Outliers_Combined_200Kb.csv')

GSBPI = read_csv('GSBPI_Colocalization_Fst_Outliers_200Kb.csv') 
  # filter(value == 'Outlier', 
  #        FST_n > 0) %>%
  # filter(row_number() %% 101 == 1) %>% 
  # filter(row_number() %% 2 == 1)

SLGBPL = read_csv('SLGBPL_Colocalization_Fst_Outliers_200Kb.csv') %>% 
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
setwd('~/Fst_sliding_window')

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

