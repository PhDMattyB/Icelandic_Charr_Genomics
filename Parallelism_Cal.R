
##################################################################
## Genomic parallelism calculation
##
## Matt Brachmann (PhDMattyB)
##
## 2022-04-04
##
##################################################################

## Load the data manipulation work horse
library(tidyverse)

## Load other packages 


## Read in the data


# GSBPI Data_Cleaning -----------------------------------------------------------


gsbpi_overlap = read_csv('GSBPI_Outlier_data.csv') %>% 
  filter(value == 'Overlapping locus') %>% 
  select(-`#Chromosome`, 
         -`Genetic distance`, 
         -`Physical position`, 
         -CHROME3, 
         -AXIS, 
         -LOADING) %>% 
  select(`Marker ID`,
         CHROME3.x,
         `Genetic distance.x`, 
         `Physical position.x`, 
         value) %>% 
  rename(MarkerID = 1, 
         AC_CHR = 2, 
         Genetic_dist = 3, 
         Physical_pos = 4, 
         Label = 5)
 
gsbpi_iso = read_csv('GSBPI_Outlier_data.csv') %>% 
  filter(value == 'Isotope locus') %>% 
  select(`Marker ID`,          
         CHROME3,
         `Genetic distance`, 
         `Physical position`, 
         value) %>% 
  rename(MarkerID = 1, 
         AC_CHR = 2, 
         Genetic_dist = 3, 
         Physical_pos = 4, 
         Label = 5)

gsbpi_morpho = read_csv('GSBPI_Outlier_data.csv') %>% 
  filter(value == 'Morphological locus') %>% 
  select(`Marker ID`,          
         CHROME3,
         `Genetic distance`, 
         `Physical position`, 
         value) %>% 
  rename(MarkerID = 1, 
         AC_CHR = 2, 
         Genetic_dist = 3, 
         Physical_pos = 4, 
         Label = 5)


bind_rows(gsbpi_overlap, 
                  gsbpi_iso, 
                  gsbpi_morpho) %>% 
  write_csv('GSBPI_RDA_outlier_clean.csv')

# SLGBPEL Data Cleaning ---------------------------------------------------

SLGBPL_overlap = read_csv('SLGBPEL_Outlier_data.csv') %>% 
  filter(value == 'Overlapping locus') %>% 
  # select(-`#Chromosome`, 
  #        -`Genetic distance`, 
  #        -`Physical position`, 
  #        -CHROME3, 
  #        -AXIS, 
  #        -LOADING) %>% 
  select(`Marker ID`,
         CHROME3,
         `Genetic distance`, 
         `Physical position`, 
         value) %>% 
  rename(MarkerID = 1, 
         AC_CHR = 2, 
         Genetic_dist = 3, 
         Physical_pos = 4, 
         Label = 5)

SLGBPL_iso = read_csv('SLGBPEL_Outlier_data.csv') %>% 
  filter(value == 'Isotope locus') %>% 
  select(`Marker ID`,          
         CHROME3,
         `Genetic distance`, 
         `Physical position`, 
         value) %>% 
  rename(MarkerID = 1, 
         AC_CHR = 2, 
         Genetic_dist = 3, 
         Physical_pos = 4, 
         Label = 5)

SLGBPL_morpho = read_csv('SLGBPEL_Outlier_data.csv') %>% 
  filter(value == 'Morphological locus') %>% 
  select(`Marker ID`,          
         CHROME3,
         `Genetic distance`, 
         `Physical position`, 
         value) %>% 
  rename(MarkerID = 1, 
         AC_CHR = 2, 
         Genetic_dist = 3, 
         Physical_pos = 4, 
         Label = 5)


bind_rows(SLGBPL_overlap, 
          SLGBPL_iso, 
          SLGBPL_morpho) %>% 
  write_csv('SLGBPL_RDA_outlier_clean.csv')

# TLGBPL Data cleaning ----------------------------------------------------

TLGBPL_overlap = read_csv('TLGBPL_Outlier_data.csv') %>% 
  filter(value == 'Overlapping locus') %>% 
  select(`Marker ID`,
         CHROME3,
         `Genetic distance`, 
         `Physical position`, 
         value) %>% 
  rename(MarkerID = 1, 
         AC_CHR = 2, 
         Genetic_dist = 3, 
         Physical_pos = 4, 
         Label = 5)

TLGBPL_iso = read_csv('TLGBPL_Outlier_data.csv') %>% 
  filter(value == 'Isotope locus') %>% 
  select(`Marker ID`,          
         CHROME3,
         `Genetic distance`, 
         `Physical position`, 
         value) %>% 
  rename(MarkerID = 1, 
         AC_CHR = 2, 
         Genetic_dist = 3, 
         Physical_pos = 4, 
         Label = 5)

TLGBPL_morpho = read_csv('TLGBPL_Outlier_data.csv') %>% 
  filter(value == 'Morphological locus') %>% 
  select(`Marker ID`,          
         CHROME3,
         `Genetic distance`, 
         `Physical position`, 
         value) %>% 
  rename(MarkerID = 1, 
         AC_CHR = 2, 
         Genetic_dist = 3, 
         Physical_pos = 4, 
         Label = 5)


bind_rows(TLGBPL_overlap, 
          TLGBPL_iso, 
          TLGBPL_morpho) %>% 
  write_csv('TLGBPL_RDA_outlier_clean.csv')



# TSBPL Data Cleaning -----------------------------------------------------

TSBPL_overlap = read_csv('TSBPL_Outlier_data.csv') %>% 
  filter(value == 'Overlapping locus') %>% 
  select(`Marker ID`,
         CHROME3,
         `Genetic distance`, 
         `Physical position`, 
         value) %>% 
  rename(MarkerID = 1, 
         AC_CHR = 2, 
         Genetic_dist = 3, 
         Physical_pos = 4, 
         Label = 5)

TSBPL_iso = read_csv('TSBPL_Outlier_data.csv') %>% 
  filter(value == 'Isotope locus') %>% 
  select(`Marker ID`,          
         CHROME3,
         `Genetic distance`, 
         `Physical position`, 
         value) %>% 
  rename(MarkerID = 1, 
         AC_CHR = 2, 
         Genetic_dist = 3, 
         Physical_pos = 4, 
         Label = 5)

TSBPL_morpho = read_csv('TSBPL_Outlier_data.csv') %>% 
  filter(value == 'Morphological locus') %>% 
  select(`Marker ID`,          
         CHROME3,
         `Genetic distance`, 
         `Physical position`, 
         value) %>% 
  rename(MarkerID = 1, 
         AC_CHR = 2, 
         Genetic_dist = 3, 
         Physical_pos = 4, 
         Label = 5)


bind_rows(TSBPL_overlap, 
          TSBPL_iso, 
          TSBPL_morpho) %>% 
  write_csv('TSBPL_RDA_outlier_clean.csv')
# Data --------------------------------------------------------------------

GSBPI = read_csv('GSBPI_RDA_outlier_clean.csv')%>% 
  filter(AC_CHR == 'AC11')
SLGBPL = read_csv('SLGBPL_RDA_outlier_clean.csv')%>% 
  filter(AC_CHR == 'AC11')
TLGBPL = read_csv('TLGBPL_RDA_outlier_clean.csv')%>% 
  filter(AC_CHR == 'AC11')
TSBPL = read_csv('TSBPL_RDA_outlier_clean.csv')%>% 
  filter(AC_CHR == 'AC11')

# GSBPI_SLGBPL = (13/427 + 13/329)*0.5
# GSBPI_TLGBPL = (13/427 + 13/310)*0.5
# GSBPI_TSBPL = (13/427 + 13/312)*0.5
# SLGBPL_TLGBPL = (13/329 + 13/310)*0.5
# SLGBPL_TSBPL = (13/329 + 13/312)*0.5
# TLGBPL_TSBPL = (13/310 + 13/310)*0.5

GSBPI_SLGBPL = inner_join(GSBPI, 
           SLGBPL)

# GSBPI_SLGBPL %>% 
#   # filter(Label == "Overlapping locus")
#   # filter(Label == "Isotope locus")
#   filter(Label == "Morphological locus")

GS_Parallelism_Cal = (56/427 + 56/329)*0.5


GSBPI_TLGBPL = inner_join(GSBPI, 
                          TLGBPL)

# GSBPI_TLGBPL %>% 
#   # filter(Label == "Overlapping locus")
#   # filter(Label == "Isotope locus")
#   filter(Label == "Morphological locus")
GT1_Parallism_Cal = (57/427 + 57/310)*0.5

GSBPI_TSBPL = inner_join(GSBPI, 
                         TSBPL)

GT2_Parallelism_Cal = (64/427 + 64/312)*0.5


SLGBPL_TLGBPL = inner_join(SLGBPL, 
                           TLGBPL)
ST1_parallelism_cal = (65/329 + 65/310)*0.5

SLGBPL_TSBPL = inner_join(SLGBPL, 
                           TSBPL)
ST2_parallelism_cal = (51/329 + 51/312)*0.5

TLGBPL_TSBPL = inner_join(TLGBPL, 
                          TSBPL)
T_parallelism_cal = (150/310 + 150/312)*0.5



# GSBPI 200Kb RDA Outlier Data cleaning -----------------------------------

GSBPI_BOTH_200kb = read_csv('GSBPI_BOTH_Colocalization_Fst_Outliers_200Kb.csv')
GSBPI_ISO_200kb = read_csv('GSBPI_ISO_Colocalization_Fst_Outliers_200Kb.csv')
GSBPI_MORPHO_200kb = read_csv('GSBPI_MORPHO_Colocalization_Fst_Outliers_200Kb.csv')

label_both = rep('Overlapping locus',
            nrow(GSBPI_BOTH_200kb)) %>% 
  as_tibble()
GSBPI_BOTH_200kb = bind_cols(GSBPI_BOTH_200kb, 
                             label_both) %>% 
  Chr_convert3() %>% 
  Mb_Conversion()

label_iso = rep('Isotope locus',
                 nrow(GSBPI_ISO_200kb)) %>% 
  as_tibble()
GSBPI_ISO_200kb = bind_cols(GSBPI_ISO_200kb, 
                             label_iso) %>% 
  Chr_convert3() %>% 
  Mb_Conversion()

label_morpho = rep('Morphological locus',
                 nrow(GSBPI_MORPHO_200kb)) %>% 
  as_tibble()
GSBPI_MORPHO_200kb = bind_cols(GSBPI_MORPHO_200kb, 
                             label_morpho) %>% 
  Chr_convert3() %>% 
  Mb_Conversion()


GSBPI_200kb_RDA = bind_rows(GSBPI_BOTH_200kb, 
                            GSBPI_ISO_200kb, 
                            GSBPI_MORPHO_200kb)

# SLGBPEL 200Kb RDA Outlier Data cleaning -----------------------------------

SLGBPEL_BOTH_200kb = read_csv('SLGBPEL_BOTH_Colocalization_Fst_Outliers_200Kb.csv')
SLGBPEL_ISO_200kb = read_csv('SLGBPEL_ISO_Colocalization_Fst_Outliers_200Kb.csv')
SLGBPEL_MORPHO_200kb = read_csv('SLGBPEL_MORPHO_Colocalization_Fst_Outliers_200Kb.csv')

label_both = rep('Overlapping locus',
                 nrow(SLGBPEL_BOTH_200kb)) %>% 
  as_tibble()
SLGBPEL_BOTH_200kb = bind_cols(SLGBPEL_BOTH_200kb, 
                             label_both) %>% 
  Chr_convert3() %>% 
  Mb_Conversion()

label_iso = rep('Isotope locus',
                nrow(SLGBPEL_ISO_200kb)) %>% 
  as_tibble()
SLGBPEL_ISO_200kb = bind_cols(SLGBPEL_ISO_200kb, 
                            label_iso) %>% 
  Chr_convert3() %>% 
  Mb_Conversion()

label_morpho = rep('Morphological locus',
                   nrow(SLGBPEL_MORPHO_200kb)) %>% 
  as_tibble()
SLGBPEL_MORPHO_200kb = bind_cols(SLGBPEL_MORPHO_200kb, 
                               label_morpho) %>% 
  Chr_convert3() %>% 
  Mb_Conversion()


SLGBPEL_200kb_RDA = bind_rows(SLGBPEL_BOTH_200kb, 
                            SLGBPEL_ISO_200kb, 
                            SLGBPEL_MORPHO_200kb)

# TLGBPL 200Kb RDA Outlier Data cleaning -----------------------------------

TLGBPL_BOTH_200kb = read_csv('TLGBPL_BOTH_Colocalization_Fst_Outliers_200Kb.csv')
TLGBPL_ISO_200kb = read_csv('TLGBPL_ISO_Colocalization_Fst_Outliers_200Kb.csv')
TLGBPL_MORPHO_200kb = read_csv('TLGBPL_MORPHO_Colocalization_Fst_Outliers_200Kb.csv')

label_both = rep('Overlapping locus',
                 nrow(TLGBPL_BOTH_200kb)) %>% 
  as_tibble()
TLGBPL_BOTH_200kb = bind_cols(TLGBPL_BOTH_200kb, 
                               label_both) %>% 
  Chr_convert3() %>% 
  Mb_Conversion()

label_iso = rep('Isotope locus',
                nrow(TLGBPL_ISO_200kb)) %>% 
  as_tibble()
TLGBPL_ISO_200kb = bind_cols(TLGBPL_ISO_200kb, 
                              label_iso) %>% 
  Chr_convert3() %>% 
  Mb_Conversion()

label_morpho = rep('Morphological locus',
                   nrow(TLGBPL_MORPHO_200kb)) %>% 
  as_tibble()
TLGBPL_MORPHO_200kb = bind_cols(TLGBPL_MORPHO_200kb, 
                                 label_morpho) %>% 
  Chr_convert3() %>% 
  Mb_Conversion()


TLGBPL_200kb_RDA = bind_rows(TLGBPL_BOTH_200kb, 
                              TLGBPL_ISO_200kb, 
                              TLGBPL_MORPHO_200kb)

# TSBPL 200Kb RDA Outlier Data cleaning -----------------------------------

TSBPL_BOTH_200kb = read_csv('TSBPL_BOTH_Colocalization_Fst_Outliers_200Kb.csv')
TSBPL_ISO_200kb = read_csv('TSBPL_ISO_Colocalization_Fst_Outliers_200Kb.csv')
TSBPL_MORPHO_200kb = read_csv('TSBPL_MORPHO_Colocalization_Fst_Outliers_200Kb.csv')

label_both = rep('Overlapping locus',
                 nrow(TSBPL_BOTH_200kb)) %>% 
  as_tibble()
TSBPL_BOTH_200kb = bind_cols(TSBPL_BOTH_200kb, 
                               label_both) %>% 
  Chr_convert3() %>% 
  Mb_Conversion()

label_iso = rep('Isotope locus',
                nrow(TSBPL_ISO_200kb)) %>% 
  as_tibble()
TSBPL_ISO_200kb = bind_cols(TSBPL_ISO_200kb, 
                              label_iso) %>% 
  Chr_convert3() %>% 
  Mb_Conversion()

label_morpho = rep('Morphological locus',
                   nrow(TSBPL_MORPHO_200kb)) %>% 
  as_tibble()
TSBPL_MORPHO_200kb = bind_cols(TSBPL_MORPHO_200kb, 
                                 label_morpho) %>% 
  Chr_convert3() %>% 
  Mb_Conversion()


TSBPL_200kb_RDA = bind_rows(TSBPL_BOTH_200kb, 
                              TSBPL_ISO_200kb, 
                              TSBPL_MORPHO_200kb)

# 200Kb parallelism -------------------------------------------------------

GSBPI_200kb = GSBPI_200kb_RDA %>% 
  group_by(AC_CHR) %>% 
  filter(AC_CHR == 'AC11')
SLGBPEL_200kb = SLGBPEL_200kb_RDA %>% 
  group_by(AC_CHR)%>% 
  filter(AC_CHR == 'AC11')
TLGBPL_200kb = TLGBPL_200kb_RDA %>% 
  group_by(AC_CHR)%>% 
  filter(AC_CHR == 'AC11')
TSBPL_200kb = TSBPL_200kb_RDA %>% 
  group_by(AC_CHR)%>% 
  filter(AC_CHR == 'AC11')


inner_join(GSBPI_200kb, 
           SLGBPEL_200kb)

inner_join(GSBPI_200kb, 
           TLGBPL_200kb)

inner_join(GSBPI_200kb, 
           TSBPL_200kb)

inner_join(SLGBPEL_200kb, 
           TLGBPL_200kb)

inner_join(SLGBPEL_200kb, 
           TSBPL_200kb)

inner_join(TLGBPL_200kb, 
           TSBPL_200kb)


GS_region_parallelism = (2/84 + 2/57)*0.5
GT1_region_parallelism = (1/84 + 1/58)*0.5
GT2_region_parallelism = (0/84 + 0/53)*0.5
ST1_region_parallelism = (2/57 + 2/58)*0.5
ST2_region_parallelism = (0/57 + 0/53)*0.5
T_region_parallelism = (4/58 + 4/53)*0.5
