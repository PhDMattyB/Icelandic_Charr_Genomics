##############################
## Overlap between morphology and isotope associated outliers
##
## Matt Brachmann (PhDMattyB)
##
## 2020-06-10
##
##############################

setwd('~/PhD/SNP Demographic modelling/Outliers_directory/RDA_Outlier_comparison')

library(tidyverse)

outlier_labels = function(Iso, Morpho){
  overlap = inner_join(Iso, 
                       Morpho, 
                       by = c('Marker ID', 
                              '#Chromosome', 
                              'Genetic distance', 
                              'Physical position', 
                              'CHROME3'))
  label1 = rep('Overlapping locus', 
               length(overlap$`Marker ID`)) %>% 
    as_tibble()
  overlap = bind_cols(overlap, 
                      label1)
  
  iso_only = anti_join(Iso, 
                       overlap, 
                       by = c('Marker ID', 
                              '#Chromosome', 
                              'Genetic distance', 
                              'Physical position', 
                              'CHROME3'))
  label2 = rep('Isotope locus', 
               length(iso_only$`Marker ID`)) %>% 
    as_tibble()
  iso_only = bind_cols(iso_only, label2)
  
  morpho_only = anti_join(Morpho,
                          overlap, 
                          by = c('Marker ID', 
                                 '#Chromosome', 
                                 'Genetic distance', 
                                 'Physical position', 
                                 'CHROME3'))
  label3 = rep('Morphological locus', 
               length(morpho_only$`Marker ID`)) %>% 
    as_tibble()
  morpho_only = bind_cols(morpho_only, 
                          label3)
  
  total = bind_rows(overlap, 
                    iso_only, 
                    morpho_only)
  return(total)
  
}


#
# test run fix issues -----------------------------------------------------
Iso = read_csv('GSBPI_RDA_outliers_mapped.csv')
# poly_iso = read_csv('Polypop_RDA_outliers_mapped_SPELCOMBO.csv')

Morpho = read_csv('GSBPI_RDA_morphology_outliers_mapped.csv')
# poly_morpho = read_csv('Polypop_RDA_morphology_outliers_total_mapped_SPELCOMBO.csv')


overlap = inner_join(Iso, 
                     Morpho, 
                     by = c('Marker ID', 
                            '#Chromosome', 
                            'Genetic distance', 
                            'Physical position', 
                            'CHROME3'))

View(overlap)

label1 = rep('Overlapping locus', 
             length(overlap$`Marker ID`)) %>% 
  as_tibble()
overlap = bind_cols(overlap, 
                    label1)

iso_only = anti_join(Iso, 
                     overlap, 
                     by = 'Marker ID')
label2 = rep('Isotope locus', 
             length(iso_only$`Marker ID`)) %>% 
  as_tibble()
iso_only = bind_cols(iso_only, label2)

morpho_only = anti_join(Morpho,
                        overlap, 
                        by = 'Marker ID')
label3 = rep('Morphological locus', 
             length(morpho_only$`Marker ID`)) %>% 
  as_tibble()
morpho_only = bind_cols(morpho_only, 
                        label3)

total = bind_rows(overlap, 
                  iso_only, 
                  morpho_only)
# Galtabol outliers  --------------------------------------------------------

g_iso = read_csv('GSBPI_RDA_outliers_mapped.csv')
# poly_iso = read_csv('Polypop_RDA_outliers_mapped_SPELCOMBO.csv')

g_morpho = read_csv('GSBPI_RDA_morphology_outliers_mapped.csv')
# poly_morpho = read_csv('Polypop_RDA_morphology_outliers_total_mapped_SPELCOMBO.csv')

g_data = outlier_labels(g_iso, g_morpho)

##check
g_data$value

# View(g_data)
# write_csv(g_data, 'GSBPI_Outlier_data2.csv')
# 
# 
# read_csv('GSBPI_Outlier_data2.csv')
# Svinavatn outliers ------------------------------------------------------

s_iso = read_csv('SLGBPEL_RDA_outliers_mapped.csv')
s_morpho = read_csv('SLGBPEL_RDA_morphology_outliers_total_mapped.csv')

s_data = outlier_labels(s_iso, s_morpho)

s_data$value

# write_csv(s_data, 
#           'SLGBPEL_Outlier_data.csv')

# TLGBPL outliers ---------------------------------------------------------

T1_iso = read_csv('Thingvallavatn_LGBPL_RDA_outliers_mapped.csv')
T1_morpho = read_csv('TLGBPL_RDA_morphology_outliers_mapped.csv')

T1_data = outlier_labels(T1_iso, 
                         T1_morpho)
# write_csv(T1_data, 
#           'TLGBPL_Outlier_data.csv')
# TSBPL outliers ---------------------------------------------------------

T2_iso = read_csv('Thingvallavatn_PLSB_RDA_outliers_mapped.csv')
T2_morpho = read_csv('TSBPL_RDA_morphology_outliers_mapped.csv')

T2_data = outlier_labels(T2_iso, 
                         T2_morpho)
# write_csv(T2_data, 
#           'TSBPL_Outlier_data.csv')

# poly outliers ---------------------------------------------------------

poly_iso = read_csv('Polypop_RDA_outliers_mapped_SPELCOMBO.csv')
poly_morpho = read_csv('Polypop_RDA_morphology_outliers_total_mapped_SPELCOMBO.csv')

poly_data = outlier_labels(poly_iso, 
                         poly_morpho)
write_csv(poly_data, 
          'Polypop_Outlier_data.csv')


# Outlier overlap populations ---------------------------------------------
g_label = rep('G: Benthic - Pelagic', 
    length(g_data$`Marker ID`)) %>% 
  as_tibble()

g_data = bind_cols(g_data, g_label) %>% 
  rename(Morph_pair = value1) %>% 
  select(`Marker ID`, 
         `#Chromosome`, 
         `Genetic distance`, 
         `Physical position`, 
         CHROME3, 
         value, 
         Morph_pair)

s_label = rep('S: Benthic - Pelagic', 
              length(s_data$`Marker ID`)) %>% 
  as_tibble()

s_data = bind_cols(s_data, s_label)%>% 
  rename(Morph_pair = value1)%>% 
  select(`Marker ID`, 
         `#Chromosome`, 
         `Genetic distance`, 
         `Physical position`, 
         CHROME3, 
         value, 
         Morph_pair)

T1_label = rep('T: Benthic 1 - Pelagic', 
              length(T1_data$`Marker ID`)) %>% 
  as_tibble()

T1_data = bind_cols(T1_data, T1_label)%>% 
  rename(Morph_pair = value1)%>% 
  select(`Marker ID`, 
         `#Chromosome`, 
         `Genetic distance`, 
         `Physical position`, 
         CHROME3, 
         value, 
         Morph_pair)

T2_label = rep('T: Benthic 2 - Pelagic', 
              length(T2_data$`Marker ID`)) %>% 
  as_tibble()

T2_data = bind_cols(T2_data, T2_label)%>% 
  rename(Morph_pair = value1)%>% 
  select(`Marker ID`, 
         `#Chromosome`, 
         `Genetic distance`, 
         `Physical position`, 
         CHROME3, 
         value, 
         Morph_pair)


# big_data = bind_rows(g_data, 
#                      s_data, 
#                      T1_data, 
#                      T2_data)

overlap = inner_join(g_data, 
                     s_data, 
                     by = c('Marker ID', 
                            '#Chromosome',
                            'Genetic distance', 
                            'Physical position', 
                            'CHROME3'))

# View(overlap)
overlap = inner_join(overlap, 
                     T1_data, 
                     by = c('Marker ID', 
                            '#Chromosome',
                            'Genetic distance', 
                            'Physical position', 
                            'CHROME3'))
overlap = inner_join(overlap, 
                     T2_data, 
                     by = c('Marker ID', 
                            '#Chromosome',
                            'Genetic distance', 
                            'Physical position', 
                            'CHROME3')) 

overlap %>% 
  select(CHROME3) %>% 
  filter(CHROME3 != 'Contigs')

overlap = overlap %>% 
  filter(CHROME3 != 'Contigs') %>% 
  select(`Marker ID`, 
         `Physical position`, 
         value.x, 
         value.y, 
         value.x.x, 
         value.y.y, 
         CHROME3) %>% 
  arrange(CHROME3)
## number overlapping in traits between all: 6/17
## 11/17 not overlapping in at least one morph pair
write_tsv(overlap, 'Overlapping_SNPs.txt')

# Graph Commmon outliers --------------------------------------------------
theme_set(theme_bw())
library(patchwork)

overlap %>% 
  filter(CHROME3 != 'Contigs') %>% 
  View()

outlier_color3 = c('#F2B90F',
                  '#D92B04',
                  '#1C6C8C')
GSBPI = overlap %>% 
  filter(CHROME3 != 'Contigs') %>% 
  ggplot()+
  geom_point(aes(x = CHROME3, 
                  y = value.x,
                  col = value.x), 
             size = 4, 
             alpha = 0.6)+
  scale_color_manual(values = outlier_color3)+
  labs(title = 'A)', 
       x = 'Chromosome')+
  guides(col = guide_legend(title = 'Outlier'))+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.ticks = element_line(size = 1), 
        panel.grid = element_blank(), 
        legend.position = 'none')

SLGBPEL = overlap %>% 
  filter(CHROME3 != 'Contigs') %>% 
  ggplot()+
  geom_point(aes(x = CHROME3, 
                 y = value.y,
                 col = value.y), 
             size = 4, 
             alpha = 0.6)+
  scale_color_manual(values = outlier_color3)+
  labs(title = 'B)', 
       x = 'Chromosome')+
  guides(col = guide_legend(title = 'Outlier'))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.ticks = element_line(size = 1), 
        panel.grid = element_blank(), 
        legend.position = 'none')

outlier_color2 = c('#F2B90F',
                   '#1C6C8C')

TLGBPL = overlap %>% 
  filter(CHROME3 != 'Contigs') %>% 
  ggplot()+
  geom_point(aes(x = CHROME3, 
                 y = value.x.x,
                 col = value.x.x), 
             size = 4, 
             alpha = 0.6)+
  scale_color_manual(values = outlier_color2)+
  labs(title = 'C)', 
       x = 'Chromosome')+
  guides(col = guide_legend(title = 'Outlier'))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.ticks = element_line(size = 1), 
        panel.grid = element_blank(), 
        legend.position = 'none')

TSBPL = overlap %>% 
  filter(CHROME3 != 'Contigs') %>% 
  ggplot()+
  geom_point(aes(x = CHROME3, 
                 y = value.y.y,
                 col = value.y.y), 
             size = 4, 
             alpha = 0.6)+
  scale_color_manual(values = outlier_color3)+
  labs(title = 'D)', 
       x = 'Chromosome')+
  guides(col = guide_legend(title = 'Outlier'))+
  theme(axis.title.y = element_blank(), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.ticks = element_line(size = 1), 
        panel.grid = element_blank(), 
        legend.position = 'none')

GSBPI/SLGBPEL/TLGBPL/TSBPL

ggsave('Overlapping_RDA_outlier_loci_alpha.tiff', 
       plot = last_plot())
