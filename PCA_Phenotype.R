##################################################################
## Body shape PCA
##
## Matt Brachmann (PhDMattyB)
##
## 2022-04-04
##
##################################################################

## Load the data manipulation work horse
library(tidyverse)


setwd('~/PhD_Genomics_Chapter3/Fst_Iceland_pops/')

# Phenotype = read_csv('RDA_Morphological_variables_PCA.csv')
# Phenotype = read_csv('RDA_Morphological_Variables_AllPops.csv')
Phenotype = read_csv('GSTVMF_Morph_Eco_Geno.csv') %>% 
  filter(LaMorph != 'S.LGB2', 
         LaMorph != 'T.PL2')
# 
Phenotype %>%
  select(LaMorph) %>%
  distinct()

Phenotype_mid = mutate(.data = Phenotype,
                         Morph2 = as.factor(case_when(
                           LaMorph == 'G.SB' ~ 'Small Benthic',
                           LaMorph == 'G.PI' ~ 'Piscivorous',
                           LaMorph == 'S.LGB1' ~ 'Large Benthic',
                           LaMorph == 'S.PL' ~ 'Planktivorous',
                           LaMorph == 'S.PI' ~ 'Planktivorous',
                           LaMorph == 'T.LGB' ~ 'Large Benthic',
                           LaMorph == 'T.SB' ~ 'Small Benthic',
                           LaMorph == 'T.PL1' ~ 'Planktivorous',
                           LaMorph == 'V.BR' ~ 'Brown',
                           LaMorph == 'V.SIL' ~ 'Silver')))


Phenotype_clean = Phenotype_mid %>% 
  filter(Lake %in% c('Galtabol', 
                     'Svinavatn', 
                     'Thingvallavatn', 
                     'Vatnshlidarvatn')) %>% 
  dplyr::select(1:4, 
         Morph2,
         LaMorph, 
         Fork.length, 
         CO_NA_PCA_1, 
         CO_NA_PCA_2, 
         NA_PCA_1, 
         NA_PCA_2)

# Phenotype_clean %>% 
#   filter(Lake == 'Svinavatn') %>% 
#   View()
# Phenotype_clean$BP2

PCA_colours = c('#467BB3', 
                '#FF1E0C', 
                '#108565', 
                '#D1BA0A')

theme_set(theme_bw())

Phenotype_mean = Phenotype_clean %>% 
  group_by(Lake, 
           Morph2,
           BP) %>% 
  summarise(mean_pc1 = mean(CO_NA_PCA_1), 
            mean_pc2 = mean(CO_NA_PCA_2)) 


PCA_BodyShape = ggplot(Phenotype_clean, 
       aes(x = CO_NA_PCA_1, 
           y = CO_NA_PCA_2))+
  geom_point(aes(col = Lake, 
                 shape = BP))+
  ## Galtabol mean line segment
  geom_segment(x = 0.0325, 
               y = 0.0003, 
               xend = 0.015, 
               yend = 0.0006, 
               col = 'black', 
               size = 2)+
  ## Svinvatn mean line segment
  geom_segment(x = -0.030, 
               y = -0.003, 
               xend = -0.031, 
               yend = 0.009, 
               col = 'black', 
               size = 2)+
  ## Thing1 mean line segment
  geom_segment(x = 0.038, 
               y = 0.006, 
               xend = 0.017, 
               yend = 0.021, 
               col = 'black', 
               size = 2)+
  ## Thing2 mean line segment
  geom_segment(x = 0.051, 
               y = 0.002, 
               xend = 0.017, 
               yend = 0.021, 
               col = 'black', 
               size = 2)+
  ## Vatn mean line segment
  geom_segment(x = -0.005, 
               y = -0.019, 
               xend = -0.001, 
               yend = -0.023, 
               col = 'black', 
               size = 2)+
  geom_point(data = Phenotype_mean, 
             aes(x = mean_pc1, 
                 y = mean_pc2, 
                 shape = BP), 
             col = 'black', 
             size = 5)+
  geom_point(data = Phenotype_mean, 
             aes(x = mean_pc1, 
                 y = mean_pc2,
                 col = Lake,
                 shape = BP),
             size = 4)+
  scale_color_manual(values = PCA_colours)+
  labs(x = 'Principal component 1 (38.8%)', 
       y = 'Principal component 2 (17.0%)')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

ggsave('PCA_BodyShape_NoAllometry.tiff',
       plot = PCA_BodyShape, 
       dpi = 'retina', 
       unit = 'cm', 
       width = 15, 
       height = 15
       )
