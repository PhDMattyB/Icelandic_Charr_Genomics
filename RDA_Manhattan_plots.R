##############################
## RDA Manhattan plots
##
## Matt Brachmann (PhDMattyB)
##
## 2020-02-05
##
##############################

setwd('~/PhD/SNP Demographic modelling/Outliers_directory')

library(patchwork)
library(rsed)
library(vegan)
library(tidyverse)
library(ggman)


theme_set(theme_bw())

# isotopic rda ------------------------------------------------------------
## SNP data
RAW_iso = read.delim("Polypop_Genotype_arranged_A.raw",
                 sep = "",
                 stringsAsFactors = F)
RAW_iso = read.delim("PolyPops_Genotype_SPELCOMBO_recodeA.raw",
                     sep = "",
                     stringsAsFactors = F)

##Environmental data; isotopic data
ENVIRONMENT = read_csv('LFMM_Enviro_SPELCOMBO.csv') %>% 
  # filter(`#FamilyID` %in% c('9', '10')) %>% 
  dplyr::select(-Sex) 
## Arrange the environmental data by familty ID
ENV_ARRANGED = ENVIRONMENT %>%  
  rename(IID = IndividualID, 
         FID = `#FamilyID`) %>% 
  arrange(IID)
## arrange the snp data by family id
RAW_CHARR_TIB = as_tibble(RAW_iso) %>% 
  arrange(IID)

## Make sure the individual ID's match up across the two data sets
ENV_ARRANGED$IID == RAW_CHARR_TIB$IID

## just the snp matrix
SNPS_iso = RAW_CHARR_TIB %>%
  dplyr::select(matches("Affx.")) %>% 
  as_tibble(RAW_CHARR_TIB)

## nas are 0. this is a consequence of the data transformation
SNPS_iso[is.na(SNPS_iso)] = 0
##data from of environmental data
CandN = ENV_ARRANGED %>% 
  dplyr::select(Carbon, Nitrogen) %>% 
  as.data.frame()

## run the rda
CHARR_RDA_iso = rda(SNPS_iso ~ Carbon + Nitrogen, data = ENV_ARRANGED,
                scale = T)

summary(eigenvals(CHARR_RDA_iso, 
                  model = "constrained"))

# AllSNPS = read_csv('Polypop_RDA_AllSNPS.csv')
KeyStone_iso = read_csv('RDA_IsoMorpho_Outliers_SNPList.csv') %>% 
  as.data.frame()
# poly_pop_iso = read_csv('Polypop_RDA_outliers_mapped.csv') %>% 
#   rename(ID = `Marker ID`)

poly_pop_iso = read_csv('Polypop_RDA_outliers_mapped_SPELCOMBO.csv') %>%
  rename(ID = `Marker ID`)


poly_outliers_iso = inner_join(poly_pop_iso, 
                           KeyStone_iso, 
                           by = 'ID')
# poly_RDA_values = read_csv('PolyPop_RDA_values.csv')

snp_scores_iso = read_csv('Poly_RDA_SNP_scores_SPELCOMBO.csv')
# test_iso = sed_substitute(snp_scores_iso$SNP, 
#                'Affx.', 
#                'Affx-')
test_iso = gsub("Affx.", 
                "Affx-", 
                snp_scores_iso$SNP)

test_iso = gsub("_.*",
                "",
                test_iso)
# test_iso = sed_substitute(test_iso, 
#                '_.*',
#                "")
test_iso = as.data.frame(test_iso) %>% 
  as_tibble() %>% 
  rename(ID = test_iso)
snp_scores_iso = bind_cols(snp_scores_iso, 
                       test_iso)
cand_scores_iso = inner_join(poly_outliers_iso, 
                         snp_scores_iso, 
                         by = 'ID')

cand_scores_iso = left_join(poly_pop_iso, 
                            snp_scores_iso, 
                            by = 'ID') %>% 
  filter(CHROME3 != 'Contigs')

plot_bp = read_csv('Polypop_RDA_BP_SPELCOMBO.csv')
biplot = read_csv('Poly_RDA_biplot.csv')
plot_snps_iso = read_csv('Polypop_RDA_SNPS_SPELCOMBO.csv')

ENV_iso = read_csv('Polypop_environmental_data.csv')

ENV_iso = mutate(.data = ENV_iso,
              `Morph type2` = as.factor(case_when(
                POP_NAME == "T.LGB" ~ "T: Benthic 1",
                POP_NAME == 'V.BR' ~ 'V: Pelagic',
                POP_NAME == 'V.SIL' ~ 'V: Benthic',
                POP_NAME == 'S.LGB' ~ 'S: Benthic',
                POP_NAME == 'S.PL'~ 'S: Pelagic 1',
                POP_NAME == 'S.PI' ~ 'S: Pelagic 2',
                POP_NAME == 'T.PL' ~ 'T: Pelagic',
                POP_NAME == 'T.SB' ~ 'T: Benthic 2',
                POP_NAME == 'G.SB' ~ 'G: Benthic',
                POP_NAME == 'G.PI' ~ 'G: Pelagic'))) 

ENV_iso = mutate(.data = ENV_iso,
                 `Morph type` = as.factor(case_when(
                   POP_NAME == "T.LGB" ~ "Benthic",
                   POP_NAME == 'V.BR' ~ 'Pelagic',
                   POP_NAME == 'V.SIL' ~ 'Benthic',
                   POP_NAME == 'S.LGB' ~ 'Benthic',
                   POP_NAME == 'S.PL'~ 'Pelagic',
                   POP_NAME == 'S.PI' ~ 'Pelagic',
                   POP_NAME == 'T.PL' ~ 'Pelagic',
                   POP_NAME == 'T.SB' ~ 'Benthic 2',
                   POP_NAME == 'G.SB' ~ 'Benthic',
                   POP_NAME == 'G.PI' ~ 'Pelagic'))) 

ENV_iso = mutate(.data = ENV_iso,
                 Population = as.factor(case_when(
                   POP_NAME == "T.LGB" ~ "Thingvallavatn",
                   POP_NAME == 'V.BR' ~ 'Vatnshlidarvatn',
                   POP_NAME == 'V.SIL' ~ 'Vatnshlidarvatn',
                   POP_NAME == 'S.LGB' ~ 'Svinavatn',
                   POP_NAME == 'S.PL'~ 'Svinavatn',
                   POP_NAME == 'S.PI' ~ 'Svinavatn',
                   POP_NAME == 'T.PL' ~ 'Thingvallavatn',
                   POP_NAME == 'T.SB' ~ 'Thingvallavatn',
                   POP_NAME == 'G.SB' ~ 'Galtabol',
                   POP_NAME == 'G.PI' ~ 'Galtabol'))) 

# Combing Svinavatn Pelagic morphs ----------------------------------------

poly_pop_iso = read_csv('Polypop_RDA_outliers_mapped_SPELCOMBO.csv') %>%
  rename(ID = `Marker ID`)

poly_pop_iso$ID


# poly_outliers_iso = inner_join(poly_pop_iso, 
#                                KeyStone_iso, 
#                                by = 'ID')
# poly_RDA_values = read_csv('PolyPop_RDA_values.csv')

snp_scores_iso = read_csv('Poly_RDA_SNP_scores_SPELCOMBO.csv') %>% 
  rename(ID = SNP)

# test_iso = gsub("Affx.", 
#                 "Affx-", 
#                 snp_scores_iso$SNP)

# test_iso = sed_substitute(snp_scores_iso$SNP, 
#                'Affx.', 
#                'Affx-')
# snp_scores_iso$ID = substr(snp_scores_iso$ID,
#                         1,
#                         regexpr("\\_", snp_scores$ID)-1)

snp_scores_iso$ID = gsub("Affx.",
                         "Affx-", 
                         snp_scores_iso$ID)
snp_scores_iso$ID = gsub("_.*",
                         "",
                         snp_scores_iso$ID)

# test_iso = sed_substitute(test_iso, 
#                '_.*',
#                "")
# test_iso = as.data.frame(test_iso) %>% 
#   as_tibble() %>% 
#   rename(ID = test_iso)
# snp_scores_iso = bind_cols(snp_scores_iso, 
#                            test_iso)
cand_scores_iso = inner_join(poly_pop_iso, 
                             snp_scores_iso, 
                             by = 'ID') %>% 
  filter(CHROME3 != 'Contigs')

cand_scores_iso = left_join(poly_pop_iso, 
                            snp_scores_iso, 
                            by = 'ID') %>% 
  filter(CHROME3 != 'Contigs')

# plot_bp = read_csv('Polypop_RDA_BP.csv')
# biplot = read_csv('Poly_RDA_biplot.csv')
plot_snps_iso = read_csv('Polypop_RDA_SNPS_SPELCOMBO.csv')

ENV_iso = read_csv('Polypop_environmental_data.csv')

ENV_iso = mutate(.data = ENV_iso,
                 `Morph type2` = as.factor(case_when(
                   POP_NAME == "T.LGB" ~ "T: Benthic 1",
                   POP_NAME == 'V.BR' ~ 'V: Pelagic',
                   POP_NAME == 'V.SIL' ~ 'V: Benthic',
                   POP_NAME == 'S.LGB' ~ 'S: Benthic',
                   POP_NAME == 'S.PL'~ 'S: Pelagic',
                   POP_NAME == 'S.PI' ~ 'S: Pelagic',
                   POP_NAME == 'T.PL' ~ 'T: Pelagic',
                   POP_NAME == 'T.SB' ~ 'T: Benthic 2',
                   POP_NAME == 'G.SB' ~ 'G: Benthic',
                   POP_NAME == 'G.PI' ~ 'G: Pelagic'))) 

ENV_iso = mutate(.data = ENV_iso,
                 `Morph type` = as.factor(case_when(
                   POP_NAME == "T.LGB" ~ "Benthic",
                   POP_NAME == 'V.BR' ~ 'Pelagic',
                   POP_NAME == 'V.SIL' ~ 'Benthic',
                   POP_NAME == 'S.LGB' ~ 'Benthic',
                   POP_NAME == 'S.PL'~ 'Pelagic',
                   POP_NAME == 'S.PI' ~ 'Pelagic',
                   POP_NAME == 'T.PL' ~ 'Pelagic',
                   POP_NAME == 'T.SB' ~ 'Benthic 2',
                   POP_NAME == 'G.SB' ~ 'Benthic',
                   POP_NAME == 'G.PI' ~ 'Pelagic'))) 

ENV_iso = mutate(.data = ENV_iso,
                 Population = as.factor(case_when(
                   POP_NAME == "T.LGB" ~ "Thingvallavatn",
                   POP_NAME == 'V.BR' ~ 'Vatnshlidarvatn',
                   POP_NAME == 'V.SIL' ~ 'Vatnshlidarvatn',
                   POP_NAME == 'S.LGB' ~ 'Svinavatn',
                   POP_NAME == 'S.PL'~ 'Svinavatn',
                   POP_NAME == 'S.PI' ~ 'Svinavatn',
                   POP_NAME == 'T.PL' ~ 'Thingvallavatn',
                   POP_NAME == 'T.SB' ~ 'Thingvallavatn',
                   POP_NAME == 'G.SB' ~ 'Galtabol',
                   POP_NAME == 'G.PI' ~ 'Galtabol'))) 

# RDA Isotope plot --------------------------------------------------------

## benthic are in blue and pelagics are in red
RDA_colors = c('#467BB3', #GSB,
               '#FF1E0C', #GPI
               '#4440B3', #slgb
               # '#FF640D', #spi 
               '#FF8E06', #spl 
               '#7D8AFF', #tlgb
               '#238CB3', #tsb
               '#FF3361', #tpl
               '#12AEB3', #VSIL
               '#FF5338') #vbr
PCA_colors = c('#467BB3',
               '#FF1E0C',
               '#108565',
               '#D1BA0A')

rda_iso = ggplot() +
  # geom_point(aes(x = CHARR_RDA$CCA$v[,1], 
  #                y = CHARR_RDA$CCA$v[,2]), 
  #            col = "gray86") +
  
  ## this geom point is for the zoomed out version
  geom_point(data = plot_snps_iso,
             aes(x = RDA1,
                 y = RDA2),
             col = '#ADA597')+
  
  ## this geom point is for the zoomed in version
  # geom_point(data = plot_snps,
  #            aes(x = RDA1,
  #                y = RDA2),
  #            col = '#C7B9A3')+
  geom_point(data = cand_scores_iso, 
             aes(x = RDA1, 
                 y = RDA2), 
             col = '#0AB33A', 
             size = 3)+
  ## These two things are only for the zoomed out rda to look
  ## at population level variation related to isotopic variation
  geom_point(data = ENV_iso,
             aes(x = RDA1,
                 y = RDA2,
                 col = Population,
                 shape = `Morph type`))+
  scale_colour_manual(values = PCA_colors)+
  
  ## THis is for the zoomed in RDA to show snp variation
  # geom_segment(aes(xend = CHARR_RDA$CCA$biplot[,1]/10,
  #                  yend = CHARR_RDA$CCA$biplot[,2]/10,
  #                  x=0,
  #                  y=0),
  #              colour="black",
  #              size=0.5,
  #              linetype=1,
  #              arrow = arrow(length = unit(0.02, "npc"))) +
  # geom_text(aes(x = 1.5*CHARR_RDA$CCA$biplot[,1]/10,
  #               y = 1.2*CHARR_RDA$CCA$biplot[,2]/10,
  #               label = colnames(ENV_ARRANGED[,3:4]))) +
  
  ## THis is for the zoomed out RDA to show var in pops
  geom_segment(aes(xend = CHARR_RDA_iso$CCA$biplot[,1],
                   yend = CHARR_RDA_iso$CCA$biplot[,2],
                   x=0,
                   y=0),
               colour="black",
               size=0.5,
               linetype=1,
               arrow = arrow(length = unit(0.02, "npc"))) +
  geom_text(aes(x = 1.5*CHARR_RDA_iso$CCA$biplot[,1],
                y = 1.2*CHARR_RDA_iso$CCA$biplot[,2],
                label = colnames(ENV_ARRANGED[,3:4]))) +
  
  labs(x = 'RDA 1 (81.7% variance explained)', 
       y = 'RDA 2 (18.3% variance explained)', 
       title = 'A)')+
  theme(#legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 15), 
        axis.ticks = element_line(size = 1), 
        plot.title = element_text(size = 15, 
                                  hjust = 0), 
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12))

rda_iso

ggsave(plot = last_plot(), 
      'RDA_Isotope_PolyPops_Full_SPELCOMBO.tiff')

ggsave(plot = last_plot(), 
       'RDA_Isotope_PolyPops_SNPS.tiff')



# Morphology RDA ----------------------------------------------------------

Genotype = read.delim("PolyPops_Genotype_SPELCOMBO_recodeA.raw",
                 sep = "",
                 stringsAsFactors = F)
Genotype = as_tibble(Genotype) %>% 
  arrange(IID)

Morphology = read_csv('RDA_Morphological_variables_SPELCOMBO.csv') %>% 
  # filter(POP_name %in% c('T.SB', 'T.PL')) %>% 
  arrange(IndividualID) %>% 
  rename(FID = `#FamilyID`, 
         IID = IndividualID) %>% 
  distinct()


Morphology = mutate(.data = Morphology,
                    Population_full = as.factor(case_when(
                      FID == "1" ~ "Thingvallavatn - Large benthic",
                      FID == '2' ~ 'Vatnshlidarvatn - Brown',
                      FID == '3' ~ 'Vatnshlidarvatn - Silver',
                      FID == '4'~ 'Svinavatn - Pelagic',
                      # FID == '5' ~ 'Svinavatn - Pelagic',
                      FID == '6' ~ 'Thingvallavatn - Planktivorous',
                      FID == '7' ~ 'Thingvallavatn - Small benthic',
                      FID == '8' ~ 'Svinavatn - Large benthic',
                      FID == '9' ~ 'Galtabol - Small benthic',
                      FID == '10' ~ 'Galtabol - Piscivorous')))

Genotype$IID == Morphology$IID

## need a genotype matrix so we will need all snps, 
## all snps start with Affx. 
SNPS = Genotype %>% 
  select(matches("Affx.")) %>% 
  as_tibble(Genotype)

SNPS[is.na(SNPS)] = 0

data_only = Morphology %>% 
  select(FL, 
         PC1, 
         PC2, 
         PC3) %>% 
  as.data.frame()

## Run the actual RDA. Remember an RDA is a PCA of a linear model. 
CHARR_RDA = rda(SNPS ~ FL + PC1 + PC2 + PC3,
                data = data_only,
                scale = T)

## Calculate the rsquared for the model. 
## This calculates the effect size that the model actually explains
## Report the adjusted Rsquared value
# RsquareAdj(CHARR_RDA)
# 
# summary(eigenvals(CHARR_RDA, model = "constrained"))

sum = summary(CHARR_RDA)

read_csv('Polypop_RDA_Morphology_outliers_SPELCOMBO.csv')


poly_morpho = read_csv('Polypop_RDA_morphology_outliers_total_mapped_SPELCOMBO.csv') %>% 
  rename(ID = `Marker ID`)
poly_outliers = inner_join(poly_morpho, 
                           KeyStone_cand_list, 
                           by = 'ID')

poly_outliers = poly_morpho
# poly_morpho = read_csv('Polypop_RDA_morphology_outliers_total_mapped.csv') %>% 
#   rename(ID = `Marker ID`)
# 
# poly_outliers = inner_join(poly_morpho, 
#                            KeyStone_cand_list, 
#                            by = 'ID')

## The significant snp on CHromosome 24 might be a false positive!!

# poly_RDA_values = read_csv('PolyPop_RDA_values.csv')

snp_scores = read_csv('Polypop_RDA_morpho_scores_total_SPELCOMBO.csv') %>% 
  rename(ID = SNP)

test = gsub("Affx.", 
            "Affx-", 
              snp_scores$ID)

test = gsub("_.*", 
            "",
            test)
# test = sed_substitute(snp_scores$SNP, 
#                         'Affx.', 
#                         'Affx-')
# test = sed_substitute(test, 
#                       '_.*',
#                       "")
test = as.data.frame(test) %>% 
  as_tibble() %>% 
  rename(ID = test)
snp_scores = bind_cols(snp_scores, 
                       test) %>% 
  dplyr::select(-ID) %>% 
  rename(ID = ID1)
cand_scores = inner_join(poly_outliers, 
                         snp_scores, 
                         by = 'ID')

cand_scores = left_join(poly_morpho, 
                        snp_scores, 
                        by = 'ID')%>% 
  filter(CHROME3 != 'Contigs')

cand_scores$RDA1

# RDA Morphology plot -----------------------------------------------------
plot_SNPS = read_csv('PolypopsRDa_SNPS_Morphology_SPELCOMBO.csv')
plot_loc = read_csv('PolypopsRDA_LOC_Morphology_SPELCOMBO.csv')
plot_bp = read_csv('PolypopsRDA_BP_Morphology_SPELCOMBO.csv')

Morphology_data = bind_cols(Morphology, 
                            plot_loc)

Morphology_data = mutate(.data = Morphology_data,
             `Morph type2` = as.factor(case_when(
               POP_name == "T.LGB" ~ "T: Benthic 1",
               POP_name == 'V.BR' ~ 'V: Pelagic',
               POP_name == 'V.SIL' ~ 'V: Benthic',
               POP_name == 'S.LGB' ~ 'S: Benthic',
               POP_name == 'S.PEL'~ 'S: Pelagic',
               # POP_name == 'S.PI' ~ 'S: Pelagic 2',
               POP_name == 'T.PL' ~ 'T: Pelagic',
               POP_name == 'T.SB' ~ 'T: Benthic 2',
               POP_name == 'G.SB' ~ 'G: Benthic',
               POP_name == 'G.PI' ~ 'G: Pelagic'))) 

Morphology_data = mutate(.data = Morphology_data,
                 `Morph type` = as.factor(case_when(
                   POP_name == "T.LGB" ~ "Benthic",
                   POP_name == 'V.BR' ~ 'Pelagic',
                   POP_name == 'V.SIL' ~ 'Benthic',
                   POP_name == 'S.LGB' ~ 'Benthic',
                   POP_name == 'S.PEL'~ 'Pelagic',
                   # POP_name == 'S.PI' ~ 'Pelagic',
                   POP_name == 'T.PL' ~ 'Pelagic',
                   POP_name == 'T.SB' ~ 'Benthic 2',
                   POP_name == 'G.SB' ~ 'Benthic',
                   POP_name == 'G.PI' ~ 'Pelagic'))) 

Morphology_data = mutate(.data = Morphology_data,
                 Population = as.factor(case_when(
                   POP_name == "T.LGB" ~ "Thingvallavatn",
                   POP_name == 'V.BR' ~ 'Vatnshlidarvatn',
                   POP_name == 'V.SIL' ~ 'Vatnshlidarvatn',
                   POP_name == 'S.LGB' ~ 'Svinavatn',
                   POP_name == 'S.PEL'~ 'Svinavatn',
                   # POP_name == 'S.PI' ~ 'Svinavatn',
                   POP_name == 'T.PL' ~ 'Thingvallavatn',
                   POP_name == 'T.SB' ~ 'Thingvallavatn',
                   POP_name == 'G.SB' ~ 'Galtabol',
                   POP_name == 'G.PI' ~ 'Galtabol'))) 

KeyStone_cand_list = read_csv('RDA_IsoMorpho_Outliers_SNPList.csv') %>% 
  as.data.frame()

## benthic are in blue and pelagics are in red
RDA_colors = c('#467BB3', #GSB,
               '#FF1E0C', #GPI
               '#4440B3', #slgb
               # '#FF640D', #spi 
               '#FF8E06', #spl 
               '#7D8AFF', #tlgb
               '#238CB3', #tsb
               '#FF3361', #tpl
               '#12AEB3', #VSIL
               '#FF5338') #vbr
PCA_colors

rda_morpho = ggplot() +
  # geom_point(aes(x = CHARR_RDA$CCA$v[,1], 
  #                y = CHARR_RDA$CCA$v[,2]), 
  #            col = "gray86") +
  
  ## this geom point is for the zoomed out version
  geom_point(data = plot_SNPS,
             aes(x = RDA1,
                 y = RDA2),
             col = '#ADA597')+

  ## this geom point is for the zoomed in version
  # geom_point(data = plot_snps,
  #            aes(x = RDA1,
  #                y = RDA2),
  #            col = '#C7B9A3')+
  geom_point(data = cand_scores,
             aes(x = RDA1, 
                 y = RDA2), 
             col = '#0AB33A', 
             size = 3)+
  ## These two things are only for the zoomed out rda to look
  ## at population level variation related to isotopic variation
  geom_point(data = Morphology_data,
             aes(x = RDA1,
                 y = RDA2,
                 col = Population,
                 shape = `Morph type`))+
  scale_colour_manual(values = PCA_colors)+
  
  ## THis is for the zoomed in RDA to show snp variation
# geom_segment(aes(xend = CHARR_RDA$CCA$biplot[,1]/10,
#                  yend = CHARR_RDA$CCA$biplot[,2]/10,
#                  x=0,
#                  y=0),
#              colour="black",
#              size=0.5,
#              linetype=1,
#              arrow = arrow(length = unit(0.02, "npc"))) +
# geom_text(aes(x = 1.5*CHARR_RDA$CCA$biplot[,1]/10,
#               y = 1.2*CHARR_RDA$CCA$biplot[,2]/10,
#               label = colnames(Morphology_data[,3:6]))) +

## THis is for the zoomed out RDA to show var in pops
geom_segment(aes(xend = CHARR_RDA$CCA$biplot[,1],
                 yend = CHARR_RDA$CCA$biplot[,2],
                 x=0,
                 y=0),
             colour="black",
             size=0.5,
             linetype=1,
             arrow = arrow(length = unit(0.02, "npc"))) +
  geom_text(aes(x = 2.0*CHARR_RDA$CCA$biplot[,1],
                y = 2.0*CHARR_RDA$CCA$biplot[,2],
                label = colnames(Morphology[,3:6]))) +
  
  labs(x = 'RDA 1 (57.3% variance explained)', 
       y = 'RDA 2 (22.2% variance explained)',
       title = 'B)')+
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 15), 
        axis.ticks = element_line(size = 1),
        plot.title = element_text(size = 15, 
                                  hjust = 0))

rda_morpho

ggsave(plot = last_plot(), 
       'RDA_Morphology_PolyPops_Full_SPELCOMBO.tiff')

ggsave(plot = last_plot(), 
       'RDA_Morphology_PolyPops_SNPS.tiff')

rda_iso/rda_morpho

ggsave(plot = last_plot(), 
       'RDA_Iso_morpho_combo_PolyPops_NoContigs_SPELCOMBO.tiff', 
       height = 9)
# RDA manhattan plot isotopes ---------------------------------------------

rda_iso_total = read_csv('Polypop_RDA_AllSNPS.csv')

rda_man_plot = ggman(rda_iso_total, 
                     chrom = 'CHROME3', 
                     pvalue = 'RDA1',
                     snp = 'MARKER_ID',
                     bp = 'PDIST',
                     logTransform = F,
                     pointSize = 2, 
                     title = 'C)',
                     xlabel = 'Chromosome',
                     ymin = -0.4, 
                     ymax = 0.4, 
                     lineColour = 'black', 
                     relative.positions = TRUE)

out_highlight = as.character(poly_outliers$ID)
light_up = ggmanHighlight(rda_man_plot, 
                          highlight = out_highlight, 
                          colour = 'Red')
label = ggmanLabel(light_up, 
                   labelDfm = KeyStone_cand_list, 
                   snp = 'ID', 
                   label = 'ID',
                   colour = 'Red')

plot = label + 
  # scale_color_manual(values = c('#FF0E0E',
  #                               '#2470FF'))+
  theme(axis.text.x =
          element_text(size  = 15,
                       angle = 90,
                       hjust = 1,
                       vjust = 1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15), 
        axis.ticks = element_line(size = 1), 
        plot.title = element_text(size = 15))

plot

ggsave(plot = last_plot(), 
       'RDA_iso_total_outliers_manplot_13.02.2020.tiff', 
       width = 8)


ggmanZoom(plot, chromosome = 'AC13')
ggmanZoom(plot, chromosome = 'AC24')
ggmanZoom(plot, chromosome = 'AC25')
ggmanZoom(plot, chromosome = 'AC33')



# rda manhattan plot morphology -------------------------------------------

rda_morpho_total = read_csv('Polypop_RDA_Morphology_allsnps__total_mapped.csv')
## morphology rda results
rda_man_plot = ggman(rda_morpho_total, 
                     chrom = 'CHROME3', 
                     pvalue = 'RDA1',
                     snp = 'MARKER_ID',
                     bp = 'PDIST',
                     logTransform = F,
                     pointSize = 2, 
                     title = 'D)',
                     xlabel = 'Chromosome',
                     ymin = -0.4, 
                     ymax = 0.4, 
                     lineColour = 'black', 
                     relative.positions = TRUE)

out_highlight = as.character(poly_outliers$ID)
light_up = ggmanHighlight(rda_man_plot, 
                          highlight = out_highlight, 
                          colour = 'Red')

lab = poly_outliers$ID %>% 
  as_tibble() %>% 
  rename(ID = value) %>% 
  as.data.frame()
label = ggmanLabel(light_up, 
                   labelDfm = lab, 
                   snp = 'ID', 
                   label = 'ID',
                   colour = 'Red')

plot2 = label + 
  # scale_color_manual(values = c('#FF0E0E',
  #                               '#2470FF'))+
  theme(axis.text.x =
          element_text(size  = 15,
                       angle = 90,
                       hjust = 1,
                       vjust = 1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15), 
        axis.ticks = element_line(size = 1), 
        plot.title = element_text(size = 15))

plot2

ggsave(plot = last_plot(), 
       'RDA_morpho_manhattan_plot_13.02.2020.tiff')

plot/plot2

ggsave(plot = last_plot(), 
       'RDA_genomic_position_13.02.2020.tiff', 
       height = 8, 
       width = 9)

plot3 = (rda_iso/rda_morpho)+(plot/plot2)

ggsave(plot = last_plot(), 
       'Fig4_RDA_genomics_13.02.2020.tiff')