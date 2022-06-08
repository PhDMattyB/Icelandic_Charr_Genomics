##################################################################
## Phenotypic trajectory analysis
##
## Matt Brachmann (PhDMattyB)
##
## 2022-04-04
##
##################################################################

setwd('~/PhD_Genomics_Chapter3/Fst_Iceland_pops/')

## Load the data manipulation work horse
library(tidyverse)
library(geomorph)
library(RRPP)

PWS_data = read_csv('GSTV_PartialWarpScores_AllometryExcluded.csv') 

Geno_indv = read_csv('GSTVMF_Morph_Eco_Geno.csv') %>% 
  filter(Lake %in% c('Galtabol', 
                     'Svinavatn', 
                     'Thingvallavatn', 
                     'Vatnshlidarvatn')) %>% 
  na.omit() %>% 
  rename(SpecimenID = Specimen.ID)

PWS_data_clean = inner_join(PWS_data, 
           Geno_indv, 
           by = 'SpecimenID') %>% 
  filter(LaMorph.x != 'S.LGB2') %>% 
  rename(Lake = Lake.x, 
         Morph = Morph.x, 
         BP = BP.x, 
         LaMorph = LaMorph.x, 
         Sex = Sex.x)
# View(PWS_data_clean)

# PWS_data_clean = mutate(.data = PWS_data_clean,
#                  BP2 = as.factor(case_when(
#                    LaMorph == 'G.SB' ~ 'Benthic',
#                    LaMorph == 'G.PI' ~ 'Pelagic',
#                    LaMorph == 'S.LGB' ~ 'Benthic',
#                    LaMorph == 'S.PL' ~ 'Pelagic',
#                    LaMorph == 'S.PI' ~ 'Pelagic',
#                    LaMorph == 'T.LGB' ~ 'Benthic',
#                    LaMorph == 'T.SB' ~ 'Benthic2',
#                    LaMorph == 'T.PL' ~ 'Pelagic',
#                    LaMorph == 'V.BR' ~ 'Pelagic',
#                    LaMorph == 'V.SIL' ~ 'Benthic')))


# PWS_data_clean = mutate(.data = PWS_data_clean,
#                         BP2 = as.factor(case_when(
#                           LaMorph == 'G.SB' ~ 'G.Benthic',
#                           LaMorph == 'G.PI' ~ 'G.Pelagic',
#                           LaMorph == 'S.LGB' ~ 'S.Benthic',
#                           LaMorph == 'S.PL' ~ 'S.Pelagic',
#                           LaMorph == 'S.PI' ~ 'S.Pelagic',
#                           LaMorph == 'T.LGB' ~ 'T.Benthic',
#                           LaMorph == 'T.SB' ~ 'T.Benthic2',
#                           LaMorph == 'T.PL' ~ 'T.Pelagic',
#                           LaMorph == 'V.BR' ~ 'V.Pelagic',
#                           LaMorph == 'V.SIL' ~ 'V.Benthic')))

PWS_data_clean = mutate(.data = PWS_data_clean,
                        Vector = as.factor(case_when(
                          LaMorph == 'G.SB' ~ 'GBP',
                          LaMorph == 'G.PI' ~ 'GBP',
                          LaMorph == 'S.LGB' ~ 'SBP',
                          LaMorph == 'S.PL' ~ 'SBP',
                          LaMorph == 'S.PI' ~ 'SBP',
                          LaMorph == 'T.LGB' ~ 'TBP',
                          LaMorph == 'T.SB' ~ 'TBP2',
                          LaMorph == 'T.PL' ~ 'TBP',
                          LaMorph == 'T.PL2' ~ 'TBP2',
                          LaMorph == 'V.BR' ~ 'VBP',
                          LaMorph == 'V.SIL' ~ 'VBP')))
# RRPP --------------------------------------------------------------------


rrpp = rrpp.data.frame(PWS_data_clean)

Bodyshape_PWS = cbind(rrpp$PW19X, 
                      rrpp$PW19Y, 
                      rrpp$PW18X, 
                      rrpp$PW18Y, 
                      rrpp$PW17X,
                      rrpp$PW17Y, 
                      rrpp$PW16X,
                      rrpp$PW16Y, 
                      rrpp$PW15X, 
                      rrpp$PW15Y, 
                      rrpp$PW14X, 
                      rrpp$PW14Y, 
                      rrpp$PW13X, 
                      rrpp$PW13Y,
                      rrpp$PW12X, 
                      rrpp$PW12Y, 
                      rrpp$PW11X, 
                      rrpp$PW11Y, 
                      rrpp$PW10X, 
                      rrpp$PW10Y, 
                      rrpp$PW9X, 
                      rrpp$PW9Y, 
                      rrpp$PW8X, 
                      rrpp$PW8Y, 
                      rrpp$PW7X, 
                      rrpp$PW7Y, 
                      rrpp$PW6X, 
                      rrpp$PW6Y, 
                      rrpp$PW5X, 
                      rrpp$PW5Y, 
                      rrpp$PW4X, 
                      rrpp$PW4Y, 
                      rrpp$PW3X, 
                      rrpp$PW3Y, 
                      rrpp$PW2X, 
                      rrpp$PW2Y, 
                      rrpp$PW1X, 
                      rrpp$PW1Y, 
                      rrpp$UNIX, 
                      rrpp$UNIY, 
                      rrpp$CS)

# Partial warps based trajectory analysis ---------------------------------

library(remotes)
# install_github("fruciano/GeometricMorphometricsMix")
library(GeometricMorphometricsMix)
library(Morpho)

PWS_PCA = prcomp(Bodyshape_PWS)

# BP_cols = c('#1d3557',
#             '#e63946')
# 
# BP_shapes = c(1, 
#               16, 
#               0, 
#               15, 
#               2, 
#               17, 
#               5, 
#               18, 
#               6, 
#               25)
summary(PWS_PCA)

PWS_PCA$x

simple_shapes = c(16, 17)
colours = c('#467BB3', 
            '#FF1E0C', 
            '#108565',
            '#2a9d8f',
            '#D1BA0A')
theme_set(theme_bw())
PCA_BP_Vectors = qplot(PWS_PCA$x[,1], 
      PWS_PCA$x[,2], 
      # col = rrpp$BP, 
      col = rrpp$Vector,
      # shape = rrpp$Vector, 
      shape = rrpp$BP,
      size = 3) +
  coord_equal()+
  # scale_color_manual(values = colours)+
  # scale_shape_manual(values = BP_shapes)+
  scale_color_manual(values = colours)+
  scale_shape_manual(values = simple_shapes)+
  labs(x = 'Benthic-pelagic trajectory 1', 
       y = 'Benthic-pelagic trajectory 2')+
  theme(panel.grid = element_blank(), 
        legend.position = 'none',
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

## take the PC1 and PC2 scores from prcomp
## then test the diffrences between the angles
## then use the function parallel_analysis


PCA_identifiers = PWS_data_clean %>% 
  dplyr::select(SpecimenID, 
                Lake, 
                Morph,
                BP, 
                LaMorph, 
                Vector)
PCA_coords = PWS_PCA$x %>% 
  as_tibble() %>% 
  dplyr::select(PC1, 
                PC2)

# PCA_coords_multi = PWS_PCA$x %>% 
#   as_tibble()

PCA_data = bind_cols(PCA_identifiers, 
                     PCA_coords)


# Test PCA angle differences ----------------------------------------------
Coords_SPLIT_Vector = split(as.data.frame(PCA_coords), 
                            list(PCA_data$Vector), 
                            drop = T)
Coords_Vector = lapply(Coords_SPLIT_Vector, 
                       colMeans)

# Coords_SPLIT_Vector = split(as.data.frame(PCA_coords_multi), 
#                             list(PCA_data$Vector), 
#                             drop = T)
# Coords_Vector = lapply(Coords_SPLIT_Vector, 
#                        colMeans)

GBP = Coords_Vector$GBP
SBP = Coords_Vector$SBP
TBP = Coords_Vector$TBP
TBP2 = Coords_Vector$TBP2
VBP = Coords_Vector$VBP

## Test for differences in vector angles
TestOfAngle(GBP, 
            SBP, 
            flip = F)

TestOfAngle(GBP, 
            TBP, 
            flip = F)

TestOfAngle(GBP, 
            TBP2, 
            flip = F)

TestOfAngle(GBP, 
            VBP, 
            flip = F)

TestOfAngle(SBP, 
            TBP, 
            flip = F)

TestOfAngle(SBP, 
            TBP2, 
            flip = F)
TestOfAngle(SBP, 
            VBP, 
            flip = F)

TestOfAngle(TBP, 
            TBP2, 
            flip = F)

TestOfAngle(TBP, 
            VBP, 
            flip = F)

TestOfAngle(TBP2, 
            VBP, 
            flip = F)


# Magnitude diffs PCA vectors ---------------------------------------------

Mag_GBP_SBP = dist_mean_boot(A = Coords_SPLIT_Vector$GBP, 
                             B = Coords_SPLIT_Vector$SBP)

Mag_GBP_TBP = dist_mean_boot(A = Coords_SPLIT_Vector$GBP, 
                             B = Coords_SPLIT_Vector$TBP)

Mag_GBP_TBP2 = dist_mean_boot(A = Coords_SPLIT_Vector$GBP, 
                              B = Coords_SPLIT_Vector$TBP2)

Mag_GBP_VBP = dist_mean_boot(A = Coords_SPLIT_Vector$GBP, 
                             B = Coords_SPLIT_Vector$VBP)

Mag_SBP_TBP = dist_mean_boot(A = Coords_SPLIT_Vector$SBP, 
                             B = Coords_SPLIT_Vector$TBP)

Mag_SBP_TBP2 = dist_mean_boot(A = Coords_SPLIT_Vector$SBP, 
                              B = Coords_SPLIT_Vector$TBP2)

Mag_SBP_VBP = dist_mean_boot(A = Coords_SPLIT_Vector$SBP, 
                             B = Coords_SPLIT_Vector$VBP)

Mag_TBP_TBP2 = dist_mean_boot(A = Coords_SPLIT_Vector$TBP, 
                              B = Coords_SPLIT_Vector$TBP2)

Mag_TBP_VBP = dist_mean_boot(A = Coords_SPLIT_Vector$TBP, 
                             B = Coords_SPLIT_Vector$VBP)

Mag_TBP2_VBP = dist_mean_boot(A = Coords_SPLIT_Vector$TBP2, 
                              B = Coords_SPLIT_Vector$VBP)



# combine graphs ----------------------------------------------------------
library(patchwork)
PCA_BP_Vectors
PCA_BodyShape

Phenotype_analyses = PCA_BodyShape|PCA_BP_Vectors

ggsave('Phenotype_trajectory_analysis_FINAL.tiff',
       plot = Phenotype_analyses, 
       dpi = 'retina', 
       unit = 'cm', 
       width = 25, 
       height = 20
)
