##############################
## RDA analysis for outlier loci - morphology
##
## Matt Brachmann (PhDMattyB)
##
## 2019-12-05
##
##############################

setwd('~/PhD/SNP Demographic modelling/Outliers_directory/Morphometric_outliers/RDA')

library(patchwork)
# library(janitor)
# library(devtools)
# library(skimr)
# library(rsed)
# library(data.table)
# library(sjPlot)
library(MASS)
library(vegan)
library(tidyverse)


theme_set(theme_bw())

# Other packages to load

# RDA outlier analysis ----------------------------------------------------
setwd('~/PhD/SNP Demographic modelling/Outliers_directory/Morphometric_outliers/RDA')

## load in the .raw file for the benthic-pelagic pair 
## we're comparing
Genotype = read.delim('SLGBPEL_Genotype_recodeA.raw', 
                      sep = "", 
                      stringsAsFactors = F)
Genotype = read.delim('PolyPops_Genotype_SPELCOMBO_recodeA.raw', 
                      sep = "", 
                      stringsAsFactors = F)

## Arrange the genotpe file by individual id
Genotype = as_tibble(Genotype) %>% 
  arrange(IID)

## load in our morphological variables, 
## select the benthic-pelagic pairs we want to compare
## arrange by individual id
## rename the columns to match the .raw file
# Morphology = read_csv('RDA_Morphological_variables.csv') %>% 
#   # filter(POP_name %in% c('T.SB', 'T.PL')) %>% 
#   arrange(IndividualID) %>% 
#   rename(FID = `#FamilyID`, 
#          IID = IndividualID) %>% 
#   distinct()
Morphology = read_csv('RDA_Morphological_variables_SPELCOMBO.csv') %>% 
  filter(POP_name %in% c('S.LGB', 'S.PEL')) %>%
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
                      # FID == '5' ~ 'Svinvatn - Piscivorous',
                      FID == '6' ~ 'Thingvallavatn - Planktivorous',
                      FID == '7' ~ 'Thingvallavatn - Small benthic',
                      FID == '8' ~ 'Svinavatn - Large benthic',
                      FID == '9' ~ 'Galtabol - Small benthic',
                      FID == '10' ~ 'Galtabol - Piscivorous')))

write_csv(Morphology, 
          'RDA_Morphological_Variables_04.04.2022.csv')

## Make sure the individual ID's match up across the two data sets
Genotype$IID == Morphology$IID

## need a genotype matrix so we will need all snps, 
## all snps start with Affx. 
SNPS = Genotype %>% 
  select(matches("Affx.")) %>% 
  as_tibble(Genotype)

SNPS[is.na(SNPS)] = 0

## calculate the percentage of missing data
(sum(is.na(SNPS))/14187)*100

## We need to get rid of na values for the analysis
## using this function to impute missing values based on the 
## most commone genotype
SNPS.imp = apply(SNPS , 
                 2,
                 function(x) replace(x, is.na(x), 
                                     as.numeric(names(which.max(table(x))))))
## make sure there is no NA's
sum(is.na(SNPS.imp))

## Need to make sure our morphological variables aren't super correlated. 
# par("mar")
# par(mar=c(2,2,2,2))
# dev.off()
# pairs.panels(Morphology[,3:6], scale = T)

## Need a data frame of only the four morphological variables
data_only = Morphology %>% 
  select(FL, 
         PC1, 
         PC2, 
         PC3) %>% 
  as.data.frame()

## Run the actual RDA. Remember an RDA is a PCA of a linear model. 
CHARR_RDA = rda(SNPS.imp ~ FL + PC1 + PC2 + PC3,
                data = data_only,
                scale = T)

## Calculate the rsquared for the model. 
## This calculates the effect size that the model actually explains
## Report the adjusted Rsquared value
RsquareAdj(CHARR_RDA)
##GSBPI = 13.57%

summary(eigenvals(CHARR_RDA, model = "constrained"))
#screeplot(CHARR_RDA)
## RDA1 appears to be explaining a shit load of the variation

## Test that carbon and nitrogen are significant in the RDA
signif.full = anova.cca(CHARR_RDA, 
                        parallel=getOption("mc.cores")) # default is permutation=999


## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## NEED to run the signifincance by each axis to see if each
## axis is actually significant. 
## This will be very computationally extensive and will lock your
## computer right up
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
signif.axis = anova.cca(CHARR_RDA,
                        by="axis",
                        parallel=getOption("mc.cores"))
# signif.axis
## Until we actually have the time to run this and run it
## not on our computer we will asssume only the first axis is significant
## As is explains almost all of the variaiton in genotype. 
## WE NEED TO DOUBLE CHECK THIS!!!!!!!!!

vif.cca(CHARR_RDA)
## All under under 5. Low multicollinearity 

sum = summary(CHARR_RDA)

plot_snps = as_tibble(sum$species)
write_csv(plot_snps, 'SLGBPEL_RDA_SNPS_Morphology.csv')
plot_SNPS = read_csv('SLGBPEL_RDA_SNPS_Morphology.csv')

plot_loc = as_tibble(sum$sites)
write_csv(plot_loc, 'SLGBPEL_RDA_LOC_Morphology.csv')
plot_loc = read_csv('SLGBPEL_RDA_LOC_Morphology.csv')

plot_bp = as.data.frame(sum$biplot)
write_csv(plot_bp, 'SLGBPEL_RDA_BP_Morphology.csv')
plot_bp = read_csv('SLGBPEL_RDA_BP_Morphology.csv')

variables = c('Fork length', 
              'Relative warp 1',
              'Relative warp 2',
              'Relative warp 3')

# plot_bp = plot_bp %>% 
#   mutate(variable = substr(variables, 1, 15))

Morphology = bind_cols(Morphology, 
                       plot_loc)

## GALTABOL ONLY !!!!
# Morphology = Morphology %>% 
#   slice(-20)
# RDA Graph explaining relationships ---------------------------------------

# rda_col = c('#FF0E0E','#2470FF')
# 
# ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ## WE NEED TO FIX THIS GRAPH. 
# ## THE LABELS FOR THE ARROWS ARE FUCKED AND I CANT FIGURE IT OUT
# ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
# 
# RDA_graph = ggplot(data = Morphology, 
#                    aes(x = RDA1, y = RDA2))+
#   geom_point(data = plot_SNPS, 
#              aes(x = RDA1, y = RDA2)) +
#   geom_point(aes(col = POP_name))+
#   scale_color_manual(values = rda_col)+
#   labs(x = 'RDA 1 (79.7% variance explained)',
#        y = 'RDA 2 (8.7% variance explained)')+
#   geom_segment(data = plot_bp, 
#                aes(x = 0, y = 0, 
#                    xend = RDA1, yend = RDA2),
#                color = '#8A8F99', 
#                arrow = arrow(length = unit(0.03, 'npc')))+
#   geom_text(data = plot_bp,
#             check_overlap = F,
#             label = c('Fork length',
#                       'Relative warp 1',
#                       'Relative warp 2',
#                       'Relative warp 3'),
#             fontface = 'bold',
#             # position = position_jitter(width = -1,
#             #                            height = 1))+
#             vjust = 1,
#             hjust = 0,
#             # nudge_y = 0.5, 
#             # nudge_x = 0.5, 
#             position = position_jitter(width = 1.5, 
#                                        height = 1))+
#   # geom_text(data = plot_bp, label = rownames(plot_bp),
#   #           vjust = 0, nudge_y = -0.45)+
#   # geom_label(data = plot_bp, label = rownames(plot_bp),
#   #          vjust = 0, nudge_y = -0.45)+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), 
#         legend.position = 'bottom')
# 
# ## Need to fix the text on the graph so that it lines up with
# ## the arrows
# RDA_graph
# # outlier_check = ENV_ARRANGED2 %>% 
# #   select(IID, 
# #          Carbon, 
# #          Nitrogen, 
# #          POP_NAME, 
# #          RDA1, 
# #          RDA2) %>% 
# #   View()
# 
# ggsave(plot = last_plot(), 
#        '.tiff', 
#        width = 20, 
#        height = 20, 
#        unit = 'cm')
# ## The ouliter has a completely normal carbon and nitrogen value. 
# ## The genotype matrix must be whats driving the outlier nature 
# ## of the individual
# 

SCORES = scores(CHARR_RDA, choices = 1,
                display = 'species')

scores2 = scores(CHARR_RDA, 
                 choices = 2, 
                 display = 'species')
#hist(SCORES[,1], main="Loadings on RDA1")

SNP_SCORES = as.data.frame(cbind(SNP = rownames(SCORES), SCORES)) %>% 
  as_tibble()

SNP_SCORES2 = as.data.frame(cbind(SNP = rownames(scores2), scores2)) %>% 
  as_tibble()

scores_total = left_join(SNP_SCORES, 
                         SNP_SCORES2, 
                         by = 'SNP') %>% 
  as.data.frame()

write_csv(SNP_SCORES, 
          'SLGBPEL_RDA_morpho_scores_total.csv')

#SNP_SCORES = read_csv('SNP_SCORES_Polypops_RDA.csv')

## RDA outliers #####
## FUNCTION FROM THE RDA VIGNETTE
outliers = function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     
  # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               
  # locus names in these tails
}

OUTLIERS1 = outliers(SCORES[,1],4)
length(OUTLIERS1)

outliers2 = outliers(scores2[,1],4)
length(outliers2)

## GSBPI has 350 outliers on the first axis for 4 stdev 
## names of the outlier loci on the first RDA axis. 

CHARR_RDA_OUT = cbind.data.frame(rep(1,
                                     times = length(OUTLIERS1)),
                                 names(OUTLIERS1),
                                 unname(OUTLIERS1))

CHARR_RDA_OUT = as_tibble(CHARR_RDA_OUT) %>% 
  rename(AXIS = 1, SNP = 2, LOADING = 3)

add_predictors = matrix(nrow=(OUTLIERS1), ncol=4)  
# 2 columns for 2 predictors
colnames(add_predictors) = c("Fork length",
                             "Relative warp 1", 
                             "Relative warp 2", 
                             "Relative warp 3")

charr_rda_out2 = cbind.data.frame(rep(1, 
                                      times = length(outliers2)),
                                  names(outliers2),
                                  unname(outliers2))

charr_rda_out2 = as_tibble(charr_rda_out2) %>% 
  rename(AXIS = 2, 
         SNP = 2, 
         LOADING = 3)
## need to trouble shoot this for loop. somethings fucky. 
## THis for loop is FUCKED 
n_cand = length(CHARR_RDA_OUT$SNP)
foo <- matrix(nrow=(n_cand), ncol=4)  
colnames(foo) <- c("Fork length",
                   "Relative warp 1", 
                   "Relative warp 2", 
                   "Relative warp 3")

test = as.data.frame(CHARR_RDA_OUT)

##THe SNP names are not the same after the _
## Get rid of the part after the _ and everything should match up

CHARR_RDA_OUT$SNP = as.character(CHARR_RDA_OUT$SNP)

# CHARR_RDA_OUT$SNP = sed_substitute(CHARR_RDA_OUT$SNP,
#                                    "_.*", 
#                                    "")
CHARR_RDA_OUT$SNP = substr(CHARR_RDA_OUT$SNP,
                           1,
                           regexpr("\\_", CHARR_RDA_OUT$SNP)-1)



SNP_imputed = as_tibble(SNPS)
SNP_name = names(SNP_imputed)
test = gsub("_.*", 
            "", 
            SNP_name)

SNP_imputed = SNP_imputed %>% 
  rename_all(funs(c(test)))


CHARR_RDA_OUT = as.data.frame(CHARR_RDA_OUT)
SNP_imputed = as.data.frame(SNP_imputed)
# head(SNP_imputed)
Morphology = as.data.frame(data_only)

for (i in 1:length(CHARR_RDA_OUT$SNP)) {
  nam <- CHARR_RDA_OUT[i,2]
  snp.gen <- SNP_imputed[,nam]
  foo[i,] <- apply(Morphology,2,function(x) cor(x,snp.gen))
}

candidates = cbind.data.frame(CHARR_RDA_OUT, foo) %>% 
  as_tibble()

## When only one axis is significant, 
## skip the step of removing the duplicate snps across the axes. 
## write this file to the directory and call it later

write_csv(candidates,
          'SLGBPEL_RDA_Morphology_outliers.csv')

## need to find the ouliers for the second axis of the RDA
n_cand = length(charr_rda_out2$SNP)
foo <- matrix(nrow=(n_cand), ncol=4)  
colnames(foo) <- c("Fork length",
                   "Relative warp 1", 
                   "Relative warp 2", 
                   "Relative warp 3")

test = as.data.frame(charr_rda_out2)

##THe SNP names are not the same after the _
## Get rid of the part after the _ and everything should match up

charr_rda_out2$SNP = as.character(charr_rda_out2$SNP)

charr_rda_out2$SNP = gsub("_.*", 
                          "",
                          charr_rda_out2$SNP)

# charr_rda_out2$SNP = sed_substitute(charr_rda_out2$SNP,
#                                    "_.*", 
#                                    "")

SNP_imputed = as_tibble(SNPS)
SNP_name = names(SNP_imputed)
test = gsub("_.*", 
            "", 
            SNP_name)

SNP_imputed = SNP_imputed %>% 
  rename_all(funs(c(test)))


charr_rda_out2 = as.data.frame(charr_rda_out2)
SNP_imputed = as.data.frame(SNP_imputed)
# head(SNP_imputed)
Morphology = as.data.frame(data_only)

for (i in 1:length(charr_rda_out2$SNP)) {
  nam <- charr_rda_out2[i,2]
  snp.gen <- SNP_imputed[,nam]
  foo[i,] <- apply(Morphology,2,function(x) cor(x,snp.gen))
}

candidates2 = cbind.data.frame(charr_rda_out2, foo) %>% 
  as_tibble()
write_csv(candidates2,
          'SLGBPEL_RDA_Morphology_outliers_axis2.csv')

# rda duplicates ----------------------------------------------------------
Outliers1 = read_csv('SLGBPEL_RDA_Morphology_outliers.csv')
Outliers2 = read_csv('SLGBPEL_RDA_Morphology_outliers_axis2.csv') %>% 
  rename(AXIS = 1)

outliers_total = left_join(Outliers1, 
                            Outliers2, 
                            by = 'SNP')

outliers_axes12 = bind_rows(Outliers1, 
                            Outliers2) 

test = intersect(Outliers1, 
                 Outliers2)


write_csv(Outliers1, 
          'SLGBPEL_RDA_Morphology_total_outliers.csv')
## RDA align to MAP ####
CHARR_RDA_OUT = read_csv('Polypop_RDA_Morphology_outliers_SPELCOMBO.csv')

# CHARR_RDA_OUT = read_csv('Polypop_RDA_Morphology_total_outliers_SPELCOMBO.csv')
CHARR_RDA_OUT = read_csv('SLGBPEL_RDA_Morphology_total_outliers.csv')

# CHARR_RDA_OUT$SNP = as.character(CHARR_RDA_OUT$SNP)
# CHARR_RDA_OUT$SNP = sed_substitute(CHARR_RDA_OUT$SNP,
#                                    "Affx.",
#                                    "Affx-")

CHARR_RDA_OUT$SNP = gsub("Affx.",
                         "Affx-",
                         CHARR_RDA_OUT$SNP)
# 
# CHARR_RDA_OUT$SNP = sed_substitute(CHARR_RDA_OUT$SNP,
#                                    "_.*",
#                                    "")

OUT_MAP = MAP3[MAP3$`Marker ID` %in% CHARR_RDA_OUT$SNP,]

OUT_RDA_MAPPED = merge(OUT_MAP, CHARR_RDA_OUT, 
                       by.x = 'Marker ID', by.y = 'SNP') %>% 
  as_tibble()


OUT_RDA_MAPPED$`Marker ID` = as.character(OUT_RDA_MAPPED$`Marker ID`)
OUT_RDA_MAPPED$`Physical position` = as.numeric(OUT_RDA_MAPPED$`Physical position`)
OUT_RDA_MAPPED$CHROME3 = as.character(OUT_RDA_MAPPED$CHROME3)

# write_tsv(OUT_RDA_MAPPED, 'Galtabol_RDA_outliers_mapped.txt',
#           col_names = T)
write_csv(OUT_RDA_MAPPED, 
          'SLGBPEL_RDA_morphology_outliers_total_mapped.csv')

RDA_MAPPED = read_csv('SLGBPEL_RDA_morphology_outliers_total_mapped.csv')

## Define outliers ####
CHROME_GROUP = rep('RDA OUTLIER', 
                   length(RDA_MAPPED$CHROME3)) %>% 
  as_tibble()

## Joins the outliers and labels them as RDA outlier loci
RDA_GROUP = bind_cols(RDA_MAPPED, CHROME_GROUP) %>% 
  as_tibble() %>% 
  rename(Morpho_out = value)

#Get RDA total map info
SNP_SCORES$SNP = gsub("Affx.", 
                      "Affx-", 
                      SNP_SCORES$SNP)
SNP_SCORES$SNP = gsub("_.*", 
                      "", 
                      SNP_SCORES$SNP)


SNP_SCORES2$SNP = as.character(SNP_SCORES2$SNP)
SNP_SCORES2$SNP = gsub("Affx.", 
                       "Affx-", 
                       SNP_SCORES2$SNP)
SNP_SCORES2$SNP = gsub("_.*", 
                       "",
                       SNP_SCORES2$SNP)
# 
# 
# SNP_SCORES2$SNP = sed_substitute(SNP_SCORES2$SNP, 
#                                 "Affx.", 
#                                 "Affx-")
# 
# SNP_SCORES2$SNP = sed_substitute(SNP_SCORES2$SNP,
#                                 "_.*", 
#                                 "")

snppy_scores = bind_cols(SNP_SCORES, 
                         SNP_SCORES2) %>% 
  select(-SNP1)


TOTAL_RDA_MAPPED = merge(MAP3, snppy_scores,
                         by.x = 'Marker ID', by.y = 'SNP') %>% 
  as_tibble()

TOTAL_RDA_MAPPED = merge(MAP3, SNP_SCORES,
                         by.x = 'Marker ID', by.y = 'SNP') %>% 
  as_tibble()


TOTAL_RDA_MAPPED = TOTAL_RDA_MAPPED %>% 
  rename(MARKER_ID = `Marker ID`, 
         PDIST = `Physical position`, 
         GDIST = `Genetic distance`)

TOTAL_RDA_MAPPED$RDA1 = as.numeric(as.character(TOTAL_RDA_MAPPED$RDA1))
TOTAL_RDA_MAPPED$RDA2 = as.numeric(as.character(TOTAL_RDA_MAPPED$RDA2))


write_csv(TOTAL_RDA_MAPPED, 
          'SLGBPEL_RDA_Morphology_allsnps_total_mapped.csv')

# RDA ggman ---------------------------------------------------------------

TOTAL_RDA_MAPPED = read_csv('GSBPI_RDA_Morphology_allsnps_mapped.csv') %>% 
  rename(MARKER_ID = 1, PDIST = 4)
RDA_Outliers = read_csv('GSBPI_RDA_morphology_outliers_mapped.csv')
KeyStone_cand_list = read_csv('RDA_IsoMorpho_Outliers_SNPList.csv') %>% 
  as.data.frame()

# str(TOTAL_RDA_MAPPED)
# 
# TOTAL_RDA_MAPPED$MARKER_ID
# sum(is.na(TOTAL_RDA_MAPPED$PDIST))
# NO_NAs = na.omit(TOTAL_RDA_MAPPED)

rda_man_plot = ggman(TOTAL_RDA_MAPPED, 
                     chrom = 'CHROME3', 
                     pvalue = 'RDA1',
                     snp = 'MARKER_ID',
                     bp = 'PDIST',
                     logTransform = F,
                     pointSize = 2, 
                     title = 'RDA - G: Benthic-Pelagic',
                     xlabel = 'Chromosome',
                     ymin = -0.4, 
                     ymax = 0.4, 
                     lineColour = 'black', 
                     relative.positions = TRUE)

out_highlight = as.character(KeyStone_cand_list$ID)
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
          element_text(size  = 8,
                       angle = 90,
                       hjust = 1,
                       vjust = 1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot

ggmanZoom(plot, chromosome = 'AC13')
ggmanZoom(plot, chromosome = 'AC24')
ggmanZoom(plot, chromosome = 'AC25')
ggmanZoom(plot, chromosome = 'AC33')

ggsave(plot = last_plot(), 
       'SLGBPL_RDA_Phenotype_Keystone_outliers.tiff', 
       height = 20,
       width = 20, 
       unit = 'cm')


rda_man_plot = ggman(TOTAL_RDA_MAPPED, 
                     chrom = 'CHROME3', 
                     pvalue = 'RDA1',
                     snp = 'MARKER_ID',
                     bp = 'PDIST',
                     logTransform = F,
                     pointSize = 2, 
                     title = 'RDA - G: Benthic-Pelagic',
                     xlabel = 'Chromosome',
                     ymin = -0.4, 
                     ymax = 0.4, 
                     lineColour = 'black', 
                     relative.positions = TRUE)
