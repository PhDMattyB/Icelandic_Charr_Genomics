##############################
## Loci under local adaptation
##
## Matt Brachmann (PhDMattyB)
##
## 2019-06-20
##
##############################

setwd('~/PhD/SNP Demographic modelling/Outliers_directory')

# library(patchwork)
# library(janitor)
# library(devtools)
# library(rsed)
# library(ggman)
# library(qvalue)
# library(sjPlot)
library(vegan)
library(psych)
library(tidyverse)

theme_set(theme_bw())

## MAP FILE #####
# update_map = read_tsv('Nov052018_plink_input_icelandic_matt.map') %>% 
#   arrange(`Marker ID`) %>% 
#   select(`#Chromosome`,
#          `Assigned-Molecule`,
#          `Marker ID`)
# write_tsv(update_map, 
#           '11.12.2019_Chromosome_linkageGroup.txt')

MAP = read_tsv('Feb202019_Poly_Plink_input.map') %>% 
  arrange(`Marker ID`)

MAP$`Marker ID` == update_map$`Marker ID`
MAP$`#Chromosome` == update_map$`#Chromosome`

# MAP %>% 
#   group_by(`#Chromosome`) %>% 
#   summarise(n()) %>% 
#   View()

sum(is.na(MAP$`Physical position`))

MAP$`#Chromosome`
chr_names = c('AC01', 
              'AC02', 
              'AC03', 
              'AC04p', 
              'AC04q.1:29', 
              'AC04q.2', 
              'AC05', 
              'AC06.1', 
              'AC06.2', 
              'AC07', 
              'AC08', 
              'AC09', 
              'AC10', 
              'AC11', 
              'AC12', 
              'AC13', 
              'AC14', 
              'AC15', 
              'AC16', 
              'AC17', 
              'AC18', 
              'AC19', 
              'AC20', 
              'AC21', 
              'AC22', 
              'AC23', 
              'AC24', 
              'AC25', 
              'AC26', 
              'AC27', 
              'AC28', 
              'AC30', 
              'AC31', 
              'AC32', 
              'AC33', 
              'AC34', 
              'AC35', 
              'AC36', 
              'AC37', 
              'contigs')

t(chr_names)
length(chr_names)

## AC 6 issues
## Arctic charr genome: AC 6 is split in to two acrocentric chromosomes
## AC6 is labeled as AC6.1 and AC6.2  
## Based on Christines stuff, AC6 is a fused metacentric chromosome
## in Icelandic Arctic charr. 

## AC4.1:29 issues
## Arctic charr genome: AC4.1q:29 is a single chromosome
## Based on Christines stuff, they are two separate chromosomes. 
## AC4.1q is everything >40,000,000 bp
## AC29 is everything <40,000,000 bp

MAP3 = mutate(.data = MAP,
              CHROME3 = as.factor(case_when(
                `#Chromosome` == '1' ~ 'AC01',
                `#Chromosome` == '2' ~ 'AC02',
                `#Chromosome` == '3' ~ 'AC03',
                `#Chromosome` == '4' ~ 'AC04p',
                `#Chromosome` == '5' ~ 'AC04q.1:29',
                `#Chromosome` == '6' ~ 'AC04q.2',
                `#Chromosome` == '7' ~ 'AC05',
                `#Chromosome` == '8' ~ 'AC06',
                `#Chromosome` == '9' ~ 'AC06',
                `#Chromosome` == '10' ~ 'AC07',
                `#Chromosome` == '11' ~ 'AC08',
                `#Chromosome` == '12' ~ 'AC09',
                `#Chromosome` == '13' ~ 'AC10',
                `#Chromosome` == '14' ~ 'AC11',
                `#Chromosome` == '15' ~ 'AC12',
                `#Chromosome` == '16' ~ 'AC13',
                `#Chromosome` == '17' ~ 'AC14',
                `#Chromosome` == '18' ~ 'AC15',
                `#Chromosome` == '19' ~ 'AC16',
                `#Chromosome` == '20' ~ 'AC17',
                `#Chromosome` == '21' ~ 'AC18',
                `#Chromosome` == '22' ~ 'AC19',
                `#Chromosome` == '23' ~ 'AC20',
                `#Chromosome` == '24' ~ 'AC21',
                `#Chromosome` == '25' ~ 'AC22',
                `#Chromosome` == '26' ~ 'AC23',
                `#Chromosome` == '27' ~ 'AC24',
                `#Chromosome` == '28' ~ 'AC25',
                `#Chromosome` == '29' ~ 'AC26',
                `#Chromosome` == '30' ~ 'AC27',
                `#Chromosome` == '31' ~ 'AC28',
                `#Chromosome` == '32' ~ 'AC30',
                `#Chromosome` == '33' ~ 'AC31',
                `#Chromosome` == '34' ~ 'AC32',
                `#Chromosome` == '35' ~ 'AC33',
                `#Chromosome` == '36' ~ 'AC34',
                `#Chromosome` == '37' ~ 'AC35',
                `#Chromosome` == '38' ~ 'AC36',
                `#Chromosome` == '39' ~ 'AC37',
                `#Chromosome` > '39' ~ 'Contigs')))

is.na(MAP3$CHROME3)
MAP3$CHROME3[is.na(MAP3$CHROME3)] = 'Contigs'
MAP3$CHROME3
#View(MAP3$CHROME3)

## Need to split chromosome AC04q.1:29 into two different chromosomes
## We split the chromosomes into AC04q.1 and AC29
AC04q.1_29 = MAP3 %>% 
  filter(CHROME3 == 'AC04q.1:29') %>% 
  rename(PP = 'Physical position')

test = mutate(.data = AC04q.1_29,
              CHROME3 = as.factor(case_when(
                PP > '40000000' ~ 'AC04q.1',
                PP < '40000000' ~ 'AC29')))

test = test %>% 
  rename(`Physical position` = PP)

MAP4 = MAP3 %>% 
  filter(CHROME3 != 'AC04q.1:29')

MAP3 = bind_rows(MAP4,
                 test) %>% 
  arrange(`Marker ID`)

#View(MAP3)

sum(MAP$`Marker ID`== MAP3$`Marker ID`)

## PED FILES #####
## Need to load the ped file with only polymorphic populations 
## Once loaded we can split the genotype file up into each populations
PED = read_tsv('Feb202019_Poly_Plink_input.ped')

# Combine S pelagic -------------------------------------------------------

S.PI = PED %>% 
  filter(`#FamilyID` == '5') %>% 
  dplyr::select(`#FamilyID`, 
                IndividualID)

value = rep(4, length(S.PI$IndividualID)) %>% 
  as.numeric() %>% 
  as_tibble()

S.PI_id = bind_cols(S.PI, value) %>% 
  select(-`#FamilyID`) %>% 
  rename(`#FamilyID` = value)

S.PI_Geno = PED %>% 
  filter(`#FamilyID` == '5') %>% 
  select(contains('Affx-'))

S.PI_data = bind_cols(S.PI_id, 
                      S.PI_Geno) %>% 
  select(`#FamilyID`, 
         IndividualID, 
         contains('Affx-'))

PED = PED %>% 
  filter(`#FamilyID` != '5')

PED = bind_rows(PED,
                S.PI_data)

write_tsv(PED, 
          'PolyPops_Genotype_SPELCOMBO.PED')
# back to analysis --------------------------------------------------------
PED = read_tsv('PolyPops_Genotype_SPELCOMBO.PED')

View(PED[1:10, 10:20])
# Gal_ped = PED %>% 
#   filter(`#FamilyID` %in% c('9', '10'))
# write_tsv(Gal_ped, 'Galtabol_Genotype.ped', col_names = T)

PED = PED %>% 
  filter(`#FamilyID` %in% c('8', '4'))

write_tsv(PED, 'SLGBPEL_Genotype.ped', col_names = T)

## We will use this ped file to delinate populations
# POPS = read_tsv('Feb202019_Poly_Plink_input.ped')
POPS = PED %>% select(`#FamilyID`, IndividualID) %>% 
  filter(`#FamilyID` %in% c('9', '10'))

# cluster = mutate(.data = POPS, 
#                  cluster = as.factor(case_when(
#                    `#FamilyID` == '9' ~ 'GSB',
#                    `#FamilyID` == '10' ~ 'GPI'
#                  )))
# 
# write_tsv(cluster, 
#           'Galtabol_pops_fst.txt', col_names = F)

POPS = mutate(.data = POPS,
              POP_name = as.factor(case_when(
                `#FamilyID` == "1" ~ "T.LGB",
                `#FamilyID` == '2' ~ 'V.BR',
                `#FamilyID` == '3' ~ 'V.SIL',
                `#FamilyID` == '4'~ 'S.PEL',
                `#FamilyID` == '5' ~ 'S.PEL',
                `#FamilyID` == '6' ~ 'T.PL',
                `#FamilyID` == '7' ~ 'T.SB',
                `#FamilyID` == '8' ~ 'S.LGB',
                `#FamilyID` == '9' ~ 'G.SB',
                `#FamilyID` == '10' ~ 'G.PI'))) 

POPS %>% group_by(POP_name) %>% 
  summarise(n = n())



## LOCAL ADAPTATION STUFF #########
## LFMM #####
## obtain the recoded data from PLINK
RECODE = read.delim2('Galtabol_recode12.ped', 
                     col.names = F, 
                     sep = "\t")

RECODE = read.delim2('Thingvallavatn_Genotype_12.ped', 
                     col.names = F, 
                     sep = "\t")


head(RECODE)


RECODE_CHARR = fread('Galtabol_recode12.ped', 
                     stringsAsFactors = F, 
                     header = F, 
                     data.table = F)
head(RECODE_CHARR)

# CHARR_POP2 = read_tsv('Vatnshlidarvatn_Genotype.ped')

ENVIRONMENT = read_csv('LFMM_Enviro_file.csv') %>% 
  filter(`#FamilyID` %in% c('4','8')) %>% 
  arrange(IndividualID) %>% 
  select(Carbon, Nitrogen)

# SLGBPI = read_csv('LFMM_Enviro_file.csv') %>% 
#   filter(`#FamilyID` %in% c('5','8')) %>% 
#   arrange(IndividualID) %>% 
#   select(IndividualID) 
# 
# write_csv(SLGBPI, 
#           'SLGBPI_enviroID.csv')
# 
# SLGBPL = ENVIRONMENT = read_csv('LFMM_Enviro_file.csv') %>% 
#   filter(`#FamilyID` %in% c('4','8')) %>% 
#   arrange(IndividualID) %>% 
#   select(IndividualID)
# write_csv(SLGBPL, 
#           'SLGBPL_enviroID.csv')

fwrite(ENVIRONMENT,'Galtabol_EnviroVars.env', 
       sep = "\t", 
       row.names = F, 
       col.names = F, 
       quote = F)

PED = read_tsv('Feb202019_Poly_Plink_input.ped') 

PED = PED %>% 
  arrange(IndividualID) %>% 
  write_tsv('Polypop_Genotype_arranged.ped')

ped2lfmm(input.file = "Svinavatn_Genotype_arranged_12.ped")

## DO NOT RUN THIS ON YOUR LAPTOP
## RUN THIS ON ONE OF THE LAB COMPUTERS
LFMM = lfmm("Vatnshlidarvatn_Genotype_12.lfmm", 
            "Vatnshlidarvatn_EnviroVars.env", 
            K = 2, 
            repetitions = 10, 
            CPU = 18, 
            project = 'new')
## carbon LFMM ####
# zs_carbon = z.scores(LFMM, K = 2, d = 1)
# zs_med_carbon = apply(zs_carbon, MARGIN = 1, median)
# lambda_carbon = median(zs_med_carbon^2)/qchisq(0.5, df = 1)
# adj_pval_carbon = pchisq(zs_med_carbon^2/lambda_carbon,
#                          df = 1, lower = F)
# hist(adj_pval_carbon, col = 'red')
# qval_carbon = qvalue(adj_pval_carbon)
# MAP3


## ran LFMM on another computer and needed to import the data
setwd('~/PhD/SNP Demographic modelling/Outliers_directory/LFMM/LFMM_data')

data = read_csv('LFMM_Svin_PLLGB_candidates_carbon_unmapped.csv')
# data = read_csv('LFMM_Galtabol_candidates_carbon_unmapped.csv')
df = bind_cols(MAP3, data)
sum(is.na(df$`Physical position`))
df = na.omit(df)
out = df[(df$qvalue <0.001),]
out
# write_csv(out, 
#           'LFMM_Galtabol_Carbon_Outlier_loci_FINAL.csv')

## THese are the OUTLIER LOCI FOR THE LFMM ANALYSIS
write_csv(out,
          'LFMM_Svin_PLLGB_Carbon_Outlier_loci_FINAL.csv')

## Need all outlier and non-outlier snps as well
TOTAL_LFMM_MAPPED = merge(MAP3, 
                          df,
                         by.x = 'Marker ID', 
                         by.y = 'Marker ID') %>% 
  as_tibble()

write_csv(TOTAL_LFMM_MAPPED, 
          'Svin_PLLGB_Carbon_AllSNPS.csv')

setwd('~/PhD/SNP Demographic modelling/Outliers_directory')

df = read_csv('Svin_PLLGB_Carbon_AllSNPS.csv') %>% 
  select(`Marker ID`, 
         `Genetic distance.x`, 
         `Physical position.x`, 
         CHROME3.x, 
         qvalue, 
         pvalue) %>% 
  rename('Genetic distance' = `Genetic distance.x`, 
         'Physical position' = `Physical position.x`, 
         Chromosome = CHROME3.x)
df_out = read_csv('LFMM_Svin_PLLGB_Carbon_Outlier_loci_FINAL.csv')

## Read in the candidate loci file for the LFMM carbon outliers
man_plot = ggman(df, 
      chrom = 'Chromosome', 
      pvalue = 'pvalue',
      snp = 'Marker ID', 
      bp = 'Physical position',
      pointSize = 1.5, 
      xlabel = 'Chromosome',
      sigLine = -log10(max(df_out$pvalue)), 
      lineColour = 'black',
      title = 'LFMM - Svinavatn - PLLGB - Carbon', 
      relative.positions = TRUE)

LFMM_carbon_man = 
  man_plot + scale_color_manual(values = c('#FF0E0E',
                                         '#2470FF', 
                                         '#1455CC'))+
  theme(axis.text.x =
          element_text(size  = 8,
                       angle = 90,
                       hjust = 1,
                       vjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  # theme(axis.text.x =
  #         element_text(size  = 8,
  #                      angle = 90,
  #                      hjust = 1,
  #                      vjust = 1))

## nitrogen LFMM ####
# zs_nitrogen = z.scores(LFMM, K = 2, d = 2)
# zs_med_nitrogen = apply(zs_nitrogen, MARGIN = 1, median)
# lambda_nitrogen = median(zs_med_nitrogen^2)/qchisq(0.5, df = 1)
# adj_pval_nitrogen = pchisq(zs_med_nitrogen^2/lambda_nitrogen,
#                            df = 1, lower = F)
# hist(adj_pval_nitrogen, col = 'red')
# #hist(adj_pval_carbon, col = 'red')
# qval_nitrogen = qvalue(adj_pval_nitrogen)
# MAP3
# 
# nitrogen_data = data.frame(cbind(MAP3, qval_nitrogen$qvalues,
#                                  qval_nitrogen$pvalues), 
#                            stringsAsFactors = F)
# 
# write.table(nitrogen_data, file = 'nitrogen_lfmm_significantSNPS',
#             quote = F, 
#             sep = "\t", 
#             row.names = F,
#             col.names = T)

## read in the new data from the rerunning of the LFMM analysis
## Read in the candidate loci file for the LFMM nitrogen outliers
nitrogen = read_csv('Svin_PLLGB_Nitrogen_AllSNPS.csv') %>% 
  select(`Marker ID`, 
         `Genetic distance.x`, 
         `Physical position.x`, 
         CHROME3.x, 
         qvalue, 
         pvalue) %>% 
  rename('Genetic distance' = `Genetic distance.x`, 
         'Physical position' = `Physical position.x`, 
         Chromosome = CHROME3.x)
nitrogen_out = read_csv('LFMM_Svin_PLLGB_Nitrogen_Outlier_loci_FINAL.csv')

man_plot = ggman(nitrogen, 
                 chrom = 'Chromosome', 
                 pvalue = 'pvalue',
                 snp = 'Marker ID', 
                 bp = 'Physical position',
                 pointSize = 1.5, 
                 xlabel = 'Chromosome',
                 sigLine = -log10(max(nitrogen_out$pvalue)), 
                 lineColour = 'black', 
                 title = 'LFMM - Svinavatn - PLLGB - Nitrogen', 
                 relative.positions = TRUE)

LFMM_nitrogen_man = man_plot + scale_color_manual(values = c('#FF0E0E',
                                                           '#2470FF', 
                                                           '#1455CC'))+
  theme(axis.text.x =
          element_text(size  = 8,
                       angle = 90,
                       hjust = 1,
                       vjust = 1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


LFMM_carbon_man/LFMM_nitrogen_man

ggsave(plot = last_plot(), 
       'LFMM_ManhattanPlot_Svinavatn_PLLGB.tiff', 
       height = 20, 
       width = 25, 
       units = 'cm')

## CN LFMM outlier overlap #####
carbon_outlier = carbon_qvalue %>% 
  arrange(`Marker ID`)

carb_label = rep('CARBON OUTLIER', 
                 length(carbon_outlier$CHROME3)) %>% 
  as_tibble()

carbon_outlier = bind_cols(carbon_outlier, 
                           carb_label)  %>% 
  rename(outlier = value) %>% 
  select(value1)

nitrogen_outlier = nitrogen_qvalue %>% 
  arrange(`Marker ID`)

nit_label = rep('NITROGEN OUTLIER', 
                 length(nitrogen_outlier$CHROME3)) %>% 
  as_tibble()

nitrogen_outlier = bind_cols(nitrogen_outlier, 
                           nit_label)  %>% 
  rename(outlier = value) 


snp_overlap = intersect(carbon_outlier$`Marker ID`,
                        nitrogen_outlier$`Marker ID`)
length(snp_overlap)
# 116 significant snps overlap between the lfmm runs for carbon
# and nitrogen stable isotope signatures

## Thingvallavatn: there are 299 snps that overlap between 
## the lfmm runs for carbon and nitrogen stable isotope signatures




LFMM_carbon = read_csv('LFMM_carbon_outliers_Svinavatn.csv')
LFMM_nitrogen = read_csv('LFMM_Nitrogen_outliers_Svinavatn.csv')
# PCA_outliers = read_tsv('') %>% 
#   na.omit()

PCA_outliers = read_csv('PCAdapt_svin_Outliers_k2_Q0.01.csv') %>% 
  na.omit()

# PCA_outliers = read.delim('Galtabol_PCA_OUTLIERS_k2.txt', 
#                           stringsAsFactors = F) %>% 
#   as.tibble()

venn_pal = c('#2776CC',
             '#CC2B1C',
             '#10B372')

#B300A2
venn.diagram(
  x = list(PCA_outliers$MARKER_ID,
           LFMM_carbon$`Marker ID`,
           LFMM_nitrogen$`Marker ID`),
  category.names = c('PCAdapt',
                     'LFMM Carbon',
                     'LFMM Nitrogen'),
  filename = 'Svinavatn_outlier_overlap.tiff',
  imagetype = 'tiff',
  # height = 400,
  # width = 400, 
  # resolution = 300,
  # compression = 'lzw',
  
  # cirlces
  lwd = 2,
  lty = 'blank',
  fill = venn_pal,
  
  # Numbers
  cex = .8,
  fontface = "bold",
  fontfamily = "sans")
  
  # Set names
  # cat.cex = 0.6,
  # cat.fontface = "bold",
  # cat.default.pos = "outer",
  # cat.pos = c(-27, 27, 135),
  # cat.dist = c(0.055, 0.055, 0.085),
  # cat.fontfamily = "sans",
  # rotation = 1)


LFMM_carbon = LFMM_carbon %>% 
  rename(MARKER_ID = `Marker ID`)
LFMM_nitrogen = LFMM_nitrogen %>% 
  rename(MARKER_ID = `Marker ID`)
PCA_outliers

# Two_method = LFMM_carbon$MARKER_ID[LFMM_carbon$MARKER_ID %in% LFMM_nitrogen$MARKER_ID]
# length(Two_method)
# Three = Two_method[Two_method %in% PCA_outliers$MARKER_ID]
# length(Three)
# 
# '%!in%' <- function(x,y)!('%in%'(x,y)) # not in
# test = LFMM_carbon$MARKER_ID[LFMM_carbon$MARKER_ID %!in% LFMM_nitrogen$MARKER_ID]
# length(test)
# test2 = test[test %!in% PCA_outliers$MARKER_ID]
# length(test2)

carb_snps_only = anti_join(LFMM_carbon,
                           LFMM_nitrogen,
                           by = 'MARKER_ID')

carb_snps_only = anti_join(carb_snps_only,
                            PCA_outliers, 
                            by = 'MARKER_ID')


nit_snps_only = anti_join(LFMM_nitrogen, 
                          LFMM_carbon, 
                          by = 'MARKER_ID')

nit_snps_only = anti_join(nit_snps_only, 
                          PCA_outliers, 
                          by = 'MARKER_ID')

three_method_overlap = inner_join(LFMM_carbon, 
                                  LFMM_nitrogen, 
                                  by = 'MARKER_ID')

three_method_overlap = inner_join(three_method_overlap, 
                                  PCA_outliers, 
                                  by = 'MARKER_ID') 

three_method_overlap %>% filter(CHROME3 != 'Contigs') %>% 
  arrange(CHROME3)

LFMM_overlap = inner_join(LFMM_carbon, 
                          LFMM_nitrogen, 
                          by = 'MARKER_ID')
LFMM_overlap = anti_join(LFMM_overlap, 
                         PCA_outliers, 
                         by = 'MARKER_ID') %>% 
  filter(CHROME3.x != 'Contigs') %>% 
  arrange(CHROME3.x) %>% 
  View()


# carb_only_label = rep('CARBON OUTLIER ONLY',
#                       length(carb_snps_only$`Marker ID`)) %>% 
#   as.tibble()
# 
# carb_snps_only = bind_cols(carb_snps_only, carb_only_label) %>% 
#   rename(CN_out = value)
# 
# nit_snps_only = anti_join(nitrogen_outlier, 
#                           carbon_outlier, 
#                           by = 'Marker ID')
# 
# nit_only_label = rep('NITROGEN OUTLIER ONLY',
#                      length(nit_snps_only$`Marker ID`)) %>% 
#   as.tibble()
# 
# nit_snps_only = bind_cols(nit_snps_only, nit_only_label) %>% 
#   rename(CN_out = value)
# 
# 
# overlap_both = semi_join(carbon_outlier, 
#                          nitrogen_outlier, 
#                          by = 'Marker ID')
# both_label = rep('CN OUTLIERS',
#                  length(overlap_both$`Marker ID`)) %>% 
#   as.tibble()
# 
# overlap_both = bind_cols(overlap_both, both_label) %>% 
#   rename(CN_out = value)
# 
# LFMM_SNP_data = bind_rows(carb_snps_only, nit_snps_only, 
#                           overlap_both)
# LFMM_SNP_data$CN_out
# 
# write_csv(LFMM_SNP_data, 
#           'LFMM_Galtabol_Outliers.csv')

# head(LFMM_SNP_data$OVERLAP)
# write.table(LFMM_SNP_data, file = 'LFMM_SNP_outliers', 
#             quote = T, 
#             sep = "\t", 
#             row.names = F, 
#             col.names = T)


## RDA ####
setwd('~/PhD/SNP Demographic modelling/Outliers_directory')
RAW = read.delim("Galtabol_recodeA.raw",
                 sep = "",
                 stringsAsFactors = F)
RAW = read.delim("PolyPops_Genotype_SPELCOMBO_recodeA.raw",
                 sep = "",
                 stringsAsFactors = F)
dim(RAW)

RAW = read.delim("SLGBPEL_Genotype_recodeA.raw", 
                  sep = "",
                 stringsAsFactors = FALSE)

# RAW = read.delim("Galtabol_recodeA.raw",
#                  sep = "",
#                  stringsAsFactors = F)

# POPS = mutate(.data = POPS,
#               POP_name = as.factor(case_when(
#                 #`#FamilyID` == "1" ~ "T.LGB",
#                 `#FamilyID` == '2' ~ 'V.BR',
#                 `#FamilyID` == '3' ~ 'V.SIL')))
#`#FamilyID` == '4'~ 'S.PL',
#`#FamilyID` == '5' ~ 'S.PI',
#`#FamilyID` == '6' ~ 'T.PL',
#`#FamilyID` == '7' ~ 'T.SB',
#`#FamilyID` == '8' ~ 'S.LGB',
#`#FamilyID` == '9' ~ 'G.SB',
#`#FamilyID` == '10' ~ 'G.PI'))) 

ENVIRONMENT = read_csv('LFMM_Enviro_file.csv') %>% 
  # filter(`#FamilyID` %in% c('9', '10')) %>% 
  select(`#FamilyID`, IndividualID, Carbon, Nitrogen) 

ENVIRONMENT = read_csv('LFMM_Enviro_SPELCOMBO.csv') %>% 
  filter(`#FamilyID` %in% c('9', '10')) %>%
  select(`#FamilyID`, IndividualID, Carbon, Nitrogen) 

## For the RDA the environmental data and the genotype files 
## need to match row for row. 
## Therefore, we need to arrange the columns by individual ID
ENV_ARRANGED = ENVIRONMENT %>%  
  rename(IID = IndividualID, 
         FID = `#FamilyID`) %>% 
  arrange(IID)

RAW_CHARR_TIB = as_tibble(RAW) %>% 
  arrange(IID)

## Make sure the individual ID's match up across the two data sets
ENV_ARRANGED$IID == RAW_CHARR_TIB$IID

SNPS = RAW_CHARR_TIB %>%
  select(matches("Affx.")) %>% 
  as_tibble(RAW_CHARR_TIB)

# View(SNPS[1:20,1:10])

SNPS[is.na(SNPS)] = 0

# (sum(is.na(SNPS))/14187)*100
sum(is.na(SNPS))
## Get rid of the na's
## This will impute the missing genotype as the average/mean genotype 
## within the dataset. 
## The way RecodeA works is that the monomorphic alleles
## are coded as both 0s or NAs. If we change the NA's to 0
## then we get around the issue of imputing data. 
## THe SNP chip gives us high quality data, there is no
## was that 68% of the data is missing. Especially when this is not
## the case in the original ped files. There is no missing data
## in the original ped files. 

SNPS.imp = apply(SNPS , 
             2,
             function(x) replace(x, is.na(x), 
                                 as.numeric(names(which.max(table(x))))))
## make sure there is no NA's
## for the Galtabol data set there are 1675 missing points
## This accounts for 11% of the overall data. 
sum(is.na(SNPS.imp))

## double check to make sure the 
## carbon and nitrogen data isn't hyper correlated
pairs.panels(ENV_ARRANGED[,3:4], scale = T)

CandN = ENV_ARRANGED %>% 
  select(Carbon, Nitrogen) %>% 
  as.data.frame()

dim(SNPS)
## Run the actual RDA. Remember an RDA is a PCA of a linear model. 
# CHARR_RDA_imp = rda(SNPS.imp ~ Carbon + Nitrogen, data = ENV_ARRANGED,
#                 scale = T)

CHARR_RDA = rda(SNPS ~ Carbon + Nitrogen, data = ENV_ARRANGED,
                    scale = T)

## Calculate the rsquared for the model. 
## This calculates the effect size that the model actually explains
## Report the adjusted Rsquared value
## adjusted rsquared for galtabol = 0.0182 or 1.82%
## Carbon and nitrogen have an effect of 1.82% on the genotype matrix.  
RsquareAdj(CHARR_RDA)
## Vatnshlidarvatn = 1.99% or 2.0%
## Thingvalalvatn = 4.5%
## Svinavatn = 6.0%
## Svinavatn PLLGB = 8.37%
## Thingvallavatn LGBPL = 6.04%
## Thingvallavatn PLSB = 5.44%

summary(eigenvals(CHARR_RDA, 
                  model = "constrained"))
## Galtabol RDA1 = 76.22% RDA2 = 23.78%
## Vatnshlidarvatn RDA1 = 68.78% RDA2 = 31.22%
## Thingvallavatn RDA1 = 84.11% RDA2 = 15.89%
## Svinavatn RDA1 = 86.63% 13.37%
## Svinavatn PLLGB RDA1 = 84.28% RDA2 = 15.72
## Thingvallavatn LGBPL RDA1 = 87.3% RDA2 = 12.7%
## Thingvallavatn PLBS RDA1 = 84.17% RDA2 = 15.83%


## Run a screeplot but we don't really have enough variables to really
## visualize it. It's a box plot with two bars. 
screeplot(CHARR_RDA)

## Test that carbon and nitrogen are significant in the RDA
signif.full = anova.cca(CHARR_RDA, 
                        parallel=getOption("mc.cores")) # default is permutation=999

## DO NOT RUN THIS AGAIN!!!!!!!!
signif.axis = anova.cca(CHARR_RDA, 
                        by="axis", 
                        parallel=getOption("mc.cores"))
signif.axis
##Galtabol:
## only the first RDA axis is actually signficant
## Vatnshlidarvatn: Only RDA1 significant
## THingvallavatn: only rda1 significant
## we will only find the outliers associated with the significant axis. 
vif.cca(CHARR_RDA)
## All under under 5. Low multicollinearity 


ENV_ARRANGED2 = mutate(.data = ENV_ARRANGED,
                       POP_NAME = as.factor(case_when(
                         FID == "1" ~ "T.LGB",
                         FID == '2' ~ 'V.BR',
                         FID == '3' ~ 'V.SIL',
                         FID == '4'~ 'S.PEL',
                         # FID == '5' ~ 'S.PI',
                         FID == '6' ~ 'T.PL',
                         FID == '7' ~ 'T.SB',
                         FID == '8' ~ 'S.LGB',
                         FID == '9' ~ 'G.SB',
                         FID == '10' ~ 'G.PI')))
#FID == '11' ~ 'Mjoavatn',
#FID == '12' ~ 'Fljotaa'))) 
# POP = POP %>% dplyr::pull(POP_name)

ENV_ARRANGED2 = mutate(.data = ENV_ARRANGED2,
                       Population = as.factor(case_when(
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

## Make sure the ID's still line up, if they don't we messed it up
ENV_ARRANGED2$IID == RAW_CHARR_TIB$IID

sum = summary(CHARR_RDA)
sum$species

plot_snps = as_tibble(sum$species)
write_csv(plot_snps, 'SLGBPEL_RDA_SNPS.csv')
plot_snps = read_csv('Galtabol_RDA_SNPS.csv')

plot_loc = as_tibble(sum$sites)
write_csv(plot_loc, 'SLGBPEL_RDA_LOC.csv')
plot_loc = read_csv('Galtabol_RDA_LOC.csv')

plot_bp = as.data.frame(sum$biplot)
write_csv(plot_bp, 'SLGBPEL_RDA_BP.csv')
plot_bp = read_csv('Galtabol_RDA_BP.csv')

ENV_ARRANGED2 = bind_cols(ENV_ARRANGED2, plot_loc)

check = ENV_ARRANGED2 %>% 
  dplyr::select(FID, IID, POP_NAME, RDA1, RDA2)  


## RDA graph #####
rda_col = c('#FF0E0E','#2470FF')
rda_col = c('#2470FF','#FF0E0E','#CC8386')
## There is one weird outlier from the small benthic population in Galtabol, 
## Not sure what's up with that. 

## Galtabol RDA1 = 76.22% RDA2 = 23.78%
## Vatnshlidarvatn RDA1 = 68.78% RDA2 = 31.22%
## Thingvallavatn RDA1 = 84.11% RDA2 = 15.89%
## Svinavatn RDA1 = 86.63% 13.37%
## Svinavatn PLLGB RDA1 = 84.28% RDA2 = 15.72

# Thingvallavatn LGBPL RDA1 = 87.3% RDA2 = 12.7%
## Thingvallavatn PLBS RDA1 = 84.17% RDA2 = 15.83%


plot2 = ggplot(data = ENV_ARRANGED2, 
       aes(x = RDA1, y = RDA2))+
  geom_point(data = plot_snps, 
             aes(x = RDA1, y = RDA2)) +
  geom_point(aes(col = POP_NAME))+
  # scale_color_manual(values = rda_col)+
  labs(x = 'RDA 1 (84.3% variance explained)',
       y = 'RDA 2 (15.7% variance explained)')+
  geom_segment(data = plot_bp, 
               aes(x = 0, y = 0, 
                   xend = RDA1, yend = RDA2),
               color = '#8A8F99', 
               arrow = arrow(length = unit(0.03, 'npc')))+
  geom_text(data = plot_bp, 
            label = c('Carbon', 'Nitrogen'), 
            vjust = -1,
            hjust = 0.5,
            nudge_y = -0.45)+
  # geom_text(data = plot_bp, label = rownames(plot_bp),
  #           vjust = 0, nudge_y = -0.45)+
  # geom_label(data = plot_bp, label = rownames(plot_bp),
  #          vjust = 0, nudge_y = -0.45)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = 'bottom')

plot2
# outlier_check = ENV_ARRANGED2 %>% 
#   select(IID, 
#          Carbon, 
#          Nitrogen, 
#          POP_NAME, 
#          RDA1, 
#          RDA2) %>% 
#   View()

ggsave(plot = last_plot(), 
       'rda_polypops.tiff', 
       width = 20, 
       height = 20, 
       unit = 'cm')
## The ouliter has a completely normal carbon and nitrogen value. 
## The genotype matrix must be whats driving the outlier nature 
## of the individual

SCORES = scores(CHARR_RDA, choices = 1,
                display = 'species')
dim(SCORES)
head(SCORES)

SCORES2 = scores(CHARR_RDA, choices = 2,
                display = 'species')
dim(SCORES2)
head(SCORES2)

SCORES = cbind(SCORES, 
               SCORES2)
head(SCORES)

hist(SCORES[,1], main="Loadings on RDA1")
hist(SCORES2[,1], main="Loadings on RDA2")


## RDA scores fot the SNPs for the first three constrained 

SNP_SCORES = as.data.frame(cbind(SNP = rownames(SCORES), SCORES)) %>% 
  as_tibble()

write_csv(SNP_SCORES, 
          'SLGBPEL_RDA_SNP_scores.csv')
## For the RDA analysis for all polymorphic populations together
## load a separate csv file with all snp scores

#SNP_SCORES = read_csv('SNP_SCORES_Polypops_RDA.csv')
## RDA outliers #####
## FUNCTION FROM THE RDA VIGNETTE
outliers = function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

OUTLIERS1 = outliers(SCORES[,1],4)
length(OUTLIERS1)
## Galtabol
## 326 outliers on the first RDA axis
## Vatnshlidarvatn: 0 outliers associated with 4 standard deviations
## Vatnshlidarvatn 0 outliers associated with 3 standard deviations
## Vatnshlidarvatn: 46 outliers associated with 2 standard deviations
## Thingvallavatn: 260 outliers
## Svinavatn 294 outliers
## Svinavatn pllgb 297 outliers
## Thingvallavatn lgbpl 264 outliers
## Thingvallavatn plsb 260 outliers
## SLGBPEL 294

## 331 outliers on RDA axis 1 for all populations
names(OUTLIERS1)
## names of the outlier loci on the first RDA axis. 

outliers2 = outliers(SCORES2[,1],4)
length(outliers2)
## 228 outliers associated with RDA axis 2 for polymorphic populations. 

CHARR_RDA_OUT = cbind.data.frame(rep(1,
                                     times = length(OUTLIERS1)),
                                 names(OUTLIERS1),
                                 unname(OUTLIERS1))
head(CHARR_RDA_OUT)
CHARR_RDA_OUT = as_tibble(CHARR_RDA_OUT) %>% 
  rename(AXIS = 1, SNP = 2, LOADING = 3)

add_predictors = matrix(nrow=(OUTLIERS1), ncol=2)  # 2 columns for 2 predictors
colnames(add_predictors) = c("Carbon","Nitrogen")

##Correlations of each candidate snp with the environmental 
## predictors. 
n_cand = length(CHARR_RDA_OUT$SNP)
foo <- matrix(nrow=(n_cand), ncol=2)  
colnames(foo) <- c("Carbon_SNP_cor","Nitrogen_SNP_cor")

test = as.data.frame(CHARR_RDA_OUT)

##THe SNP names are not the same after the _
## Get rid of the part after the _ and everything should match up

CHARR_RDA_OUT$SNP = as.character(CHARR_RDA_OUT$SNP)

# head(CHARR_RDA_OUT$SNP)
# gsub("\\_*","",CHARR_RDA_OUT$SNP)
CHARR_RDA_OUT$SNP = substr(CHARR_RDA_OUT$SNP, 1, regexpr("\\_", CHARR_RDA_OUT$SNP)-1)

# CHARR_RDA_OUT$SNP = sed_substitute(CHARR_RDA_OUT$SNP,
#                                    "_.*", 
#                                    "")
# SNP_imputed = as_tibble(SNPS.imp)
# SNP_name = names(SNP_imputed)
# test = gsub("_.*", 
#             "", 
#             SNP_name)
# 
# SNP_imputed = SNP_imputed %>% 
#   rename_all(funs(c(test)))
# 
# CHARR_RDA_OUT = as.data.frame(CHARR_RDA_OUT)
# SNP_imputed = as.data.frame(SNP_imputed)
# # head(SNP_imputed)
# CandN = as.data.frame(CandN)

SNP = as_tibble(SNPS)
SNP_name = names(SNP)

test = substr(SNP_name, 
              1, 
              regexpr("\\_", SNP_name)-1)

# test = gsub("_.*", 
#             "", 
#             SNP_name)

# test = sed_substitute(test, 
#                "Affx.", 
#                "Affx-")

SNP = SNP %>% 
  rename_all(funs(c(test)))

CHARR_RDA_OUT = as.data.frame(CHARR_RDA_OUT)
SNP = as.data.frame(SNP)
# head(SNP_imputed)
CandN = as.data.frame(CandN)


for (i in 1:length(CHARR_RDA_OUT$SNP)) {
  nam <- CHARR_RDA_OUT[i,2]
  snp.gen <- SNP[,nam]
  foo[i,] <- apply(CandN,2,function(x) cor(x,snp.gen))
}

candidates = cbind.data.frame(CHARR_RDA_OUT, foo) %>% 
  as_tibble()
dim(CHARR_RDA_OUT)
head(CHARR_RDA_OUT)
dim(foo)
head(foo)
## When only one axis is significant, 
## skip the step of removing the duplicate snps across the axes. 
## write this file to the directory and call it later
View(candidates)

# write_tsv(candidates, 
#           'Svinavatn_PLLGB_RDA_Outliers.txt', 
#           col_names = TRUE)

write_csv(candidates,
          'GSBPI_RDA_outliers.csv')

# CHARR_RDA_OUT = cbind.data.frame(rep(2,
#                                      times = length(OUT2)),
#                                  names(OUT2),
#                                  unname(OUT2))
# head(CHARR_RDA_OUT)
# CHARR_RDA_OUT2 = as_tibble(CHARR_RDA_OUT) %>% 
#   rename(AXIS = 1,
#          SNP = 2,
#          LOADING = 3)

# CHARR_RDA_OUT = rbind(CHARR_RDA_OUT1, 
#                       CHARR_RDA_OUT2)
# length(CHARR_RDA_OUT$SNP[duplicated(CHARR_RDA_OUT$SNP)])
##Galtabol: 44 duplicates across the 2 RDA axes

## RDA align to MAP ####
CHARR_RDA_OUT = read_csv('GSBPI_RDA_Outliers.csv')
CHARR_RDA_OUT = read_csv('Svinavatn_PLLGB_RDA_Outliers.csv')
CHARR_RDA_OUT = read_csv('polypop_RDA_outliers.csv')
CHARR_RDA_OUT = read_csv('polypop_RDA_outliers_SPELCOMBO.csv')
CHARR_RDA_OUT = read_csv('SLGBPEL_RDA_outliers.csv')

# CHARR_RDA_OUT$SNP = as.character(CHARR_RDA_OUT$SNP)

# CHARR_RDA_OUT$SNP = substr(CHARR_RDA_OUT$SNP, 
#               1, 
#               regexpr("\\-", CHARR_RDA_OUT$SNP)-1)

CHARR_RDA_OUT$SNP = gsub("Affx.", "Affx-", CHARR_RDA_OUT$SNP)
# CHARR_RDA_OUT$SNP = sed_substitute(CHARR_RDA_OUT$SNP,
#                                    "Affx.",
#                                    "Affx-")
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
          'GSBPI_RDA_outliers_mapped.csv')

OUT_RDA_MAPPED = read_csv('GSBPI_RDA_outliers_mapped.csv')

## Define outliers ####
CHROME_GROUP = rep('RDA OUTLIER', 
                   length(OUT_RDA_MAPPED$CHROME3)) %>% 
  as_tibble()

## Joins the outliers and labels them as RDA outlier loci
OUT_RDA_GROUP = bind_cols(OUT_RDA_MAPPED, CHROME_GROUP) %>% 
  as_tibble() %>% 
  rename(CN_out = value)

#Get RDA total map info
SNP_SCORES$SNP = as.character(SNP_SCORES$SNP)
SNP_SCORES$SNP = gsub("Affx.", 
                      "Affx-", 
                      SNP_SCORES$SNP)
# SNP_SCORES$SNP = sed_substitute(SNP_SCORES$SNP, 
#                                 "Affx.", 
#                                 "Affx-")

SNP_SCORES$SNP = substr(SNP_SCORES$SNP,
              1,
              regexpr("\\_", SNP_SCORES$SNP)-1)

# SNP_SCORES$SNP = sed_substitute(SNP_SCORES$SNP,
#                                 "_.*", 
#                                 "")


TOTAL_RDA_MAPPED = merge(MAP3, SNP_SCORES,
                         by.x = 'Marker ID', by.y = 'SNP') %>% 
  as_tibble()

TOTAL_RDA_MAPPED = TOTAL_RDA_MAPPED %>% 
  rename(MARKER_ID = `Marker ID`, 
         PDIST = `Physical position`, 
         GDIST = `Genetic distance`)

TOTAL_RDA_MAPPED$RDA1 = as.numeric(as.character(TOTAL_RDA_MAPPED$RDA1))

write_csv(TOTAL_RDA_MAPPED, 
          'GSBPI_RDA_AllSNPS.csv')

gal = read_csv('GSBPI_RDA_AllSNPS.csv')

RDA_values = CHARR_RDA$CCA$v %>% 
  as_tibble()
write_csv(RDA_values, 
          'GSBPI_RDA_values.csv')
## RDA ggMAN #####
TOTAL_RDA_MAPPED = read_csv('Svinavatn_PLLGB_AllSNPS.csv') %>% 
  rename(MARKER_ID = 1, PDIST = 4)
RDA_Outliers = read_csv('Svinavatn_PLLGB_RDA_outliers_mapped.csv')
KeyStone_cand_list = read_csv('RDA_IsoMorpho_Outliers_SNPList.csv') %>% 
  as.data.frame()

poly_pop = read_csv('Polypop_RDA_outliers_mapped.csv') %>% 
  rename(ID = `Marker ID`)

inner_join(poly_pop, KeyStone_cand_list, by = 'ID')

# str(TOTAL_RDA_MAPPED)
# 
TOTAL_RDA_MAPPED$MARKER_ID
sum(is.na(TOTAL_RDA_MAPPED$PDIST))
NO_NAs = na.omit(TOTAL_RDA_MAPPED)

rda_man_plot = ggman(TOTAL_RDA_MAPPED, 
                     chrom = 'CHROME3', 
                     pvalue = 'RDA1',
                     snp = 'MARKER_ID',
                     bp = 'PDIST',
                     logTransform = F,
                     pointSize = 2, 
                     title = 'RDA - Isotope - SLGBPL',
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

plot2/plot

ggsave(plot = last_plot(), 
       'SLGBPL_RDA_Isotope_Keystone_outliers.tiff', 
       height = 20,
       width = 20, 
       unit = 'cm')

LFMM_carbon_man/LFMM_nitrogen_man/plot2/plot
plot2/(LFMM_carbon_man|LFMM_nitrogen_man|plot)


ggsave(plot = last_plot(), 
       'Svinavatn_PLLGB_Outlier_stacked.tiff', 
       height = 25, 
       width = 25, 
       units = 'cm')

## LFMM RDA Overlap New ######
# LFMM = read_csv('LFMM_Galtabol_outliers.csv') %>%
#   select(-outlier,
#          -nitrogen_only_outlier,
#          -cn_shared_outliers,
#          -carbon_only_outlier)

setwd('~/PhD/SNP Demographic modelling/Outliers_directory')
LFMM_carbon = read_csv('LFMM_Thing_LGBPL_Carbon_Outlier_loci_FINAL.csv')

## NO RDA OUTLIERS IN VATNSHLIDARVATN!!!!!!!
RDA = read_csv('Thingvallavatn_LGBPL_RDA_outliers_mapped.csv')

## overlap between rda and carbon outliers
LFMMCarbon_RDA = intersect(RDA$`Marker ID`,
                    LFMM_carbon$`Marker ID`)
length(LFMMCarbon_RDA)

LFMMCarbon_RDA = inner_join(LFMM_carbon, 
          RDA, 
          by = 'Marker ID') %>% 
  select(`Marker ID`, 
         CHROME3.x) %>% 
  arrange(CHROME3.x) %>% 
  rename(Chromosome = CHROME3.x)


label = rep('LFMM (Carbon) - RDA overlap', 
            length(LFMMCarbon_RDA$`Marker ID`)) %>% 
  as_tibble()

LFMMCarbon_RDA = bind_cols(LFMMCarbon_RDA, 
                           label) %>% 
  rename(Methods = value)

## overlap between rda and nitrogen outliers
LFMM_Nitrogen = read_csv('LFMM_Thing_LGBPL_Nitrogen_Outlier_loci_FINAL.csv')

LFMM_Nitrogen_RDA = intersect(RDA$`Marker ID`,
                    LFMM_Nitrogen$`Marker ID`)
length(LFMM_Nitrogen_RDA)

LFMM_Nitrogen_RDA = inner_join(LFMM_Nitrogen, 
                            RDA, 
                            by = 'Marker ID') %>% 
  select(`Marker ID`, 
         CHROME3.x) %>% 
  arrange(CHROME3.x) %>% 
  rename(Chromosome = CHROME3.x)


label = rep('LFMM (Nitrogen) - RDA overlap', 
            length(LFMM_Nitrogen_RDA$`Marker ID`)) %>% 
  as_tibble()

LFMM_Nitrogen_RDA = bind_cols(LFMM_Nitrogen_RDA, 
                           label) %>% 
  rename(Methods = value)


outliers = bind_rows(LFMMCarbon_RDA, 
                     LFMM_Nitrogen_RDA) %>% 
  arrange(Chromosome)


# inner_join(LFMM_carbon, 
#            LFMM_Nitrogen, 
#            by = 'Marker ID')

# ## LFMM CN outlier overlap 
# 
LFMM_CN = inner_join(LFMM_carbon,
                     LFMM_Nitrogen,
                     by = 'Marker ID')

## Vatnshlidarvatn only!!!!!!
LFMM_CN = inner_join(LFMM_carbon,
                     LFMM_Nitrogen,
                     by = 'Marker ID') %>%
  select(`Marker ID`,
         CHROME3.x) %>%
  arrange(CHROME3.x) %>%
  rename(Chromosome = CHROME3.x)

View(LFMM_CN)
label = rep('LFMM Carbon & Nitrogen overlap',
            length(LFMM_CN$`Marker ID`)) %>%
  as_tibble()

LFMM_CN = bind_cols(LFMM_CN,
                       label) %>%
  rename(Methods = value)

label = rep('Thingvallavatn - PLSB',
            length(LFMM_CN$`Marker ID`)) %>%
  as_tibble()
Vatnshlidarvatn_outliers = bind_cols(LFMM_CN,
                                     label) %>%
  rename(Population = value) %>%
  select(Population,
         `Marker ID`,
         Chromosome,
         Methods) %>%
  distinct()
# 

## Other populations that are not Vantshlidarvatn !!!!
LFMMCN_RDA = inner_join(LFMM_CN,
                              RDA,
                              by = 'Marker ID') %>%
  select(`Marker ID`,
         CHROME3.x) %>%
  arrange(CHROME3.x) %>%
  rename(Chromosome = CHROME3.x)


label = rep('LFMM (C & N) - RDA overlap',
            length(LFMMCN_RDA$`Marker ID`)) %>%
  as_tibble()

LFMMCN_RDA = bind_cols(LFMMCN_RDA,
                             label) %>%
  rename(Methods = value)

outliers = bind_rows(outliers, 
                              LFMMCN_RDA) %>% 
  arrange(Chromosome)

label = rep('Thingvallavatn - LGBPL', 
            length(outliers$`Marker ID`)) %>% 
  as_tibble()



## Other populations that are not Vatnshlidarvatn!!!!!!
Thingvallavatn_LGBPL_outliers = bind_cols(outliers, 
                                         label) %>% 
  rename(Population = value) %>% 
  select(Population, 
         `Marker ID`, 
         Chromosome, 
         Methods) %>% 
  distinct()

Thingvallavatn_PLSB_outliers = bind_cols(outliers, 
                                    label) %>% 
  rename(Population = value) %>% 
  select(Population, 
         `Marker ID`, 
         Chromosome, 
         Methods) %>% 
  distinct()

Thingvallavatn_outliers = bind_cols(outliers, 
                              label) %>% 
  rename(Population = value) %>% 
  select(Population, 
         `Marker ID`, 
         Chromosome, 
         Methods) %>% 
  distinct()

Svinavatn_outliers = bind_cols(outliers, 
                              label) %>% 
  rename(Population = value) %>% 
  select(Population, 
         `Marker ID`, 
         Chromosome, 
         Methods) %>% 
  distinct()

Galtabol_outliers = bind_cols(outliers, 
                              label) %>% 
  rename(Population = value) %>% 
  select(Population, 
         `Marker ID`, 
         Chromosome, 
         Methods) %>% 
  distinct()

View(Thingvallavatn_outliers)

tab_df(Thingvallavatn_LGBPL_outliers, 
       file = 'Thingvallavatn_LGBPL_RDA_LFMMCarbon_overlap.doc')


G_S_pop = bind_rows(Galtabol_outliers,
                    Svinavatn_outliers)

G_S_T = bind_rows(G_S_pop,
                  Thingvallavatn_outliers)
G_S_T_V = bind_rows(G_S_T,
                    Vatnshlidarvatn_outliers)

tab_df(G_S_T_V, 
       file = 'FINAL_Overlap_outlier_detection_methods_AllPops.doc')
## VennDiagram ####
PCA_outliers = read_csv('PolyPopns_PCAdapt_Outliers_k9_Q0.01.csv') %>% 
  na.omit()
LFMM_Nitrogen = read_csv('LFMM_Vatnshlidarvatn_Nitrogen_Outlier_loci_FINAL.csv')
LFMM_carbon = read_csv('LFMM_Vatnshlidarvatn_Carbon_Outlier_loci_FINAL.csv')
RDA = read_csv('Thingvallavatn_RDA_outliers_mapped.csv')


# PCA_outliers = read.delim('Thingvallavatn_PCA_OUTLIERS_k2.txt', 
#                           stringsAsFactors = F) %>% 
#   as.tibble()
venn_pal = c('#2776CC',
             '#CC2B1C')

venn_pal = c('#2776CC',
             '#CC2B1C',
             '#10B372')

venn_pal = c('#2776CC',
             '#CC2B1C',
             '#10B372',
             '#B300A2')

library(VennDiagram)
#B300A2
venn.diagram(
  x = list(#RDA$`Marker ID`,
           #PCA_outliers$MARKER_ID),
           LFMM_carbon$`Marker ID`,
           LFMM_Nitrogen$`Marker ID`),
  category.names = c(#'RDA',
                     #'PCAdapt',
                     'LFMM Carbon',
                     'LFMM Nitrogen'),
  filename = 'Vatnshlidarvatn_VennDiagram_outlier_loci.tiff',
  imagetype = 'tiff',
  # height = 400,
  # width = 400, 
  # resolution = 300,
  # compression = 'lzw',
  
  # cirlces
  lwd = 2,
  lty = 'blank',
  fill = venn_pal,
  
  # Numbers
  cex = .8,
  fontface = "bold",
  fontfamily = "sans")
  
  # Set names
  # cat.cex = 0.6,
  # cat.fontface = "bold",
  # cat.default.pos = "outer",
  # cat.pos = c(-27, 27, 135),
  # cat.dist = c(0.055, 0.055, 0.085),
  # cat.fontfamily = "sans",
  # rotation = 1)



test = inner_join(RDA, 
           PCA_outliers, 
           by = 'MARKER_ID') %>% 
  select(MARKER_ID)

write_csv(test, 
          'Polypops_RDA_PCAdapt_overlap.csv')
## OUtlier comparison####

PCA_outliers = read.delim('Galtabol_PCA_OUTLIERS_k2.txt', 
                          stringsAsFactors = F) %>% 
  as.tibble()


RDA_OUTLIER_NODUP = read_csv('Galtabol_RDA_OUTLIERS_NODUP_MARKERS.csv') %>% 
  arrange(SNP) %>% 
  rename(MARKER_ID = SNP)
RDA_OUTLIER_NODUP$MARKER_ID = sed_substitute(RDA_OUTLIER_NODUP$MARKER_ID, 
                                             "Affx.", 
                                             "Affx-")
RDA_OUTLIER_NODUP$MARKER_ID = sed_substitute(RDA_OUTLIER_NODUP$MARKER_ID,
                                             "_.*", 
                                             "")
LFMM_data = read.delim('LFMM_SNP_outliers', 
                       stringsAsFactors = F) %>% 
  as.tibble() %>% 
  select(-value, -value2)


## RDA PCA OVERLAP ####
RDA_PCA = intersect(RDA_OUTLIER_NODUP$MARKER_ID, 
                    PCA_outliers$MARKER_ID)
length(RDA_PCA)

PCA_RDA = intersect(PCA_outliers$MARKER_ID,
                    RDA_OUTLIER_NODUP$MARKER_ID)
length(PCA_RDA)

## LFMM PCA OVERLAP #####
CN_together = LFMM_data %>% 
  filter(OVERLAP == 'CN OUTLIERS')

CN_PCA = intersect(CN_together$Marker.ID, PCA_outliers$MARKER_ID)
length(CN_PCA)

C_only = LFMM_data %>% 
  filter(OVERLAP == 'CARBON OUTLIER ONLY')

C_PCA = intersect(C_only$Marker.ID, PCA_outliers$MARKER_ID)
length(C_PCA)

N_only = LFMM_data %>% 
  filter(OVERLAP == 'NITROGEN OUTLIER ONLY')

N_PCA = intersect(N_only$Marker.ID, PCA_outliers$MARKER_ID)
length(N_PCA)


## LFMM RDA OVERLAP ####
CN_RDA = intersect(CN_together$Marker.ID, RDA_OUTLIER_NODUP$MARKER_ID)
length(CN_RDA)

C_RDA = intersect(C_only$Marker.ID, RDA_OUTLIER_NODUP$MARKER_ID)
length(C_RDA)

N_RDA = intersect(N_only$Marker.ID, RDA_OUTLIER_NODUP$MARKER_ID)
length(N_RDA)
