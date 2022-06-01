##############################
## Comparing model likelihoods
##
## Matt Brachmann (PhDMattyB)
##
## 2019-10-01
##
##############################

setwd('C:/Users/Matt Brachmann/Documents/Fsc26/Galtabol/Model_lhoods')
setwd('C:/Users/Matt Brachmann/Documents/Fsc26/Vatnshlidarvatn/Model_lhoods')
setwd('C:/Users/Matt Brachmann/Documents/Fsc26/Thing_files/model_lhoods')
setwd('C:/Users/Matt Brachmann/Documents/Fsc26/Svin_files/')

library(tidyverse)
library(patchwork)
# Make lhoods data frams --------------------------------------------------

model = read_tsv('S_BP_CIM.lhoods', 
                 col_names = FALSE) %>%
  select(-X1)

## Make a label for each model that we are comparing
label = rep('S_CIM', length(model$X2)) %>%
  as_tibble()

## Combine the columns from the original data frame
## with the label we made for the model
model = bind_cols(model, label) %>%
  dplyr::rename(OG_model = value)

## Write the new data frame to a file and re-read in the data above
write_tsv(model, 'S_CIM_colname.lhoods')


# Read in newly created data from above -----------------------------------


## Read in the data sets we have created. 
IM = read_tsv("G_IM3_100000.lhoods") %>% 
  dplyr::rename(lhoods = 1, 
                OG_model = 2)

SC = read_tsv("G_SC.lhoods") %>% 
  dplyr::rename(lhoods = 1, 
                OG_model = 2)
# IM_MT = read_tsv("T_IM_MT_colname.LHOODS") %>% 
#   dplyr::rename(lhoods = 1, 
                # OG_model = 2)
# SC_MT = read_tsv("T_SC_MT_colname.LHOODS") %>% 
#   dplyr::rename(lhoods = 1, 
#                 OG_model = 2)

CIM = read_tsv('G_IM_Change.lhoods') %>% 
  dplyr::rename(lhoods = 1, 
                OG_model = 2)

# IM_change = read_tsv('GIM_Change.lhoods') 
# SC = read_tsv('GSC_colname.lhoods') %>% 
#   rename(lhoods = 1, 
#          OG_model = 2)

## This section allows you to create the data frames 
## that we will be using to graph


## Combine the data from all of the models that we have run
models = bind_rows(IM,
                   CIM, 
                   SC)
models = bind_rows(IM, SC, IM_MT, SC_MT)

## Create a population label as we will be dealing with models 
## from multiple different populations
popn_label = rep('Galtabol', 
                 length(models$lhoods)) %>% 
  as_tibble()

## Make the data set we will be using and order the data properly
models = bind_cols(models, 
                   popn_label) %>% 
  rename(Population = value) %>% 
  dplyr::select(Population, 
                OG_model, 
                lhoods)
## Create a new column with the full model name
models = mutate(.data = models, 
                Model = as.factor(case_when(
                  OG_model == 'SC' ~ 'Secondary contact',
                  OG_model == 'IM' ~ 'Isolation by migration',
                  OG_model == "CIM" ~ 'Change in migration')))

models = mutate(.data = models, 
                Model = as.factor(case_when(
                  OG_model == 'T_IM' ~ 'Isotlation by migration',
                  OG_model == 'T_SC' ~ 'Secondary contact',
                  OG_model == "T_IM_MT" ~ 'Isolation by migration (MT)', 
                  OG_model == 'T_SC_MT' ~ 'Secondary contact (MT)')))

# View(models)
write_tsv(models, 'Svinavatn_model_lhoods.txt')
# Likelihood visualization ------------------------------------------------
setwd('~/PhD/SNP Demographic modelling/Fsc26')

Galtabol = read_tsv('Galtabol_model_lhoods.txt')
Vatnshlidarvatn = read_tsv('Vatnshlidarvatn_model_lhoods.txt')
Thingvallavatn = read_tsv('Thingvallavatn_model_lhoods.txt')
Svinavatn = read_tsv('Svinavatn_model_lhoods.txt')

models = bind_rows(Galtabol,
                   Svinavatn,
                   Thingvallavatn,
                   Vatnshlidarvatn)

# View(models)
models = mutate(.data = models, 
                New_model = as.factor(case_when(
                  OG_model == 'IM' ~ '2IM',
                  OG_model == 'GIM_Change' ~ 'CIM',
                  OG_model == "GSC" ~ '2SC', 
                  OG_model == 'V_IM' ~ '2IM',
                  OG_model == 'V_IM_Change' ~ 'CIM',
                  OG_model == 'V_SC' ~ '2SC',
                  OG_model == 'T_IM' ~ '3IM', 
                  OG_model == 'T_SC' ~ '3SC', 
                  OG_model == 'T_IM_MT' ~ '3IM_MT', 
                  OG_model == 'T_SC_MT' ~ '3SC_MT', 
                  OG_model == 'S_IM' ~ '2IM',
                  OG_model == 'S_SC' ~ '2SC', 
                  OG_model == 'S_CIM' ~ 'CIM')))

## Colour palettes for the graph
mod_pal = c("#258CE6", "#597499", '#FF6114', '#FF1C08')
mod_pal = c('#258CE6', '#FF1C08')

## order the models so the graph makes sense
## all of the IM models first and then the secondary contact models
models = models %>%
  mutate(New_model = fct_relevel(New_model,
                             '2IM', 
                             'CIM',
                             '3IM_MT',
                             '2SC', 
                             '3SC_MT'))
# View(models)

## Actual graph using the data and palettes that we have created
theme_set(theme_bw())

Galtabol_graph = models %>% 
  filter(Population == 'Galtabol') %>% 
  ggplot(aes(x = New_model, 
             y = lhoods))+
  geom_boxplot(aes(fill = Model), 
               col = 'Black')+
  scale_fill_manual(values = mod_pal)+
  labs(y = 'Maximum likelihood', 
       title = 'A)')+
  facet_grid(~Population)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 45,
        #                            vjust = 0.9,
        #                            hjust = 0.9),
        axis.title.x = element_blank(),
        legend.position = 'none',
        strip.background = element_rect(fill = 'light grey'),
        strip.text = element_text(face = 'bold',
                                  size = 12))

Svinavatn_graph = models %>% 
  filter(Population == 'Svinavatn') %>% 
  ggplot(aes(x = New_model, 
             y = lhoods))+
  geom_boxplot(aes(fill = Model), 
               col = 'Black')+
  scale_fill_manual(values = mod_pal)+
  labs(y = 'Maximum likelihood', 
       title = 'B)')+
  facet_grid(~Population)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 45,
        #                            vjust = 0.9,
        #                            hjust = 0.9),
        axis.title.x = element_blank(),
        legend.position = 'none',
        strip.background = element_rect(fill = 'light grey'),
        strip.text = element_text(face = 'bold',
                                  size = 12))

Thingvallavatn_graph = models %>% 
  filter(Population == 'Thingvallavatn') %>% 
  ggplot(aes(x = New_model, 
             y = lhoods))+
  geom_boxplot(aes(fill = Model), 
               col = 'Black')+
  scale_fill_manual(values = mod_pal)+
  labs(y = 'Maximum likelihood', 
       title = 'C)')+
  facet_grid(~Population)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45,
                                   vjust = 0.9,
                                   hjust = 0.9),
        axis.title.x = element_blank(),
        legend.position = 'none',
        strip.background = element_rect(fill = 'light grey'),
        strip.text = element_text(face = 'bold',
                                  size = 12))

Vatnshlidarvatn_graph = models %>% 
  filter(Population == 'Vatnshlidarvatn') %>% 
  ggplot(aes(x = New_model, 
             y = lhoods))+
  geom_boxplot(aes(fill = Model), 
               col = 'Black')+
  scale_fill_manual(values = mod_pal)+
  labs(y = 'Maximum likelihood', 
       title = 'D)')+
  facet_grid(~Population)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45,
                                   vjust = 0.9,
                                   hjust = 0.9),
        axis.title.x = element_blank(),
        legend.position = 'none',
        strip.background = element_rect(fill = 'light grey'),
        strip.text = element_text(face = 'bold',
                                  size = 12))

combo = (Galtabol_graph + Svinavatn_graph)/(Thingvallavatn_graph + Vatnshlidarvatn_graph)

## Save graph to a file 
ggsave('Fsc26_model_output.tiff', 
       plot = combo, 
       width = 20, 
       height = 20, 
       units = 'cm')

## Read in the data that has the AIC values for each model 
## so that we can compare them
aic = read_tsv("Model_AIC_values.txt") %>% 
  dplyr::select(-DeltaL) %>% 
  rename(DeltaL = 3)




# Galtabol aic ------------------------------------------------------------
setwd('C:/Users/Matt Brachmann/Documents/Fsc26/Galtabol/Model_lhoods')
# setwd('C:/Users/Matt Brachmann/Documents/Fsc26/Vatnshlidarvatn/Model_lhoods')
# setwd('C:/Users/Matt Brachmann/Documents/Fsc26/Thing_files/model_lhoods')
# setwd('C:/Users/Matt Brachmann/Documents/Fsc26/Svin_files/')

GIM_aic = read_tsv('G_IM3_100000.aic')
im_lab = rep('Isolation by migration', 
             length(GIM_aic$AIC)) %>% 
  as_tibble()
GIM_aic = bind_cols(GIM_aic, 
                    im_lab)

GIMC_aic = read_tsv('G_IM_Change.aic')
imc_lab = rep('Change in migration', 
              length(GIMC_aic$AIC)) %>% 
  as_tibble()
GIMC_aic = bind_cols(GIMC_aic, 
                     imc_lab)

GSC_aic = read_tsv('GSC2.aic')
sc_lab = rep('Secondary contact', 
             length(GSC_aic$AIC)) %>% 
  as_tibble()
GSC_aic = bind_cols(GSC_aic, 
                    sc_lab)

AIC = bind_rows(GIM_aic, 
                GIMC_aic, 
                GSC_aic) %>% 
  rename(model = value)

AIC_lab = rep('Galtabol', 
              length(models$model)) %>% 
  as_tibble()

AIC = bind_cols(models, 
                mod_lab) %>% 
  rename(Population = value)

AIC$AIC

AIC %>% 
  write_csv('Galtabol_Model_AIC.csv')
# Galtabol plot ----------------------------------------------------------------
galtabol = models %>% 
  filter(Population == 'Galtabol')

galtabol = mutate(.data = galtabol, 
                AIC = as.factor(case_when(
                  New_model == '2SC' ~ 'AIC = 10655.2',
                  New_model == '2IM' ~ 'AIC = 10779.2',
                  New_model == 'CIM' ~ 'AIC = 9891.8')))

galtabol = mutate(.data = galtabol, 
                  AIC = as.factor(case_when(
                    New_model == '2SC' ~ 'Delta-AIC = 763.4',
                    New_model == '2IM' ~ 'Delta-AIC = 887.4',
                    New_model == 'CIM' ~ 'Delta-AIC = 0')))

prot2 = distinct(galtabol, New_model) %>%
  arrange(New_model)

prot2$yloc = max(galtabol$lhoods) + 45
prot2$label = c("AIC = 10779.2",
                "AIC = 9891.8",
                "AIC = 10655.2")

prot3 = distinct(galtabol, New_model) %>%
  arrange(New_model)

prot3$yloc = max(galtabol$lhoods) + 30
prot3$label = c("Delta-AIC = 887.4",
                "Delta-AIC = 0",
                "Delta-AIC = 763.4")

Gal_model = ggplot(data = galtabol,
       aes(x = New_model,
           y = lhoods))+
  geom_boxplot(aes(fill = New_model), 
               col = 'Black')+
  geom_label(aes(label = AIC))+
  geom_text(data = prot2,
            aes(y = yloc,
                label = label),
            position = position_dodge(width = 0.75),
            size = 5)+
  geom_text(data = prot3,
            aes(y = yloc,
                label = label),
            position = position_dodge(width = 0.75),
            size = 5)+
  scale_fill_manual(values = mod_pal)+
  labs(y = 'Maximum likelihood', 
       title = 'A)')+
  # facet_grid(~Population)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, 
                                   hjust = 0.9, 
                                   size = 15),
        axis.title.x = element_blank(), 
        legend.position = 'none', 
        strip.background = element_rect(fill = 'light grey'), 
        strip.text = element_text(face = 'bold', 
                                  size = 12))
Gal

## Save graph to a file 
ggsave('Galtabol_fsc26.tiff', 
       plot = Gal_model)





# svinavatn aic -----------------------------------------------------------

setwd('C:/Users/Matt Brachmann/Documents/Fsc26/Svin_files/')
SIM_aic = read_tsv('')
im_lab = rep('Isolation by migration', 
             length(GIM_aic$AIC)) %>% 
  as_tibble()
GIM_aic = bind_cols(GIM_aic, 
                    im_lab)

GIMC_aic = read_tsv('G_IM_Change.aic')
imc_lab = rep('Change in migration', 
              length(GIMC_aic$AIC)) %>% 
  as_tibble()
GIMC_aic = bind_cols(GIMC_aic, 
                     imc_lab)

GSC_aic = read_tsv('GSC2.aic')
sc_lab = rep('Secondary contact', 
             length(GSC_aic$AIC)) %>% 
  as_tibble()
GSC_aic = bind_cols(GSC_aic, 
                    sc_lab)

AIC = bind_rows(GIM_aic, 
                GIMC_aic, 
                GSC_aic) %>% 
  rename(model = value)

AIC_lab = rep('Galtabol', 
              length(models$model)) %>% 
  as_tibble()

AIC = bind_cols(models, 
                mod_lab) %>% 
  rename(Population = value)

AIC$AIC

AIC %>% 
  write_csv('Galtabol_Model_AIC.csv')

# Vatnshlidarvatn AIC -----------------------------------------------------
setwd('~/PhD/SNP Demographic modelling/Fsc26/Vatnshlidarvatn')

VIM_aic = read_tsv('V_IM.aic')
im_lab = rep('Isolation by migration', 
             length(VIM_aic$AIC)) %>% 
  as_tibble()
VIM_aic = bind_cols(VIM_aic, 
                    im_lab)

VIMC_aic = read_tsv('V_IM_Change.aic')
imc_lab = rep('Change in migration', 
              length(VIMC_aic$AIC)) %>% 
  as_tibble()
VIMC_aic = bind_cols(VIMC_aic, 
                     imc_lab)

VSC_aic = read_tsv('V_SC.aic')
sc_lab = rep('Secondary contact', 
             length(VSC_aic$AIC)) %>% 
  as_tibble()
VSC_aic = bind_cols(VSC_aic, 
                    sc_lab)

AIC = bind_rows(VIM_aic, 
                VIMC_aic, 
                VSC_aic) %>% 
  rename(model = value)

AIC_lab = rep('Vatnshlidarvatn', 
              length(AIC$model)) %>% 
  as_tibble()

AIC = bind_cols(AIC, 
                AIC_lab) %>% 
  rename(Population = value)

AIC %>% 
  write_csv('Vatnshlidarvatn_model_AIC.csv')
# Vatnshlidarvatn ---------------------------------------------------------
vatnshlidarvatn = models %>% 
  filter(Population == 'Vatnshlidarvatn')

vatn2 = distinct(vatnshlidarvatn, New_model) %>%
  arrange(New_model)

vatn2$yloc = max(vatnshlidarvatn$lhoods) + 5
vatn2$label = c("AIC = 1440.8", 
                "AIC = 1450.7",
                "AIC = 1445.1")

vatn3 = distinct(vatnshlidarvatn, New_model) %>%
  arrange(New_model)

vatn3$yloc = max(vatnshlidarvatn$lhoods) + 4.5
vatn3$label = c("Delta-AIC = 0", 
                "Delta-AIC = 9.9",
                "Delta-AIC = 4.3")


vatn_model = ggplot(data = vatnshlidarvatn,
       aes(x = New_model,
           y = lhoods))+
  geom_boxplot(aes(fill = New_model), 
               col = 'Black')+
  geom_text(data = vatn2,
            aes(y = yloc,
                label = label),
            position = position_dodge(width = 0.75),
            size = 5)+
  geom_text(data = vatn3,
            aes(y = yloc,
                label = label),
            position = position_dodge(width = 0.75),
            size = 5)+
  scale_fill_manual(values = mod_pal)+
  labs(y = 'Maximum likelihood', 
       title = 'D)')+
  # facet_grid(~Population)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        axis.text.x = element_text(size = 15, 
                                   angle = 45, 
                                   vjust = 0.9, 
                                   hjust = 0.9),
        axis.title.x = element_blank(), 
        legend.position = 'none', 
        strip.background = element_rect(fill = 'light grey'), 
        strip.text = element_text(face = 'bold', 
                                  size = 12))

vatn_model
## Save graph to a file 
ggsave('Vatnshlidarvatn_fsc26_AIClables.tiff', 
       plot = vatn_model)



# Thingvallavatn AIC ------------------------------------------------------
Thingvallavatn = models %>% 
  filter(Population == 'Thingvallavatn')


# Thingvallavatn = read_tsv('Thingvallavatn_model_lhoods.txt') 
# View(Thingvallavatn)
# thing_aic = Thingvallavatn %>% 
#   distinct(New_model) %>% 
#   arrange(New_model)

thing_aic = distinct(Thingvallavatn, New_model) %>%
  arrange(New_model)

thing_aic$yloc = max(Thingvallavatn$lhoods) + 100
thing_aic$label = c("AIC = 5254.7",
                "AIC = 5329.6",
                "AIC = 3080.2",
                "AIC = 5358.5")
thing_aic = thing_aic %>% 
  rename(AIC = label) %>% 
  select(-yloc)

thing_aic2 = distinct(Thingvallavatn, New_model) %>%
  arrange(New_model)

thing_aic2$yloc = max(Thingvallavatn$lhoods) + 60
thing_aic2$label = c("Delta-AIC = 2174.5",
                    "Delta-AIC = 2249.4",
                    "Delta-AIC = 0",
                    "Delta-AIC = 2278.3")
thing_aic2 = thing_aic2 %>% 
  rename(DeltaL = label)

inner_join(thing_aic, 
           thing_aic2)

Thing_models = ggplot(data = Thingvallavatn,
       aes(x = New_model,
           y = lhoods))+
  geom_boxplot(aes(fill = New_model), 
               col = 'Black')+
  geom_text(data = thing_aic,
            aes(y = yloc,
                label = label),
            position = position_dodge(width = 0.75),
            size = 5)+
  geom_text(data = thing_aic2,
            aes(y = yloc,
                label = label),
            position = position_dodge(width = 0.75),
            size = 5)+
  scale_fill_manual(values = mod_pal)+
  labs(y = 'Maximum likelihood', 
       title = 'C)')+
  ylim(-1200, -500)+
  # facet_grid(~Population)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        axis.text.x = element_text(size = 15, 
                                   angle = 45, 
                                   vjust = 0.9, 
                                   hjust = 0.9),
        axis.title.x = element_blank(), 
        legend.position = 'none', 
        strip.background = element_rect(fill = 'light grey'), 
        strip.text = element_text(face = 'bold', 
                                  size = 12))

Thing_models
## Save graph to a file 
ggsave('Thingvallavatn_fsc26_AIClables.tiff', 
       plot = Thing_models)




# combo save --------------------------------------------------------------
library(patchwork)

