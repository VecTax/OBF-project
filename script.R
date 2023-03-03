rm(list = ls())

##############################################################################################################
## Packages  --------------------------------------------------

library(stringr)
library(doBy)
library(readxl)
library(tidyverse)
library(dplyr)
library(car)
library(lme4)
library(lmerTest)
library(glmm)
library(jtools)
library(broom)
library(broom.mixed)
library(DHARMa)
library(glmmTMB)
library(reshape2)
library(lsmeans)
library(ggeffects)
library(FactoMineR)
library(factoextra)

##############################################################################################################
## Data import  --------------------------------------------------

nomenclature = read.csv("OFB_total_2007a2020_taxons_cites_retenus_test.csv", sep = ";")
releves = read.csv("se_ofb_pa.csv", encoding = "latin1", sep = ";", dec = ",")
stations_suppr = read.csv("Tab_suivi_CSE_2016_2020_majfin_20210318.csv")
dates_correction = read.table("releves_double_dateOK.txt", sep = "\t", header = T)
variables_09_simp = read.csv("Tab_infoeco_CSE_2009_2013_majOKsaufSol_20210318.csv", sep = ";", fileEncoding="latin1")
variables_16_simp = read.csv("Tab_infoeco_CSE_2016_2020_maj_20210318.csv", sep = ";", fileEncoding="latin1")
baseflor = read.csv2("Liste_Bourgogne_v20201203_Taxref12.csv", sep = ";", skip = 1, fileEncoding="latin1")
baseflor_correction = read.csv("indices_inconnus(Recuperation automatique).csv", sep = ";", dec =",")
listsp_prot = read.csv("liste espece protegee.csv", sep = ";")

##############################################################################################################
## Station's information formatting  --------------------------------------------------

infos_stations = releves[1:15, ]
infos_stations = subset(infos_stations, select = -c(X, X.1, X.2))
infos_stations = as.data.frame(t(infos_stations))

## Renaming
names(infos_stations) = c("id_station", "commune", "code_INSEE", "num_perso_releve", "date", "observateur", "id_releve",
                          "commentaire", "recouvrement_arbore", "recouvrement_arbustif", "recouvrement_herbace",
                          "recouvrement_muscinal", "recouvrement_tot", "nb_taxons", "aucun_taxon")
infos_stations = infos_stations[-1, ]

## Numerical transformation
infos_stations$recouvrement_arbore = as.numeric(infos_stations$recouvrement_arbore)
infos_stations$recouvrement_arbustif = as.numeric(infos_stations$recouvrement_arbustif)
infos_stations$recouvrement_herbace = as.numeric(infos_stations$recouvrement_herbace)
infos_stations$recouvrement_muscinal = as.numeric(infos_stations$recouvrement_muscinal)
infos_stations$recouvrement_tot = as.numeric(infos_stations$recouvrement_tot)
infos_stations$nb_taxons = as.numeric(infos_stations$nb_taxons)
infos_stations$station = row.names(infos_stations)
row.names(infos_stations) = seq_len(nrow(infos_stations))

## Keep the unique number of station
infos_stations$station = gsub("^(.*).\\d+$", "\\1", perl = T, infos_stations$station)

## Wrongly recorded dates
dates_correction$Date2 = format(as.Date(dates_correction$Date2, "%d/%m/%Y"), "%m/%d/%Y")
dates_correction$date.1 = format(as.Date(dates_correction$date.1, "%d/%m/%Y"), "%m/%d/%Y")
for (i in dates_correction$ID.releve2.1) {
  infos_stations$date[infos_stations$id_releve == i] = dates_correction$Date2[dates_correction$ID.releve2.1 == i]
}

## Adding date and campaign variables
infos_stations$annee = as.numeric(gsub("^.*/(\\d{4})$", "\\1", perl = T, infos_stations$date))
infos_stations$campagne = ifelse(infos_stations$annee %in% 2009:2013, "2009-2013", "2016-2020")

## Season variable
infos_stations$newdate = as.Date(infos_stations$date, "%m/%d/%Y")
infos_stations$mois = strftime(infos_stations$newdate, "%m")
src0 = list(c("04","05","06"), c("07","08","09","10","11"))
tgt0 = list("printemps","ete")
infos_stations$saison_rel = recodeVar(infos_stations$mois, src=src0, tgt=tgt0)

##############################################################################################################
## Formatting survey file --------------------------------------------------

releves_flore = releves[- (1:15), ]
names(releves_flore)[5:ncol(releves_flore)] = releves[7, 5:ncol(releves)]
releves_flore = renameCol(releves_flore, c("Numero.perso.Station", "X", "X.1", "X.2"),
                          c("especes_non_ok", "cd_nom_terrain", "strate", "indigenat"))

## Switching to long format
releves_flore = reshape(releves_flore,
                        varying = 5:ncol(releves_flore),
                        idvar = 1:4, direction = "long", timevar = "id_releve",
                        v.names = "abondance",
                        times = names(releves_flore)[5:ncol(releves_flore)],
                        new.row.names = 1:(nrow(releves_flore)
                                           * (ncol(releves_flore) - 4)))

## Deleting arables species
releves_flore$indigenat = as.factor(releves_flore$indigenat)
releves_flore = releves_flore %>% 
  filter(!(indigenat %in% c("Cultivee", "Plantation certaine", "Plantation probable")))

## Adding dates and station informations
releves_flore = merge(releves_flore,
                      infos_stations[, c("id_station", "date", "id_releve", "campagne", "station", "saison_rel")],
                      by = "id_releve", all.x = T)

## Wrongly recorded occurences
releves_flore$abondance[releves_flore$station == "EM95_A" & releves_flore$cd_nom_terrain == 81457 &
                          releves_flore$date == "5/28/2020"]  = "r"
releves_flore$abondance[releves_flore$station == "DY73_F" & releves_flore$cd_nom_terrain == 92302 &
                          releves_flore$date == "4/28/2020"] = "+"

## Add year
releves_flore$annee = as.numeric(gsub("^.*/(\\d{4})$", "\\1", perl = T, releves_flore$date))

## Add presence data
src1 = list(c(""), c("i", "r", "+", "1", "2", "3", "4", "5"))
tgt1 = list(0, 1)
releves_flore$presence = recodeVar(releves_flore$abondance, src=src1, tgt=tgt1)
releves_flore$presence = as.numeric(releves_flore$presence)

## Formatting
variables = c("saison_rel", "campagne", "station","especes_non_ok")
releves_flore[variables] = lapply(releves_flore[variables], factor)

#### Adding permanent nomenclature  
## Erase duplicate data
nomenclature = unique(nomenclature)
releves_flore$cd_nom_terrain = as.character(releves_flore$cd_nom_terrain)


releves_flore = merge(releves_flore, nomenclature, by.x = "cd_nom_terrain", by.y = "CD_OFB_BH", all.x = T)
releves_flore = renameCol(releves_flore, c("Nom_OFB_OK"),
                          c("especes_ok"))

#### Deleting unprospected stations
stations_suppr = stations_suppr[, c("Cd_CSE", "Comparaison.des.P12")]
liste_suppr = stations_suppr$Cd_CSE[stations_suppr$Comparaison.des.P12 == 0]

## Capital letters formatting
liste_suppr = toupper(liste_suppr)
releves_flore$station = toupper(releves_flore$station)
variables_09_simp$Cd_CSE = toupper(variables_09_simp$Cd_CSE)
variables_16_simp$Cd_CSE = toupper(variables_16_simp$Cd_CSE)

releves_flore = releves_flore[!releves_flore$station %in% liste_suppr, ]
variables_09_simp = variables_09_simp[!variables_09_simp$Cd_CSE %in% liste_suppr, ]
variables_16_simp = variables_09_simp[!variables_09_simp$Cd_CSE %in% liste_suppr, ]

#### Stations environmental data
variables_09 = variables_09_simp[, c("Cd_CSE", "Topo", "Indice_topo","Pente_categories","Pente_.valeur..",
                                     "Exposition", "Altitude", "Humus.contexte", "Nom_ensnat", "Cd_maille2_.EnsNat",
                                     "niveau.1...contextes.d.usages.des.sols", "niveau.2...principaux.usages.du.sol","niveau.3...usages.du.sol",
                                     "Nom_regnat", "Cd_regnat", "Usage.a.la.date.d.observation.p1.et.commentaires","Usage.a.la.date.d.observation.p2.et.commentaires")]

variables_16 = variables_16_simp[, c("Cd_CSE", "Topo", "Indice_topo", "Pente_categories", "Pente_.valeur..",
                                     "Exposition", "Altitude", "Humus.contexte", "Nom_ensnat", "Cd_maille2_.EnsNat",
                                     "niveau.1...contextes.d.usages.des.sols", "niveau.2...principaux.usages.du.sol","niveau.3...usages.du.sol",
                                     "Nom_regnat", "Cd_regnat", "Usage.a.la.date.d.observation.p1.et.commentaires","Usage.a.la.date.d.observation.p2.et.commentaires")]

variables_09 = renameCol(variables_09,
                         c("Cd_CSE", "Topo", "Indice_topo", "Pente_categories", "Pente_.valeur..",
                           "Exposition", "Altitude", "Humus.contexte", "Nom_ensnat", "Cd_maille2_.EnsNat",
                           "niveau.1...contextes.d.usages.des.sols", "niveau.2...principaux.usages.du.sol","niveau.3...usages.du.sol",
                           "Nom_regnat", "Cd_regnat", "Usage.a.la.date.d.observation.p1.et.commentaires","Usage.a.la.date.d.observation.p2.et.commentaires"),
                         c("station", "topo", "indice_topo", "pente_categorie", "pente_valeur", "exposition",
                           "altitude", "humus_contexte", "nom_ensnat", "maille", "occ_sol_n1", "occ_sol_n2","occ_sol_n3",
                           "nom_regnat", "cd_regnat", "usage_p1", "usage_p2"))
variables_16 = renameCol(variables_16,
                         c("Cd_CSE", "Topo", "Indice_topo", "Pente_categories", "Pente_.valeur..",
                           "Exposition", "Altitude", "Humus.contexte", "Nom_ensnat", "Cd_maille2_.EnsNat",
                           "niveau.1...contextes.d.usages.des.sols", "niveau.2...principaux.usages.du.sol","niveau.3...usages.du.sol",
                           "Nom_regnat", "Cd_regnat", "Usage.a.la.date.d.observation.p1.et.commentaires","Usage.a.la.date.d.observation.p2.et.commentaires"),
                         c("station", "topo", "indice_topo", "pente_categorie", "pente_valeur", "exposition",
                           "altitude", "humus_contexte", "nom_ensnat", "maille", "occ_sol_n1", "occ_sol_n2","occ_sol_n3",
                           "nom_regnat", "cd_regnat", "usage_p1", "usage_p2"))

## Create a column only for station name
variables_09$maille = str_split_fixed(variables_09$station,"_",2)[,1] 
variables_16$maille = str_split_fixed(variables_16$station,"_",2)[,1] 

## Merging informations into survey file
releves_flore_09 = merge(releves_flore[releves_flore$annee <= 2013,], variables_09[,c("station","occ_sol_n1","occ_sol_n2","occ_sol_n3")], by = "station", all = FALSE)
releves_flore_16 = merge(releves_flore[releves_flore$annee >= 2016,], variables_16[,c("station","occ_sol_n1","occ_sol_n2","occ_sol_n3")], by = "station", all = FALSE)
releves_flore = rbind(releves_flore_09, releves_flore_16)
releves_flore[c("occ_sol_n1","occ_sol_n2","occ_sol_n3","especes_ok")] = lapply(releves_flore[c("occ_sol_n1","occ_sol_n2","occ_sol_n3","especes_ok")], factor)
releves_flore$maille = str_split_fixed(releves_flore$station,"_",2)[,1]

## Add new soil nomenclature
src2 = list(c("ACA"), c("ACF","ACJ","ACP","NTN","NTT"),c("ADB","HBA","HBB","HBI","HLA","HLB","HLF","NAS"),
            c("APA","APF","APL","APP"),c("FAA","FAR"),c("FNA","FNR"),c("AMC","AMP","FMF","FMI"))
tgt2 = list(1,2,3,4,5,6,7)
releves_flore$occ_sol_mix = recodeVar(releves_flore$occ_sol_n3, src=src2, tgt=tgt2)
releves_flore$occ_sol_mix = as.factor(releves_flore$occ_sol_mix)

##############################################################################################################
## Formatting traits data --------------------------------------------------

## Selecting useful variables
variables_int = c("CD_nom.Taxref.12","Nom.valide.Taxref.12","Nom.simple","Nom.commun","FAMILLE",
                  "Hyb","Envahissante","LRR.Bourgogne","LR.France..2019.","LR.Europe..2012.","Messicoles.PNA",
                  "Messicoles.Bourgogne","Stat..1","Stat..2", colnames(baseflor[,58:74]))
baseflor = baseflor[,variables_int]
baseflor = within(baseflor, rm("CATMINAT"))

## Traits file correction
baseflor_correction = baseflor_correction[,1:9]
baseflor_correction[c("Humid_Ell_tbphyto", "pH_Ell_tbphyto", "Trophie_Ell_tbphyto", "Lum_Ell_tbphyto", 
                      "Temp_Ell_tbphyto", "Kont_Ell_tbphyto")] = lapply(baseflor_correction[c("Humid_Ell_tbphyto", "pH_Ell_tbphyto", "Trophie_Ell_tbphyto", "Lum_Ell_tbphyto", 
                                                                                              "Temp_Ell_tbphyto", "Kont_Ell_tbphyto")], as.numeric)

## Formatting
variables2 = c("Hyb","Envahissante","LRR.Bourgogne","LR.France..2019.","LR.Europe..2012.", "Messicoles.PNA",
               "Messicoles.Bourgogne", "pollinisation","Stat..1","Stat..2")
baseflor[variables2] = lapply(baseflor[variables2], factor)

variables3 = c("Humid_Ell_tbphyto", "pH_Ell_tbphyto", "Trophie_Ell_tbphyto", "Lum_Ell_tbphyto", 
               "Temp_Ell_tbphyto", "Kont_Ell_tbphyto", "Salin_Ell_tbphyto")
baseflor[variables3] = lapply(baseflor[variables3], as.numeric)

infos_stations$mois = strftime(infos_stations$newdate, "%m")

## Protection status
nom_prot = colnames(listsp_prot[, 2:12])
listsp_prot[nom_prot] = lapply(listsp_prot[nom_prot], factor)

##Taxonomic groups gathering list
list_sp_tot = releves_flore %>% 
  select(especes_ok, especes_non_ok) %>% 
  distinct(especes_non_ok, .keep_all = T) 

list_sp_tot$especes_ok = as.character(list_sp_tot$especes_ok)
list_sp_tot$especes_non_ok = as.character(list_sp_tot$especes_non_ok)

list_sp_tax = list_sp_tot %>% 
  mutate(keep = case_when(especes_non_ok != especes_ok ~ 1, TRUE ~ 0)) %>% 
  filter(keep == 1) %>% 
  mutate(keep_gr =  stringr::str_split(especes_ok, " ", n= 3)) %>% 
  filter(especes_ok != "Taraxacum F.H.Wigg.") %>% 
  mutate(keep_gr = map_chr(keep_gr, 3)) %>% 
  filter(keep_gr == "Gr.") %>% 
  select(especes_non_ok, especes_ok)

##############################################################################################################
## Species presence probability over time trends --------------------------------------------------

## Grouping occurrences with the same naturality / stratum / name without correction
releves_flore_sp = releves_flore %>% 
  filter(!(strate == "A")) %>% 
  arrange(id_releve, especes_ok, desc(strate)) %>% 
  distinct(especes_ok, id_releve, annee, .keep_all = T) %>% 
  group_by(especes_ok, annee, id_releve) %>% 
  distinct(.keep_all = T)

## Model function
pres_model = function(df) {
  glmer(presence ~ saison_rel + scale(annee) + (1|maille),
        data = df, family = "binomial"(link = "logit"), control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))
}

## Species minimal occurrence limit
listsp_occ = releves_flore_sp %>% 
  group_by(especes_ok) %>% 
  filter(presence == 1) %>% 
  count() %>% 
  filter(n >= 50)

listsp_occ = listsp_occ$especes_ok

## Nesting format
releves_nested = releves_flore_sp %>% 
  group_by(especes_ok, CD_OFB_OK) %>%
  filter(especes_ok %in% listsp_occ) %>% 
  nest()

## Model data extraction (model's and anova's parameters)
releves_nested = releves_nested %>% 
  mutate(model. = map(data, pres_model), tidy_out = map(model., broom.mixed::tidy), 
         anov = map(model., Anova), anovresu = map(anov, broom::tidy)) %>% 
  unnest(tidy_out, .drop = TRUE) %>% 
  unnest(anovresu, .drop = TRUE)

## Over-dispersion test and extraction of these data
releves_nested = releves_nested %>% 
  mutate(testdisper = map(model., testDispersion, plot = FALSE, alternative = "greater", type = "DHARMa"), output_dispersion = map(testdisper, broom::tidy)) %>% 
  unnest(output_dispersion, .drop = TRUE)

## Check for any model errors
releves_nested = releves_nested %>% 
  mutate(singularity = map(model., isSingular, tol = 1e-4)) %>% 
  filter(p.value != 0)

## P-value correction
releves_pvalue.correct = releves_nested %>% 
  group_by(especes_ok) %>% 
  summarise(sum1 = sum(p.value1), sum2 = sum(p.value2)) %>% 
  arrange(sum1)

k = length(unique(releves_nested$especes_ok))
releves_nested$p.value1.correct = releves_nested$p.value1
releves_nested$p.value2.correct = releves_nested$p.value2
for (j in releves_pvalue.correct$especes_ok) {
  releves_nested$p.value1.correct[releves_nested$especes_ok == j] = (releves_nested$p.value1[releves_nested$especes_ok == j])*k
  releves_nested$p.value2.correct[releves_nested$especes_ok == j] = (releves_nested$p.value2[releves_nested$especes_ok == j])*k
  k = k-1
}

## Student test for average trend difference 
t.test(estim_density$estimate, mu = 0)

## Native / introduced species 
estim_density_stat = merge(estim_density, baseflor[,c("CD_nom.Taxref.12", "pollinisation", "floraison","Stat..1", "Stat..2" , variables3)],
                           by.x = "CD_OFB_OK", by.y = "CD_nom.Taxref.12", all.x = T)
estim_density_stat = merge(estim_density_stat, unique(releves_flore[,c("CD_OFB_OK","indigenat")]),
                           by = "CD_OFB_OK", all.x = T)
## Student test for average trend difference between native and introduced species
t.test(estim_density_stat$estimate[estim_density_stat$Stat..1 == "Ind."], estim_density_stat$estimate[estim_density_stat$Stat..1 == "Nat. (E.)"])

## Protected species
estim_density_prot = merge(estim_density, listsp_prot,
                           by.x = "CD_OFB_OK", by.y = "CD_REF", all.x = T)
estim_density_prot %>% filter_at(nom_prot, any_vars(!is.na(.)))

##############################################################################################################
## Species diversity --------------------------------------------------

############ Species richness
releves_rich_sais = releves_flore[releves_flore$presence == 1,]
releves_rich_sais$annee = as.numeric(as.character(releves_rich_sais$annee))

## Count of richness
releves_rich_sais_maille = releves_rich_sais %>% 
  group_by(annee,saison_rel, maille, station, id_releve, especes_ok, occ_sol_n1, occ_sol_n2, occ_sol_mix) %>% 
  count() %>% 
  group_by(annee,saison_rel, maille, station, id_releve, occ_sol_n1, occ_sol_n2, occ_sol_mix) %>% 
  count() 

## Formating
releves_rich_sais_maille = renameCol(releves_rich_sais_maille, c("n"),
                                     c("richesse"))
releves_rich_sais_maille$occ_sol_n1 = as.factor(releves_rich_sais_maille$occ_sol_n1)
releves_rich_sais_maille$occ_sol_n2 = as.factor(releves_rich_sais_maille$occ_sol_n2)
releves_rich_sais_maille$occ_sol_mix = as.factor(releves_rich_sais_maille$occ_sol_mix)

## Model
mx = glmmTMB(richesse ~ saison_rel + scale(annee)*occ_sol_mix + (1|maille/station), 
             family = "poisson", data = releves_rich_sais_maille, contrasts = list(occ_sol_mix = "contr.sum"))
summary(mx)
Anova(mx)

## Over-dispersion test
testDispersion(mx, plot = F, type = "DHARMa", alternative = "greater")

## Post-hoc test for the interaction parameter
leastsquare = lsmeans(mx, pairwise ~ occ_sol_mix:scale(annee)) 
leastsquare$lsmeans
leastsquare$contrasts

############ SHannon and Pielou indexes

## Calculation of indexes
releves_shannon = releves_ellen %>% 
  mutate(shannon_sp = (abondance_num.norm)*log(abondance_num.norm)) %>% 
  group_by(id_releve, annee, saison_rel, occ_sol_mix, station, maille) %>% 
  summarise(richesse_brute = n(), shannon = -sum(shannon_sp), maxshannon = log(richesse_brute), eveness = shannon/maxshannon)

## Shannon index model
mshannon = lmer(shannon ~ saison_rel + scale(annee)*occ_sol_mix + (1|maille/station),
                data = releves_shannon, contrasts = list(occ_sol_mix = "contr.sum"))

testDispersion(mshannon, alternative = "greater")
Anova(mshannon)
summary(mshannon)

ggemmeans(mshannon, terms = c("annee","occ_sol_mix"), type = "fe") # Model predictions

## Pielou index model
mpielou = lmer(eveness ~ saison_rel + scale(annee)*occ_sol_mix + (1|maille/station),
               data = releves_shannon, contrasts = list(occ_sol_mix = "contr.sum"))

testDispersion(mpielou, alternative = "greater")
Anova(mpielou)
summary(mpielou)

ggemmeans(mpielou, terms = c("annee","occ_sol_mix"), type = "fe")

## Post-hoc test for the interaction parameter
leastsquare = lsmeans(mshannon, pairwise ~ occ_sol_mix:scale(annee), pbkrtest.limit = 5185)
leastsquare$lsmeans
leastsquare$contrasts

##############################################################################################################
## Traits Ellenberg analyses --------------------------------------------------

## Selection of present species
releves_ellen = releves_flore[releves_flore$presence == 1,]

## Adding traits variables
releves_ellen = merge(releves_ellen, baseflor[,c("CD_nom.Taxref.12", "pollinisation", "floraison","Stat..1", "Stat..2" , variables3)], by.x = "CD_OFB_OK", by.y = "CD_nom.Taxref.12", all.x = T)

## Adding blank information
for (i in baseflor_correction$CD_OFB_OK) {
  releves_ellen[releves_ellen$CD_OFB_OK == i, "Humid_Ell_tbphyto"] = baseflor_correction[baseflor_correction$CD_OFB_OK == i, "Humid_Ell_tbphyto"]
  releves_ellen[releves_ellen$CD_OFB_OK == i, "pH_Ell_tbphyto"] = baseflor_correction[baseflor_correction$CD_OFB_OK == i, "pH_Ell_tbphyto"]
  releves_ellen[releves_ellen$CD_OFB_OK == i, "Trophie_Ell_tbphyto"] = baseflor_correction[baseflor_correction$CD_OFB_OK == i, "Trophie_Ell_tbphyto"]
  releves_ellen[releves_ellen$CD_OFB_OK == i, "Lum_Ell_tbphyto"] = baseflor_correction[baseflor_correction$CD_OFB_OK == i, "Lum_Ell_tbphyto"]
  releves_ellen[releves_ellen$CD_OFB_OK == i, "Temp_Ell_tbphyto"] = baseflor_correction[baseflor_correction$CD_OFB_OK == i, "Temp_Ell_tbphyto"]
  releves_ellen[releves_ellen$CD_OFB_OK == i, "Kont_Ell_tbphyto"] = baseflor_correction[baseflor_correction$CD_OFB_OK == i, "Kont_Ell_tbphyto"]
}

releves_ellen$id_releve = as.factor(releves_ellen$id_releve)

## Numerical transformation of abondance
releves_ellen$abondance_num = as.numeric(recodeVar(releves_ellen$abondance, c("i", "r", "+", "1", "2", "3", "4", "5"),
                                                   c(0.00025, 0.00045, 0.0045, 0.03, 0.15, 0.375, 0.625, 0.875)))

variables4 = paste0(variables3,"_abond_weighted")

## Sum of recovery rates for occurrences
releves_ellen = releves_ellen %>% 
  group_by(especes_ok, id_releve, indigenat, strate) %>% 
  mutate(abondance_num = sum(abondance_num)) %>% 
  distinct(especes_ok, id_releve, indigenat, strate, .keep_all = T)

releves_ellen = releves_ellen %>% 
  group_by(especes_ok, id_releve, indigenat) %>% 
  mutate(abondance_num = sum(abondance_num)) %>% 
  distinct(especes_ok, id_releve, indigenat, .keep_all = T)

releves_ellen = releves_ellen %>% 
  group_by(especes_ok, id_releve) %>% 
  mutate(abondance_num = sum(abondance_num)) %>% 
  distinct(especes_ok, id_releve, .keep_all = T)

## Recovery rates normalisation
normalize = function(x, na.rm = TRUE) {
  return(x/sum(x))
}

releves_ellen = releves_ellen %>% 
  group_by(id_releve) %>% 
  mutate(abondance_num.norm = normalize(abondance_num))

releves_ellen = releves_ellen %>% 
  filter(abondance_num.norm != 0)

## Deleting species without any trait information
pasindice = releves_ellen %>%
  group_by(especes_ok) %>% 
  distinct(especes_ok, .keep_all = T) %>%
  select(especes_ok, all_of(variables3)) %>%
  select(-Salin_Ell_tbphyto) %>%
  filter_all(any_vars(is.na(.)))

## Calculation of average EIV mean for each survey
mean_ellen = releves_ellen %>% 
  mutate_at(variables3, list(abond_weighted = ~(.*abondance_num.norm))) %>%   
  group_by(id_releve, maille, station, annee, saison_rel, occ_sol_mix) %>% 
  summarise_at(variables4, list(mean_EIV = ~(sum(., na.rm = TRUE)/sum(abondance_num.norm[!is.na(.)])))) 

## PCA combined with diversity results
mean_ellen2 = merge(mean_ellen, mean_polli[,c("id_releve","indice")], by = "id_releve")
ellen_rich = merge(mean_ellen2, releves_rich_sais_maille, by = "id_releve") %>% 
  drop_na() %>% 
  select(occ_sol_mix.x, annee.x, richesse, paste0(variables4,"_mean_EIV"), indice)
ellen_rich = ellen_rich[, -10]

pca_elrich = PCA(ellen_rich[,-1], graph = F)

ACP1 = fviz_pca_biplot(pca_elrich, axes = c(1,2), label ="var", habillage = ellen_rich$occ_sol_mix.x, 
                       addEllipses = T, col.var = "black", labelsize = 5, ellipse.level=0.95, invisible = NA)

## Long format
mean_ellen = melt(mean_ellen, measure.vars = paste0(variables4,"_mean_EIV"), value.name = "EIV_value", variable.name = "EIV_type")

## Deleting surveys without any information
mean_ellen = mean_ellen %>% 
  filter(!is.na(EIV_value)) 

## Nesting format
ellen_nested_sum = mean_ellen %>% 
  group_by(EIV_type) %>%
  nest()

## Model
ellen_mod_sum = function(df) {
  lmer(EIV_value ~ saison_rel + scale(annee)*occ_sol_mix + (1|maille/station),
       data = df, contrasts = list(occ_sol_mix = "contr.sum"))
}

ellen_mod_sum_glmm = function(df) {
  glmmTMB(EIV_value ~ saison_rel + scale(annee)*occ_sol_mix + (1|maille/station),
          data = df, contrasts = list(occ_sol_mix = "contr.sum"))
} ## Convergence anomaly

ellen_nested_sum = ellen_nested_sum %>% 
  mutate(model. = case_when(EIV_type != "Humid_Ell_tbphyto_abond_weighted_mean_EIV" ~ map(data, ellen_mod_sum), 
                            EIV_type == "Humid_Ell_tbphyto_abond_weighted_mean_EIV"~ map(data, ellen_mod_sum_glmm)))

## Anova / Over-dispersion tests / P-value correction
ellen_nested_sum = ellen_nested_sum %>% 
  mutate(tidy_out = map(model., broom.mixed::tidy), resid = map(model., broom.mixed::augment), 
         anov = map(model., Anova), anovresu = map(anov, broom::tidy), 
         testdisper = map(model., testDispersion, plot = F, type = "DHARMa", alternative = "greater"), 
         disperout = map(testdisper, broom.mixed::tidy)) %>% 
  unnest(disperout, .drop = T) %>% 
  mutate(p.value.corr.testdis = (p.value*6)) %>% 
  rename(p.value.testdis = p.value)

ellen_nested_sum = ellen_nested_sum %>% 
  unnest(tidy_out, .drop = T) %>%
  rename(p.value.param = p.value) %>% 
  mutate(p.value.param.corr = p.value.param*6) %>% 
  unnest(anovresu, .drop = T) %>%  
  mutate(p.value.correct.aov = (p.value*6)) %>% 
  rename(p.value.aov = p.value)