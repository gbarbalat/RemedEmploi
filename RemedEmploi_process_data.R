rm(list=ls())

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## TITRE DU PROJET : RemedEmploi 
##
## Investigateur coordonnateur : Nicolas Franck + Sophie CERVELLO
## N°ANSM : 2019-A00124-53
## Acquisition des data : juillet 2025
##                                                                  
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Rapport monito de visite de Monitoring de clôture 
# C : 86 + 24 + 30 = 140
# E : 49 + 9 + 5 = 63
# F : 37 + 11 + 21 = 69
# 
# ESAT Lyon - 01
# C Nombre de patients inclus 86
# D Nombre de patients en cours d’étude (au moment de la visite) 0
# E Nombre de patients ayant terminé l’étude conformément au protocole 49
# F Nombre de patients sortis prématurément de l’étude (décès, perdu de vue, EIG…) 37
# 
# PASSAGE PRO ALLONNE - 02
# C Nombre de patients inclus 24
# D Nombre de patients en cours 4
# E Nombre de patients ayant terminé l’étude 9
# F Nombre de patients sortis prématurément d'essai 11: c("001-D-J", "003-N-S", "009-G-A", "010-D-M", "011-M-N", "014-W-J", "015-M-C", "016-N-J", "018-P-C", "019-C-V", "024-B-E")

# MESSIDOR Bobigny - 03
# C Nombre de patients inclus 30
# D Nombre de patients en cours 4
# E Nombre de patients ayant terminé l’étude 5
# F Nombre de patients sortis prématurément d'essai 21

#setwd("X:\\RemedEmploi\\Work")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                        Import packages and libraries                      ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require('tidyverse')
require('gtsummary')
require('sjmisc')
library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Import and preprocessing data                     ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# IMPORT DES DATA (fusion des .sav fait par Julien)
df.full <- read.csv2("data_fusion.csv")

# cleaning
df.full %>%
  rename_with(~ stringr::str_remove(.x, "(_E1_C1)|(_E2_C2)|(_E7_C7)|(_E5_C5)")) %>%
  rename_with(~ stringr::str_replace_all(string = .x, pattern = "_E3_C3",replacement = "_t0" )) %>%
  rename_with(~ stringr::str_replace_all(string = .x, pattern = "_E4_C4",replacement = "_t1" )) %>%
  rename_with(~ stringr::str_replace_all(string = .x, pattern = "_E6_C6",replacement = "_t2" )) %>%
    mutate(ProtocolID = str_remove(ProtocolID, "PHRCN_RemedEmploi - "),
         across(where(is.numeric), ~ na_if(., 8888)),
         across(where(is.character), ~ na_if(., "8888"))) %>% # !! 8888 en NA
  rename_with(~ stringr::str_remove(.x, "^(RE_)|(IR_)")) %>%
  rename(CENTRE = ProtocolID,
         AGE = CALC_AGE) %>% 
  mutate(ASSIDUITE_SOIN = as.numeric(ASSIDUITE_SOIN)) %>% 
  # set_variable_labels(
  #   CENTRE = "Centre", 
  #   MINI_DIAG1 = "Diagnostic PRINCIPAL de trouble mental sévère",
  #   MINI_DIAG2 = "Diagnostics secondaires de TMS") %>% 
  #drop_unused_value_labels() %>% 
  to_label() -> df.full

# conversion des dates
df.full %>% 
  mutate(MOIS_M1=as.character(MOIS_M1)) %>%
  mutate(across(c(contains("DATE"),"EMBAUCHE","MOIS_M1"),as.Date, format="%d/%m/%Y")) -> df.full
str(df.full$DATE_CONSENTEMENT)
str(df.full$MOIS_M1)

df.full %>% 
  select(all_of(contains(c("HEURES_PREVUES_","HEURES_EFFECTUEES_","DUREE_ABSENCE_")))) -> tmp


# Total Primary Outcome
# faut recalculer (car sinon les 8888 compte comme 8888h..)
df.full %>%
  select(-PREVUES_TOTAL,-EFFECTUEES_TOTAL) -> df.full

# 10 patients  ont stoppés l'étude puis retour (décision avec EG : on ne garde que les _BIS)
idpatbis <- grep(pattern = "_BIS$", x = df.full$StudySubjectID, value = T)
idpatwithoutbis <- str_remove(idpatbis,"_BIS")
BDD_BIS <- df.full %>%
  filter(StudySubjectID %in% idpatbis | StudySubjectID %in% idpatwithoutbis)

#WHAT TO DO with 10? ----
df.full %>%
  #filter(StudySubjectID %notin% idpatwithoutbis) %>%
  filter(!StudySubjectID %in% idpatwithoutbis) %>%
  #filter(!StudySubjectID %in% idpatbis) %>%
  mutate(StudySubjectID = str_remove(StudySubjectID, "_BIS$"))  -> df.full

# Liste rando
df.rando <- read.csv2("listeRando_RemedEmploi_fusion.csv") %>%
  mutate(across(c("Inclusion", "Naissance", "Entree_ESAT", "Rando"), ~ as.Date(.x,format="%d/%m/%Y")),
         "AGE_INCLUS" = as.numeric(Inclusion - Naissance) / 365.25
  ) -> df.rando

# pb des caractères invisibles
clean_invisible <- function(df) {
  df |> mutate(across(
    where(is.character),
    ~ str_replace_all(.x, "[[:cntrl:]]|\\u00A0|\\u200B|\\u200C|\\u200D|\\uFEFF", "")
  ))
}

df.full  <- clean_invisible(df.full)
df.rando <- clean_invisible(df.rando)

setdiff(df.full$StudySubjectID,df.rando$StudySubjectID) # aucun id dans rehabase qui n'est pas dans la liste inclusion -> OK !!
setdiff(df.rando$StudySubjectID,df.full$StudySubjectID) # inverse -> OK !!

  

# FIN : fusion avec la liste rando
df.full %>% 
  left_join(df.rando, by="StudySubjectID") %>% 
  relocate(names(df.rando), .after = StudySubjectID) -> df.full

# traitement des dates d'inclusion pour les 10 qui sont revenus dans l'étude 
#(-> à la date de la rando et pas de sa première inclusion)
df.full %>% 
  filter(StudySubjectID %in% idpatwithoutbis) %>% 
  select(StudySubjectID,Inclusion,Rando,EMBAUCHE)

df.full$Inclusion[df.full$StudySubjectID %in% idpatwithoutbis] <- df.full$Rando[df.full$StudySubjectID %in% idpatwithoutbis] 
df.full$Inclusion[df.full$StudySubjectID == "03016MA"] <- "2021-04-30" # pas de date de rando pour lui, on lui met la date de sa première venue

# calcul du délai de l'embauche à l'ESAT
# A NOTER : date EMBAUCHE 1970-12-05... erreur,
#j'ai pris la date Entree_ESAT du fichier liste rando et supprimé la date du 1970-12-05 dans EMBAUCHE 
#(03003AM, 03002AS, 03005SR,03007EW,03008JS )
df.full$EMBAUCHE[which(df.full$EMBAUCHE == "1970-12-05")] <- NA

df.full %>% 
  select(StudySubjectID,EMBAUCHE,Entree_ESAT) %>% 
  mutate("diff_ESAT"=abs(EMBAUCHE-Entree_ESAT)) %>% 
  filter(diff_ESAT>=90) # diff sup à 3 mois, seul 2 patients (03019DA,03019DA), ils ne sont pas randomisés

df.full <- df.full %>% 
  mutate("DIFF_ESAT_INCLUS" = ifelse(!is.na(EMBAUCHE), EMBAUCHE - Inclusion, Entree_ESAT - Inclusion), .after = AGE_INCLUS,
         "DIFF_ESAT_INCLUS_crit" = ifelse(DIFF_ESAT_INCLUS < (-30.42*18),"non","oui"),
         "IAG_crit" = ifelse(APT_GEN < 70,"non","oui"),
         #"ASSIDUITE_SOIN_crit" = ifelse(ASSIDUITE_SOIN < 7,"oui","non"),
         "ASSIDUITE_SOIN_crit" = ifelse(ASSIDUITE_SOIN < 7,"non","oui")
         ) 

# vérif DIFF_ESAT_INCLUS
df.full %>%
  filter(StudySubjectID == "03014NM") %>%
  select(StudySubjectID ,Inclusion,EMBAUCHE,Entree_ESAT)

# sans le traitement des _BIS pour vérifier les dates des 10 réinclus
df.full %>%
  filter(StudySubjectID %in% c(idpatbis,idpatwithoutbis)) %>%
  arrange(StudySubjectID) %>%
  select(StudySubjectID, DATE_CONSENTEMENT,Inclusion,Rando, Entree_ESAT, EMBAUCHE,DATE_SORTIE, HEURES_PREVUES_1,HEURES_EFFECTUEES_1)

# ESAT_Inclusion > 18 months
df.full %>%
  filter(DIFF_ESAT_INCLUS_crit == "non") %>%
  select(StudySubjectID, DIFF_ESAT_INCLUS, RATIO_1, SERS_SC_TOT_t1, SERS_SC_TOT_t2)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Non rando et équilibre ----
# A noter : 03023JA sortie avant rando et avec bras (ARRET maladie de LONGUE durée (Sep23); 
# A sortir de la randomisation; remplacé par 03028SJ dans le Gp1 RCS), j'ai supprimé le bras dans rando_fusion
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df.full %>% 
  filter(!is.na(BRAS)) %>% 
  tbl_summary(include = BRAS)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Vérif des critère d'inclusion ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df.full <- df.full %>% 
  filter(!is.na(BRAS)) %>%
  filter(BRAS!="")

df.full %>% 
  filter(DIFF_ESAT_INCLUS>0)

df.full %>%
  tbl_summary(include = c(CENTRE,AGE_INCLUS, IAG_crit, MINI_DIAG1, DIFF_ESAT_INCLUS, DIFF_ESAT_INCLUS_crit), 
              type = all_continuous() ~ "continuous2", statistic = list(all_continuous() ~ c(
    "{min} - {max}",
    "{mean} ({sd})"
  )))


# differentiels de date
important_dates <- c("StudySubjectID", "Rando", "CONSENTEMENT_MEDIC",  "MEDIC_t0", "ENPSY_MEDIC",  
                     "M1_Rando" ,"M1_MEDIC" , "M1", "M1_t0" , "t1_t0","t1_Rando","t2_t0","t2_Rando",
                     "M1_t1" , "t2_t1")
library(lubridate)

df.full <- df.full %>%
  mutate(
    # interval() creates the span, time_length(unit = "month") converts it accurately
    CONSENTEMENT_MEDIC = time_length(interval(MEDIC_DATE,DATE_CONSENTEMENT), "month")%>% round(),
    ENPSY_MEDIC        = time_length(interval(MEDIC_DATE, ENPSY_DATE), "month") %>% round(),
    MEDIC_t0           = time_length(interval(ISMI_DATE_t0,MEDIC_DATE), "month") %>% round(), 
    M1_Rando           = time_length(interval(Rando,MOIS_M1), "month") %>% round(),
    M1_MEDIC           = time_length(interval(MEDIC_DATE, MOIS_M1), "month") %>% round(),
    M1                 = HEURES_EFFECTUEES_1,
    
    M1_t0              = time_length(interval(ISMI_DATE_t0, MOIS_M1), "month") %>% round(),
    
    
    t1_t0              = time_length(interval(ISMI_DATE_t0, ISMI_DATE_t1), "month") %>% round(),
    t1_Rando              = time_length(interval(Rando, ISMI_DATE_t1), "month") %>% round(),
    
    t2_t0              = time_length(interval(ISMI_DATE_t0, ISMI_DATE_t2), "month") %>% round(),
    t2_Rando              = time_length(interval(Rando, ISMI_DATE_t2), "month") %>% round(),
    
    
    M1_t1              = time_length(interval(ISMI_DATE_t1, MOIS_M1), "month") %>% round(),
    t2_t1              = time_length(interval(ISMI_DATE_t1, ISMI_DATE_t2), "month") %>% round()
  ) 

tmp <- df.full %>%
  select(all_of(important_dates)) -> tmp
head(tmp)
summary(tmp)
purrr::map(tmp, table, useNA="always")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Vérif de l'ordonnancement ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Eval ESAT post ttt
df.full %>% 
  #mutate(MOIS_M1_rec = paste0("01/", MOIS_M1) |> as.Date(format = "%d/%m/%Y")  ) %>% 
  mutate(MOIS_M1_rec = MOIS_M1) %>% 
  
  filter(MOIS_M1_rec<SERS_DATE_t0) %>% 
  select(StudySubjectID,SERS_DATE_t0, SERS_DATE_t1, SERS_DATE_t2, MOIS_M1, SORTIE) # %>%   pull(StudySubjectID)




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Descriptif de la population ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# on ne garde que les scores totaux ou sous-dimensions
df.full %>% 
  select(!starts_with(c("SERS_DATE","SERS_Q",
                        "ISMI_DATE","ISMI_Q",
                        "STORI_DATE","STORI_Q",
                        "STAIA_DATE","STAIA_Q",
                        "STAIB_DATE","STAIB_Q",
                        "ENPSY_DATE",
                        "CALC_AGE_","SEXE_","NIVETUD_"))) -> df.full

# Pré inclusion ----

varnames.preinclusion <- c("MINI_DIAG1", "TSA", "MATR_NB", "MATR_STDR", 
                           "CUBES_NB", "CUBES_STDR", 
                           "PUZZ_NB", "PUZZ_STDR", "SIMIL_NB", "SIMIL_STDR", "INFO_NB", "INFO_STDR", "VOCAB_NB", "VOCAB_STDR", "IRP", "ICV", "APT_GEN")
table(df.full$MINI_DIAG1)
df.full$MINI_DIAG1 <- as.factor(df.full$MINI_DIAG1)
df.full %>% 
  select(BRAS,CENTRE,AGE,varnames.preinclusion) %>% 
  tbl_summary(by = BRAS) %>% 
  add_overall() 

# Évaluation 
varnames.eval <- c("SERS_POS", "SERS_NEG", "SERS_SC_TOT", "ISMI_ALIEN", "ISMI_APPROB_STEREO", "ISMI_EXP_DISCRIM", "ISMI_RET_SOC", "ISMI_INTERN_STIGMA", "ISMI_RESIST_STIGMA", "ISMI_SC_TOT", "STORI_MORATOIRE", "STORI_CONSCIENCE", "STORI_PREPARATION", "STORI_RECONSTRUCTION", "STORI_CROISSANCE", "STORI_ST_RETABLISSEMENT", "STAIA_TOT", "STAIB_TOT", "ACSO", "ACSO_H", "MASC_NB", "MASC_SIGMA", "AIHQ_HB_NB", "AIHQ_HB_SIGMA", "AIHQ_ATTRIB_RESP_NB", "AIHQ_ATTRIB_RESP_SIGMA", "AIHQ_AB_NB", "AIHQ_AB_SIGMA", "PERSO_FLU_PONC", "PERSO_INTER_LIB", "PERSO_INTER_INDIC", "PERSO_INTER_TOTAL", "PERSO_CONV_SOC", "TREF_BR_DEGOUT_NB", "TREF_BR_DEGOUT_SIGMA", "TREF_BR_MEPRIS_NB", "TREF_BR_MEPRIS_SIGMA", "TREF_BR_JOIE_NB", "TREF_BR_JOIE_SIGMA", "TREF_BR_PEUR_NB", "TREF_BR_PEUR_SIGMA", "TREF_BR_TRISTESSE_NB", "TREF_BR_TRISTESSE_SIGMA", "TREF_BR_COLERE_NB", "TREF_BR_COLERE_SIGMA", "TREF_BR_TOTAL_NB", "TREF_BR_TOTAL_SIGMA", "TREF_SD_DEGOUT_NB", "TREF_SD_DEGOUT_SIGMA", "TREF_SD_MEPRIS_NB", "TREF_SD_MEPRIS_SIGMA", "TREF_SD_JOIE_NB", "TREF_SD_JOIE_SIGMA", "TREF_SD_PEUR_NB", "TREF_SD_PEUR_SIGMA", "TREF_SD_TRISTESSE_NB", "TREF_SD_TRISTESSE_SIGMA", "TREF_SD_COLERE_NB", "TREF_SD_COLERE_SIGMA", "TREF_SD_TOTAL_NB", "TREF_SD_TOTAL_SIGMA")
varnames.eval.quanti <- c("PERSO_INTER_INDIC", "PERSO_CONV_SOC")

# Evaluation T0 ----
df.full %>% 
  select(paste0(varnames.eval,"_t0")) %>% 
  tbl_summary(type = paste0(varnames.eval.quanti,"_t0") ~  "continuous")

# Évaluation T1 ----
df.full %>% 
  select(BRAS,paste0(varnames.eval,"_t1")) %>% 
  tbl_summary(by = BRAS, type = c(paste0(varnames.eval.quanti,"_t1"),starts_with("TREF_")) ~  "continuous")

# Évaluation T2 ----
df.full %>% 
  select(paste0(varnames.eval,"_t2")) %>% 
  tbl_summary(type = c(paste0(varnames.eval.quanti,"_t2"),starts_with("TREF_")) ~  "continuous")


# Sortie d'étude ----
# Assiduité aux séances du groupe de soins : !! Attention!! La participation à moins de 60% du contenu du groupe (7/12) entraîne une Sortie d’Étude prématurée !
df.full %>% 
  select("SORTIE", "ETAPE", "MOTIF", "MOTIF_OTH","ASSIDUITE_SOIN_crit") %>% 
  tbl_summary()
df.full %>%
  select(SORTIE, SERS_SC_TOT_t0, SERS_SC_TOT_t1, MOIS_M1, M1_MEDIC, M1_t0, HEURES_PREVUES_1, SERS_SC_TOT_t2) -> tmp

# vérif équilibre des sorties /bras
df.full %>% 
  #filter(SORTIE == 0 | ASSIDUITE_SOIN_crit == "oui") %>% 
  #filter(SORTIE == 1) %>% 
  filter(ASSIDUITE_SOIN_crit == "oui") %>% 
  
  pull(StudySubjectID) -> idpat.sortie

df.full$SORTIE_TOT <- df.full$StudySubjectID %in% idpat.sortie
chisq.test(table(df.full$SORTIE_TOT,df.full$BRAS))


# Assiduite des soins/Heures travaillees ----
df.full %>% 
  select("ASSIDUITE_SOIN",
         all_of(contains(c("PREVUES_","EFFECTUEES_", "ABSENCE_")))
         )%>% 
  tbl_summary(type = ASSIDUITE_SOIN ~  "categorical") 


#Prepare anal ----
library(lmtp)
# lmtp avec analyse for each time point as in Hoffman tutorial and then same graph
MI <- FALSE
run_ITT_MI <- FALSE
ClaCos_anal <- FALSE
m <- 30 #for MI even if MI is FALSE just one MI anal (of the outcome) will be performed
maxit <- 15
outcome <- "RATIO_TOT"
all_outcomes_CJP <- c( "RATIO_TOT", "RATIO_H_T_1", "RATIO_H_T_2", "RATIO_H_T_3", "RATIO_H_T_4", "RATIO_H_T_5", "RATIO_H_T_6")
all_outcomes_CJS <- c("SERS_SC_TOT", "ISMI_SC_TOT", "STORI_CROISSANCE", 
                   "STAIA_TOT", "STAIB_TOT")
if (ClaCos_anal) {
all_outcomes_CJS <- c( "ACSO", "ACSO_H", 
                    "MASC_NB", 
                    "AIHQ_HB_NB","AIHQ_ATTRIB_RESP_NB", "AIHQ_AB_NB",  
                    "PERSO_FLU_PONC", "PERSO_INTER_LIB", "PERSO_INTER_INDIC", "PERSO_INTER_TOTAL", "PERSO_CONV_SOC", 
                    "TREF_BR_TOTAL_NB","TREF_SD_TOTAL_NB")
}
baseline_CJP <- c("Sexe", "Niv_IAG", "AGE_INCLUS", "CENTRE", "dx", "Study_Year", #"Rando","Study_Wave","Study_Year"
              #"ASSIDUITE_SOIN",
              "ISMI_SC_TOT_t0"#"SERS_SC_TOT_t0", "ISMI_SC_TOT_t0", "STORI_CROISSANCE_t0", "STAIA_TOT_t0", "STAIB_TOT_t0"
              ) 
baseline_CJS <- c("Sexe", "Niv_IAG", "AGE_INCLUS", "CENTRE", "dx", "Study_Year") #"Rando","Study_Wave","Study_Year"
              #"ASSIDUITE_SOIN",
              #"STORI_CROISSANCE_t0"#"SERS_SC_TOT_t0", "ISMI_SC_TOT_t0", "STORI_CROISSANCE_t0", "STAIA_TOT_t0", "STAIB_TOT_t0"
mediator <- "ASSIDUITE_SOIN"#for crumble use only

learners_trt <- c("SL.mean", "SL.glm", "SL.ranger", "SL.bayesglm", "SL.glmnet", "SL.earth", "SL.xgboost")#, "SL.xgboost", "SL.bayesglm", "SL.earth", "SL.glmnet"
learners_trt <- c("SL.mean", "SL.glm","SL.glmnet") #"SL.mean", "SL.earth"
#correcting for "Bad Luck" Randomization
#In small to medium RCTs, randomization doesn't always work perfectly.
#By chance, your treatment group might end up slightly older or sicker than your control group.
#Including those variables as covariates "levels the playing field" mathematically, 
#ensuring the treatment effect you see isn't just a result of these baseline imbalances.

learners_outcome <- c("SL.mean", "SL.glm", "SL.ranger", "SL.bayesglm", "SL.glmnet", "SL.earth", "SL.xgboost")#, "SL.xgboost", "SL.bayesglm", "SL.earth", "SL.glmnet"
learners_outcome <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.earth")#, "SL.xgboost", "SL.bayesglm", "SL.earth", "SL.glmnet"
nfolds <- 10

df.full <- df.full %>%
  mutate(dx=case_when(MINI_DIAG1 %in% c(1,2,3,5) ~ "Mood",
                      MINI_DIAG1 %in% c(6,7,8,9,10,16) ~ "Anx",
                      MINI_DIAG1 %in% c(13,0,17) ~ "SCZ"
                      #MINI_DIAG1==13 ~ "SCZ"
                      )
         ) %>%
    rowwise() %>%
    mutate("PREVUES_TOTAL" = ifelse(all(is.na(c_across(starts_with("HEURES_PREVUES_")))),NA,sum(c_across(starts_with("HEURES_PREVUES_")), na.rm = T)),
           "EFFECTUEES_TOTAL" =  ifelse(all(is.na(c_across(starts_with("HEURES_EFFECTUEES_")))),NA, sum(c_across(starts_with("HEURES_EFFECTUEES_")), na.rm = T)),
           "ABSENCE_TOTAL" =   ifelse(all(is.na(c_across(starts_with("DUREE_ABSENCE_")))),NA,sum(c_across(starts_with("DUREE_ABSENCE_")), na.rm = T)), .before = "MOIS_M1") %>%
  ungroup %>%
  
  #mutate(RATIO_TOT = ifelse(PREVUES_TOTAL==0,NA,min(1,round(EFFECTUEES_TOTAL/PREVUES_TOTAL,2)))) %>% # on tronque à 1 si > 1
  mutate(RATIO_TOT = round(EFFECTUEES_TOTAL/PREVUES_TOTAL, 2), .before = "MOIS_M1") %>% # on tronque à 1 si > 1
  
  mutate(
    RATIO_H_T_1 = ifelse(HEURES_PREVUES_1 == 0, NA, pmin(1, round(HEURES_EFFECTUEES_1 / HEURES_PREVUES_1, 2))),
    RATIO_H_T_2 = ifelse(HEURES_PREVUES_2 == 0, NA, pmin(1, round(HEURES_EFFECTUEES_2 / HEURES_PREVUES_2, 2))),
    RATIO_H_T_3 = ifelse(HEURES_PREVUES_3 == 0, NA, pmin(1, round(HEURES_EFFECTUEES_3 / HEURES_PREVUES_3, 2))),
    RATIO_H_T_4 = ifelse(HEURES_PREVUES_4 == 0, NA, pmin(1, round(HEURES_EFFECTUEES_4 / HEURES_PREVUES_4, 2))),
    RATIO_H_T_5 = ifelse(HEURES_PREVUES_5 == 0, NA, pmin(1, round(HEURES_EFFECTUEES_5 / HEURES_PREVUES_5, 2))),
    RATIO_H_T_6 = ifelse(HEURES_PREVUES_6 == 0, NA, pmin(1, round(HEURES_EFFECTUEES_6 / HEURES_PREVUES_6, 2))),

  ) %>%
  mutate(BRAS=case_when(BRAS=="I" ~ 0,
                        BRAS=="RCS" ~ 1)
         ) %>%

  mutate(Study_Wave = case_when(
    Rando < as.Date("2021-01-01") ~ "Late_2020",
    
    Rando >= as.Date("2021-01-01") & Rando < as.Date("2021-09-01") ~ "Early_2021",
    Rando >= as.Date("2021-09-01") & Rando < as.Date("2022-01-01") ~ "Late_2021",
    
    Rando >= as.Date("2022-01-01") & Rando < as.Date("2022-09-01") ~ "Early_2022",
    Rando >= as.Date("2022-09-01") & Rando < as.Date("2023-01-01") ~ "Late_2022",
    
    Rando >= as.Date("2023-01-01") & Rando < as.Date("2023-09-01") ~ "Early_2023",
    Rando >= as.Date("2023-09-01") & Rando < as.Date("2024-01-01") ~ "Late_2023",
    
    Rando >= as.Date("2024-01-01") & Rando < as.Date("2024-09-01") ~ "Early_2024",
    Rando >= as.Date("2024-09-01") & Rando < as.Date("2025-01-01") ~ "Late_2024"
  )) %>%
  mutate(Study_Year = case_when(
    Rando < as.Date("2021-01-01") ~ "Late_2020",
    
    Rando >= as.Date("2021-01-01") & Rando < as.Date("2022-01-01") ~ "2021",
    
    Rando >= as.Date("2022-01-01") & Rando < as.Date("2023-01-01") ~ "2022",
    
    Rando >= as.Date("2023-01-01") & Rando < as.Date("2024-01-01") ~ "2023",
    
    Rando >= as.Date("2024-01-01") & Rando < as.Date("2025-01-01") ~ "2024",
    
  ))
df.full$Sexe <- as.factor(df.full$Sexe)
df.full$CENTRE <- as.factor(df.full$CENTRE)
df.full$dx <- as.factor(df.full$dx)
df.full$Study_Wave <- as.factor(df.full$Study_Wave)
df.full$Rando <- as.factor(df.full$Rando)
df.full$Study_Year <- as.factor(df.full$Study_Year)
df.full$ASSIDUITE_SOIN_crit <- as.factor(df.full$ASSIDUITE_SOIN_crit)

table(df.full$dx, useNA = "always")
table(df.full$Study_Wave, useNA = "always")
table(df.full$Rando, useNA = "always")
table(df.full$Study_Year, useNA = "always")

# 1. Identify the periods you have (e.g., 1 to 6)
n_periods <- 6 

# 2. Iterate through each period
for (i in 1:n_periods) {
  
  # Column names for the current iteration
  mois_col  <- paste0("MOIS_ANNEE_", i)
  ratio_col <- paste0("RATIO_", i)
  
  # Get the unique month values present in this specific MOIS_ANNEE_i column
  # (Ignoring NAs)
  unique_months <- unique(df.full[[mois_col]])
  unique_months <- unique_months[!is.na(unique_months)]
  
  for (unique_months_idx in unique_months) {
    # Create the new column name
    new_col_name <- paste0("RATIO_H_T_", unique_months_idx)
    
    # If the column doesn't exist yet, initialize with NA
    if (!(new_col_name %in% names(df.full))) {
      df.full[[new_col_name]] <- NA_real_
    }
    
    # Logic: If MOIS_ANNEE_i equals 'm', put RATIO_i into 'RATIO_H_T_m'
    mask <- !is.na(df.full[[mois_col]]) & df.full[[mois_col]] == unique_months_idx
    df.full[mask, new_col_name] <- df.full[mask, ratio_col]
  }
}

#Impute missing data
if (
df.full %>% 
  select(all_of(c(baseline_CJP,paste0(all_outcomes_CJS, "_t0")))) %>% 
  anyNA()
) { 
  # https://journals.sagepub.com/doi/full/10.1177/0962280216683570 
  # For missing data in a baseline covariate, 
  # simpler approaches such as mean imputation and the missing indicator method 
  # can outperform MI. 
  df.full_imputed <- df.full
  # 1. Loop through each variable in your baseline vector
  for (var in c(baseline_CJP,paste0(all_outcomes_CJS, "_t0"))) {
    
    # Check if the column actually has missing data
    if (any(is.na(df.full_imputed[[var]]))) {
      
      # Define the name for the new indicator variable
      indicator_name <- paste0(var, "_NA")
      
      # 2. Create the missing indicator (1 if NA, 0 otherwise)
      df.full_imputed[[indicator_name]] <- ifelse(is.na(df.full_imputed[[var]]), 1, 0)
      
      # 3. Fill the original column with its mean
      # Note: we use na.rm = TRUE to calculate the mean of existing values
      var_mean <- mean(df.full_imputed[[var]], na.rm = TRUE)
      df.full_imputed[[var]][is.na(df.full_imputed[[var]])] <- var_mean
      
      message(paste("Processed:", var, "- Created:", indicator_name))
    }
  }
  
  #even if MI is FALSE still perform MI
  if (run_ITT_MI) {
    split_data <- split(df.full, df.full$BRAS)
    imputed_list <- lapply(split_data, function(data) {
    mice::mice(
        data %>% 
          select(
            all_of(baseline_CJP), 
            starts_with(all_outcomes_CJP),
            starts_with(all_outcomes_CJS),
        
            all_of(mediator)
            ),
          m = m,               
          maxit = maxit,
          seed = 123)
        })
    completed_data_list <- lapply(imputed_list, function(imputed_data) {
      lapply(1:m, function(m_idx) complete(imputed_data, m_idx))
      })
    pooled_imputations <- lapply(1:m, function(m_idx) {
      rbind(completed_data_list[["0"]][[m_idx]] %>% mutate(BRAS=0),
            completed_data_list[["1"]][[m_idx]] %>% mutate(BRAS=1)
            )
    })
    df.full_imputed_MI <- map(pooled_imputations, ~{
      .x %>%
      arrange(df.full$StudySubjectID)
    }) 
  }
  
  }
baseline_CJP_MI <- baseline_CJP
baseline_CJP <- c(baseline_CJP, "STORI_CROISSANCE_t0_NA")
baseline_CJP <- c(baseline_CJP, "ISMI_SC_TOT_t0_NA")


run_tmle <- function(df, outcome, baseline, cens, MI=FALSE, anal="ITT", CompCovar=FALSE) {
  
  df_mdf <- df %>%
    mutate(SORTIE=df.full$SORTIE,
           cens = case_when(is.na(df.full[[outcome]]) ~ 0,
                            !is.na(df.full[[outcome]]) ~ 1)
    )
  
  if (anal=="PP") {
    df_mdf <- df_mdf %>%
      mutate(cens=case_when(cens==0 ~ 0,
                            cens==1 & SORTIE==0 ~ 1,
                            cens==1 & SORTIE==1 ~ 0)
             )
  }
    
  if (is.null(cens) & MI==FALSE ) {df_mdf <- df_mdf %>% filter(cens==1)} 
  
  #rmv missing covariate
  #use of df.full (not imputed)
  if(CompCovar==TRUE) {
    df_mdf <- df_mdf %>% 
      drop_na(all_of(baseline))
  }
  
  set.seed(123)
  tmle_control <- lmtp::lmtp_tmle(data = df_mdf,
                                  trt = "BRAS",
                                  outcome = outcome, 
                                  baseline = baseline,
                                  outcome_type = "continuous",
                                  cens=cens, 
                                  shift = static_binary_off, 
                                  folds=nfolds,
                                  learners_trt=learners_trt,
                                  learners_outcome=learners_outcome)
  set.seed(123)
  tmle_trt <- lmtp::lmtp_tmle(data = df_mdf,
                              trt = "BRAS",
                              outcome = outcome, 
                              baseline = baseline,
                              outcome_type = "continuous",
                              cens=cens, 
                              shift = static_binary_on,
                              folds=nfolds,
                              learners_trt=learners_trt,
                              learners_outcome=learners_outcome)
  #lmtp_contrast(tmle_trt, ref=tmle_control) %>% print
  
  tmle_contrast <- lmtp_contrast(tmle_trt, ref=tmle_control)
  
  #moderation analysis
  df_mdf_eif <- df_mdf %>% cbind(eif=tmle_contrast$eifs)
  # Identify which variables in baseline vector contain "_NA"
  vars_to_exclude <- grep("_NA|Study_Year|CENTRE", baseline, value = TRUE)  # Subtract those from the original baseline list
  moderators <- setdiff(baseline, vars_to_exclude)
  #print(moderators)
  formula <- formula(paste0("eif ~ ", paste0(moderators, collapse = "+")))
  effect_mdf <- lm(formula, df_mdf_eif)
  
  return(list(df_mdf_eif=df_mdf_eif,
              tmle_control=tmle_control,
              tmle_trt=tmle_trt,
              tmle_contrast=tmle_contrast,
              effect_mdf=effect_mdf
              )
         )
}

extract_Rubin <- function(which_list) {
  # pool the results using Rubin's rule 
  # 1. Extract thetas from the 'vals' sub-list of each of the 5 elements
  thetas <- sapply(which_list, function(x) x$tmle_contrast$vals$theta)
  
  # 2. Extract standard errors from the same 'vals' sub-list
  # Squaring them to provide the variance (standard for MI results)
  variances <- sapply(which_list, function(x) (x$tmle_contrast$vals$std.error)^2)
  
  # 3. Create the results_miint object
  results_miint <- data.frame(
    tmle_est_imp = thetas,
    tmle_var_imp = variances
  )
  
  # Verify the result
  #print(results_miint)
  
  tmle_est <- mean(results_miint$tmle_est_imp)
  var_inter <- (1/(m-1))*sum((results_miint$tmle_est_imp-tmle_est)^2)
  var_intra <- mean(results_miint$tmle_var_imp)
  tmle_var <- var_intra + (1 + 1/m)*var_inter
  r<-(1+1/m)*var_inter/var_intra
  vvv<-(m-1)*(1+1/r)^2 
  tal<-qt(0.025,df=vvv,lower.tail=F) 
  tmle_lCI <- tmle_est - tal*sqrt(tmle_var)
  tmle_uCI <- tmle_est + tal*sqrt(tmle_var) 
  tmle_pval<-2*(1-pt(q=abs(tmle_est/sqrt(tmle_var)),df=vvv))  
  TMLEres_MI_int <-data.frame(cbind(tmle_est,tmle_var,tmle_lCI,tmle_uCI,tmle_pval)) 
  
  return(TMLEres_MI_int)
}

# Vf CJP ----

# vérification EFFECTUEES_TOTAL <= PREVUES_TOTAL
df.full %>%
  filter(EFFECTUEES_TOTAL > PREVUES_TOTAL) %>%
  select(StudySubjectID,EFFECTUEES_TOTAL,PREVUES_TOTAL)

# vérification EFFECTUEES_TOTAL+ABSENCE_TOTAL ~= PREVUES_TOTAL
df.full %>%
  filter(abs((EFFECTUEES_TOTAL  + ABSENCE_TOTAL) - PREVUES_TOTAL) >= 5) %>% # marge de 5h..
  select(StudySubjectID,EFFECTUEES_TOTAL,ABSENCE_TOTAL,PREVUES_TOTAL)

# Nombre de sorties d'étude avec pas assez de séances
df.full %>% 
  filter(SORTIE == 1 & ASSIDUITE_SOIN_crit == "non") %>% 
  select("StudySubjectID","SORTIE", "ETAPE", "MOTIF", "MOTIF_OTH","ASSIDUITE_SOIN_crit") %>% 
  pull(StudySubjectID)

# sorties d'étude et critère assiduité des séances
df.full %>% 
  filter(SORTIE == 1 & ASSIDUITE_SOIN_crit == "oui") %>% 
  select("StudySubjectID","SORTIE", "ETAPE", "MOTIF", "MOTIF_OTH","ASSIDUITE_SOIN_crit") %>% 
  pull(StudySubjectID)


#Worst case scenario
# mis à 0 pour les sortie d'étude
# df.full %>%
#   mutate(
#     RATIO_H_T = ifelse(is.na(RATIO_H_T) & StudySubjectID %in% idpat.sortie, 0, RATIO_H_T)
#   ) -> df.full

#Vf CJS ----

# 1. Evaluer l’impact sur la cognition sociale, l’estime de soi, l’auto-stigmatisation et le rétablissement de l’entrainement de la cognition sociale pratiqué dans l’environnement professionnel ;

df.full %>% 
  select(!matches("SIGMA")) -> df.full

varnames.eval <- c("SERS_POS", "SERS_NEG", "SERS_SC_TOT", 
                   "ISMI_ALIEN", "ISMI_APPROB_STEREO", "ISMI_EXP_DISCRIM", "ISMI_RET_SOC", "ISMI_INTERN_STIGMA", "ISMI_RESIST_STIGMA", "ISMI_SC_TOT",
                   "STORI_MORATOIRE", "STORI_CONSCIENCE", "STORI_PREPARATION", "STORI_RECONSTRUCTION", "STORI_CROISSANCE", "STORI_ST_RETABLISSEMENT", 
                   "STAIA_TOT", "STAIB_TOT",
                   "ACSO", "ACSO_H", 
                   "MASC_NB", "MASC_SIGMA", 
                   "AIHQ_HB_NB", "AIHQ_HB_SIGMA", "AIHQ_ATTRIB_RESP_NB", "AIHQ_ATTRIB_RESP_SIGMA", "AIHQ_AB_NB", "AIHQ_AB_SIGMA", 
                   "PERSO_FLU_PONC", "PERSO_INTER_LIB", "PERSO_INTER_INDIC", "PERSO_INTER_TOTAL", "PERSO_CONV_SOC", 
                   "TREF_BR_DEGOUT_NB", "TREF_BR_DEGOUT_SIGMA", "TREF_BR_MEPRIS_NB", "TREF_BR_MEPRIS_SIGMA", "TREF_BR_JOIE_NB", "TREF_BR_JOIE_SIGMA", "TREF_BR_PEUR_NB", "TREF_BR_PEUR_SIGMA", "TREF_BR_TRISTESSE_NB", "TREF_BR_TRISTESSE_SIGMA", "TREF_BR_COLERE_NB", "TREF_BR_COLERE_SIGMA", "TREF_BR_TOTAL_NB", "TREF_BR_TOTAL_SIGMA", "TREF_SD_DEGOUT_NB", "TREF_SD_DEGOUT_SIGMA", "TREF_SD_MEPRIS_NB", "TREF_SD_MEPRIS_SIGMA", "TREF_SD_JOIE_NB", "TREF_SD_JOIE_SIGMA", "TREF_SD_PEUR_NB", "TREF_SD_PEUR_SIGMA", "TREF_SD_TRISTESSE_NB", "TREF_SD_TRISTESSE_SIGMA", "TREF_SD_COLERE_NB", "TREF_SD_COLERE_SIGMA", "TREF_SD_TOTAL_NB", "TREF_SD_TOTAL_SIGMA")

varnames.eval <- c( "SERS_SC_TOT", "ISMI_SC_TOT", "STORI_CROISSANCE", 
                    "STAIA_TOT", "STAIB_TOT",
                    "ACSO", "ACSO_H", 
                    "MASC_NB", 
                    "AIHQ_HB_NB","AIHQ_ATTRIB_RESP_NB", "AIHQ_AB_NB",  
                    "PERSO_FLU_PONC", "PERSO_INTER_LIB", "PERSO_INTER_INDIC", "PERSO_INTER_TOTAL", "PERSO_CONV_SOC", 
                    "TREF_BR_TOTAL_NB","TREF_SD_TOTAL_NB")

# Conversion en long seulement le time
df.full %>%
  select(c(StudySubjectID, BRAS, starts_with(varnames.eval), all_of(baseline_CJS))) %>%
  pivot_longer(
    cols = -c(StudySubjectID, BRAS, all_of(baseline_CJS)),
    names_to = c(".value", "time"),
    names_pattern = "(.*)_t(\\d+)"
  ) |>
  mutate(
    time = factor(time, levels = c("0", "1", "2")),
    BRAS = factor(BRAS)
  ) %>%
  arrange(StudySubjectID, time) -> df.full.time.long

# Conversion en long le time et le type de variable
df.full.time.long  %>% 
  mutate(across(-c(StudySubjectID, BRAS, time), as.character)) %>%
  pivot_longer(
    cols = -c(StudySubjectID,BRAS,time),
    names_to = "variable"
  ) %>% 
  mutate(variable = factor(variable, levels = varnames.eval)
  ) -> df.full.timeVariable.long

# plot pour chaque echelle/test
plot.echelle <- function(varnames.plot, varnames.prefix, plot.save = F){
  
  plot_df <- df.full.time.long %>%
    select(StudySubjectID, BRAS, time, varnames.plot) %>%
    pivot_longer(
      cols =  starts_with(varnames.prefix),
      names_to = "measure",
      values_to = "score"
    ) 
  ggplot(plot_df, aes(
    x = time,
    y = score,
    color = BRAS,
    linetype = BRAS,
    group = interaction(measure, BRAS)
  )) +
    stat_summary(fun = mean, geom = "line", size = 1.2) +
    stat_summary(fun = mean, geom = "point", size = 2.5) +
    scale_linetype_manual(values =   c("0" = "dashed", "1" = "solid") ) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      title = paste("Évolution des scores",varnames.prefix,"par temps et par bras"),
      x = "",
      y = "Score moyen",
      color = "Dimension:",
      linetype = "Bras:"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.key.width = unit(1.2, "cm"),    #  élargit les segments
      legend.text = element_text(size = 10, family = "serif"),
      legend.title = element_text(size = 11, face = "bold", family = "serif")
    ) -> plot.res
  
  if (!plot.save) {
    return(plot.res)
  } else {
    ggsave(
      filename = paste0("evol_", varnames.prefix, ".png"),
      plot = plot.res,
      width = 16, height = 10, dpi = 300
    )
  }
}

#analysis
df.full.time.long <- df.full.time.long %>%
  mutate(time = as.factor(time))
anal_lmer <- function(varnames.plot) {
  
  # 2. Build the formula parts
  # Fixed effects: Interaction (BRAS * time) + all baseline variables
  fixed_parts <- c("BRAS * time", baseline_CJS)
  
  # 3. Create the formula object
  # reformulate(termlabels, response)
  base_formula <- reformulate(termlabels = fixed_parts, 
                              response = varnames.plot)
  
  # 4. Add the random effect component
  final_formula <- as.formula(paste(format(base_formula), " (1 | StudySubjectID)"))
  
  # 5. Run the model using lmerTest for p-values
  model_result <- lmerTest::lmer(final_formula, data = df.full.time.long)
  
  # 6. View results
  summary(model_result) %>% print
  
}

# SERS
varnames.prefix <- "SERS"
varnames.plot <-  grep(paste0("^",varnames.prefix), names(df.full.time.long), value = TRUE, perl = TRUE)
plot.echelle(varnames.plot = varnames.plot, varnames.prefix = varnames.prefix, plot.save = F)
anal_lmer(varnames.plot = varnames.plot)

# ISMI
varnames.prefix <- "ISMI"
varnames.plot <-  grep(paste0("^",varnames.prefix), names(df.full.time.long), value = TRUE, perl = TRUE)
plot.echelle(varnames.plot = varnames.plot, varnames.prefix = varnames.prefix, plot.save = F)
anal_lmer(varnames.plot = varnames.plot)

# STORI
varnames.prefix <- "STORI"
varnames.plot <-  grep(paste0("^",varnames.prefix), names(df.full.time.long), value = TRUE, perl = TRUE)
plot.echelle(varnames.plot = varnames.plot, varnames.prefix = varnames.prefix, plot.save = F)
anal_lmer(varnames.plot = varnames.plot)

# STAI
varnames.prefix <- "STAIA"
varnames.plot <-  grep(paste0("^",varnames.prefix), names(df.full.time.long), value = TRUE, perl = TRUE)
plot.echelle(varnames.plot = varnames.plot, varnames.prefix = varnames.prefix, plot.save = F)
anal_lmer(varnames.plot = varnames.plot)


# ACSO
varnames.prefix <- "ACSO"
varnames.plot <-  grep(paste0("^",varnames.prefix), names(df.full.time.long), value = TRUE, perl = TRUE)
plot.echelle(varnames.plot = varnames.plot, varnames.prefix = varnames.prefix, plot.save = F)
anal_lmer(varnames.plot = varnames.plot)

# MASC
varnames.prefix <- "MASC"
varnames.plot <-  grep(paste0("^",varnames.prefix), names(df.full.time.long), value = TRUE, perl = TRUE)
plot.echelle(varnames.plot = varnames.plot, varnames.prefix = varnames.prefix, plot.save = F)
anal_lmer(varnames.plot = varnames.plot)

# AIHQ
varnames.prefix <- "AIHQ"
varnames.plot <-  grep(paste0("^",varnames.prefix), names(df.full.time.long), value = TRUE, perl = TRUE)
plot.echelle(varnames.plot = varnames.plot, varnames.prefix = varnames.prefix, plot.save = F)

# PERSO
varnames.prefix <- "PERSO"
varnames.plot <-  grep(paste0("^",varnames.prefix), names(df.full.time.long), value = TRUE, perl = TRUE)
plot.echelle(varnames.plot = varnames.plot, varnames.prefix = varnames.prefix, plot.save = F)

# TREF
varnames.prefix <- "TREF_BR"
#varnames.plot <-  grep("^TREF_BR(?!.*TOTAL)", names(df.full.time.long), value = TRUE, perl = TRUE)
varnames.plot <-  grep(paste0("^",varnames.prefix), names(df.full.time.long), value = TRUE, perl = TRUE)
plot.echelle(varnames.plot = varnames.plot, varnames.prefix = varnames.prefix, plot.save = F)
anal_lmer(varnames.plot = varnames.plot)

# ITT simple ----
df.full %>%
  # select only the variables you want to summarize
  select(BRAS, baseline_CJS, all_of(all_outcomes_CJP) , all_of(starts_with(all_outcomes_CJS))) %>% 
  tbl_summary(
    by = BRAS,
    # Use 'continuous' type for ratios to get Mean (SD) or Median (IQR)
    type = all_continuous() ~ "continuous",
    statistic = all_continuous() ~ "{mean} ({sd})",
    digits = all_continuous() ~ 2,
    missing = "always" # Optional: hides the 'Unknown' row if you have NAs
  ) %>%
  add_p() %>% # Optional: adds a p-value to compare ratios between groups
  bold_labels()

# PP simple ----
#PP = ITT - SORTIE
df.full %>% 
  filter(SORTIE==0) %>%
  # select only the variables you want to summarize
  select(BRAS, baseline_CJS, all_of(all_outcomes_CJP) , all_of(starts_with(all_outcomes_CJS))) %>% 
  tbl_summary(
    by = BRAS,
    # Use 'continuous' type for ratios to get Mean (SD) or Median (IQR)
    type = all_continuous() ~ "continuous",
    statistic = all_continuous() ~ "{mean} ({sd})",
    digits = all_continuous() ~ 2,
    missing = "always" # Optional: hides the 'Unknown' row if you have NAs
  ) %>%
  add_p() %>% # Optional: adds a p-value to compare ratios between groups
  bold_labels()


#Compare ITT-CC and missing ----

# 1. Create the missingness indicator
df_compare <- df.full %>%
  mutate(
    ratio_status = ifelse(is.na(RATIO_TOT), "Missing TOTAL_RATIO", "Present TOTAL_RATIO")
  ) %>%
  # 2. Select only the variables you want to summarize + the grouping variable
  select(
    ratio_status, 
    baseline_CJP_MI, 
    baseline_CJS, 
    all_outcomes_CJP, 
    paste0(all_outcomes_CJS,"_t0"), 
    BRAS
  )

# 3. Create the comparison table
comparison_table <- df_compare %>%
  tbl_summary(
    by = ratio_status, # Split the table by missingness status
    missing = "always"
  ) %>%
  add_p() %>%                # Add p-values to see if missingness is non-random
  add_overall() %>%          # Optional: add a column for the whole sample
  bold_labels() %>%
  italicize_levels()

# 4. View the table
comparison_table


#Compare CC and total sample ----
# 1. Define the variables to summarize (to keep code clean)
vars_to_summarize <- c(
  baseline_CJP_MI, 
  baseline_CJS, 
  all_outcomes_CJP, 
  paste0(all_outcomes_CJS, "_t0"), 
  "BRAS"
)

# 2. Create Table A: The Full Sample (The Reference)
tab_full <- df.full %>%
  select(all_of(vars_to_summarize)) %>%
  tbl_summary(
    missing = "always"
  ) %>%
  modify_header(label = "**Variable**", stat_0 = "**Full Sample (N={n})**")

# 3. Create Table B: Only the 'Present' Ratio group
tab_present <- df.full %>%
  filter(!is.na(RATIO_TOT)) %>% # Keep only present
  select(all_of(vars_to_summarize)) %>%
  tbl_summary(
    missing = "always"
  ) %>%
  modify_header(stat_0 = "**Present Ratio (N={n})**")

# 4. Merge them into one table
# This displays them side-by-side for direct comparison
comparison_table <- tbl_merge(
  tbls = list(tab_full, tab_present),
  tab_spanner = c("**Reference**", "**Sub-group**")
) %>%
  bold_labels() %>%
  italicize_levels()

# 5. View the table
comparison_table

stop("next step: lmtp analysis")

#Compare PP-CC and missing ----

# 1. Create the missingness indicator
df_compare <- df.full %>%
  # 2. Select only the variables you want to summarize + the grouping variable
  select(
    SORTIE, 
    baseline_CJP_MI, 
    baseline_CJS, 
    all_outcomes_CJP, 
    paste0(all_outcomes_CJS,"_t0"), 
    BRAS
  )

# 3. Create the comparison table
comparison_table <- df_compare %>%
  tbl_summary(
    by = SORTIE, # Split the table by missingness status
    missing = "always"
  ) %>%
  add_p() %>%                # Add p-values to see if missingness is non-random
  add_overall() %>%          # Optional: add a column for the whole sample
  bold_labels() %>%
  italicize_levels()

# 4. View the table
comparison_table


