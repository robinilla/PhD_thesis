# ---------------------------------------------
#  PhD Thesis: Large-scale monitoring of wild mammal abundance: 
#  modeling frameworks for structured and unstructured data
#  Chapter 3. Modeling frameworks for abundance patterns of European wild mammals
#  Author: Sonia Illanas
#  Institution: Institute of Game and Wildlife Research
#  Date of last modification: 13/02/2026
# ---------------------------------------------
## R version 4.5.2
## tidyverse version: 2.0.0
## sf version: 1.0-21
## viridis version: 0.6.5
## MASS version: 7.3-65
## glmmTMB version: 1.1.13
## Hmisc version: 5.2-4
## ggpubr version: 0.6.1
## patchwork version: 1.3.2
## ggh4x version: 0.3.1

library(tidyverse)
library(sf)
library(viridis)
library(MASS)
library(glmmTMB)
library(Hmisc)

# rm(list=ls())
cal_data <- function(pred, obs, g = 9){
  bins <-cut2(pred, g = g) # defining bins (percentiles) on the predicted
  tibble(
    pred_mean = tapply(pred, bins, mean), # calculating mean values for predicted
    obs_mean = tapply(obs, bins, mean), # calculating mean values for observed
    n = tapply (obs, bins, length) )
}

#Charge widespread species  data

# BADGER
# --------------------
data<-st_read(dsn="data_20260207.gpkg", layer="data_spatial_Badger")
load("densityModel_Badger_glmmmTMB_Europe.Rdata")

set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]

#Calibration plots - validation data
data_val<-data_val %>% mutate(Pred_model=predict(Badger_DensityModel, data_val, type="response"), Pred_species=Pred_model/10000) #%>% rename(Bioregion=BR_todasr1)

data_pos<-data_val %>% filter(dens>0)
cpwBadger<-bind_rows(
  cal_data(pred = data_val$Pred_species, obs=data_val$dens) %>% mutate(Type = "All data")#,
  # cal_data(pred = data_pos$Pred_species, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
) %>% 
  mutate(species="f) Badger")


#2.3.1.) Calibration plots for each bioregion
data_val_BR1<-data_val %>% filter(Bioregion==1)
# data_val_BR2<-data_val %>% filter(Bioregion==2)
data_val_BR3<-data_val %>% filter(Bioregion==3)
data_val_BR4<-data_val %>% filter(Bioregion==4)

cpwBadger.br<-
bind_rows(cal_data(pred = data_val_BR1$Pred_species, obs=data_val_BR1$dens) %>% mutate(Type = "All data", Bioregion="Northern bioregion"),
          # cal_data(pred = data_val_BR2$Pred_species, obs=data_val_BR2$dens) %>% mutate(Type = "All data", Bioregion="Southern bioregion"),
          cal_data(pred = data_val_BR3$Pred_species, obs=data_val_BR3$dens) %>% mutate(Type = "All data", Bioregion="Eastern bioregion"),
          cal_data(pred = data_val_BR4$Pred_species, obs=data_val_BR4$dens) %>% mutate(Type = "All data", Bioregion="Western bioregion")
  ) %>% 
  mutate(species="f) Badger")

rm(list = setdiff(ls(), c("cal_data", "cpwBadger.br", "cpwBadger")))
# --------------------

# RED FOX
# --------------------
data<-st_read(dsn="data_20260207.gpkg", layer="data_spatial_RedFox")
load("densityModel_RedFox_GLM_xy_without0_20221220.Rdata")
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]

data_val<-data_val %>% mutate(Pred_model=predict(RedFox_DensityModel, data_val, type="response"), Pred_species=Pred_model/10000) #%>% rename(Bioregion=BR_todasr1)

data_pos<-data_val %>% filter(dens>0)
cpwFox<-bind_rows(
  cal_data(pred = data_val$Pred_species, obs=data_val$dens) %>% mutate(Type = "All data")#,
  # cal_data(pred = data_pos$Pred_species, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
) %>%  mutate(species = "e) Red fox")


#2.3.1.) Calibration plots for each bioregion
data_val_BR1<-data_val %>% filter(Bioregion==1)
data_val_BR2<-data_val %>% filter(Bioregion==2)
data_val_BR3<-data_val %>% filter(Bioregion==3)
data_val_BR4<-data_val %>% filter(Bioregion==4)

cpwFox.br<-
  bind_rows(cal_data(pred = data_val_BR1$Pred_species, obs=data_val_BR1$dens) %>% mutate(Type = "All data", Bioregion="Northern bioregion"),
            cal_data(pred = data_val_BR2$Pred_species, obs=data_val_BR2$dens) %>% mutate(Type = "All data", Bioregion="Southern bioregion"),
            cal_data(pred = data_val_BR3$Pred_species, obs=data_val_BR3$dens) %>% mutate(Type = "All data", Bioregion="Eastern bioregion"),
            cal_data(pred = data_val_BR4$Pred_species, obs=data_val_BR4$dens) %>% mutate(Type = "All data", Bioregion="Western bioregion")
            ) %>% 
  mutate(species = "e) Red fox")

rm(list = setdiff(ls(), c("cal_data", "cpwBadger.br", "cpwBadger", "cpwFox", "cpwFox.br")))

# --------------------

# WILD BOAR
# --------------------
data<-st_read(dsn="F:/IREC/PhD/D/Docs/Sonia/0_PhD/7_DocTesis/Chapter_ENETWILD/GIS/data_20260207.gpkg", layer="data_spatial_WildBoar")
df_br1<-data %>% filter(Bioregion==1); st_geometry(df_br1)<-NULL
df_br2<-data %>% filter(Bioregion==2); st_geometry(df_br2)<-NULL
df_br3<-data %>% filter(Bioregion==3); st_geometry(df_br3)<-NULL
df_br4<-data %>% filter(Bioregion==4); st_geometry(df_br4)<-NULL

load(file="F:/IREC/PhD/D/Docs/Sonia/0_PhD/7_DocTesis/Chapter_ENETWILD/GIS/ModelResults/densityModel_WildBoar_NorthernBioregion_Europe.Rdata")
load(file="F:/IREC/PhD/D/Docs/Sonia/0_PhD/7_DocTesis/Chapter_ENETWILD/GIS/ModelResults/densityModel_WildBoar_SouthernBioregion_Europe.Rdata")
load(file="F:/IREC/PhD/D/Docs/Sonia/0_PhD/7_DocTesis/Chapter_ENETWILD/GIS/ModelResults/densityModel_WildBoar_EasternBioregion_Europe.Rdata")
load(file="F:/IREC/PhD/D/Docs/Sonia/0_PhD/7_DocTesis/Chapter_ENETWILD/GIS/ModelResults/densityModel_WildBoar_WesternBioregion_Europe.Rdata")

#1.1 Creamos data set de calibrado y de validacion (80% / 20%) para cada Bioregion
set.seed(500)
s_model_BR1 <- sample(nrow(df_br1), nrow(df_br1)*0.8) ; data_tra_BR1 <- df_br1[s_model_BR1,]; data_val_BR1 <- df_br1[-s_model_BR1,]
s_model_BR2 <- sample(nrow(df_br2), nrow(df_br2)*0.8) ; data_tra_BR2 <- df_br2[s_model_BR2,]; data_val_BR2 <- df_br2[-s_model_BR2,]
s_model_BR3 <- sample(nrow(df_br3), nrow(df_br3)*0.8) ; data_tra_BR3 <- df_br3[s_model_BR3,]; data_val_BR3 <- df_br3[-s_model_BR3,]
s_model_BR4 <- sample(nrow(df_br4), nrow(df_br4)*0.8) ; data_tra_BR4 <- df_br4[s_model_BR4,]; data_val_BR4 <- df_br4[-s_model_BR4,]

#Calibration plots  validation data
data_val_BR1<-data_val_BR1 %>% mutate(Pred_model=predict(WildBoar_NorthernBioregion_DensityModel, data_val_BR1, type="response"), Pred_Red=Pred_model/10000)
data_val_BR2<-data_val_BR2 %>% mutate(Pred_model=predict(WildBoar_SouthernBioregion_DensityModel, data_val_BR2, type="response"), Pred_Red=Pred_model/10000)
data_val_BR3<-data_val_BR3 %>% mutate(Pred_model=predict(WildBoar_EasternBioregion_DensityModel, data_val_BR3, type="response"), Pred_Red=Pred_model/10000)
data_val_BR4<-data_val_BR4 %>% mutate(Pred_model=predict(WildBoar_WesternBioregion_DensityModel, data_val_BR4, type="response"), Pred_Red=Pred_model/10000)

cpwWB.br<-
  bind_rows(cal_data(pred = data_val_BR1$Pred_Red, obs=data_val_BR1$dens) %>% mutate(Type = "All data", Bioregion="Northern bioregion"),
            cal_data(pred = data_val_BR2$Pred_Red, obs=data_val_BR2$dens) %>% mutate(Type = "All data", Bioregion="Southern bioregion"),
            cal_data(pred = data_val_BR3$Pred_Red, obs=data_val_BR3$dens) %>% mutate(Type = "All data", Bioregion="Eastern bioregion"),
            cal_data(pred = data_val_BR4$Pred_Red, obs=data_val_BR4$dens) %>% mutate(Type = "All data", Bioregion="Western bioregion")#,
            ) %>%
  mutate(species = "a) Wild boar")

rm(list = setdiff(ls(), c("cal_data", "cpwBadger.br", "cpwBadger", "cpwFox", "cpwFox.br", "cpwWB.br")))
# --------------------

# RED DEER
# --------------------
data<-st_read(dsn="data_20260207.gpkg", 
              layer="data_spatial_RedDeer")
load(file="densityModel_RedDeer_Europe.Rdata")

set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]

#Calibration plots  validation data
data_val<-data_val %>% mutate(Pred_model=predict(RedDeer_DensityModel, data_val, type="response"), Pred_Red=Pred_model/10000)

data_pos<-data_val %>% filter(dens>0)
cpwRedDeer<-bind_rows(
  cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data")#,
  # cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
) %>%  mutate(species="c) Red deer")

#2.3.1.) Calibration plots for each bioregion
data_val_BR1<-data_val %>% filter(Bioregion==1)
data_val_BR2<-data_val %>% filter(Bioregion==2)
data_val_BR3<-data_val %>% filter(Bioregion==3)
data_val_BR4<-data_val %>% filter(Bioregion==4)

cpwRedDeer.br<-
  bind_rows(cal_data(pred = data_val_BR1$Pred_Red, obs=data_val_BR1$dens) %>% mutate(Type = "All data", Bioregion="Northern bioregion"),
            cal_data(pred = data_val_BR2$Pred_Red, obs=data_val_BR2$dens) %>% mutate(Type = "All data", Bioregion="Southern bioregion"),
            cal_data(pred = data_val_BR3$Pred_Red, obs=data_val_BR3$dens) %>% mutate(Type = "All data", Bioregion="Eastern bioregion"),
            cal_data(pred = data_val_BR4$Pred_Red, obs=data_val_BR4$dens) %>% mutate(Type = "All data", Bioregion="Western bioregion")
            ) %>% 
  mutate(species="c) Red deer")

rm(list = setdiff(ls(), c("cal_data", "cpwBadger.br", "cpwBadger", "cpwFox", "cpwFox.br", "cpwWB.br", "cpwRedDeer", "cpwRedDeer.br")))
# --------------------

# ROE DEER
# --------------------
data<-st_read(dsn="data_20260207.gpkg", layer="data_spatial_RoeDeer")
load(file="densityModel_RoeDeer_Europe.Rdata")
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]

#Calibration plots  validation data
data_val<-data_val %>% mutate(Pred_model=predict(RoeDeer_DensityModel, data_val, type="response"), Pred_Red=Pred_model/10000)

data_pos<-data_val %>% filter(dens>0)
cpwRoe<-bind_rows(
  cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data")#,
  # cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
) %>% mutate(species="b) Roe deer")

#2.3.1.) Calibration plots for each bioregion
data_val_BR1<-data_val %>% filter(Bioregion==1)
data_val_BR2<-data_val %>% filter(Bioregion==2)
data_val_BR3<-data_val %>% filter(Bioregion==3)
data_val_BR4<-data_val %>% filter(Bioregion==4)

cpwRoe.br<-
  bind_rows(cal_data(pred = data_val_BR1$Pred_Red, obs=data_val_BR1$dens) %>% mutate(Type = "All data", Bioregion="Northern bioregion"),
            cal_data(pred = data_val_BR2$Pred_Red, obs=data_val_BR2$dens) %>% mutate(Type = "All data", Bioregion="Southern bioregion"),
            cal_data(pred = data_val_BR3$Pred_Red, obs=data_val_BR3$dens) %>% mutate(Type = "All data", Bioregion="Eastern bioregion"),
            cal_data(pred = data_val_BR4$Pred_Red, obs=data_val_BR4$dens) %>% mutate(Type = "All data", Bioregion="Western bioregion")
            ) %>% 
  mutate(species="b) Roe deer")

rm(list = setdiff(ls(), c("cal_data", "cpwBadger.br", "cpwBadger", "cpwFox", "cpwFox.br", "cpwWB.br", "cpwRedDeer", "cpwRedDeer.br",
                          "cpwRoe", "cpwRoe.br")))
# --------------------

# FALLOW DEER
# --------------------
data<-st_read(dsn="data_20260207.gpkg", 
              layer="data_spatial_FallowDeer")
load(file="densityModel_Fd_Europe.Rdata")
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]

#Calibration plots  validation data
data_val<-data_val %>% mutate(Pred_model=predict(Fd_DensityModel, data_val, type="response"), Pred_Red=Pred_model/10000)

data_pos<-data_val %>% filter(dens>0)
cpwFd<-bind_rows(
  cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data")#,
  # cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
) %>% mutate(species="d) Fallow deer")

#2.3.1.) Calibration plots for each bioregion
data_val_BR1<-data_val %>% filter(Bioregion==1)
data_val_BR2<-data_val %>% filter(Bioregion==2)
data_val_BR3<-data_val %>% filter(Bioregion==3)
data_val_BR4<-data_val %>% filter(Bioregion==4)

cpwFd.br<-
  bind_rows(cal_data(pred = data_val_BR1$Pred_Red, obs=data_val_BR1$dens) %>% mutate(Type = "All data", Bioregion="Northern bioregion"),
            cal_data(pred = data_val_BR2$Pred_Red, obs=data_val_BR2$dens) %>% mutate(Type = "All data", Bioregion="Southern bioregion"),
            cal_data(pred = data_val_BR3$Pred_Red, obs=data_val_BR3$dens) %>% mutate(Type = "All data", Bioregion="Eastern bioregion"),
            cal_data(pred = data_val_BR4$Pred_Red, obs=data_val_BR4$dens) %>% mutate(Type = "All data", Bioregion="Western bioregion")
            ) %>% 
  mutate(species="d) Fallow deer")


rm(list = setdiff(ls(), c("cal_data", "cpwBadger.br", "cpwBadger", "cpwFox", "cpwFox.br", "cpwWB.br", "cpwRedDeer", "cpwRedDeer.br",
                          "cpwRoe", "cpwRoe.br", "cpwFd", "cpwFd.br")))
# --------------------


#Charge constrained species  data

# ALPINE IBEX
# # --------------------
# 
# data<-st_read(dsn = "F:/IREC/PhD/D/Docs/Sonia/0_PhD/7_DocTesis/Chapter_ENETWILD/GIS/data_20260207.gpkg", layer = "data_spatial_AlpineIbex")
# load("F:/IREC/PhD/D/Docs/Sonia/0_PhD/7_DocTesis/Chapter_ENETWILD/GIS/ModelResults/ConstrainedSpecies/densityStepModel_Alpineibex_Europe.Rdata") 
# 
# s_model <- sample(nrow(data), nrow(data)*0.8) ;  data_tra <- data[s_model,]; data_val <- data[-s_model,]
# 
# #2.3) Calibration plots  validation data
# data_val<-data_val %>% 
#   mutate(Pred_model = as.integer(predict(fitbe, data_val,type = "response"))) %>%
#   mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))
# 
# data_pos<-data_val %>% filter(dens>0)
# cpcAlpineIbex<-bind_rows(
#   #cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
#   cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
# ) %>% mutate(species="g) Alpine ibex")
# 
# rm(list = setdiff(ls(), c("cpwBadger.br", "cpwBadger", "cpwFox", "cpwFox.br", "cpwWB.br", "cpwRedDeer", "cpwRedDeer.br", "cpwRoe", "cpwRoe.br", "cpwFd", "cpwFd.br",
#                           "cpcAlpineIbex")))
# # --------------------

# AOUDAD
# --------------------
data<-st_read(dsn = "data_20260207.gpkg", layer = "data_spatial_Aoudad")
load("densityStepModel_Aoudad_Europe.Rdata") 
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ;  data_tra <- data[s_model,]; data_val <- data[-s_model,]

#2.3) Calibration plots  validation data
data_val<-data_val %>%
  mutate(Pred_model = as.integer(predict(Aoudad_DensityModel, data_val,type = "response"))) %>%
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

data_pos<-data_val %>% filter(dens>0)
cpcAoudad<-bind_rows(
  # cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
) %>% 
  mutate(species="k) Aoudad")

rm(list = setdiff(ls(), c("cal_data", "cpwBadger.br", "cpwBadger", "cpwFox", "cpwFox.br", "cpwWB.br", "cpwRedDeer", "cpwRedDeer.br", "cpwRoe", "cpwRoe.br", "cpwFd", "cpwFd.br",
                          "cpcAlpineIbex","cpcAoudad")))
# --------------------

# CHAMOIS
# --------------------
data<-st_read(dsn = "data_20260207.gpkg", layer = "data_spatial_NorthChamois")
load("densityStepModel_NorthChamois_Europe.Rdata") 
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]

#2.3) Calibration plots  validation data
data_val<-data_val %>% 
  mutate(Pred_model = as.integer(predict(NorthChamois_DensityModel, data_val,type = "response"))) %>%
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

data_pos<-data_val %>% filter(dens>0)
cpcChamois<-bind_rows(
  # cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
) %>%  mutate(species="i) Northern Chamois")

rm(list = setdiff(ls(), c("cal_data","cpwBadger.br", "cpwBadger", "cpwFox", "cpwFox.br", "cpwWB.br", "cpwRedDeer", "cpwRedDeer.br", "cpwRoe", "cpwRoe.br", "cpwFd", "cpwFd.br",
                          "cpcAlpineIbex","cpcAoudad", "cpcChamois")))
# --------------------

# GOAT
# --------------------
data<-st_read(dsn = "data_20260207.gpkg", layer = "data_spatial_IbWildGoat")
load("densityStepModel_IbWildGoat_Europe.Rdata") 
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]

#2.3) Calibration plots  validation data
data_val<-data_val %>% 
  mutate(Pred_model = as.integer(predict(IbWildGoat_DensityModel, data_val,type = "response"))) %>%
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

data_pos<-data_val %>% filter(dens>0)
cpcGoat<-bind_rows(
  # cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
) %>% mutate(species="g) Iberian wild goat")

rm(list = setdiff(ls(), c("cal_data", "cpwBadger.br", "cpwBadger", "cpwFox", "cpwFox.br", "cpwWB.br", "cpwRedDeer", "cpwRedDeer.br", "cpwRoe", "cpwRoe.br", "cpwFd", "cpwFd.br",
                          "cpcAlpineIbex","cpcAoudad", "cpcChamois", "cpcGoat")))
# --------------------

# ISARD
# --------------------
data<-st_read(dsn = "data_20260207.gpkg", layer = "data_spatial_PyrChamois")
load("densityStepModel_PyrChamois_Europe.Rdata") 
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]

#2.3) Calibration plots  validation data
data_val<-data_val %>% 
  mutate(Pred_model = as.integer(predict(PyrChamois_DensityModel, data_val,type = "response"))) %>%
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

data_pos<-data_val %>% filter(dens>0)
cpcIsard<-bind_rows(
  # cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
) %>% mutate(species="h) Pyrenean chamois")

rm(list = setdiff(ls(), c("cal_data", "cpwBadger.br", "cpwBadger", "cpwFox", "cpwFox.br", "cpwWB.br", "cpwRedDeer", "cpwRedDeer.br", "cpwRoe", "cpwRoe.br", "cpwFd", "cpwFd.br",
                          "cpcAlpineIbex","cpcAoudad", "cpcChamois", "cpcGoat", "cpcIsard")))
# --------------------

# MOOSE
# --------------------
data<-st_read(dsn = "data_20260207.gpkg", layer = "data_spatial_Moose")
load("densityStepModel_Moose_Europe.Rdata") 
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]

#2.3) Calibration plots  validation data
data_val<-data_val %>% 
  mutate(Pred_model = as.integer(predict(Moose_DensityModel, data_val,type = "response"))) %>%
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

data_pos<-data_val %>% filter(dens>0)
cpcMoose<-bind_rows(
  # cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
) %>% mutate(species="j) Moose")

rm(list = setdiff(ls(), c("cal_data", "cpwBadger.br", "cpwBadger", "cpwFox", "cpwFox.br", "cpwWB.br", "cpwRedDeer", "cpwRedDeer.br", "cpwRoe", "cpwRoe.br", "cpwFd", "cpwFd.br",
                          "cpcAlpineIbex","cpcAoudad", "cpcChamois", "cpcGoat", "cpcIsard", "cpcMoose")))
# --------------------

# MOUFLON
# --------------------
data<-st_read(dsn = "data_20260207.gpkg", layer = "data_spatial_Mouflon")
load("densityStepModel_Mouflon_Europe.Rdata") 
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]

#2.3) Calibration plots  validation data
data_val<-data_val %>% 
  mutate(Pred_model = as.integer(predict(Mouflon_DensityModel, data_val,type = "response"))) %>%
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

data_pos<-data_val %>% filter(dens>0)
cpcMouflon<-bind_rows(
  # cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
) %>% mutate(species="l) Mouflon")

rm(list = setdiff(ls(), c("cal_data", "cpwBadger.br", "cpwBadger", "cpwFox", "cpwFox.br", "cpwWB.br", "cpwRedDeer", "cpwRedDeer.br", "cpwRoe", "cpwRoe.br", "cpwFd", "cpwFd.br",
                          "cpcAlpineIbex","cpcAoudad", "cpcChamois", "cpcGoat", "cpcIsard", "cpcMoose", "cpcMouflon")))
# --------------------

# # REINDEER
# # --------------------
# data<-st_read(dsn = "data_20260207.gpkg", layer = "data_spatial_Reindeer")
# load("densityStepModel_Reindeer_poisson_Europe.Rdata") 
# 
# s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]
# 
# #2.3) Calibration plots  validation data
# data_val<-data_val %>% 
#   mutate(Pred_model = as.integer(predict(Reindeer_DensityModel, data_val,type = "response"))) %>%
#   mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))
# 
# data_pos<-data_val %>% filter(dens>0)
# cpcReindeer<-bind_rows(
#   #cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data")#,
#   cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
# ) %>%  mutate(species="n) Reindeer")

# rm(list = setdiff(ls(), c("cal_data", "cpwBadger.br", "cpwBadger", "cpwFox", "cpwFox.br", "cpwWB.br", "cpwRedDeer", "cpwRedDeer.br", "cpwRoe", "cpwRoe.br", "cpwFd", "cpwFd.br",
#                           "cpcAlpineIbex","cpcAoudad", "cpcChamois", "cpcGoat", "cpcIsard", "cpcMoose", "cpcMouflon", 
#                           "cpcReindeer")))
# # --------------------

# SIKA DEER
# --------------------
data<-st_read(dsn = "data_20260207.gpkg", layer = "data_spatial_SikaDeer")
load("densityStepModel_SikaDeer_Europe.Rdata") 
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]

#2.3) Calibration plots  validation data
data_val<-data_val %>% 
  mutate(Pred_model = as.integer(predict(SikaDeer_DensityModel , data_val,type = "response"))) %>%
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

data_pos<-data_val %>% filter(dens>0)
cpcSika<-bind_rows(
  # cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
) %>%  mutate(species="m) Sika deer")

rm(list = setdiff(ls(), c("cpwBadger.br", "cpwBadger", "cpwFox", "cpwFox.br", "cpwWB.br", "cpwRedDeer", "cpwRedDeer.br", "cpwRoe", "cpwRoe.br", "cpwFd", "cpwFd.br",
                          "cpcAlpineIbex","cpcAoudad", "cpcChamois", "cpcGoat", "cpcIsard", "cpcMoose", "cpcMouflon", 
                          "cpcReindeer", "cpcSika")))
# --------------------


library(ggpubr)
library(patchwork)
library(ggh4x)
# -----------------------------
# a<-
# ggplot(cpwWB.br, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=2)+ #scale_colour_grey()+ 
#   scale_color_manual(values=c("#7C4D87", "#FAEE62"))+
#   scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="#62B9AD", #"#AAE280", "#62B9AD", "#FAEE62"
#               lty=1, cex=0.8, 
#               alpha=0.8)+
#   theme_bw()+
#   ylab("Observed HDs")+xlab("Predicted HDs")+
#   theme(plot.title = element_text(hjust = 0.5),
#         strip.background = element_rect(fill = scales::alpha("#62B9AD"), color = "black"),
#         strip.text = element_text(color = "white", face = "bold"),
#         
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(color = "grey90"), 
#         legend.position = "none",
#         legend.title = element_text(face = "bold"),
#         legend.key = element_blank()) +
#   facet_wrap(~Bioregion, ncol=4, nrow=1, scales="free")
# a
# 
# ggsave(filename = "F:/IREC/PhD/D/Docs/Sonia/0_PhD/7_DocTesis/Chapter_ENETWILD/CalibrationPlot_WB.png", # filename of the plot
#        plot = a,           # plot to save
#        width = 18, height = 5,   # dimensions
#        units = "cm",              # dimension units
#        dpi = 300)
# 
# b2<-
#   ggplot(cpwRoe.br, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=2)+ #scale_colour_grey()+ 
#   scale_color_manual(values=c("#7C4D87", "#FAEE62"))+
#   scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="#62B9AD", #"#AAE280", "#62B9AD", "#FAEE62"
#               lty=1, cex=0.8, 
#               alpha=0.8)+
#   theme_bw()+
#   ggtitle("b) Roe deer")+
#   ylab("Observed HDs")+xlab("Predicted HDs")+
#   theme(plot.title = element_text(hjust = 0.5),
#         strip.background = element_rect(fill = scales::alpha("#62B9AD", 1), color = "black"),
#         strip.text = element_text(color = "white", face = "bold"),
#         axis.title.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(color = "grey90"), 
#         legend.position = "bottom",
#         legend.title = element_text(face = "bold"),
#         legend.key = element_blank(),
#         title = element_text(face = "bold"))+
#   facet_wrap(~Bioregion, ncol=4, nrow=1, scales="free")
# 
# c2<-
#   ggplot(cpwRedDeer.br, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=2)+ #scale_colour_grey()+ 
#   scale_color_manual(values=c("#7C4D87", "#FAEE62"))+
#   scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="#62B9AD", #"#AAE280", "#62B9AD", "#FAEE62"
#               lty=1, cex=0.8, 
#               alpha=0.8)+
#   theme_bw()+
#   ggtitle("c) Red deer")+
#   ylab("Observed HDs")+xlab("Predicted HDs")+
#   theme(plot.title = element_text(hjust = 0.5),
#         strip.background = element_rect(fill = scales::alpha("#62B9AD", 1), color = "black"),
#         strip.text = element_text(color = "white", face = "bold"),
#         axis.title.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(color = "grey90"), 
#         legend.position = "bottom",
#         legend.title = element_text(face = "bold"),
#         legend.key = element_blank(),
#         title = element_text(face = "bold"))+
#   facet_wrap(~Bioregion, ncol=4, nrow=1, scales="free")
# 
# d2<-  
#   ggplot(
#     cpwFd.br
#     , aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=2)+ #scale_colour_grey()+ 
#   scale_color_manual(values=c("#7C4D87", "#FAEE62"))+
#   scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="#62B9AD", #"#AAE280", "#62B9AD", "#FAEE62"
#               lty=1, cex=0.8, 
#               alpha=0.8)+
#   theme_bw()+
#   ggtitle("d) Fallow deer")+
#   ylab("Observed HDs")+xlab("Predicted HDs")+
#   theme(plot.title = element_text(hjust = 0.5),
#         strip.background = element_rect(fill = scales::alpha("#62B9AD", 1), color = "black"),
#         strip.text = element_text(color = "white", face = "bold"),
#         axis.title.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(color = "grey90"), 
#         legend.position = "bottom",
#         legend.title = element_text(face = "bold"),
#         legend.key = element_blank(),
#         title = element_text(face = "bold")
#   ) +
#   facet_wrap(~Bioregion, ncol=4, nrow=1, scales="free")
# 
# 
# 
# wrap_plots(a, b2, c2, d2, ncol = 1) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")




# 
# b1<-
#   ggplot(cpwRoe, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#     geom_point(size=2)+ #scale_colour_grey()+ 
#     scale_color_manual(values=c("#7C4D87", "#FAEE62"))+
#     scale_shape_manual(values=c(16, 17))+
#     geom_abline(col="#62B9AD", #"#AAE280", "#62B9AD", "#FAEE62"
#                 lty=1, cex=0.8, 
#                 alpha=0.8)+
#     theme_bw()+
#     ylab("Observed HDs")+xlab("Predicted HDs")+
#     theme(plot.title = element_text(hjust = 0.5),
#           strip.background = element_rect(fill = scales::alpha("#62B9AD"), color = "black"),
#           strip.text = element_text(color = "white", face = "bold"),
#           
#           panel.grid.minor = element_blank(),
#           panel.grid.major = element_line(color = "grey90"), 
#           legend.position = "bottom",
#           legend.title = element_text(face = "bold"),
#           legend.key = element_blank())+
#     facet_wrap(~species, ncol=1, nrow=1, scales="free")
# 
# 
# c1<-
#   ggplot(cpwRedDeer, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#     geom_point(size=2)+ #scale_colour_grey()+ 
#     scale_color_manual(values=c("#7C4D87", "#FAEE62"))+
#     scale_shape_manual(values=c(16, 17))+
#     geom_abline(col="#62B9AD", #"#AAE280", "#62B9AD", "#FAEE62"
#                 lty=1, cex=0.8, 
#                 alpha=0.8)+
#     theme_bw()+
#     ylab("Observed HDs")+xlab("Predicted HDs")+
#     theme(plot.title = element_text(hjust = 0.5),
#           strip.background = element_rect(fill = scales::alpha("#62B9AD"), color = "black"),
#           strip.text = element_text(color = "white", face = "bold"),
#           
#           panel.grid.minor = element_blank(),
#           panel.grid.major = element_line(color = "grey90"), 
#           legend.position = "bottom",
#           legend.title = element_text(face = "bold"),
#           legend.key = element_blank())+
#     facet_wrap(~species, ncol=1, nrow=1, scales="free")
# 
# 
# d1<-
#   ggplot(cpwFd, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#     geom_point(size=2)+ #scale_colour_grey()+ 
#     scale_color_manual(values=c("#7C4D87", "#FAEE62"))+
#     scale_shape_manual(values=c(16, 17))+
#     geom_abline(col="#62B9AD", #"#AAE280", "#62B9AD", "#FAEE62"
#                 lty=1, cex=0.8, 
#                 alpha=0.8)+
#     theme_bw()+
#     ylab("Observed HDs")+xlab("Predicted HDs")+
#     theme(plot.title = element_text(hjust = 0.5),
#           strip.background = element_rect(fill = scales::alpha("#62B9AD"), color = "black"),
#           strip.text = element_text(color = "white", face = "bold"),
#           
#           panel.grid.minor = element_blank(),
#           panel.grid.major = element_line(color = "grey90"), 
#           legend.position = "bottom",
#           legend.title = element_text(face = "bold"),
#           legend.key = element_blank())+
#     facet_wrap(~species, ncol=1, nrow=1, scales="free")
# 
# p1<-wrap_plots(b1, c1, d1, ncol = 3) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")
# p1
# ggsave(filename = "F:/IREC/PhD/D/Docs/Sonia/0_PhD/7_DocTesis/Chapter_ENETWILD/CalibrationPlot_Wid1_v2.png", # filename of the plot
#        plot = p1,           # plot to save
#        width = 20, height = 24,   # dimensions
#        units = "cm",              # dimension units
#        dpi = 300)

# 
# 
# ggarrange(
#   ggplot(
#     cpwFox
#     , aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#     geom_point(size=2)+ #scale_colour_grey()+ 
#     scale_color_manual(values=c("#7C4D87", "#FAEE62"))+
#     scale_shape_manual(values=c(16, 17))+
#     geom_abline(col="#62B9AD", #"#AAE280", "#62B9AD", "#FAEE62"
#                 lty=1, cex=0.8, 
#                 alpha=0.8)+
#     theme_bw()+
#     ylab("Observed HDs")+xlab("Predicted HDs")+
#     theme(plot.title = element_text(hjust = 0.5),
#           strip.background = element_rect(fill = scales::alpha("#62B9AD"), color = "black"),
#           strip.text = element_text(color = "white", face = "bold"),
#           
#           panel.grid.minor = element_blank(),
#           panel.grid.major = element_line(color = "grey90"), 
#           legend.position = "bottom",
#           legend.title = element_text(face = "bold"),
#           legend.key = element_blank()
#     )+
#     facet_wrap(~species, ncol=1, nrow=1, scales="free")
#   , 
#   
#   ggplot(
#     cpwFox.br
#     , aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#     geom_point(size=2)+ #scale_colour_grey()+ 
#     scale_color_manual(values=c("#7C4D87", "#FAEE62"))+
#     scale_shape_manual(values=c(16, 17))+
#     geom_abline(col="#62B9AD", #"#AAE280", "#62B9AD", "#FAEE62"
#                 lty=1, cex=0.8, 
#                 alpha=0.8)+
#     theme_bw()+
#     ylab("Observed HDs")+xlab("Predicted HDs")+
#     theme(plot.title = element_text(hjust = 0.5),
#           strip.background = element_rect(fill = scales::alpha("#62B9AD"), color = "black"),
#           strip.text = element_text(color = "white", face = "bold"),
#           
#           panel.grid.minor = element_blank(),
#           panel.grid.major = element_line(color = "grey90"), 
#           legend.position = "bottom",
#           legend.title = element_text(face = "bold"),
#           legend.key = element_blank()) +
#     facet_wrap(~Bioregion, ncol=2, nrow=2, scales="free"), 
#   common.legend = TRUE, legend = "bottom"
# )
# 
# 
# ggarrange(
# ggplot(
#     cpwBadger
#          , aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=2)+ #scale_colour_grey()+ 
#   scale_color_manual(values=c("#7C4D87", "#FAEE62"))+
#   scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="#62B9AD", #"#AAE280", "#62B9AD", "#FAEE62"
#               lty=1, cex=0.8, 
#               alpha=0.8)+
#   theme_bw()+
#   ylab("Observed HDs")+xlab("Predicted HDs")+
#   theme(plot.title = element_text(hjust = 0.5),
#         strip.background = element_rect(fill = scales::alpha("#62B9AD"), color = "black"),
#         strip.text = element_text(color = "white", face = "bold"),
#         
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(color = "grey90"), 
#         legend.position = "bottom",
#         legend.title = element_text(face = "bold"),
#         legend.key = element_blank()
#   )+
#   facet_wrap(~species, ncol=1, nrow=1, scales="free")
# , 
# 
# ggplot(
#   cpwBadger.br
#   , aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=2)+ #scale_colour_grey()+ 
#   scale_color_manual(values=c("#7C4D87", "#FAEE62"))+
#   scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="#62B9AD", #"#AAE280", "#62B9AD", "#FAEE62"
#               lty=1, cex=0.8, 
#               alpha=0.8)+
#   theme_bw()+
#   ylab("Observed HDs")+xlab("Predicted HDs")+
#   theme(plot.title = element_text(hjust = 0.5),
#         strip.background = element_rect(fill = scales::alpha("#62B9AD"), color = "black"),
#         strip.text = element_text(color = "white", face = "bold"),
#         
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(color = "grey90"), 
#         legend.position = "bottom",
#         legend.title = element_text(face = "bold"),
#         legend.key = element_blank()) +
#   facet_wrap(~Bioregion, ncol=2, nrow=2, scales="free"), 
# common.legend = TRUE, legend = "bottom"
# )

# 
# ggsave(filename = "F:/IREC/PhD/D/Docs/Sonia/0_PhD/7_DocTesis/Chapter_ENETWILD/CalibrationPlot_Widespread.png", # filename of the plot
#        plot = cw,           # plot to save
#        width = 5, height = 18,   # dimensions
#        units = "cm",              # dimension units
#        dpi = 300)
# ------------------------

# -----------------------
# cpw.br<-rbind(cpwWB.br, 
#               cpwRoe.br, 
#               cpwRedDeer.br,
#               cpwFd.br, 
#               cpwFox.br, 
#               cpwBadger.br)
# 
# all<-ggplot(cpw.br, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=2)+ #scale_colour_grey()+ 
#   scale_color_manual(values=c("#7C4D87", "#FAEE62"))+
#   scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="#62B9AD", #"#AAE280", "#62B9AD", "#FAEE62"
#               lty=1, cex=0.8, 
#               alpha=0.8)+
#   theme_bw()+
#   ylab("Observed HDs")+xlab("Predicted HDs")+
#   theme(plot.title = element_text(hjust = 0.5),
#         strip.background = element_rect(fill = scales::alpha("#62B9AD", 1), color = "black"),
#         strip.text = element_text(color = "white", face = "bold"),
#         
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(color = "grey90"), 
#         legend.position = "bottom",
#         legend.title = element_blank(),
#         legend.key = element_blank(), 
#         title = element_text(face = "bold"))+
#   facet_wrap2(species~Bioregion, 
#              scales="free", 
#              ncol=4, 
#              strip.position = "left", 
#              strip = strip_nested())
# all

# ggsave(filename = "F:/IREC/PhD/D/Docs/Sonia/0_PhD/7_DocTesis/Chapter_ENETWILD/CalibrationPlot_WidespreadBRAll_v2.png", # filename of the plot
#        plot = all,           # plot to save
#        width = 18, height = 24,   # dimensions
#        units = "cm",              # dimension units
#        dpi = 300)
# -----------------------

plot.wbr <- function(data, label = "title", scales = "free") {
  ggplot(data, aes(x = pred_mean, y = obs_mean, col = Type, shape = Type)) +
    geom_point(size = 2) +
    scale_color_manual(values = c("#7C4D87", "#FAEE62")) +
    scale_shape_manual(values = c(16, 17)) +
    geom_abline(
      col = "#62B9AD",
      lty = 1,
      alpha = 0.8
    ) +
    theme_bw() +
    ylab("Observed HDs") +
    xlab("Predicted HDs") +
    
    facet_wrap(~ Bioregion, scales = scales, ncol = 4) +
    
    ## ---- TEXTO ÚNICO A LA DERECHA ----
  labs(tag = label) +
    
    theme(
      strip.background = element_rect(
        fill = "#62B9AD",
        color = "black"
      ),
      strip.text = element_text(
        color = "white",
        face = "bold"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size=13),
      
      axis.title = element_text(size=14),
      strip.text.y = element_text(size=14),
      strip.text.x = element_text(size=14),
      
      ## ---- CONTROL DEL TEXTO ----
      plot.tag.position = c(1.03, 0.5),   # derecha, centrado vertical
      plot.tag = element_text(
        angle = -90,
        face = "bold",
        size = 17,
        color = "#62B9AD"
      ),
      
      ## margen derecho para que no se corte
      plot.margin = margin(0.2, 35, 0.2, 0.2)
    )
}

cpwBadger.br$Bioregion <- factor(
  cpwBadger.br$Bioregion,
  levels = c("Eastern bioregion", "Northern bioregion", "Southern bioregion", "Western bioregion")
)

cpwBadger.br2 <- cpwBadger.br %>% complete(species, Bioregion)
cpwBadger.br2$fake_col <- cpwBadger.br2$Bioregion
# plot.wbr(cpwBadger.br2)

all<-ggarrange(
  plot.wbr(cpwWB.br, label = "a) Wild boar", scales = "fixed") + xlab(NULL),
  plot.wbr(cpwRoe.br, label = "b) Roe deer", scales = "fixed") + xlab(NULL),
  plot.wbr(cpwRedDeer.br, label = "c) Red deer", scales = "fixed")+ xlab(NULL),
  plot.wbr(cpwFd.br, label = "d) Fallow deer", scales = "fixed")  + xlab(NULL),
  plot.wbr(cpwFox.br, label = "e) Red fox", scales = "fixed") + xlab(NULL),
  plot.wbr(cpwBadger.br2, label="f) European badger", scales = "fixed"),
  nrow=6, common.legend = TRUE, legend = "none", 
  align = "hv"
)
# all
ggsave(filename = "F:/IREC/PhD/D/Docs/Sonia/0_PhD/7_DocTesis/Chapter_ENETWILD/CalibrationPlot_WidespreadBRAll_v6.3.png", # filename of the plot
       plot = all,           # plot to save
       width = 28, height = 40,   # dimensions
       units = "cm",              # dimension units
       dpi = 300)


objs <- ls()
a<-objs[grepl("cpw", objs) & !grepl("br", objs)]
list.w<-list()
for ( i in 1:length(a)){
  df.w<-get(a[i])
  df.w<-df.w %>% mutate(name=a[i])
  list.w[[i]]<-df.w
}
list.w<-do.call(rbind, list.w)

# CONSTRAINED SPECIES
list.c<-list()
for ( i in 1:length(ls(pattern="cpc"))){
  df.w<-get(ls(pattern="cpc")[i])
  df.w<-df.w %>% mutate(name=ls(pattern="cpc")[i])
  list.c[[i]]<-df.w
}
list.c<-do.call(rbind, list.c)

listE<-rbind(list.w, list.c) %>% 
          mutate(species = ifelse(species=="j) Goat", "g) Iberian wild goat", 
                                         ifelse(species=="k) Isard", "h) Pyrenean chamois", 
                                                ifelse(species=="l) Moose", "j) Moose", species))))

# -------------
# cp<-
# ggplot(list.c, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=2)+ #scale_colour_grey()+ 
#   scale_color_manual(values=c("#7C4D87", "#FAEE62"))+
#   scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="#62B9AD", #"#AAE280", "#62B9AD", "#FAEE62"
#               lty=1, cex=0.8, 
#               alpha=0.8)+
#   theme_bw()+
#   ylab("Observed HDs")+xlab("Predicted HDs")+
#   theme(plot.title = element_text(hjust = 0.5),
#     strip.background = element_rect(fill = scales::alpha("#62B9AD"), color = "black"),
#     strip.text = element_text(color = "white", face = "bold"),
#     
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(color = "grey90"), 
#     legend.position = "bottom",
#     legend.title = element_text(face = "bold"),
#     legend.key = element_blank()
#   )+
#   facet_wrap(~species, ncol=3, nrow=3, scales="free")
# cp
# ggsave(filename = "F:/IREC/PhD/D/Docs/Sonia/0_PhD/7_DocTesis/Chapter_ENETWILD/CalibrationPlot_Constrained.png", # filename of the plot
#        plot = cp,           # plot to save
#        width = 18, height = 18,   # dimensions
#        units = "cm",              # dimension units
#        dpi = 300)
# -------------

allE<-ggplot(listE, aes(x = pred_mean, y = obs_mean, col = Type, shape = Type))+
  geom_point(size = 2)+ #scale_colour_grey()+ 
  scale_color_manual(values = c("#7C4D87", "#FAEE62"))+
  scale_shape_manual(values = c(16, 17))+
  geom_abline(col = "#62B9AD", #"#AAE280", "#62B9AD", "#FAEE62"
              lty = 1, 
              cex = 0.8, 
              alpha = 0.8)+
  theme_bw()+
  ylab("Observed HDs") + xlab("Predicted HDs") +
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill = scales::alpha("#62B9AD", 1), color = "black"),
        strip.text = element_text(color = "white", face = "bold"),
        axis.title = element_text(size=12, face ="plain"),
        strip.text.y = element_text(size=14),
        strip.text.x = element_text(size=14),
        
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90"), 
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_blank(), 
        legend.text = element_text(size=13),
        title = element_text(face = "bold"))+
  facet_wrap( ~ species, scales = "fixed", ncol=3)
allE
# 
ggsave(filename = "CalibrationPlot_AllE_3.6.png", # filename of the plot
       plot = allE,           # plot to save
       width = 20, height = 27,   # dimensions
       units = "cm",              # dimension units
       dpi = 300)
