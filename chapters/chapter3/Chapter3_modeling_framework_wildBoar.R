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

library(tidyverse)
library(sf)
#rm(list=ls())

# 0. Load wild boar data (dont be afraid for the name of objects (Red))
# ----------------------------------------------------------------------
# try not to use this same name for other objects, as if something fails you may need 
# to reload this data and it takes a while.... 
# If so, and something is wrong you may reuse this object form line 11, and save time! :)
spatial<-st_read(dsn="WB_WLDM_merged_v1.5.gpkg", layer='WBEuropeanStandarWLDM_AllMerged_LongFormat_variables')

# IMPORTANT WARNING1: 
# For removing Bayern federal estate it's needed that all rows have a locality not NA value. 
# Therefore, those ones that are na we'll record their locationID+country there
spatialID<-spatial %>% 
                mutate(locality=ifelse(is.na(locality), paste(locationID, country), locality)) %>% 
                dplyr::filter(locality!='Bayern')

# IMPORTANT WARNING2: 
# be careful with the locationID name. Some countries have only numeric codes 
# could be repeated among them and that would be a problem. To solve that, 
# I added the country code before the locatinID. Here there are some countries 
# which locationID was numeric, check if they are still numeric or have changed,
# or if other countries have change it to numeric!
spatialID<-spatialID %>% mutate(locationID=
                                        ifelse(country=="Norway", paste("NO_", locality, sep="_"), 
                                               ifelse(country=="Austria", paste("AT_", locality, sep="_"), 
                                                      ifelse(country=="Belarus", paste("BLR_", locality, sep="_"), 
                                                             ifelse(country=="Finland", paste("FI_", locality, sep="_"), 
                                                                    ifelse(country=="Lithuania", paste("LT_", locality, sep="_"), 
                                                                           ifelse(country=="Latvia", paste("LV", locality, sep="_"), 
                                                                                  ifelse(country=="Belgium", paste("BE", locality, sep="_"), locationID))))))#)
                                        ), 
                                      country=ifelse(country=="ES" | country=="Spain", "Spain", country)) %>% 
                    # CHECK IF YOU NEED ANY COLUMN MORE BUT AT FIRST GLANCE NO NEEDED
                    dplyr::select(locationID, country, locality, dataTime, harvTot, area_km2, NUT, starts_with("lc"), ends_with("_mean"), BR_todasr1)

spatialID1<-spatialID %>% 
  dplyr::select(BR_todasr1, area_km2, country, dataTime, harvTot, locationID) %>% 
  mutate(dataTime=as.numeric(substr(dataTime, 1, 4))) %>%  
  mutate(harvTot=ifelse(harvTot>0, harvTot, 0)) # where we have NA in harvest complete with a 0  

# transform from long towide
spatialID1<-spatialID1 %>% 
  # important: transform to table and remove geometry column, it can be geom or geometry depending if you load the data from a shapefile or geopackage
  as_tibble() %>% dplyr::select(-geom) %>% 
  filter(dataTime>2014) %>%  # filter years that we're going to use for transforming data (long-wide)
  # retain Sweden finest data 
  filter(country!="Sweden" | dataTime!=2012) %>%  
  filter(country!="Sweden" | dataTime!=2013) %>%
  filter(country!="Sweden" | dataTime!=2014) %>%
  filter(country!="Sweden" | dataTime!=2015) %>%
  filter(country!="Sweden" | dataTime!=2016) %>%
  filter(country!="Sweden" | dataTime!=2017) %>%
  filter(country!="Sweden" | dataTime!=2018) %>%
  mutate(dataTime = str_c("harvTot", dataTime)) %>% 
  pivot_wider(names_from = dataTime, values_from = c("harvTot"), names_sep="", values_fn=sum)  # transform from long to wide

unique<-spatialID %>% mutate(locationID_country=paste(locationID, country, sep="_"))

# remove duplicated data by locationID_country  
# --> expecting to have just one administrative unit per country and remove duplicated polygons
unique<-unique[!duplicated(unique$locationID_country), ]

# Join the non duplicated polygon layer with wide table by locationID and country
# remove columns from each part that may cause confusion
spatial_<- unique %>% dplyr::select(-harvTot, -dataTime) %>%  
  inner_join(spatial1 %>% dplyr::select(-c(BR_todasr1, area_km2)), 
             by=c("locationID", "country")) %>%
  rename(alt=alt_mean, Eucmean=EucalyptusSpp_mean, hfp=hfp_mean, snow=snow_mean,  sun=sun_mean) %>% 
  mutate(NUTS=as.factor(ifelse(NUT=="comuni" | NUT=="Concelho" | NUT=="District" | NUT=="Kommuner" | NUT=="Municipality", 1, ifelse(NUT=="NUTS0" |NUT=="NUTS1" |NUT=="NUTS2" |NUT=="NUTS3", 2, 0))),
         locationID=as.factor(locationID),
         locality=as.factor(locality), 
         country=as.factor(country), 
         Eucmean = replace_na(Eucmean, 0)) 

# calculate area (km2)
# this step can take a few
spatial_<-spatial_ %>% mutate(area_km2=st_area(spatial_)/1000000) ; attributes(spatial_$area_km2) = NULL

spatial_ <-spatial_ %>% 
  # calculate harvest Maximum between the years that we have data. 
  # Be aware if new data is introduced you may add
  # some more for example Red_spatial_$harvTot2022...etc.
  mutate(harvMax=pmax(spatial_$harvTot2015, spatial_$harvTot2016, spatial_$harvTot2017, spatial_$harvTot2018, spatial_$harvTot2019, spatial_$harvTot2020, spatial_$harvTot2021, na.rm = TRUE)) %>% 
  mutate(harvMax=ifelse(is.na(harvMax), 0, harvMax)) %>% 
  # calculate density
  mutate(dens=round(harvMax/area_km2, 2)) %>% 
  mutate(dens_r_trans=dens*10000) %>% 
  dplyr::select(locationID, locality, country, NUTS, #species,
                matches("excel"), matches("shp"),
                matches("bio"), grow_mean, alt, Eucmean, hfp, snow, sun, area_km2, BR_todasr1, #Bioregion,
                matches("lc"), 
                # if you have added new columns for calculate harvMax remember 
                # to select them here following the example above harvTot2022...
                harvTot2015, harvTot2016, harvTot2017, harvTot2018, harvTot2019, harvTot2020, harvTot2021, harvMax, dens, dens_r_trans, suit_mean)

# calculate coordinates
# this step can take a few too
spatial_<- spatial_ %>% mutate(x_cen=st_coordinates(st_centroid(spatial_))[,1], y_cen=st_coordinates(st_centroid(spatial_))[,2])

# it's nice if you save it in a geopackage, and have a look that your data 
# it's okay or anything is missing, etc.

df_<-spatial_ %>%
  mutate(lc10=lc_10_sum/lc_10_count*100, lc11=lc_11_sum/lc_11_count*100, lc12=lc_12_sum/lc_12_count*100, lc20=lc_20_sum/lc_20_count*100, lc30=lc_30_sum/lc_30_count*100, lc40=lc_40_sum/lc_40_count*100, lc50=lc_50_sum/lc_50_count*100, lc60=lc_60_sum/lc_60_count*100, lc61=lc_61_sum/lc_61_count*100, lc62=lc_62_sum/lc_62_count*100, lc70=lc_70_sum/lc_70_count*100, lc71=lc_71_sum/lc_71_count*100, lc72=lc_72_sum/lc_72_count*100, lc80=lc_80_sum/lc_80_count*100, lc81=lc_81_sum/lc_81_count*100, lc82=lc_82_sum/lc_82_count*100, lc90=lc_90_sum/lc_90_count*100, lc100=lc_100_sum/lc_100_count*100, lc110=lc_110_sum/lc_110_count*100,lc120=lc_120_sum/lc_120_count*100,lc121=lc_121_sum/lc_121_count*100,lc122=lc_122_sum/lc_122_count*100,lc130=lc_130_sum/lc_130_count*100,lc140=lc_140_sum/lc_140_count*100,lc150=lc_150_sum/lc_150_count*100,lc152=lc_152_sum/lc_152_count*100,lc153=lc_153_sum/lc_153_count*100,lc160=lc_160_sum/lc_160_count*100,lc170=lc_170_sum/lc_170_count*100,lc180=lc_180_sum/lc_180_count*100,lc190=lc_190_sum/lc_190_count*100,lc200=lc_200_sum/lc_200_count*100,lc201=lc_201_sum/lc_201_count*100,lc202=lc_202_sum/lc_202_count*100,lc210=lc_210_sum/lc_210_count*100,lc220=lc_220_sum/lc_220_count*100) %>% 
  mutate(Eucmean=ifelse(Eucmean<0, 0.01, Eucmean)) %>% #Como para wild boar transformamos la variable de Eucmean para eliminar los valores negativos por valores de 0 muy peque?os
  dplyr::select(!matches("_cou") & !matches("_sum")) %>% 
  # quito variables que no vaya a usar en el modelo 
  dplyr::select(-c( #excel_file,excel_path, shp_file, shp_path, #species, #me cargo variables que no necesito para el modelo
    grow_mean, #la variable grow tiene datos raros, asi que la voy a eliminar de nuestro set de datos
    lc50, lc62, lc72, lc81, lc82,lc121, lc170 #, #Me cargo variables lc que en el summary son todo 0
  )) %>% 
  filter(dens<50) #tres poligonos de >Spain con areas muy muy muy muy small

df_<-df_[!is.na(df_$dens),] #Remove if there are NA values from RV
colnames(df_)<-gsub("_mean", "", colnames(Red_))


#2) Modeling phase
# --------------------------------
#2.1) scale ALL of your covariates: check after that the mean summary value will be 0 for all of them
# select the names of the covariates that you want to scale
# check it its place 23:77 ALL covariates: bio, lc, x_cen, y_cen,  etc. 
# dont overwrite Red_ object. You will need mean and sd of the covariates for 
# standardizing the 10kmgrid
data<-df_ %>% dplyr::select(locationID, locality, country, NUTS, BR_todasr1, area_km2, matches("harv"), 
                             dens, dens_r_trans,
                             matches("bio"), alt, snow, sun, Eucmean, hfp, 
                             matches("lc"), 
                             x_cen, y_cen) %>% rename(Bioregion=BR_todasr1)

#remove NA values from the dataset
summary(data)
data<-data[!is.na(data$bio_01),] #Remove NA values from the dataset
data<-data[!is.na(data$snow),] #Remove NA values from the dataset
data<-data[!is.na(data$hfp),] #Remove NA values from the dataset
data<-data[!is.na(data$alt),] #Remove NA values from the dataset
data<-data[!is.na(data$lc10),] #Remove NA values from the dataset
# check again if there is still some row with NA values
summary(data) 

vars_std <- c(colnames(data %>% as_tibble())[c(17:71)])
data<-data  %>% mutate(across(all_of(vars_std),~ as.numeric(scale(.))))
summary(data) #mean and sd must to be 0


# 2.0.1 split the dataset by each bioregion
# -------------------------------------------
df_br1<-data %>% filter(Bioregion==1); st_geometry(df_br1)<-NULL
df_br2<-data %>% filter(Bioregion==2); st_geometry(df_br2)<-NULL
df_br3<-data %>% filter(Bioregion==3); st_geometry(df_br3)<-NULL
df_br4<-data %>% filter(Bioregion==4); st_geometry(df_br4)<-NULL

# 2.0.2 Creamos data set de calibrado y de validacion (80% / 20%) para cada Bioregion
# ----------------------------------------------------------------------------------
set.seed(500)
s_model_BR1 <- sample(nrow(df_br1), nrow(df_br1)*0.8) ; data_tra_BR1 <- df_br1[s_model_BR1,]; data_val_BR1 <- df_br1[-s_model_BR1,]
s_model_BR2 <- sample(nrow(df_br2), nrow(df_br2)*0.8) ; data_tra_BR2 <- df_br2[s_model_BR2,]; data_val_BR2 <- df_br2[-s_model_BR2,]
s_model_BR3 <- sample(nrow(df_br3), nrow(df_br3)*0.8) ; data_tra_BR3 <- df_br3[s_model_BR3,]; data_val_BR3 <- df_br3[-s_model_BR3,]
s_model_BR4 <- sample(nrow(df_br4), nrow(df_br4)*0.8) ; data_tra_BR4 <- df_br4[s_model_BR4,]; data_val_BR4 <- df_br4[-s_model_BR4,]


#2.1 VIF <2
# seleccionamos variables para el vif, para eso miramos como es la distribucion de las variables
# ----------------------------------------------------------------------------------
hist(data_tra_BR3$Eucmean)
hist(data_tra_BR3$alt) 
hist(data_tra_BR3$snow)

hist(data_tra_BR3$lc10)
hist(data_tra_BR3$lc11)
hist(data_tra_BR3$lc12) #
hist(data_tra_BR3$lc20) #
hist(data_tra_BR3$lc30)
hist(data_tra_BR3$lc40) 
hist(data_tra_BR3$lc60) 
hist(data_tra_BR3$lc70)
hist(data_tra_BR3$lc90) 
hist(data_tra_BR3$lc100)
hist(data_tra_BR3$lc120) 
hist(data_tra_BR3$lc130)

hist(data_tra_BR3$lc61) #
hist(data_tra_BR3$lc71) #
hist(data_tra_BR3$lc80) #
hist(data_tra_BR3$lc110) #
hist(data_tra_BR3$lc122) #
hist(data_tra_BR3$lc140) #
hist(data_tra_BR3$lc150) #
hist(data_tra_BR3$lc152) #
hist(data_tra_BR3$lc153) #
hist(data_tra_BR3$lc160) #
hist(data_tra_BR3$lc180) #
hist(data_tra_BR3$lc190) #
hist(data_tra_BR3$lc200) #
hist(data_tra_BR3$lc201) #
hist(data_tra_BR3$lc202) #
hist(data_tra_BR3$lc210) #
hist(data_tra_BR3$lc220) #


#2.1 VIF <2
variables.br1<-df_br1 %>% dplyr::select(matches("bio_"), alt, snow, #Eucmean,
                                        lc10, lc11, lc30, lc40, lc60, lc70, lc90, #lc100, #lc120, #lc130, 
                                        lc180, #lc190, 
                                        lc210)#, x_cen, y_cen)
variables.br2<-df_br2 %>% dplyr::select(matches("bio_"), alt, snow, Eucmean,lc10, lc11, lc12, lc30, lc40, lc60, lc70, #lc90, 
                                        lc100, lc120, lc130)#, x_cen, y_cen)
variables.br3<-df_br3 %>% dplyr::select(matches("bio_"), alt, snow, #Eucmean,
                                        lc10, lc11, lc30, lc40, lc60, lc70, lc90, 
                                        lc100, #lc120, 
                                        lc130)#, x_cen, y_cen)
variables.br4<-df_br4 %>% dplyr::select(matches("bio_"), alt, snow, Eucmean,lc10, lc11, lc30, lc40, lc60, lc70, lc90, 
                                        lc100, lc130, lc190#, x_cen, y_cen
                                        )

#2.1 hacemos el vif
resultado.vif1<-usdm::vifstep(variables.br1, th=2)
resultado.vif2<-usdm::vifstep(variables.br2, th=2)
resultado.vif3<-usdm::vifstep(variables.br3, th=2)
resultado.vif4<-usdm::vifstep(variables.br4, th=2)


#2.1) Density Models: binomial negativa necesita que sea un numero entero, por eso se ha transformado la variable densidad *10000
# ----------------------------------------------------------------------------------
library(MASS)
#2.1) Density Models
empty.model1<-glm.nb(dens_r_trans ~ 1, data=data_tra_BR1 %>% filter(dens_r_trans>0), link=log)
empty.model2<-glm.nb(dens_r_trans ~ 1, data=data_tra_BR2 %>% filter(dens_r_trans>0), link=log)
empty.model3<-glm.nb(dens_r_trans ~ 1, data=data_tra_BR3 %>% filter(dens_r_trans>0), link=log)
empty.model4<-glm.nb(dens_r_trans ~ 1, data=data_tra_BR4 %>% filter(dens_r_trans>0), link=log)

resultado.vif1
resultado.vif2
resultado.vif3
resultado.vif4

variables_BR11<-formula(dens_r_trans ~ 
                          bio_02 + bio_15 + bio_18 + 
                          snow + 
                          # lc10 + 
                          lc11 + lc30 + lc40 + lc60 + lc90 + lc180 + lc210)
variables_BR22<-formula(dens_r_trans ~ 
                          bio_03 + bio_04 + bio_08 + bio_09 + 
                          # alt + 
                          snow + Eucmean+
                          lc10 + lc12 + lc30 + lc40 + lc60 + lc70 + lc100 + lc120 + lc130)
variables_BR33<-formula(dens_r_trans ~
                          bio_08 + bio_09 + bio_13 + bio_15 +
                          # alt + 
                          lc10 + lc30 + lc40 + lc60 + lc70 + lc90 + lc100 + lc130)
variables_BR44<-formula(dens_r_trans ~
                          bio_02 + bio_08 + bio_15 + bio_18 +
                          Eucmean + alt +
                          lc10 + lc11 + lc30 + lc40 + lc70 + lc90 + lc100 + lc130 + lc190)

DensityModel11 <- stepAIC(empty.model1, trace=FALSE, direction="both", scope=variables_BR11); summary(DensityModel11)
DensityModel22 <- stepAIC(empty.model2, trace=FALSE, direction="both", scope=variables_BR22); summary(DensityModel22)
DensityModel33 <- stepAIC(empty.model3, trace=FALSE, direction="both", scope=variables_BR33); summary(DensityModel33)
DensityModel44 <- stepAIC(empty.model4, trace=FALSE, direction="both", scope=variables_BR44); summary(DensityModel44)

WildBoar_NorthernBioregion_DensityModel<-DensityModel11; save("WildBoar_NorthernBioregion_DensityModel", file="./ModelResults/densityModel_WildBoar_NorthernBioregion_Europe.Rdata")
WildBoar_SouthernBioregion_DensityModel<-DensityModel22; save("WildBoar_SouthernBioregion_DensityModel", file="./ModelResults/densityModel_WildBoar_SouthernBioregion_Europe.Rdata")
WildBoar_EasternBioregion_DensityModel<-DensityModel33; save("WildBoar_EasternBioregion_DensityModel",   file="./ModelResults/densityModel_WildBoar_EasternBioregion_Europe.Rdata")
WildBoar_WesternBioregion_DensityModel<-DensityModel44; save("WildBoar_WesternBioregion_DensityModel",   file="./ModelResults/densityModel_WildBoar_WesternBioregion_Europe.Rdata")
st_write(data, 
         dsn = "data_20260207.gpkg",
         layer = "data_spatial_WildBoar_BR", 
         driver = "GPKG", append=TRUE)


#2.3. Calibration plots  validation data
# ----------------------------------------------------------------------------------
data_val_BR1<-data_val_BR1 %>% mutate(Pred_model=predict(DensityModel11, data_val_BR1, type="response"), Pred_Red=Pred_model/10000)
data_val_BR2<-data_val_BR2 %>% mutate(Pred_model=predict(DensityModel22, data_val_BR2, type="response"), Pred_Red=Pred_model/10000)
data_val_BR3<-data_val_BR3 %>% mutate(Pred_model=predict(DensityModel33, data_val_BR3, type="response"), Pred_Red=Pred_model/10000)
data_val_BR4<-data_val_BR4 %>% mutate(Pred_model=predict(DensityModel44, data_val_BR4, type="response"), Pred_Red=Pred_model/10000)



library(Hmisc)
s_class_BRtodas1 <- cut2(data_val_BR1$Pred_Red, g=9) # defining bins (percentiles) on the predicted
s_mean_BRtodas1 <- as.matrix(tapply(data_val_BR1$Pred_Red, s_class_BRtodas1, mean)) # calculating mean values for predicted
y_mean_BRtodas1 <- as.matrix(tapply(data_val_BR1$dens, s_class_BRtodas1, mean)) # calculating mean values for observed
corbr1<-cor.test(data_val_BR1$Pred_Red, data_val_BR1$dens)
data_val_BR1w0<-data_val_BR1 %>%filter(dens_r_trans>0)
s_class_BRtodas1w0<-cut2(data_val_BR1w0$Pred_Red, g=9)
s_mean_BRtodas1w0 <- as.matrix(tapply(data_val_BR1w0$Pred_Red, s_class_BRtodas1w0,mean)) # calculating mean values for predicted
y_mean_BRtodas1w0 <- as.matrix(tapply(data_val_BR1w0$dens, s_class_BRtodas1w0,mean)) # calculating mean values for observed

s_class_BRtodas2 <- cut2(data_val_BR2$Pred_Red, g=9) # defining bins (percentiles) on the predicted
s_mean_BRtodas2 <- as.matrix(tapply(data_val_BR2$Pred_Red, s_class_BRtodas2,mean)) # calculating mean values for predicted
y_mean_BRtodas2 <- as.matrix(tapply(data_val_BR2$dens, s_class_BRtodas2,mean)) # calculating mean values for observed
corbr2<-cor.test(data_val_BR2$Pred_Red, data_val_BR2$dens)
data_val_BR2w0<-data_val_BR2 %>%filter(dens_r_trans>0)
s_class_BRtodas2w0<-cut2(data_val_BR2w0$Pred_Red, g=9)
s_mean_BRtodas2w0 <- as.matrix(tapply(data_val_BR2w0$Pred_Red, s_class_BRtodas2w0,mean)) # calculating mean values for predicted
y_mean_BRtodas2w0 <- as.matrix(tapply(data_val_BR2w0$dens, s_class_BRtodas2w0,mean)) # calculating mean values for observed

s_class_BRtodas3 <- cut2(data_val_BR3$Pred_Red, g=9) # defining bins (percentiles) on the predicted
s_mean_BRtodas3 <- as.matrix(tapply(data_val_BR3$Pred_Red, s_class_BRtodas3,mean)) # calculating mean values for predicted
y_mean_BRtodas3 <- as.matrix(tapply(data_val_BR3$dens, s_class_BRtodas3,mean)) # calculating mean values for observed
corbr3<-cor.test(data_val_BR3$Pred_Red, data_val_BR3$dens)
data_val_BR3w0<-data_val_BR3 %>%filter(dens_r_trans>0)
s_class_BRtodas3w0<-cut2(data_val_BR3w0$Pred_Red, g=9)
s_mean_BRtodas3w0 <- as.matrix(tapply(data_val_BR3w0$Pred_Red, s_class_BRtodas3w0,mean)) # calculating mean values for predicted
y_mean_BRtodas3w0 <- as.matrix(tapply(data_val_BR3w0$dens, s_class_BRtodas3w0,mean)) # calculating mean values for observed

s_class_BRtodas4 <- cut2(data_val_BR4$Pred_Red, g=9) # defining bins (percentiles) on the predicted
s_mean_BRtodas4 <- as.matrix(tapply(data_val_BR4$Pred_Red, s_class_BRtodas4, mean)) # calculating mean values for predicted
y_mean_BRtodas4 <- as.matrix(tapply(data_val_BR4$dens, s_class_BRtodas4, mean)) # calculating mean values for observed
corbr4<-cor.test(data_val_BR4$Pred_Red, data_val_BR4$dens)
data_val_BR4w0<-data_val_BR4 %>%filter(dens_r_trans>0)
s_class_BRtodas4w0<-cut2(data_val_BR4w0$Pred_Red, g=9)
s_mean_BRtodas4w0 <- as.matrix(tapply(data_val_BR4w0$Pred_Red, s_class_BRtodas4w0,mean)) # calculating mean values for predicted
y_mean_BRtodas4w0 <- as.matrix(tapply(data_val_BR4w0$dens, s_class_BRtodas4w0,mean)) # calculating mean values for observed

cpbr<-rbind(data.frame(s_mean=s_mean_BRtodas1,   y_mean=y_mean_BRtodas1, Type=as.factor("All data"), Bioregion="Northern bioregion"),
            data.frame(s_mean=s_mean_BRtodas1w0, y_mean=y_mean_BRtodas1w0, Type=as.factor("Data without 0"), Bioregion="Northern bioregion"),
            data.frame(s_mean=s_mean_BRtodas2,   y_mean=y_mean_BRtodas2, Type=as.factor("All data"), Bioregion="Southern bioregion"),
            data.frame(s_mean=s_mean_BRtodas2w0, y_mean=y_mean_BRtodas2w0, Type=as.factor("Data without 0"), Bioregion="Southern bioregion"),
            data.frame(s_mean= s_mean_BRtodas3, y_mean=y_mean_BRtodas3, Type=as.factor("All data"), Bioregion="Eastern bioregion"),
            data.frame(s_mean=s_mean_BRtodas3w0, y_mean=y_mean_BRtodas3w0, Type=as.factor("Data without 0"), Bioregion="Eastern bioregion"),
            data.frame(s_mean=s_mean_BRtodas4,   y_mean=y_mean_BRtodas4, Type=as.factor("All data"), Bioregion="Western bioregion"),
            data.frame(s_mean=s_mean_BRtodas4w0, y_mean=y_mean_BRtodas4w0, Type=as.factor("Data without 0"), Bioregion="Western bioregion"))

ggplot(cpbr, aes(x=s_mean, y=y_mean, shape=Type, col=Type))+
  geom_point(size=4)+ scale_colour_grey()+ scale_shape_manual(values=c(16, 17))+
  geom_abline(col="red", lty=2, cex=0.7)+
  theme_bw()+
  xlim(0,max(cpbr$s_mean, max(na.omit(cpbr$y_mean))))+
  ylim(0,max(max(na.omit(cpbr$y_mean)), cpbr$s_mean))+
  ylab("Observed hunted bag density")+xlab("Predicted")+
  ggtitle(paste("Calibration plot Wild boar"
  ))+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~Bioregion)#, scales = "free")




#3. Spatial layer: 10x10km grid
# ----------------------------------------------------------------------------------
grid_10km<-st_read("grid_10km_wb.shp") %>% dplyr::select(-BR_T10) %>% 
  rename(NUTS=NUT, Eucmean=Eucalyptus) %>% 
  mutate(NUTS=as.factor(NUTS), 
         Bioregion=as.factor(Bioregion),
         lc_50=as.numeric(lc_50),
         lc_62=as.numeric(lc_62),
         lc_72=as.numeric(lc_72), 
         lc_122=as.numeric(lc_122), 
         lc_153=as.numeric(lc_153), 
         Eucmean = replace_na(Eucmean, 0), Eucmean = ifelse(Eucmean<0, 0,Eucmean))
grid_10km<- grid_10km %>% mutate(x_cen=st_coordinates(st_centroid(grid_10km))[,1], y_cen=st_coordinates(st_centroid(grid_10km))[,2])

# check covariates have same name as modeling covariates
colnames(grid_10km)[32:65]<-gsub("lc_", "lc", colnames(grid_10km)[32:65])
colnames(grid_10km)<-gsub("_mean", "", colnames(grid_10km))
colnames(grid_10km)[9:17]<-gsub("bio_", "bio_0", colnames(grid_10km)[9:17])

Red_  %>% dplyr::select(locationID, locality, country, NUTS, BR_todasr1, area_km2, matches("harv"), 
                        dens, dens_r_trans,
                        matches("bio"), alt, snow, sun, Eucmean, hfp, 
                        matches("lc"), 
                        x_cen, y_cen) %>% 
  rename(Bioregion=BR_todasr1) %>% 
  as_tibble() %>% dplyr::select(17:71) %>% 
  colnames() 

means<-Red_ %>% dplyr::select(locationID, locality, country, NUTS, BR_todasr1, area_km2, matches("harv"),  dens, dens_r_trans, matches("bio"), alt, snow, sun, Eucmean, hfp,  matches("lc"),  x_cen, y_cen) %>% 
                 rename(Bioregion=BR_todasr1) %>% 
                 as_tibble() %>% dplyr::select(17:71) %>% 
                 summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
sdss<-Red_ %>% dplyr::select(locationID, locality, country, NUTS, BR_todasr1, area_km2, matches("harv"),  dens, dens_r_trans, matches("bio"), alt, snow, sun, Eucmean, hfp,  matches("lc"),  x_cen, y_cen) %>% 
                rename(Bioregion=BR_todasr1) %>% 
                as_tibble() %>% dplyr::select(17:71) %>% 
                summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

# standarize grid covariates with territorial unit means and sd!
grid_10km<-grid_10km %>% mutate(across(all_of(colnames(means)), ~ . - means[[cur_column()]])) %>% 
  mutate(across(all_of(colnames(means)), ~ . / sdss[[cur_column()]]))  


library(modEvA) #MESS
# ---------------------------------------------------------------------------------- 
summary(DensityModel11) #Solo ponemos las variables que selecciona el modelo que son sobre las que vamos a hacer la prediccion
data_<-data %>% filter(Bioregion==1) %>% as.data.frame() %>% dplyr::select("bio_02",  "bio_18",
                                                                           "lc11", "lc40", "lc60", "lc90",  "lc180", "lc210", "snow")

#10km
grid_10km_<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km))) %>% filter(Bioregion==1)
grid10km_df<-grid_10km_ %>% as.data.frame()  %>% dplyr::select("id_buena",
                                                               "bio_02",  "bio_18","lc11", "lc40", "lc60", "lc90",  "lc180", "lc210", "snow")
mess_10km<-MESS(V=data_, P=grid10km_df, id.col=1)

# grid_10km_<-grid_10km_[c(-69)]; names(grid_10km_)
mess_10km_<-mess_10km[c(-2:-10)] ; names(mess_10km_)
grid10km_mess1<-grid_10km_ %>% filter(Bioregion==1) %>%  full_join(mess_10km_, by="id_buena")


summary(DensityModel22) #Solo ponemos las variables que selecciona el modelo que son sobre las que vamos a hacer la prediccion
data_<-data %>% filter(Bioregion==2) %>% as.data.frame() %>% dplyr::select("bio_03", "bio_04", "bio_08", "bio_09", "lc10", "lc12", "lc30", "lc40",  "lc60",  "lc70","lc100", "lc120", "lc130", "snow", "Eucmean")
#10km
grid_10km_<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km))) %>% filter(Bioregion==2)
grid10km_df<-grid_10km_ %>% as.data.frame()  %>% dplyr::select("id_buena",
                                                               "bio_03", "bio_04", "bio_08", "bio_09", "lc10", "lc12", "lc30", "lc40",  "lc60",  "lc70","lc100", "lc120", "lc130", "snow", "Eucmean")
mess_10km<-MESS(V=data_, P=grid10km_df, id.col=1)

# grid_10km_<-grid_10km_[c(-1:-69)]
mess_10km_<-mess_10km[c(-2:-16)] ; names(mess_10km_)
grid10km_mess2<-grid_10km_ %>% filter(Bioregion==2) %>%  full_join(mess_10km_, by="id_buena")


summary(DensityModel33) #Solo ponemos las variables que selecciona el modelo que son sobre las que vamos a hacer la prediccion
data_<-data %>% filter(Bioregion==3) %>% as.data.frame() %>% dplyr::select("bio_09", "bio_13", "bio_15",  "lc10",  "lc60", "lc70", "lc90", "lc100", "lc130")
#10km
grid_10km_<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km))) %>% filter(Bioregion==3)
grid10km_df<-grid_10km_ %>% as.data.frame()  %>% dplyr::select("id_buena", 
                                                               "bio_09", "bio_13", "bio_15",  "lc10",  "lc60", "lc70", "lc90", "lc100", "lc130")
mess_10km<-MESS(V=data_, P=grid10km_df, id.col=1)

# grid_10km_<-grid_10km_[c(-1:-69)]
mess_10km_<-mess_10km[c(-2:-10)] ; names(mess_10km_)
grid10km_mess3<-grid_10km_ %>% filter(Bioregion==3) %>%  full_join(mess_10km_, by="id_buena")


summary(DensityModel44) #Solo ponemos las variables que selecciona el modelo que son sobre las que vamos a hacer la prediccion
data_<-data %>% filter(Bioregion==4) %>% as.data.frame() %>% dplyr::select("bio_02", "bio_15", "bio_18", "lc10", "lc11", "lc30", "lc70", "lc100", "lc130", "lc190", "alt", "Eucmean")
#10km
grid_10km_<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km))) %>% filter(Bioregion==4)
grid10km_df<-grid_10km_ %>% as.data.frame()  %>% dplyr::select("id_buena",
                                                               "bio_02", "bio_15", "bio_18", "lc10", "lc11", "lc30", "lc70", "lc100", "lc130", "lc190", "alt", "Eucmean")
mess_10km<-MESS(V=data_, P=grid10km_df, id.col=1)

# grid_10km_<-grid_10km_[c(-1:-69)]
mess_10km_<-mess_10km[c(-2:-13)] ; names(mess_10km_)
grid10km_mess4<-grid_10km_ %>% filter(Bioregion==4) %>%  full_join(mess_10km_, by="id_buena")


# 3.2. save predictions in a same spatial layer
# ----------------------------------------------------------------------------------
pred10km<-rbind(grid10km_mess1 %>% mutate(Pred_model=predict(DensityModel11, grid10km_mess1, type="response")) %>% mutate(Pred=Pred_model/10000),
                grid10km_mess2 %>% mutate(Pred_model=predict(DensityModel22, grid10km_mess2, type="response")) %>% mutate(Pred=Pred_model/10000), 
                grid10km_mess3 %>% mutate(Pred_model=predict(DensityModel33, grid10km_mess3, type="response")) %>% mutate(Pred=Pred_model/10000), 
                grid10km_mess4 %>% mutate(Pred_model=predict(DensityModel44, grid10km_mess4, type="response")) %>% mutate(Pred=Pred_model/10000))

st_write(predRedDeer10km,
         dsn="predicciones_20260203.gpkg",
         layer="DensityModel_grid10km_WildBoar_BR_final4_MESS_BR",
         driver = "GPKG",append=TRUE)

grid10km_mess1<-wb.layer %>% filter(Bioregion==1)
grid10km_mess2<-wb.layer %>% filter(Bioregion==2)
grid10km_mess3<-wb.layer %>% filter(Bioregion==3)
grid10km_mess4<-wb.layer %>% filter(Bioregion==4)

pred1<-predict(WildBoar_NorthernBioregion_DensityModel, grid10km_mess1, type="link", se.fit=TRUE)
pred2<-predict(WildBoar_SouthernBioregion_DensityModel, grid10km_mess2, type="link", se.fit=TRUE)
pred3<-predict(WildBoar_EasternBioregion_DensityModel, grid10km_mess3, type="link", se.fit=TRUE)
pred4<-predict(WildBoar_WesternBioregion_DensityModel, grid10km_mess4, type="link", se.fit=TRUE)

z <- qnorm(0.975)

predRedDeer10km<-rbind( grid10km_mess1 %>% mutate( 
  # pred_model = exp(pred1$fit)/10000,
  pred_lower = exp(pred1$fit - z * pred1$se.fit)/10000,
  pred_upper = exp(pred1$fit + z * pred1$se.fit)/10000),
  
  grid10km_mess2 %>% mutate( 
    # pred_model = exp(pred2$fit)/10000,
    pred_lower = exp(pred2$fit - z * pred2$se.fit)/10000,
    pred_upper = exp(pred2$fit + z * pred2$se.fit)/10000),
  
  grid10km_mess3 %>% mutate( 
    # pred_model = exp(pred3$fit)/10000,
    pred_lower = exp(pred3$fit - z * pred3$se.fit)/10000,
    pred_upper = exp(pred3$fit + z * pred3$se.fit)/10000),
  
  grid10km_mess4 %>% mutate( 
    # pred_model = exp(pred4$fit)/10000,
    pred_lower = exp(pred4$fit - z * pred4$se.fit)/10000,
    pred_upper = exp(pred4$fit + z * pred4$se.fit)/10000)
)

st_write(predRedDeer10km,
         dsn="predicciones_20260203.gpkg",
         layer="DensityModel_grid10km_WildBoar_BR_final4_MESS_BR_IC",
         driver = "GPKG", append=TRUE)
