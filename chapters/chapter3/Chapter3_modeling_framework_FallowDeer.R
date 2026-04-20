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
#Charge Red deer data
spatial<-st_read(dsn="wildRuminants_v5.gpkg", layer='Fd')

# IMPORTANT WARNING1: 
# For removing Bayern federal estate it's needed that all rows have a locality not NA value. 
# Therefore, those ones that are na we'll record their locationID+country there
spatialID<-spatial %>% 
  mutate(locality=ifelse(is.na(locality), paste(locationID, country), locality)) %>% 
  dplyr::filter(locality!='Bayern') 

spatialID<-spatialID %>% mutate(locationID=ifelse(country=="PL", paste("PL_H_", locality, sep="_"), 
                                                        ifelse(country=="Estonia", paste("EE_H_", locality, sep="_"), 
                                                               ifelse(country=="Norway", paste("NO_", locality, sep="_"), locationID))), 
                                      country=ifelse(country=="ES" | country=="Spain", "Spain", country))

spatial1<-spatialID %>% 
  dplyr::select(starts_with("lc"), ends_with("_mean"), BR_todasr1, area_km2, country, dataTime, harvTot, locationID) %>% 
  mutate(dataTime=as.numeric(substr(dataTime, 1, 4))) %>% 
  mutate(harvTot=ifelse(harvTot>0, harvTot, 0)) 
spatial1<-spatial1 %>% 
  as_tibble() %>%
  filter(dataTime>2014) %>%  
  filter(country!="Sweden" | dataTime!=2012) %>%  
  filter(country!="Sweden" | dataTime!=2013) %>%
  filter(country!="Sweden" | dataTime!=2014) %>%
  filter(country!="Sweden" | dataTime!=2015) %>%
  filter(country!="Sweden" | dataTime!=2016) %>%
  filter(country!="Sweden" | dataTime!=2017) %>%
  filter(country!="Sweden" | dataTime!=2018) %>%
  mutate(dataTime = str_c("harvTot", dataTime)) %>%
  pivot_wider(names_from = dataTime, values_from = c("harvTot"), names_sep="", values_fn=sum)  # long to wide

unique<-spatialID %>% mutate(locationID_country=paste(locationID, country, sep="_"))
unique<-unique[!duplicated(unique$locationID_country), ]

spatial_<- unique %>% dplyr::select(-harvTot, -species, -dataTime) %>%  
  inner_join(spatial1 %>% dplyr::select(!matches("bio") & !matches("lc") & !matches("mean"), -c(BR_todasr1, geom, area_km2)), by=c("locationID", "country")) %>%
  rename(alt=alt_mean, Eucmean=EucalyptusSpp_mean, hfp=hfp_mean, snow=snow_mean,  sun=sun_mean) %>% 
  mutate(NUTS=as.factor(ifelse(NUTS=="comuni" | NUTS=="Concelho" | NUTS=="District" | NUTS=="Kommuner" | NUTS=="Municipality", 1, ifelse(NUTS=="NUTS0" |NUTS=="NUTS1" |NUTS=="NUTS2" |NUTS=="NUTS3", 2, 0))),
         locationID=as.factor(locationID),locality=as.factor(locality), country=as.factor(country), 
         Eucmean = replace_na(Eucmean, 0)) 
spatial_<-spatial_ %>% mutate(area_km2=st_area(spatial_)/1000000) ; attributes(spatial_$area_km2) = NULL

spatial_ <-spatial_ %>% 
  mutate(harvMax=pmax(spatial_$harvTot2015, spatial_$harvTot2016, spatial_$harvTot2017, spatial_$harvTot2018, spatial_$harvTot2019, spatial_$harvTot2020, spatial_$harvTot2021, na.rm = TRUE)) %>% 
  mutate(harvMax=ifelse(is.na(harvMax), 0, harvMax)) %>% 
  mutate(dens=round(harvMax/area_km2, 2)) %>% 
  mutate(dens_r_trans=dens*10000) %>% 
  dplyr::select(locationID, locality, country, NUTS, #species,
                matches("excel"), matches("shp"),
                matches("bio"), grow_mean, alt, Eucmean, hfp, snow, sun, area_km2, BR_todasr1, #Bioregion,
                matches("lc"), 
                harvTot2015,harvTot2016, harvTot2017,harvTot2018, harvTot2019,harvTot2020, harvTot2021, harvMax, dens, dens_r_trans, suit_mean)

spatial_<- spatial_ %>% mutate(x_cen=st_coordinates(st_centroid(spatial_))[,1], y_cen=st_coordinates(st_centroid(spatial_))[,2])  


df_<-spatial_ %>%
  mutate(lc10=lc_10_sum/lc_10_count*100, lc11=lc_11_sum/lc_11_count*100, lc12=lc_12_sum/lc_12_count*100, lc20=lc_20_sum/lc_20_count*100, lc30=lc_30_sum/lc_30_count*100, lc40=lc_40_sum/lc_40_count*100, lc50=lc_50_sum/lc_50_count*100, lc60=lc_60_sum/lc_60_count*100, lc61=lc_61_sum/lc_61_count*100, lc62=lc_62_sum/lc_62_count*100, lc70=lc_70_sum/lc_70_count*100, lc71=lc_71_sum/lc_71_count*100, lc72=lc_72_sum/lc_72_count*100, lc80=lc_80_sum/lc_80_count*100, lc81=lc_81_sum/lc_81_count*100, lc82=lc_82_sum/lc_82_count*100, lc90=lc_90_sum/lc_90_count*100, lc100=lc_100_sum/lc_100_count*100, lc110=lc_110_sum/lc_110_count*100,lc120=lc_120_sum/lc_120_count*100,lc121=lc_121_sum/lc_121_count*100,lc122=lc_122_sum/lc_122_count*100,lc130=lc_130_sum/lc_130_count*100,lc140=lc_140_sum/lc_140_count*100,lc150=lc_150_sum/lc_150_count*100,lc152=lc_152_sum/lc_152_count*100,lc153=lc_153_sum/lc_153_count*100,lc160=lc_160_sum/lc_160_count*100,lc170=lc_170_sum/lc_170_count*100,lc180=lc_180_sum/lc_180_count*100,lc190=lc_190_sum/lc_190_count*100,lc200=lc_200_sum/lc_200_count*100,lc201=lc_201_sum/lc_201_count*100,lc202=lc_202_sum/lc_202_count*100,lc210=lc_210_sum/lc_210_count*100,lc220=lc_220_sum/lc_220_count*100) %>% 
  mutate(Eucmean = ifelse(Eucmean<0, 0.01, Eucmean)) %>% 
  dplyr::select(!matches("_cou") & !matches("_sum")) %>% 
  # remove columns that are not covariates of the model or have weird values
  dplyr::select(-c(excel_file,excel_path, shp_file, shp_path, 
                   grow_mean, #weird values
                   lc50, lc62, lc72, lc81, lc82,lc121, lc170#, # all values are 0
  )) %>% #quito las bioregiones porque no las vamos a usar para el modelo
  filter(dens<50) # remove silver polygons 

df_<-df_[!is.na(df_$dens),] #Remove NA values from RV: it does not remove any value

colnames(df_)<-gsub("_mean", "", colnames(df_))

df_<-df_ %>% dplyr::select(locationID, locality, country, NUTS, BR_todasr1, area_km2, matches("harv"), 
                             dens, dens_r_trans,
                             matches("bio"), alt, snow, sun, Eucmean, hfp, 
                             matches("lc"), 
                             x_cen, y_cen) %>% 
  dplyr::select(-c("individualCount...39", "individualCount...51")) %>% 
  rename(Bioregion=BR_todasr1)

#remove NA values from the dataset
summary(df_)
df_<-df_[!is.na(df_$bio_01),] #Remove NA values from the Red_set
df_<-df_[!is.na(df_$snow),] #Remove NA values from the Red_set
df_<-df_[!is.na(df_$hfp),] #Remove NA values from the Red_set
df_<-df_[!is.na(df_$alt),] #Remove NA values from the Red_set
df_<-df_[!is.na(df_$lc10),] #Remove NA values from the Red_set
summary(df_)

# Modeling phase
# 1) data
# scale all variables.
vars_std <- c(colnames(df_ %>% as_tibble())[c(17:71)])
data<-df_  %>% mutate(across(all_of(vars_std),~ as.numeric(scale(.))))

#1.1 Create calibration and validation data set (80% / 20%)
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]

#2.0 VIF <2
#2.1 select covariates for the VIF
hist(data_tra$Eucmean)
hist(data_tra$alt) 
hist(data_tra$snow)

hist(data_tra$lc10)
hist(data_tra$lc11)
hist(data_tra$lc12) #
hist(data_tra$lc20) #
hist(data_tra$lc30)
hist(data_tra$lc40) 
hist(data_tra$lc60) 
hist(data_tra$lc70)
hist(data_tra$lc90) 
hist(data_tra$lc100)
hist(data_tra$lc120) 
hist(data_tra$lc130)

hist(data_tra$lc61) #
hist(data_tra$lc71) #
hist(data_tra$lc80) #
hist(data_tra$lc110) #
hist(data_tra$lc122) #
hist(data_tra$lc140) #
hist(data_tra$lc150) #
hist(data_tra$lc152) #
hist(data_tra$lc153) #
hist(data_tra$lc160) #
hist(data_tra$lc180) #
hist(data_tra$lc190) #
hist(data_tra$lc200) #
hist(data_tra$lc201) #
hist(data_tra$lc202) #
hist(data_tra$lc210) #
hist(data_tra$lc220) #

#2.1 vif
variables<-data_tra %>% dplyr::select(matches("bio_"), alt, snow, Eucmean,lc10, lc11, lc30, lc40, lc60, lc70, #lc90, 
                                      lc100, lc120, 
                                      lc130, x_scale, y_scale
) ; st_geometry(variables)<-NULL
resultado.vif0<-usdm::vifstep(variables, th=2)

#1.1) Density Models: negative binomial need integer -> transform density *100
library(MASS)
empty.model0<-glm.nb(dens_r_trans ~ 1, data=data_tra %>% filter(dens_r_trans > 0), link=log)
resultado.vif0
variables0<-formula(dens_r_trans ~ x_scale + y_scale +
                      bio_07_mean +
                      alt + snow + Eucmean +
                      lc10 + lc30 + lc40 + lc60 + lc70 + lc100 + lc120 + lc130)


DensityModel0_7 <- stepAIC(empty.model0, trace=FALSE, direction="both", scope=variables0); summary(DensityModel0_7)
Fd_DensityModel<-DensityModel0_7
save("Fd_DensityModel", file="densityModel_Fd_Europe.Rdata")


#Predictions to Spatial layer: 10x10km grid
grid_10km<-st_read("grid_10km_wb.shp") %>% dplyr::select(-Bioregion) %>% 
  rename(Bioregion=BR_T10, NUTS=NUT, Eucmean=Eucalyptus) %>% 
  mutate(NUTS=as.factor(NUTS), 
         Bioregion=as.factor(Bioregion),
         lc_50=as.numeric(lc_50),
         lc_62=as.numeric(lc_62),
         lc_72=as.numeric(lc_72), 
         lc_122=as.numeric(lc_122), 
         lc_153=as.numeric(lc_153), 
         Eucmean = replace_na(Eucmean, 0),
         Eucmean = ifelse(Eucmean<0, 0, Eucmean)) %>% 
  mutate(x_cen=st_coordinates(st_centroid(geometry))[,1], y_cen=st_coordinates(st_centroid(geometry))[,2])

colnames(grid_10km)[32:65]<-gsub("lc_", "lc", colnames(grid_10km)[32:65])
colnames(grid_10km)<-gsub("_mean", "", colnames(grid_10km))
colnames(grid_10km)[9:17]<-gsub("bio_", "bio_0", colnames(grid_10km)[9:17])

means<-df_ %>% as_tibble() %>% dplyr::select(17:71) %>% 
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
sdss<-df_ %>% as_tibble() %>% dplyr::select(17:71) %>% 
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

# standarize grid covariates with territorial unit means and sd!
grid_10km<-grid_10km %>% mutate(across(all_of(colnames(means)), ~ . - means[[cur_column()]])) %>% 
  mutate(across(all_of(colnames(means)), ~ . / sdss[[cur_column()]])) 

library(modEvA) #MESS
summary(DensityModel0_7) #Solo ponemos las variables que selecciona el modelo que son sobre las que vamos a hacer la prediccion

data_<-data %>% as.data.frame() %>% dplyr::select("x_cen", "y_cen",
                                                  "lc60", "lc30", "lc40", "Eucmean", "lc100", "lc70", "alt", "snow", "lc130")
#10km
grid_10km<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km)))
grid10km_df<-grid_10km %>% as.data.frame()  %>% dplyr::select("id_buena", 
                                                              "x_cen", "y_cen",
                                                              "lc60", "lc30", "lc40", "Eucmean", "lc100", "lc70", "alt", "snow", "lc130")
mess_10km<-MESS(V=data_, P=grid10km_df, id.col=1)

grid_10km_<-grid_10km[c(-1:-7)]
mess_10km_<-mess_10km[c(-2:-10)] ; names(mess_10km_)
grid10km_mess<-grid_10km %>%  full_join(mess_10km_, by="id_buena")

pred10km<-grid10km_mess %>% mutate(Pred_model=predict(DensityModel0_7, grid_10km, type="response")) %>% mutate(Pred=Pred_model/10000)
st_write(pred10km, 
         dsn="predicciones_20260203.gpkg",
         layer="DensityModel_grid10km_Fd_final_xy_MESS", driver = "GPKG",append=TRUE)
