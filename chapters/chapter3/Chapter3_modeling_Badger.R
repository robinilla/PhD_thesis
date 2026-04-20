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

library(tidyverse)
library(sf)
library(viridis)
# rm(list=ls())

#Charge Badger data
data_spatial<-st_read(dsn="Europe_Carnivores_WLDM_compiled_v2.0.gpkg", layer='Badger')
data_spatial<-data_spatial %>% mutate(BadgerSuit=round(BadgerSuit, 7), area_km2=round(area_km2, 2), alt_mean=round(alt_mean, 2), grow_mean=round(grow_mean, 4), snow_mean=round(snow_mean, 5), hfp_mean=round(hfp_mean, 3), sun_mean=round(sun_mean, 3), euc_mean=round(euc_mean, 7)) %>% 
  mutate(locationID=ifelse(country=="Norway", paste("NO_", locality, sep=""), locationID))  %>% 
  mutate(locationID_cntry=paste(locationID, country, sep="_")) %>% 
  mutate(locality=ifelse(is.na(locality), paste(locationID, country), locality))
  
data_spatialID<-data_spatial %>% dplyr::filter(locality!='Bayern')

# IMPORTANT WARNING1: 
# For removing Bayern federal estate it's needed that all rows have a locality not NA value. 
data_spatial1<-data_spatial %>% dplyr::filter(locality!='Bayern') %>% 
                  mutate(dataTime=as.numeric(substr(dataTime, 1, 4)), harvTot=ifelse(harvTot>0, harvTot, 0),
                         lc10=lc_10_sum/lc_10_coun*100, lc11=lc_11_sum/lc_11_coun*100, lc12=lc_12_sum/lc_12_coun*100, lc20=lc_20_sum/lc_20_coun*100, lc30=lc_30_sum/lc_30_coun*100, lc40=lc_40_sum/lc_40_coun*100, lc50=lc_50_sum/lc_50_coun*100, lc60=lc_60_sum/lc_60_coun*100, lc61=lc_61_sum/lc_61_coun*100, lc62=lc_62_sum/lc_62_coun*100, lc70=lc_70_sum/lc_70_coun*100, lc71=lc_71_sum/lc_71_coun*100, lc72=lc_72_sum/lc_72_coun*100, lc80=lc_80_sum/lc_80_coun*100, lc81=lc_81_sum/lc_81_coun*100, lc82=lc_82_sum/lc_82_coun*100, lc90=lc_90_sum/lc_90_coun*100, lc100=lc_100_sum/lc_100_cou*100, lc110=lc_110_sum/lc_110_cou*100,lc120=lc_120_sum/lc_120_cou*100,lc121=lc_121_sum/lc_121_cou*100,lc122=lc_122_sum/lc_122_cou*100,lc130=lc_130_sum/lc_130_cou*100,lc140=lc_140_sum/lc_140_cou*100,lc150=lc_150_sum/lc_150_cou*100,lc152=lc_152_sum/lc_152_cou*100,lc153=lc_153_sum/lc_153_cou*100,lc160=lc_160_sum/lc_160_cou*100,lc170=lc_170_sum/lc_170_cou*100,lc180=lc_180_sum/lc_180_cou*100,lc190=lc_190_sum/lc_190_cou*100,lc200=lc_200_sum/lc_200_cou*100,lc201=lc_201_sum/lc_201_cou*100,lc202=lc_202_sum/lc_202_cou*100,lc210=lc_210_sum/lc_210_cou*100,lc220=lc_220_sum/lc_220_cou*100) %>% 
                  dplyr::select(locationID_cntry, locationID, country, dataTime, harvTot, area_km2, starts_with("lc"), ends_with("_mean"), BR_todasr1, BadgerSuit)


data_spatial1<-data_spatial1 %>% as_tibble() %>%
  filter(dataTime>2011) %>%  #filtro los a?os que vamos a usar para pasar de formato long a wide
  filter(country!="Sweden" | dataTime!=2012) %>%  #quito los datos de Sweden que tienen peor resolucion espacial
  filter(country!="Sweden" | dataTime!=2013) %>%
  filter(country!="Sweden" | dataTime!=2014) %>%
  filter(country!="Sweden" | dataTime!=2015) %>%
  filter(country!="Sweden" | dataTime!=2016) %>%
  filter(country!="Sweden" | dataTime!=2017) %>%
  filter(country!="Sweden" | dataTime!=2018) %>%
  filter(country!="Poland" | dataTime!=2021) %>% 
  filter(country!="Hungary" | dataTime!=2014) %>% 
  filter(country!="Hungary" | dataTime!=2015) %>% 
  filter(country!="Hungary" | dataTime!=2016) %>% 
  mutate(dataTime = str_c("harvTot", dataTime)) %>% 
  #dplyr::select(-geometry) %>% 
  dplyr::select(-geom) %>% 
  pivot_wider(names_from = dataTime, values_from = c("harvTot"), names_sep="", values_fn=sum)  # paso de long a wide


data_spatial1[duplicated(data_spatial1$locationID_cntry),] %>% 
  mutate(country=as.factor(country)) %>% dplyr::select(country) %>% summary()


data_unique<-data_spatialID %>% mutate(dataTime=as.numeric(substr(dataTime, 1, 4)))
data_unique<-data_unique %>% 
  filter(dataTime>2011) %>%   # filter years that we're going to use for transforming data (long-wide)
  filter(country!="Sweden" | dataTime!=2012) %>%    # retain Sweden finest data 
  filter(country!="Sweden" | dataTime!=2013) %>%
  filter(country!="Sweden" | dataTime!=2014) %>%
  filter(country!="Sweden" | dataTime!=2015) %>%
  filter(country!="Sweden" | dataTime!=2016) %>%
  filter(country!="Sweden" | dataTime!=2017) %>%
  filter(country!="Sweden" | dataTime!=2018) %>%
  filter(country!="Poland" | dataTime!=2021) %>% 
  filter(country!="Hungary" | dataTime!=2014) %>% 
  filter(country!="Hungary" | dataTime!=2015) %>% 
  filter(country!="Hungary" | dataTime!=2016)

# remove duplicated data by locationID_country  
# --> expecting to have just one administrative unit per country and remove duplicated polygons
data_unique<-data_unique %>% mutate(locationID_cntry=paste(locationID, country, sep="_"))
data_unique<-data_unique[!duplicated(data_unique$locationID_cntry), ]

# Join the non duplicated polygon layer with wide table by locationID and country
# remove columns from each part that may cause confusion
data_spatial_<- data_unique %>% dplyr::select(-harvTot, -species, -dataTime) %>%  
  inner_join(data_spatial1 %>% dplyr::select(!matches("bio") & !matches("lc") & !matches("mean"), -c(BR_todasr1, area_km2, BadgerSuit, locationID_cntry)), by=c("locationID", "country")) %>%
  rename(alt=alt_mean, Eucmean=euc_mean, hfp=hfp_mean, snow=snow_mean,  sun=sun_mean) %>% 
  mutate(NUT=as.factor(ifelse(NUT=="HG"  |NUT=="Hunting grounds"  |NUT=="Hunting ground", 0, 
                              ifelse(NUT=="NUT0" |NUT=="NUT1" |NUT=="NUTS1" |NUT=="NUT2" |NUT=="NUT3" |NUT=="NUTS3" |NUT=="Country", 2, 1))),
         locationID=as.factor(locationID),locality=as.factor(locality), country=as.factor(country), 
         Eucmean = replace_na(Eucmean, 0)) 

# calculate area (km2)
data_spatial_<-data_spatial_ %>% mutate(area_km2=st_area(data_spatial_)/1000000) ; attributes(data_spatial_$area_km2) = NULL

data_spatial_ <-data_spatial_ %>% 
  # calculate harvest Maximum between the years that we have data. 
  # Be aware if new data is introduced you may add
  # some more for example Red_spatial_$harvTot2022...etc.
  mutate(harvMax=pmax(harvTot2012, harvTot2013, harvTot2014, harvTot2015, harvTot2016, harvTot2017, harvTot2018, harvTot2019, harvTot2020, harvTot2021, na.rm = TRUE)) %>% 
  mutate(harvMax=ifelse(is.na(harvMax), 0, harvMax)) %>% 
  # calculate density
  mutate(dens=round(harvMax/area_km2, 2)) %>% 
  mutate(dens_r_trans=dens*10000) %>% 
  dplyr::select(locationID, locality, country, NUT, #species,
                matches("excel"), matches("shp"),
                matches("bio"), grow_mean, alt, Eucmean, hfp, snow, sun, area_km2, BR_todasr1, #Bioregion, matches("suit") ,
                matches("lc"), BadgerSuit,
                harvTot2012, harvTot2013, harvTot2014, harvTot2015,harvTot2016, harvTot2017,harvTot2018, harvTot2019,harvTot2020, harvTot2021, harvMax, dens, dens_r_trans)
# calculate coordinates
data_spatial_<- data_spatial_ %>% 
  mutate(x_cen=st_coordinates(st_centroid(geom))[,1], y_cen=st_coordinates(st_centroid(geom))[,2])
#st_write(data_spatial_, dsn="H:/IREC_Sonia/2_WIP/Report_202210/2_WIP/2_GIS/europe_carnivores/predicciones_badger_20221215.gpkg", layer="data_spatial_v1.1", driver = "GPKG",append=TRUE)


summary(data_spatial_ )
data_<-data_spatial_ %>%
  mutate(lc10=lc_10_sum/lc_10_coun*100, lc11=lc_11_sum/lc_11_coun*100, lc12=lc_12_sum/lc_12_coun*100, lc20=lc_20_sum/lc_20_coun*100, lc30=lc_30_sum/lc_30_coun*100, lc40=lc_40_sum/lc_40_coun*100, lc50=lc_50_sum/lc_50_coun*100, lc60=lc_60_sum/lc_60_coun*100, lc61=lc_61_sum/lc_61_coun*100, lc62=lc_62_sum/lc_62_coun*100, lc70=lc_70_sum/lc_70_coun*100, lc71=lc_71_sum/lc_71_coun*100, lc72=lc_72_sum/lc_72_coun*100, lc80=lc_80_sum/lc_80_coun*100, lc81=lc_81_sum/lc_81_coun*100, lc82=lc_82_sum/lc_82_coun*100, lc90=lc_90_sum/lc_90_coun*100, lc100=lc_100_sum/lc_100_cou*100, lc110=lc_110_sum/lc_110_cou*100,lc120=lc_120_sum/lc_120_cou*100,lc121=lc_121_sum/lc_121_cou*100,lc122=lc_122_sum/lc_122_cou*100,lc130=lc_130_sum/lc_130_cou*100,lc140=lc_140_sum/lc_140_cou*100,lc150=lc_150_sum/lc_150_cou*100,lc152=lc_152_sum/lc_152_cou*100,lc153=lc_153_sum/lc_153_cou*100,lc160=lc_160_sum/lc_160_cou*100,lc170=lc_170_sum/lc_170_cou*100,lc180=lc_180_sum/lc_180_cou*100,lc190=lc_190_sum/lc_190_cou*100,lc200=lc_200_sum/lc_200_cou*100,lc201=lc_201_sum/lc_201_cou*100,lc202=lc_202_sum/lc_202_cou*100,lc210=lc_210_sum/lc_210_cou*100,lc220=lc_220_sum/lc_220_cou*100) %>% 
  dplyr::select(!matches("_cou") & !matches("_sum")) %>% 
  # remove columns that are not covariates of the model or have weird values
  dplyr::select(-c(excel_file,excel_path, shp_file, shp_path, 
                   grow_mean, #weird values
                   lc50, lc62, lc72, lc81, lc82,lc121, lc170 # all values are 0
  )) %>% 
  filter(dens<50) # remove silver polygons 


data_$Eucmean<-ifelse(data_$Eucmean<0, 0.01, data_$Eucmean) #remove negative values from Eucmean
data_<-data_[!is.na(data_$dens),] #Remove NA values from RV: it does not remove any value

# colnames(data_)
data_<-data_ %>% dplyr::select(locationID, locality, country, NUT, BR_todasr1, Bioregion, area_km2, matches("harv"), dens, dens_r_trans, IndividualCount, BadgerSuit,
                matches("bio"), alt, snow, sun, Eucmean, hfp, 
                matches("lc"), 
                x_cen, y_cen)

#Modeling phase
#1) data
#scale all variables.
vars_std <- c(colnames(data_ %>% as_tibble())[c(23:77)])
data_<-data_ %>% filter(country!="North Macedonia") %>%  mutate(country=factor(country)) 
#remove NA values from the dataset
summary(data)
data<-data[!is.na(data$bio_1),] #Remove NA values from the dataset
data<-data[!is.na(data$snow),] #Remove NA values from the dataset
data<-data[!is.na(data$hfp),] #Remove NA values from the dataset
data<-data[!is.na(data$alt),] #Remove NA values from the dataset
data<-data[!is.na(data$lc10),] #Remove NA values from the dataset

data<-data_  %>% mutate(across(all_of(vars_std),~ as.numeric(scale(.))))



#1.1 Create calibration and validation data set (80% / 20%)
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]

#2.0 VIF <2
#2.1 select covariates for the VIF
hist(data_tra$Eucmean)
hist(data_tra$alt) 
hist(data_tra$snow)
hist(data_tra$hfp)
hist(data_tra$sun)

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
hist(data_tra$lc120) #
hist(data_tra$lc130)

hist(data_tra$lc61) #
hist(data_tra$lc71) #
hist(data_tra$lc80) #
hist(data_tra$lc110) #
hist(data_tra$lc120) #
hist(data_tra$lc122) #
hist(data_tra$lc140) #
hist(data_tra$lc150) #
hist(data_tra$lc152) #
# hist(data_tra$lc153) #
hist(data_tra$lc160) #
hist(data_tra$lc180) #
hist(data_tra$lc190) 
hist(data_tra$lc200) #
hist(data_tra$lc201) #
hist(data_tra$lc202) #
hist(data_tra$lc210) #
hist(data_tra$lc220) #

#2.1 vif
variables<-data_tra %>% dplyr::select(matches("bio_"), alt, snow, sun, hfp,lc10, lc11,
                                      lc30, lc40, lc60, lc70, lc90, lc100, lc130, lc190 #, x_scale, y_scale
) ; st_geometry(variables)<-NULL
resultado.vif0<-usdm::vifstep(variables, th=2)
resultado.vif0
#2.2) Density Models:  negative binomial needs an integer value, so that multiply density variable  *10000
library(MASS)
library(glmmTMB)
data.f<-data_tra %>% filter(dens_r_trans>0)
DensityModel_notMK_without0_random_TMB = glmmTMB(dens_r_trans ~ x_cen+ y_cen+
                                                   bio_3+bio_4+bio_15+alt+lc10+lc30+lc40+lc60+lc70+lc90+lc100+lc130+lc190+ (1 | country), 
                                                 family = nbinom2, 
                                                 data=data.f, REML=F)
summary(DensityModel_notMK_without0_random_TMB); #plot(effects::allEffects(DensityModel0_7))

Badger_DensityModel<-DensityModel_notMK_without0_random_TMB
save(Badger_DensityModel, file="densityModel_Badger_glmmmTMB_Europe.Rdata")

#Calibration plots - validation data
data_val<-data_val %>% mutate(Pred_model=predict(Badger_DensityModel, data_val, type="response"), Pred_species=Pred_model/10000)

library(Hmisc)
s_class <- cut2(data_val$Pred_species, g=9) # defining bins (percentiles) on the predicted
s_mean <- as.matrix(tapply(data_val$Pred_species, s_class, mean)) # calculating mean values for predicted
y_mean <- as.matrix(tapply(data_val$dens, s_class, mean)) # calculating mean values for observed
corbr1<-cor.test(data_val$Pred_species, data_val$dens)
data_val_w0<-data_val %>%filter(dens_r_trans>0)
s_class_w0 <- cut2(data_val_w0$Pred_species, g=9) # defining bins (percentiles) on the predicted
s_mean_w0 <- as.matrix(tapply(data_val_w0$Pred_species, s_class_w0, mean)) # calculating mean values for predicted
y_mean_w0 <- as.matrix(tapply(data_val_w0$dens, s_class_w0, mean)) # calculating mean values for observed
cpE<-rbind(data.frame(s_mean=s_mean,   y_mean=y_mean, Type=as.factor("All data")),
           data.frame(s_mean=s_mean_w0, y_mean=y_mean_w0, Type=as.factor("Data without 0"))
)

cpE1xy<-data.frame(s_mean=s_mean,   y_mean=y_mean, Type=as.factor("All data"), model=as.factor("DensityModelRandom1xy_without0TMB"))

ggplot(cpE, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
  geom_point(size=4)+ scale_colour_grey()+ scale_shape_manual(values=c(16, 17))+
  geom_abline(col="red", lty=2, cex=0.5)+
  theme_bw()+
  ylab(bquote("Observed hunting bag/km "^2))+xlab(bquote("Predicted hunting bag/km"^2 *" with GLMM"))+
  ggtitle(paste("Hunted badger density"))+
  theme(plot.title = element_text(hjust = 0.5))#+
facet_wrap(~model)


#2.3.1.) Calibration plots for each bioregion
data_val_BR1<-data_val %>% filter(Bioregion==1)
data_val_BR2<-data_val %>% filter(Bioregion==2)
data_val_BR3<-data_val %>% filter(Bioregion==3)
data_val_BR4<-data_val %>% filter(Bioregion==4)

s_class_BRtodas1 <- cut2(data_val_BR1$Pred_species, g=9) # defining bins (percentiles) on the predicted
s_mean_BRtodas1 <- as.matrix(tapply(data_val_BR1$Pred_species, s_class_BRtodas1, mean)) # calculating mean values for predicted
y_mean_BRtodas1 <- as.matrix(tapply(data_val_BR1$dens, s_class_BRtodas1, mean)) # calculating mean values for observed
corbr1<-cor.test(data_val_BR1$Pred_species, data_val_BR1$dens)
data_val_BR1w0<-data_val_BR1 %>%filter(dens_r_trans>0)
s_class_BRtodas1w0<-cut2(data_val_BR1w0$Pred_species, g=9)
s_mean_BRtodas1w0 <- as.matrix(tapply(data_val_BR1w0$Pred_species, s_class_BRtodas1w0,mean)) # calculating mean values for predicted
y_mean_BRtodas1w0 <- as.matrix(tapply(data_val_BR1w0$dens, s_class_BRtodas1w0,mean)) # calculating mean values for observed

# s_class_BRtodas2 <- cut2(data_val_BR2$Pred_species, g=9) # defining bins (percentiles) on the predicted
# s_mean_BRtodas2 <- as.matrix(tapply(data_val_BR2$Pred_species, s_class_BRtodas2,mean)) # calculating mean values for predicted
# y_mean_BRtodas2 <- as.matrix(tapply(data_val_BR2$dens, s_class_BRtodas2,mean)) # calculating mean values for observed
# corbr2<-cor.test(data_val_BR2$Pred_species, data_val_BR2$dens)
# data_val_BR2w0<-data_val_BR2 %>%filter(dens_r_trans>0)
# s_class_BRtodas2w0<-cut2(data_val_BR2w0$Pred_species, g=9)
# s_mean_BRtodas2w0 <- as.matrix(tapply(data_val_BR2w0$Pred_species, s_class_BRtodas2w0,mean)) # calculating mean values for predicted
# y_mean_BRtodas2w0 <- as.matrix(tapply(data_val_BR2w0$dens, s_class_BRtodas2w0,mean)) # calculating mean values for observed
# 
s_class_BRtodas3 <- cut2(data_val_BR3$Pred_species, g=9) # defining bins (percentiles) on the predicted
s_mean_BRtodas3 <- as.matrix(tapply(data_val_BR3$Pred_species, s_class_BRtodas3,mean)) # calculating mean values for predicted
y_mean_BRtodas3 <- as.matrix(tapply(data_val_BR3$dens, s_class_BRtodas3,mean)) # calculating mean values for observed
corbr3<-cor.test(data_val_BR3$Pred_species, data_val_BR3$dens)
data_val_BR3w0<-data_val_BR3 %>%filter(dens_r_trans>0)
s_class_BRtodas3w0<-cut2(data_val_BR3w0$Pred_species, g=9)
s_mean_BRtodas3w0 <- as.matrix(tapply(data_val_BR3w0$Pred_species, s_class_BRtodas3w0,mean)) # calculating mean values for predicted
y_mean_BRtodas3w0 <- as.matrix(tapply(data_val_BR3w0$dens, s_class_BRtodas3w0,mean)) # calculating mean values for observed

s_class_BRtodas4 <- cut2(data_val_BR4$Pred_species, g=9) # defining bins (percentiles) on the predicted
s_mean_BRtodas4 <- as.matrix(tapply(data_val_BR4$Pred_species, s_class_BRtodas4, mean)) # calculating mean values for predicted
y_mean_BRtodas4 <- as.matrix(tapply(data_val_BR4$dens, s_class_BRtodas4, mean)) # calculating mean values for observed
corbr4<-cor.test(data_val_BR4$Pred_species, data_val_BR4$dens)
data_val_BR4w0<-data_val_BR4 %>%filter(dens_r_trans>0)
s_class_BRtodas4w0<-cut2(data_val_BR4w0$Pred_species, g=9)
s_mean_BRtodas4w0 <- as.matrix(tapply(data_val_BR4w0$Pred_species, s_class_BRtodas4w0,mean)) # calculating mean values for predicted
y_mean_BRtodas4w0 <- as.matrix(tapply(data_val_BR4w0$dens, s_class_BRtodas4w0,mean)) # calculating mean values for observed

cpbr<-rbind(data.frame(s_mean=s_mean_BRtodas1,   y_mean=y_mean_BRtodas1, Type=as.factor("All data"), Bioregion="Northern bioregion"),
            #data.frame(s_mean=s_mean_BRtodas2,   y_mean=y_mean_BRtodas2, Type=as.factor("All data"), Bioregion="Southern bioregion"),
            data.frame(s_mean= s_mean_BRtodas3, y_mean=y_mean_BRtodas3, Type=as.factor("All data"), Bioregion="Eastern bioregion"),
            data.frame(s_mean=s_mean_BRtodas4,   y_mean=y_mean_BRtodas4, Type=as.factor("All data"), Bioregion="Western bioregion")
)

ggplot(cpbr, aes(x=s_mean, y=y_mean, #shape=Type, col=Type
))+
  geom_point(size=4)+ scale_colour_grey()+ scale_shape_manual(values=c(16, 17))+
  geom_abline(col="red", lty=2, cex=0.7)+
  theme_bw()+
  xlim(0,max(cpbr$s_mean, max(na.omit(cpbr$y_mean))))+
  ylim(0,max(max(na.omit(cpbr$y_mean)), cpbr$s_mean))+
  ylab("Observed hunted bag density")+xlab("Predicted")+
  ggtitle(paste("Calibration plot badger"))+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~Bioregion, ncol=2)#, scales = "free")



#Spatial layer: 10x10km grid
grid_10km<-st_read("grid10_ToPredict_constrained.gpkg",  
                   layer="grid_10km_badger") %>% dplyr::select(-Bioregion) %>% 
  rename(Bioregion=BR_T10, NUTS=NUT, Eucmean=Eucalyptus) %>% 
  mutate(NUTS=as.factor(NUTS), 
         Bioregion=as.factor(Bioregion),lc_50=as.numeric(lc_50), lc_62=as.numeric(lc_62), lc_72=as.numeric(lc_72), lc_122=as.numeric(lc_122), lc_153=as.numeric(lc_153), 
         Eucmean = replace_na(Eucmean, 0), Eucmean = ifelse(Eucmean<0, 0,Eucmean)) %>% 
  mutate(x_cen=st_coordinates(st_centroid(geom))[,1], y_cen=st_coordinates(st_centroid(geom))[,2])
colnames(grid_10km)[32:65]<-gsub("lc_", "lc", colnames(grid_10km)[32:65])

data %>% colnames()
means<-data_  %>% as_tibble() %>% 
  dplyr::select(23:77) %>%          # selecciona columnas
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
sdss<-data_  %>% as_tibble() %>% 
  dplyr::select(23:77) %>%          # selecciona columnas
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))


library(modEvA) #MESS
summary(Badger_DensityModel) #Select the variables that the model select because are the ones that we're gonna use for projecting the model
data2<-data %>% as.data.frame() %>% dplyr::select("x_cen", "y_cen", "bio_3", "bio_4", "bio_15", "alt", "lc10", "lc30", "lc40", "lc60", "lc70", "lc90",  "lc100", "lc130", "lc190")

#10km
grid_10km<-grid_10km %>% mutate(across(all_of(colnames(means)), ~ . - means[[cur_column()]])) %>% 
                         mutate(across(all_of(colnames(means)), ~ . / sdss[[cur_column()]]))

grid_10km<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km)))
grid10km_df<-grid_10km %>% as.data.frame()  %>% dplyr::select("id_buena", 
                                                              "x_cen", "y_cen", "bio_3", "bio_4", "bio_15", "alt", "lc10", "lc30", "lc40", "lc60", "lc70", "lc90",  "lc100", "lc130", "lc190")
mess_10km<-MESS(V=data2, P=grid10km_df, id.col=1)
grid_10km_<-grid_10km[c(-1:-74)]; names(grid_10km_)
mess_10km_<-mess_10km[c(-2:-16)] ; names(mess_10km_)

grid10km_mess<-grid_10km %>%  full_join(mess_10km_, by="id_buena")

#betas<-fixef(Badger_DensityModel)
#ranef(Badger_DensityModel)
grid10km_mess<-grid10km_mess %>% mutate(country ='Switzerland')

predBadger10km<-grid10km_mess %>% 
  mutate(Pred_model=predict(Badger_DensityModel, grid10km_mess, type="response")) %>% 
  mutate(Pred_species=Pred_model/10000)

st_write(predBadger10km, 
         dsn="predicciones_20260203.gpkg",
         layer="DensityModel_grid10km_Badger_GLMM_xy_noMK_with0_20221220_CH", 
         driver = "GPKG",append=TRUE)
