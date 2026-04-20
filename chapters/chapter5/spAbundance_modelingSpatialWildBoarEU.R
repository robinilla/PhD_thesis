# ---------------------------------------------
#  PhD Thesis: Large-scale monitoring of wild mammal abundance: 
#  modeling frameworks for structured and unstructured data
#  Chapter 5. Accounting for spatial dependence to smooth regional patterns of wild boar abundance across Europe
#  Author: Sonia Illanas
#  Institution: Institute of Game and Wildlife Research
#  Date of last modification: 13/02/2026
# ---------------------------------------------
## R version 4.5.1
## sf version: 1.0-21
## tidyverse version: 2.0.0
## spAbundance version: 0.2.2
## coda version: 0.19-4.1
## MCMCvis version: 0.16.3

library(sf)
library(tidyverse)
library(spAbundance)
library(coda)

# rm(list=ls()) # remove all objects of your global environment

# 1. Load your data
# -----------------------------
data<-st_read("data/dataModeling.gpkg", "wildBoar_data") %>% 
  rename(Eucmean=EucalyptusSpp_mean, hfp=hfp_v3_mean, area=area_km2) %>% 
  mutate(Eucmean=ifelse(is.na(Eucmean), 0, ifelse(Eucmean<0, 0.01, Eucmean)), 
         snow_mean=ifelse(is.na(snow_mean), 0, snow_mean))

# we transform the units of the crs from meters to km because of phi: SO IMPORTANT"!!!
data<-st_transform(data, 
                   st_crs("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=km +no_defs"))

data <- data %>% mutate(Var1=round(st_coordinates(st_point_on_surface(geom))[,1], digits=6),
                        Var2=round(st_coordinates(st_point_on_surface(geom))[,2], digits=6))

colSums(is.na(st_drop_geometry(data))) > 0 # it says which cols has NA values 

# remove any NA values that there is still in the dataset
data<-data[!is.na(data$grow2_mean),]        # Remove NA values from the dataset
data<-data[!is.na(data$elev_mean),]         # Remove NA values from the dataset
data<-data[!is.na(data$hfp),]               # Remove NA values from the dataset
data<-data[!is.infinite(data$hMaxWb),]      # Remove -Inf values from the RV

colnames(data)<-gsub("_mean", "", colnames(data))
colnames(data)[39:74]<-gsub("_", "", colnames(data)[39:74])


# remove duplicated coordinates, that is, there are overlapping polygons
# we need to remove them for making the model. Duplicated locations will
# cause problems (it won't allow to run the model)
data<-data %>% distinct(Var1, Var2, .keep_all = TRUE)
data$country %>% unique() %>% sort()


## 1.1. Select the variables you want to model
# the below variables have been previously selected with a VIF<2 between the
# following ones: select(matches("bio_"), elev, snow, Eucmean, lc10, lc11, lc12, 
# lc30, lc40, lc60, lc70, lc90, lc100, lc120, lc130, lc180, lc190)
data.bm<-data %>% tibble() %>%
                  select(hMaxWbZero,       # response variable
                         locationID2,      #id code
                         area,             # area surface at square kilometers
                         # bioreg,         # bioregion
                         Var1,             # x_centroid
                         Var2,             # y_centroid
                         bio_01:bio_19,    # wordclim variables
                         elev,
                         snow, 
                         Eucmean, 
                         hfp,              # human foot print
                         lc10,             # land cover variables
                         lc11, 
                         lc12, 
                         lc30, 
                         lc40, 
                         lc60, 
                         lc70, 
                         lc90, 
                         lc100, 
                         lc120, 
                         lc130, 
                         lc180, 
                         lc190)

colnames(data.bm)[1]<-"counts" # rename the species' name as counts
                               # "counts" is the response variable

## 1.2. Remove densities >30 indiv. harvested 
# Select the type of dataset you want to analyse (full vs. reduced)
data.bm<-data.bm[!(data.bm$counts/data.bm$area>30),] # for any case: Remove densities >30 indiv. harvested 
# for the full dataset do NOT remove the following comment
# for the reduced dataset (>0 individuals harvested) remove the comment below
# data.bm<-data.bm[data.bm$counts>0,] # remove data which VR is 0.

## 1.3. Standarise your covariates
#       calculate mean and sd values for standarizing covs
covs_st<-list() ; covs_st_values<-list()
for (i in 1:length(colnames(data.bm)[6:41])){
  covariate<-colnames(data.bm)[6:41][i]
  mean<-data.bm[[covariate]] %>% mean(na.rm=TRUE)
  sd<-data.bm[[covariate]] %>% sd(na.rm=TRUE)
  cov_name<-paste(covariate, sep="_")
  covs_st_values[[i]]<-tibble(covariate=cov_name, 
                              mean=mean,
                              sd=sd)
  covs_st[[i]]<-tibble(!!cov_name:=((data.bm[[covariate]]-mean)/sd))
}
# covs_st_values<-do.call(rbind, covs_st_values)
covs_st<-do.call(cbind, covs_st) %>% as_tibble() ; rm(i, covariate, cov_name, mean, sd)

# data with standardized covariates
x <-cbind(data.bm[1:5], covs_st) %>% as_tibble() %>%  #join the data with the standardized covs
  mutate(bio15_2=bio_15^2, 
         elev_2 =elev^2, 
         area_log=log(area))

x$area_log_st<-((x$area_log - mean(x$area_log))/(sd(x$area_log)))


# 2. Make your analysis
#---------------------
# define model parameters
model.family<-'NB' # Negative Binomial

n.batchs <- 10 
batchs.length <- 10000 
n.burns <-  10000 
n.thins <-  10  
n.chains <- 3


# Spatial Random Effects
# -----------------------------
coords<-cbind(x$Var1, x$Var2) 
dist.coords<-dist(coords)

max.dist <- max(dist.coords); max.dist # maximum km from side to side in Europe 
mean.dist<-mean(dist.coords); mean.dist
# We think it is better to be more realistic. We are going to reduce the maximum
# distance to 200 km. So if there were some variable that we are not considering
# in the model but that could have any possible effect in our response variable, 
# we think that neither management decisions nor other variables that are not 
# being considered would over pass that distance. 
max.dist<-200
mean.dist<-70  

data.Zero<-list(y=x$counts,
                covs= x %>% select(bio_02, 
                                   elev, lc40,
                                   Eucmean, lc10, lc30, lc90, lc100, lc130,
                                   area_log_st),
                coords=coords
)

# check again that covariates to be used have a VIF<2
usdm::vifstep(data.frame(data.Zero$covs), th=2)

# MCMC parameters
inits <- list(beta = 0, 
              kappa = 0.5, 
              sigma.sq = 0.5,
              phi = 4/mean.dist,
              w = rep(0, length(data.Zero$y))
)

priors <- list(beta.normal = list(mean = 0, var = 1),
               phi.unif = c(10/max.dist, 2/10),
               sigma.sq.ig = c(2, 1),
               kappa.unif = c(0, 100)
)

tunings <- list(phi = 0.5, kappa = 0.5, beta = 0.1, beta.star = 0.1, w = 0.5)

out_EU<- spAbund(formula = ~  bio_02 + elev + lc40 + Eucmean + lc10 + lc30 + lc90 + lc100 + lc130 + area_log_st,
                 data = data.Zero,
                 family = model.family,
                   
                 inits = inits,
                 priors = priors,
                   
                 n.batch = n.batchs,
                 batch.length = batchs.length,
                 tuning = tunings,
                   
                 cov.model = 'exponential',
                   
                 NNGP = TRUE,
                 n.neighbors = 15, 
                 search.type = 'cb',
                   
                 accept.rate = 0.43,
                 n.omp.threads = 7,
                 n.report = 200,
                 verbose = TRUE,
                 n.burn = n.burns,
                 n.thin = n.thins,
                 n.chains = n.chains)
summary(out_EU)

# save model result
save(out_EU, file="out_EU.Rdata") 



# 3. make your predictions
# -----------------------------
# 3.0. thinning the posterior draws
# -----------------------------
# we got many posterior samples, and we're not able to spread the uncertainty 
# in predictions (not enough RAM memory and R crashes). Therefore, we've to 
# reduce (thin) randomly the posterior we've done it as follows:
# 
# # # define thinning parameter: how many draws we keep
# load(file="out_EU.Rdata") # load model results in case you don't have in your global environment
# 
# uncomment lines 211 - 239 for predictions
# thin_by <- 10   # thin 27 000 to 2 700 samples
#  
# # thin matrices of draws
# thin_matrix <- function(x, thin_by) {
#   if (is.null(x)) return(NULL)
#   keep_idx <- seq(1, nrow(x), by = thin_by)
#   x[keep_idx, , drop = FALSE]
# }
# # 
# # # thin arrays
# thin_array <- function(x, thin_by) {
#   if (is.null(x)) return(NULL)
#   keep_idx <- seq(1, dim(x)[1], by = thin_by)
#   x[keep_idx,,, drop = FALSE]
# }
# # 
# # # apply the above functions to spAbund object
# # # adjust depending on the parameters of your models
# # # not all models have the same fields 
# # # be aware that some parameters are matrices and others arrays 
# # # Use the function that correspond in each case
# out_thin<-out_EU
# out_thin$beta.samples <- thin_matrix(out_thin$beta.samples, thin_by)
# out_thin$w.samples    <- thin_matrix(out_thin$w.samples, thin_by)
# out_thin$alpha.samples <- thin_matrix(out_thin$theta.samples, thin_by)
# out_thin$mu.samples   <- thin_array(out_thin$mu.samples, thin_by)
# out_thin$kappa.samples   <- thin_matrix(out_thin$kappa.samples, thin_by)
# out_thin$n.thin<-10
# out_thin$n.post<-(12000-3000)/10
# summary(out_thin)
# -----------------------------


# -----------------------------
# 3.2 grid: 10 km x 10 km
# -----------------------------
# repeat the same steps as for 2 km x 2 km grid predictions
grid10km<-st_read("grid_10km_wb.shp")
grid10km<-st_transform(grid10km, 
                       st_crs("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=km +no_defs"))

colnames(grid10km)[32:65]<-gsub("_", "", colnames(grid10km)[32:65])
colnames(grid10km)[68]<-"area"
colnames(grid10km)[30]<-"Eucmean"
colnames(grid10km)[10]<-"bio_02"
colnames(grid10km)[29]<-"elev"

grid10km<- grid10km %>% 
            dplyr::select(bio_02, elev, Eucmean, lc10, lc30, lc40, lc90, lc100, lc130, area, Bioregion) %>% 
            mutate(Eucmean = replace_na(Eucmean, 0),
                   Eucmean = ifelse(Eucmean<0, 0, Eucmean),
                   Var1=st_coordinates(st_centroid(geometry))[,1], Var2=st_coordinates(st_centroid(geometry))[,2], 
                   id_buena=row_number())
summary(grid10km)

# 3.2.1 scale your variables with the values of the data used for modeling
grid10km_scale<-grid10km %>% mutate(
  bio_02=((bio_02-covs_st_values[2,]$mean)/covs_st_values[2,]$sd), 
  elev = ((elev-covs_st_values[20,]$mean)/covs_st_values[20,]$sd), 
  Eucmean=((Eucmean-covs_st_values[22,]$mean)/covs_st_values[22,]$sd), 
  lc10 = ((lc10/100-covs_st_values[24,]$mean)/covs_st_values[24,]$sd),
  lc30 = ((lc30/100-covs_st_values[27,]$mean)/covs_st_values[27,]$sd), 
  lc40 = ((lc40/100-covs_st_values[28,]$mean)/covs_st_values[28,]$sd), 
  lc90 = ((lc90/100-covs_st_values[31,]$mean)/covs_st_values[31,]$sd), 
  lc100 =((lc100/100-covs_st_values[32,]$mean)/covs_st_values[32,]$sd), 
  lc130 =((lc130/100-covs_st_values[34,]$mean)/covs_st_values[34,]$sd), 
  area   = area, 
  area_log=log(area),
  Var1 = Var1,  # x_centroid in km
  Var2 = Var2,   # y_centroid in km
  area_log_st = (area_log - mean(x$area_log))/(sd(x$area_log))
)

gc()

# 3.2.2. remove those cells which covariates have NA values
# not needed for our 10 x 10 km prediction layer
# summary(grid10km_scale)

# 3.2.3. format variables for predicting
coords.grid<-cbind(grid10km_scale$Var1, grid10km_scale$Var2)
coords.0<-as.matrix(coords.grid); rm(coords.grid)

X.0 <- cbind(1, grid10km_scale$bio_02, grid10km_scale$elev, grid10km_scale$lc40, grid10km_scale$Eucmean, grid10km_scale$lc10,
             grid10km_scale$lc30, grid10km_scale$lc90, grid10km_scale$lc100, grid10km_scale$lc130,
             grid10km_scale$area_log_st)
colnames(X.0) <- c('(Intercept)',  "bio_02", 'elev', 'lc40', 
                   'Eucmean', 'lc10', 'lc30', 'lc90', 'lc100', 'lc130',
                   'area_log_st')

# 3.2.4 make the predictions 
#       with the thinned draws of the posterior
out.m4b.thin <- predict(out4b_thin, X.0, coords.0, ignore.RE = FALSE, n.omp.threads = 7) #change it according to your computer's capacity
mu.0.quants.m4b.grid <- apply(out.m4b.thin$mu.0.samples, 2, quantile, c(0.025, 0.5, 0.975), na.rm=T)

grid10km_scale<-grid10km_scale %>% 
                    mutate("spA.median" = mu.0.quants.m4b.grid[2,], 
                           "spA.lowCI" = mu.0.quants.m4b.grid[1,],
                           "spA.highCI" = mu.0.quants.m4b.grid[3,]
                    )

# 3.2.5 make environmental only predictions (without spaial random effect)
library(MCMCvis)
summ<-MCMCsummary(out4b_EU$beta.samples, round=4) %>% rownames_to_column()
grid10km_scale<-
  grid10km_scale %>% 
  mutate("spA.median.e"=
           exp(summ[1,]$mean +
                 summ[2,]$mean * bio_02 +
                 summ[3,]$mean * elev + 
                 summ[4,]$mean * lc40 +
                 summ[5,]$mean * Eucmean +
                 summ[6,]$mean * lc10 +
                 summ[7,]$mean * lc30 +
                 summ[8,]$mean * lc90 +
                 summ[9,]$mean * lc100 +
                 summ[10,]$mean * lc130 +
                 summ[11,]$mean * area_log_st), 
         
         "spA.lowCI.e"=
           exp(summ[1,]$`2.5%` +
                 summ[2,]$`2.5%` * bio_02 +
                 summ[3,]$`2.5%` * elev + 
                 summ[4,]$`2.5%` * lc40 +
                 summ[5,]$`2.5%` * Eucmean +
                 summ[6,]$`2.5%` * lc10 +
                 summ[7,]$`2.5%` * lc30 +
                 summ[8,]$`2.5%` * lc90 +
                 summ[9,]$`2.5%` * lc100 +
                 summ[10,]$`2.5%` * lc130 +
                 summ[11,]$`2.5%` * area_log_st), 
         
         "spA.highCI.e"=
           exp(summ[1,]$`97.5%` +
                 summ[2,]$`97.5%` * bio_02 +
                 summ[3,]$`97.5%` * elev + 
                 summ[4,]$`97.5%` * lc40 +
                 summ[5,]$`97.5%` * Eucmean +
                 summ[6,]$`97.5%` * lc10 +
                 summ[7,]$`97.5%` * lc30 +
                 summ[8,]$`97.5%` * lc90 +
                 summ[9,]$`97.5%` * lc100 +
                 summ[10,]$`97.5%` * lc130 +
                 summ[11,]$`97.5%` * area_log_st)
  )

# 3.2.6 save results 
st_write(grid10km_scale, 
         dsn = "outSpeciesPredictions_2025December.gpkg", 
         layer="out_EU_grid10km_environmental")
