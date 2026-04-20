library(tidyverse)
library(sf)
# rm(list=ls())
setwd("/home/vant/Documentos/PhD/HB_chapter/")

library(pscl)
stepwise_zeroinfl_alternate <- function(data, response,
                                        count_fixed = character(0),
                                        count_candidates,
                                        zi_fixed = character(0),
                                        zi_candidates,
                                        start_in = c("zi","count"),
                                        dist = "negbin",
                                        criterion = c("AIC","BIC"),
                                        min_improve = 0,   # usa 2 para ΔAIC>2
                                        max_steps = 50,
                                        trace = FALSE) {
  
  start_in  <- match.arg(start_in)
  criterion <- match.arg(criterion)
  
  # ---------- helpers ----------
  rhs <- function(fixed, vars) {
    v <- c(fixed, vars)
    if (length(v) == 0) "1" else paste(v, collapse = " + ")
  }
  
  make_formula <- function(count_vars, zi_vars) {
    as.formula(
      paste(response, "~",
            rhs(count_fixed, count_vars),
            "|",
            rhs(zi_fixed, zi_vars))
    )
  }
  
  fit_one <- function(form) {
    warn <- character(0)
    ok <- TRUE
    m <- tryCatch(
      withCallingHandlers(
        zeroinfl(form, data = data, dist = dist, trace = trace),
        warning = function(w) {
          warn <<- c(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) { ok <<- FALSE; NULL }
    )
    
    list(
      model     = m,
      ok        = ok,
      converged = if (ok) isTRUE(m$converged) else FALSE,
      AIC       = if (ok) AIC(m) else NA_real_,
      BIC       = if (ok) BIC(m) else NA_real_,
      logLik    = if (ok) as.numeric(logLik(m)) else NA_real_,
      warnings  = paste(unique(warn), collapse = " | ")
    )
  }
  
  score <- function(fit)
    if (criterion == "AIC") fit$AIC else fit$BIC
  
  # ---------- modelo inicial ----------
  current_count <- character(0)
  current_zi    <- character(0)
  
  current_form <- make_formula(current_count, current_zi)
  current_fit  <- fit_one(current_form)
  
  if (!current_fit$ok || !current_fit$converged || !is.finite(score(current_fit)))
    stop("El modelo inicial no es válido.")
  
  history <- list()
  all_tried <- list()
  
  mode  <- start_in
  steps <- 0
  
  repeat {
    
    steps <- steps + 1
    if (steps > max_steps) break
    
    # ---------- generar candidatos ----------
    if (mode == "zi") {
      pool <- setdiff(zi_candidates, current_zi)
      if (length(pool) == 0) break
      
      candidates <- lapply(pool, function(v) {
        zi_try <- c(current_zi, v)
        form   <- make_formula(current_count, zi_try)
        fit    <- fit_one(form)
        list(add = v, part = "zi",
             count = current_count, zi = zi_try,
             formula = form, fit = fit)
      })
      
    } else { # count
      pool <- setdiff(count_candidates, current_count)
      if (length(pool) == 0) break
      
      candidates <- lapply(pool, function(v) {
        count_try <- c(current_count, v)
        form      <- make_formula(count_try, current_zi)
        fit       <- fit_one(form)
        list(add = v, part = "count",
             count = count_try, zi = current_zi,
             formula = form, fit = fit)
      })
    }
    
    # ---------- tabla del step ----------
    step_df <- do.call(rbind, lapply(candidates, function(x) {
      data.frame(
        step = steps,
        part = x$part,
        added = x$add,
        count_terms = if (length(x$count)==0) "(none)" else paste(x$count, collapse=" + "),
        zi_terms    = if (length(x$zi)==0) "(none)" else paste(x$zi, collapse=" + "),
        converged = x$fit$converged,
        AIC = x$fit$AIC,
        BIC = x$fit$BIC,
        logLik = x$fit$logLik,
        warnings = x$fit$warnings,
        stringsAsFactors = FALSE
      )
    }))
    
    # ---------- FILTRADO CRÍTICO ----------
    valid <- step_df$converged & is.finite(step_df[[criterion]])
    
    if (!any(valid)) {
      # no hay modelos válidos → parar
      break
    }
    
    step_df_valid <- step_df[valid, ]
    step_df_valid <- step_df_valid[order(step_df_valid[[criterion]]), ]
    
    best_row <- step_df_valid[1, , drop = FALSE]
    
    # recuperar el fit correcto
    idx <- which(sapply(candidates, function(x)
      x$add == best_row$added && x$part == best_row$part))[1]
    
    best_fit <- candidates[[idx]]$fit
    
    # ---------- comprobar mejora ----------
    current_score <- score(current_fit)
    best_score    <- score(best_fit)
    
    improve <- current_score - best_score
    
    step_df$winner <- FALSE
    step_df$winner[valid][which.min(step_df_valid[[criterion]])] <- TRUE
    
    history[[steps]] <- list(
      mode = mode,
      current_formula = current_form,
      current_score = current_score,
      step_table = step_df,
      improve = improve
    )
    
    all_tried[[steps]] <- candidates
    
    if (!is.finite(best_score) || improve < min_improve) {
      # no mejora → parar
      break
    }
    
    # ---------- aceptar ganador ----------
    if (mode == "zi") {
      current_zi <- c(current_zi, best_row$added)
    } else {
      current_count <- c(current_count, best_row$added)
    }
    
    current_form <- make_formula(current_count, current_zi)
    current_fit  <- best_fit
    
    # alternar componente
    mode <- if (mode == "zi") "count" else "zi"
  }
  
  # ---------- salida ----------
  all_steps_df <- do.call(rbind, lapply(history, `[[`, "step_table"))
  
  list(
    best_model = current_fit$model,
    best_formula = current_form,
    best_count_terms = c(count_fixed, current_count),
    best_zi_terms = c(zi_fixed, current_zi),
    history = history,
    tried = all_tried,
    table = all_steps_df
  )
}

library(Hmisc)
cal_data <- function(pred, obs, g = 9){
  bins <-cut2(pred, g = g) # defining bins (percentiles) on the predicted
  tibble(
    pred_mean = tapply(pred, bins, mean), # calculating mean values for predicted
    obs_mean = tapply(obs, bins, mean), # calculating mean values for observed
    n = tapply (obs, bins, length) )
}

#Charge constrained species  data
# Alpine ibex #
# --------------------
# rm(list=ls())
data_sp<-st_read(dsn="scripts/data/predicciones_20220811.gpkg",
                 layer='Ibex_spatial_v5')

sp<-data_sp %>%
  mutate(lc10=lc_10_sum/lc_10_count*100, lc11=lc_11_sum/lc_11_count*100, lc12=lc_12_sum/lc_12_count*100, lc20=lc_20_sum/lc_20_count*100, lc30=lc_30_sum/lc_30_count*100, lc40=lc_40_sum/lc_40_count*100, lc50=lc_50_sum/lc_50_count*100, lc60=lc_60_sum/lc_60_count*100, lc61=lc_61_sum/lc_61_count*100, lc62=lc_62_sum/lc_62_count*100, lc70=lc_70_sum/lc_70_count*100, lc71=lc_71_sum/lc_71_count*100, lc72=lc_72_sum/lc_72_count*100, lc80=lc_80_sum/lc_80_count*100, lc81=lc_81_sum/lc_81_count*100, lc82=lc_82_sum/lc_82_count*100, lc90=lc_90_sum/lc_90_count*100, lc100=lc_100_sum/lc_100_count*100, lc110=lc_110_sum/lc_110_count*100,lc120=lc_120_sum/lc_120_count*100,lc121=lc_121_sum/lc_121_count*100,lc122=lc_122_sum/lc_122_count*100,lc130=lc_130_sum/lc_130_count*100,lc140=lc_140_sum/lc_140_count*100,lc150=lc_150_sum/lc_150_count*100,lc152=lc_152_sum/lc_152_count*100,lc153=lc_153_sum/lc_153_count*100,lc160=lc_160_sum/lc_160_count*100,lc170=lc_170_sum/lc_170_count*100,lc180=lc_180_sum/lc_180_count*100,lc190=lc_190_sum/lc_190_count*100,lc200=lc_200_sum/lc_200_count*100,lc201=lc_201_sum/lc_201_count*100,lc202=lc_202_sum/lc_202_count*100,lc210=lc_210_sum/lc_210_count*100,lc220=lc_220_sum/lc_220_count*100) %>%
  dplyr::select(!matches("_cou") & !matches("_sum")) %>%
  #quito variables que no vaya a usar en el modelo
  dplyr::select(-c(excel_file,excel_path, shp_file, shp_path, #species, #me cargo variables que no necesito para el modelo
                   grow_mean, #la variable grow tiene datos raros, as? que la voy a eliminar de nuestro set de datos
                   lc50, lc62, lc72, lc81, lc82,lc121, lc170,#, #Me cargo variables lc que en el summary son todo 0
                   #BR_todasr1
                   "individualCount...39", "individualCount...51"
  )) %>% #quito las bioregiones porque no las vamos a usar para el modelo
  filter(dens<50) %>%  #tres poligonos de >Spain con areas muy muy muy muy small
  # filter(dens>0 & dens<50)
  select(locationID, locality, country, NUTS, BR_todasr1, area_km2, matches("harv"), dens, dens_r_trans,
         matches("bio"), alt, snow, sun, Eucmean, hfp,
         matches("lc"),
         x_cen, y_cen)

sp$Eucmean<-ifelse(sp$Eucmean<0, 0.01, sp$Eucmean) #Como para wild boar transformamos la variable de Eucmean para eliminar los valores negativos por valores de0 muy peque?os
sp<-sp[!is.na(sp$dens),] #Remove NA values from RV: it does not remove any value

#Modeling phase
#1) data
colnames(sp)[16:34]<-gsub("_mean", "", colnames(sp)[16:34])
vars_std <- c(colnames(sp %>% as_tibble())[16:70])

#remove NA values from the dataset
sp<-sp[!is.na(sp$bio_01),] #Remove NA values from the spset
sp<-sp[!is.na(sp$snow),] #Remove NA values from the spset
sp<-sp[!is.na(sp$hfp),] #Remove NA values from the spset
sp<-sp[!is.na(sp$alt),] #Remove NA values from the spset
sp<-sp[!is.na(sp$lc10),] #Remove NA values from the spset

data<-sp  %>% mutate(across(all_of(vars_std),~ as.numeric(scale(.))))
summary(data)

#1.1 Creamos data set de calibrado y de validacion (80% / 20%)
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]; st_geometry(data_tra)<-NULL; st_geometry(data_val)<-NULL

#2.0 VIF <2
#2.1 seleccionamos variables para el vif, para eso miramos como es la distribucion de las variables
# hist(data_tra$Eucmean)
# hist(data_tra$alt)
# hist(data_tra$snow)
#
# hist(data_tra$lc10)
# hist(data_tra$lc11)
# hist(data_tra$lc12) #
# hist(data_tra$lc20) #
# hist(data_tra$lc30)
# hist(data_tra$lc40)
# hist(data_tra$lc60)
# hist(data_tra$lc70)
# hist(data_tra$lc90)
# hist(data_tra$lc100)
# hist(data_tra$lc120)
# hist(data_tra$lc130)
#
# hist(data_tra$lc61) #
# hist(data_tra$lc71) #
# hist(data_tra$lc80) #
# hist(data_tra$lc110) #
# hist(data_tra$lc122) #
# hist(data_tra$lc140) #
# hist(data_tra$lc150) #
# hist(data_tra$lc152) #
# hist(data_tra$lc153) #
# hist(data_tra$lc160) #
# hist(data_tra$lc180) #
# hist(data_tra$lc190) #
# hist(data_tra$lc200) #
# hist(data_tra$lc201) #
# hist(data_tra$lc202) #
# hist(data_tra$lc210) #
# hist(data_tra$lc220) #

#2.1 hacemos el vif
variables<-data_tra %>% dplyr::select(matches("bio_"), alt, snow, Eucmean,lc10, lc11, lc12, lc30, lc40, lc60, lc70, lc90,
                                      lc100,
                                      lc130, x_cen, y_cen)
resultado.vif0<-usdm::vifstep(variables, th=2); resultado.vif0


#2.2) Density Models: binomial negativa necesita que sea un numero entero, por eso se ha transformado la variable densidad *10000
#ZINB
count_vars <- c("bio_02","bio_03", "bio_15","lc10","lc30","lc40","lc60","lc70","lc100")
zi_vars <- count_vars

res <- stepwise_zeroinfl_alternate(
  data = data_tra,
  response = "dens_r_trans",
  count_fixed = c("x_cen","y_cen"),              # SIEMPRE dentro del count
  count_candidates = count_vars,
  zi_fixed = c("x_cen","y_cen"),                 # o c("x_cen","y_cen") si también los quieres fijos en ZI
  zi_candidates = zi_vars,
  start_in = "zi",                               # como tu gráfico: primero inflado, luego conteo
  dist = "negbin",
  criterion = "AIC",
  min_improve = 0                                # pon 2 si quieres “solo si mejora >2”
)

res$best_formula
summary(res$best_model)

AlpineIbex_DensityModel<-res$best_model
# save(AlpineIbex_DensityModel, file="scripts/ModelResults/densityStepModel_Alpineibex_Europe.Rdata")

#2.3) Calibration plots  validation data
data_val<-data_val %>%
  mutate(Pred_model = as.integer(predict(AlpineIbex_DensityModel, data_val,type = "response"))) %>% #predecir zeroinfl
  #predecir hurdle
  # mutate(p_hat = predict(mh, data_val, type = "zero"),
  #        mu_hat=predict(mh, data_val, type="count"),
  #        # y_hat=p_pos*mu_pos) %>%
  #        y_hat = predict(mh, data_val, type="response")) %>%
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

data_pos<-data_val %>% filter(dens>0)
cal_tot<-bind_rows(
  cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
)

ggplot(cal_tot, aes(x = pred_mean, y = obs_mean)) +
  geom_point(aes(size = n), shape = 21, fill = "grey70") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
  facet_wrap(~Type, scales = "free") +
  scale_size_continuous(name = "N per bin", range = c(2, 6)) +
  labs(
    x = "Predicted mean",
    y = "Observed mean"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )
# s_class <- cut2(data_val$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean <- as.matrix(tapply(data_val$Pred_Red, s_class, mean)) # calculating mean values for predicted
# y_mean <- as.matrix(tapply(data_val$dens, s_class, mean)) # calculating mean values for observed
# # corbr1<-cor.test(data_val$Pred_Red, data_val$dens)
# data_val_w0<-data_val %>%filter(dens_r_trans>0)
# s_class_w0 <- cut2(data_val_w0$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean_w0 <- as.matrix(tapply(data_val_w0$Pred_Red, s_class_w0, mean)) # calculating mean values for predicted
# y_mean_w0 <- as.matrix(tapply(data_val_w0$dens, s_class_w0, mean)) # calculating mean values for observed

# cpE<-rbind(data.frame(s_mean=s_mean,   y_mean=y_mean, Type=as.factor("All data")),
#            data.frame(s_mean=s_mean_w0, y_mean=y_mean_w0, Type=as.factor("Data without 0")))
# ggplot(cpE, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=4)+ scale_colour_grey()+ scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="red", lty=2, cex=0.5)+
#   theme_bw()+
#   # xlim(0,max(s_mean_w0, max(na.omit(y_mean_w0))))+
#   # ylim(0,max(max(na.omit(y_mean_w0)), s_mean_w0))+
#   ylab("Observed Alpine ibex hunted bag density")+xlab("Predicted")+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# #3.0) Spatial layer: 10x10km grid
grid_10km<-st_read(dsn="scripts/data/data_constrained_sp_20220902.gpkg", layer='grid_AlpineIbex') %>% dplyr::select(-Bioregion) %>%
  rename(Bioregion=BR_T10, NUTS=NUT, Eucmean=Eucalyptus) %>%
  mutate(NUTS=as.factor(NUTS),
         Bioregion=as.factor(Bioregion),
         lc_50=as.numeric(lc_50),
         lc_62=as.numeric(lc_62),
         lc_72=as.numeric(lc_72),
         lc_122=as.numeric(lc_122),
         lc_153=as.numeric(lc_153),
         Eucmean = replace_na(Eucmean, 0))
grid_10km<-grid_10km %>%
  mutate(Eucmean = ifelse(grid_10km$Eucmean<0, 0, Eucmean)) %>%
  mutate(x_cen=st_coordinates(st_centroid(grid_10km))[,1], y_cen=st_coordinates(st_centroid(grid_10km))[,2])
colnames(grid_10km)[32:65]<-gsub("lc_", "lc", colnames(grid_10km)[32:65])
colnames(grid_10km)[9:17]<-gsub("_", "_0", colnames(grid_10km)[9:17])

data %>% colnames()
means<-sp %>% as_tibble() %>%
  select(16:70) %>%          # selecciona columnas
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
sdss<-sp %>% as_tibble() %>%
  select(16:70) %>%          # selecciona columnas
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

library(modEvA) #MESS
summary(AlpineIbex_DensityModel) #Solo ponemos las variables que selecciona el modelo que son sobre las que vamos a hacer la prediccion
data_<-data %>% as.data.frame() %>% dplyr::select(c(x_cen, y_cen, lc30))

# #10km
grid_10km<-grid_10km %>% mutate(across(all_of(colnames(means)), ~ . - means[[cur_column()]])) %>%
  mutate(across(all_of(colnames(means)), ~ . / sdss[[cur_column()]]))

grid_10km<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km)))
grid10km_df<-grid_10km %>% as.data.frame()  %>% dplyr::select("id_buena",
                                                              c(x_cen, y_cen, lc30))

mess_10km<-MESS(V=data_, P=grid10km_df, id.col=1)

grid_10km_<-grid_10km[c(-1:-73)]
mess_10km_<-mess_10km[c(-2:-4)]; names(mess_10km_)

grid10km_mess<-grid_10km %>%  full_join(mess_10km_, by="id_buena")

predRedDeer10km<-grid10km_mess %>% mutate(Pred_model = as.integer(predict(AlpineIbex_DensityModel, grid_10km, type = "response"))) %>% mutate(Pred_Red = Pred_model/10000)
st_write(predRedDeer10km,
         dsn="scripts/ModelResults/predictions/predicciones_20260203.gpkg",
         layer="ZINB_stepModel_grid10km_AlpineIbex", driver = "GPKG",append=TRUE)
# # --------------------


# Aoudad #
# --------------------
data_sp<-st_read(dsn="scripts/data/predicciones_20220811.gpkg",
                 layer='Aoudad_spatial_v5')

sp<-data_sp %>%
  mutate(lc10=lc_10_sum/lc_10_count*100, lc11=lc_11_sum/lc_11_count*100, lc12=lc_12_sum/lc_12_count*100, lc20=lc_20_sum/lc_20_count*100, lc30=lc_30_sum/lc_30_count*100, lc40=lc_40_sum/lc_40_count*100, lc50=lc_50_sum/lc_50_count*100, lc60=lc_60_sum/lc_60_count*100, lc61=lc_61_sum/lc_61_count*100, lc62=lc_62_sum/lc_62_count*100, lc70=lc_70_sum/lc_70_count*100, lc71=lc_71_sum/lc_71_count*100, lc72=lc_72_sum/lc_72_count*100, lc80=lc_80_sum/lc_80_count*100, lc81=lc_81_sum/lc_81_count*100, lc82=lc_82_sum/lc_82_count*100, lc90=lc_90_sum/lc_90_count*100, lc100=lc_100_sum/lc_100_count*100, lc110=lc_110_sum/lc_110_count*100,lc120=lc_120_sum/lc_120_count*100,lc121=lc_121_sum/lc_121_count*100,lc122=lc_122_sum/lc_122_count*100,lc130=lc_130_sum/lc_130_count*100,lc140=lc_140_sum/lc_140_count*100,lc150=lc_150_sum/lc_150_count*100,lc152=lc_152_sum/lc_152_count*100,lc153=lc_153_sum/lc_153_count*100,lc160=lc_160_sum/lc_160_count*100,lc170=lc_170_sum/lc_170_count*100,lc180=lc_180_sum/lc_180_count*100,lc190=lc_190_sum/lc_190_count*100,lc200=lc_200_sum/lc_200_count*100,lc201=lc_201_sum/lc_201_count*100,lc202=lc_202_sum/lc_202_count*100,lc210=lc_210_sum/lc_210_count*100,lc220=lc_220_sum/lc_220_count*100) %>%
  dplyr::select(!matches("_cou") & !matches("_sum")) %>%
  mutate(Eucmean = replace_na(Eucmean, 0), 
         Eucmean = ifelse(Eucmean<0, 0.01, Eucmean)) %>% #Como para wild boar transformamos la variable de Eucmean para eliminar los valores negativos por valores de0 muy peque?os
  #quito variables que no vaya a usar en el modelo
  dplyr::select(-c(excel_file,excel_path, shp_file, shp_path, #species, #me cargo variables que no necesito para el modelo
                   grow_mean, #la variable grow tiene datos raros, as? que la voy a eliminar de nuestro set de datos
                   lc50, lc62, lc72, lc81, lc82,lc121, lc170,#, #Me cargo variables lc que en el summary son todo 0
                   #BR_todasr1
                   "individualCount...39", "individualCount...51"
  )) %>% #quito las bioregiones porque no las vamos a usar para el modelo
  filter(dens<50) %>%  #tres poligonos de >Spain con areas muy muy muy muy small
  # filter(dens>0 & dens<50)
  select(locationID, locality, country, NUTS, BR_todasr1, area_km2, matches("harv"), dens, dens_r_trans,
         matches("bio"), alt, snow, sun, Eucmean, hfp,
         matches("lc"),
         x_cen, y_cen)
sp<-sp[!is.na(sp$dens),] #Remove NA values from RV: it does not remove any value

#Modeling phase
#1) data
colnames(sp)[15:33]<-gsub("_mean", "", colnames(sp)[15:33])
vars_std <- c(colnames(sp %>% as_tibble())[15:69])

#remove NA values from the dataset
sp<-sp[!is.na(sp$bio_01),] #Remove NA values from the spset
sp<-sp[!is.na(sp$snow),] #Remove NA values from the spset
sp<-sp[!is.na(sp$hfp),] #Remove NA values from the spset
sp<-sp[!is.na(sp$alt),] #Remove NA values from the spset
sp<-sp[!is.na(sp$lc10),] #Remove NA values from the spset

data<-sp  %>% mutate(across(all_of(vars_std),~ as.numeric(scale(.))))
summary(data)

#1.1 Creamos data set de calibrado y de validacion (80% / 20%)
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]; st_geometry(data_tra)<-NULL; st_geometry(data_val)<-NULL

#2.0 VIF <2
#2.1 seleccionamos variables para el vif, para eso miramos como es la distribucion de las variables
# hist(data_tra$Eucmean)
# hist(data_tra$alt)
# hist(data_tra$snow)
#
# hist(data_tra$lc10)
# hist(data_tra$lc11)
# hist(data_tra$lc12) #
# hist(data_tra$lc20) #
# hist(data_tra$lc30)
# hist(data_tra$lc40)
# hist(data_tra$lc60)
# hist(data_tra$lc70)
# hist(data_tra$lc90)
# hist(data_tra$lc100)
# hist(data_tra$lc120)
# hist(data_tra$lc130)
#
# hist(data_tra$lc61) #
# hist(data_tra$lc71) #
# hist(data_tra$lc80) #
# hist(data_tra$lc110) #
# hist(data_tra$lc122) #
# hist(data_tra$lc140) #
# hist(data_tra$lc150) #
# hist(data_tra$lc152) #
# hist(data_tra$lc153) #
# hist(data_tra$lc160) #
# hist(data_tra$lc180) #
# hist(data_tra$lc190) #
# hist(data_tra$lc200) #
# hist(data_tra$lc201) #
# hist(data_tra$lc202) #
# hist(data_tra$lc210) #
# hist(data_tra$lc220) #

#2.1 hacemos el vif
variables<-data_tra %>% dplyr::select(matches("bio_"), alt, snow, Eucmean,lc10, lc11, lc12, lc30, lc40, lc60, lc70, lc90,
                                      lc100, #lc120,
                                      lc130, x_cen, y_cen)
resultado.vif0<-usdm::vifstep(variables, th=2); resultado.vif0


#2.2) Density Models: binomial negativa necesita que sea un numero entero, por eso se ha transformado la variable densidad *10000
#ZINB
# 1) ajustamos modelo con ZI solo intercepto  + Count variables
count_vars <- c("bio_02","snow","lc10","lc30","lc40","lc60","lc70","lc100","lc130","Eucmean")
zi_vars <- count_vars
res <- stepwise_zeroinfl_alternate(
  data = data_tra,
  response = "dens_r_trans",
  count_fixed = c("x_cen","y_cen"),              # SIEMPRE dentro del count
  count_candidates = count_vars,
  zi_fixed = c("x_cen","y_cen"),                 # o c("x_cen","y_cen") si también los quieres fijos en ZI
  zi_candidates = zi_vars,
  start_in = "zi",                               # como tu gráfico: primero inflado, luego conteo
  dist = "negbin",
  criterion = "AIC",
  min_improve = 0                                # pon 2 si quieres “solo si mejora >2”
)

res$best_model
Aoudad_DensityModel<-res$best_model
# save("Aoudad_DensityModel", file="scripts/ModelResults/densityStepModel_Aoudad_Europe.Rdata")

#2.3) Calibration plots  validation data
data_val<-data_val %>%
  mutate(Pred_model = as.integer(predict(Aoudad_DensityModel, data_val,type = "response"))) %>% #predecir zeroinfl
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

# calibration plot
data_pos<-data_val %>% filter(dens>0)
cal_tot<-bind_rows(
  cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
)

ggplot(cal_tot, aes(x = pred_mean, y = obs_mean)) +
  geom_point(aes(size = n), shape = 21, fill = "grey70") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
  facet_wrap(~Type, scales = "free") +
  scale_size_continuous(name = "N per bin", range = c(2, 6)) +
  labs(
    x = "Predicted mean",
    y = "Observed mean"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

# s_class <- cut2(data_val$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean <- as.matrix(tapply(data_val$Pred_Red, s_class, mean)) # calculating mean values for predicted
# y_mean <- as.matrix(tapply(data_val$dens, s_class, mean)) # calculating mean values for observed
# # corbr1<-cor.test(data_val$Pred_Red, data_val$dens)
# data_val_w0<-data_val %>%filter(dens_r_trans>0)
# s_class_w0 <- cut2(data_val_w0$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean_w0 <- as.matrix(tapply(data_val_w0$Pred_Red, s_class_w0, mean)) # calculating mean values for predicted
# y_mean_w0 <- as.matrix(tapply(data_val_w0$dens, s_class_w0, mean)) # calculating mean values for observed
#
# cpE<-rbind(data.frame(s_mean=s_mean,   y_mean=y_mean, Type=as.factor("All data")),
#            data.frame(s_mean=s_mean_w0, y_mean=y_mean_w0, Type=as.factor("Data without 0")))
# ggplot(cpE, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=4)+ scale_colour_grey()+ scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="red", lty=2, cex=0.5)+
#   theme_bw()+
#   # xlim(0,max(s_mean_w0, max(na.omit(y_mean_w0))))+
#   # ylim(0,max(max(na.omit(y_mean_w0)), s_mean_w0))+
#   ylab("Observed Alpine ibex hunted bag density")+xlab("Predicted")+
#   theme(plot.title = element_text(hjust = 0.5))


#3.0) Spatial layer: 10x10km grid
grid_10km<-st_read(dsn="scripts/data/data_constrained_sp_20220902.gpkg",
                   layer='grid_Aoudad') %>% dplyr::select(-Bioregion) %>%
  rename(Bioregion=BR_T10, NUTS=NUT, Eucmean=Eucalyptus) %>%
  mutate(NUTS=as.factor(NUTS),
         Bioregion=as.factor(Bioregion),
         lc_50=as.numeric(lc_50),
         lc_62=as.numeric(lc_62),
         lc_72=as.numeric(lc_72),
         lc_122=as.numeric(lc_122),
         lc_153=as.numeric(lc_153),
         Eucmean = replace_na(Eucmean, 0))
grid_10km<-grid_10km %>%
  mutate(Eucmean = ifelse(grid_10km$Eucmean<0, 0, Eucmean)) %>%
  mutate(x_cen=st_coordinates(st_centroid(grid_10km))[,1], y_cen=st_coordinates(st_centroid(grid_10km))[,2])
colnames(grid_10km)[32:65]<-gsub("lc_", "lc", colnames(grid_10km)[32:65])
colnames(grid_10km)[9:17]<-gsub("_", "_0", colnames(grid_10km)[9:17])

data %>% colnames()
means<-sp %>% as_tibble() %>%
  select(15:69) %>%          # selecciona columnas
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
sdss<-sp %>% as_tibble() %>%
  select(15:69) %>%          # selecciona columnas
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

library(modEvA) #MESS
summary(Aoudad_DensityModel) #Solo ponemos las variables que selecciona el modelo que son sobre las que vamos a hacer la prediccion
data_<-data %>% as.data.frame() %>% dplyr::select(c(x_cen, y_cen, lc30))

#10km
# model_vars<-c(x_cen, y_cen, bio_03, bio_10, bio_15, lc10, lc30, lc40, lc60, lc70, lc90, lc100, lc130, Eucmean)
grid_10km<-grid_10km %>% mutate(across(all_of(colnames(means)), ~ . - means[[cur_column()]])) %>%
  mutate(across(all_of(colnames(means)), ~ . / sdss[[cur_column()]]))

grid_10km<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km)))
grid10km_df<-grid_10km %>% as.data.frame()  %>% dplyr::select("id_buena",
                                                              c(x_cen, y_cen, lc30))
mess_10km<-MESS(V=data_, P=grid10km_df, id.col=1)

grid_10km_<-grid_10km[c(-1:-73)]
mess_10km_<-mess_10km[c(-2:-4)]

grid10km_mess<-grid_10km %>%  full_join(mess_10km_, by="id_buena")

predRedDeer10km<-grid10km_mess %>% mutate(Pred_model = as.integer(predict(Aoudad_DensityModel, grid_10km, type = "response"))) %>% mutate(Pred_Red = Pred_model/10000)
st_write(predRedDeer10km,
         dsn="scripts/ModelResults/predictions/predicciones_20260203.gpkg",
         layer="ZINB_stepModel_grid10km_Aoudad", driver = "GPKG",append=TRUE)
# # --------------------

# Iberian wild goat #
# --------------------
# rm(list=ls())
data_sp<-st_read(dsn="scripts/data/predicciones_20220811.gpkg", 
                 layer='Goat_spatial_v5')

sp<-data_sp %>%
  mutate(lc10=lc_10_sum/lc_10_count*100, lc11=lc_11_sum/lc_11_count*100, lc12=lc_12_sum/lc_12_count*100, lc20=lc_20_sum/lc_20_count*100, lc30=lc_30_sum/lc_30_count*100, lc40=lc_40_sum/lc_40_count*100, lc50=lc_50_sum/lc_50_count*100, lc60=lc_60_sum/lc_60_count*100, lc61=lc_61_sum/lc_61_count*100, lc62=lc_62_sum/lc_62_count*100, lc70=lc_70_sum/lc_70_count*100, lc71=lc_71_sum/lc_71_count*100, lc72=lc_72_sum/lc_72_count*100, lc80=lc_80_sum/lc_80_count*100, lc81=lc_81_sum/lc_81_count*100, lc82=lc_82_sum/lc_82_count*100, lc90=lc_90_sum/lc_90_count*100, lc100=lc_100_sum/lc_100_count*100, lc110=lc_110_sum/lc_110_count*100,lc120=lc_120_sum/lc_120_count*100,lc121=lc_121_sum/lc_121_count*100,lc122=lc_122_sum/lc_122_count*100,lc130=lc_130_sum/lc_130_count*100,lc140=lc_140_sum/lc_140_count*100,lc150=lc_150_sum/lc_150_count*100,lc152=lc_152_sum/lc_152_count*100,lc153=lc_153_sum/lc_153_count*100,lc160=lc_160_sum/lc_160_count*100,lc170=lc_170_sum/lc_170_count*100,lc180=lc_180_sum/lc_180_count*100,lc190=lc_190_sum/lc_190_count*100,lc200=lc_200_sum/lc_200_count*100,lc201=lc_201_sum/lc_201_count*100,lc202=lc_202_sum/lc_202_count*100,lc210=lc_210_sum/lc_210_count*100,lc220=lc_220_sum/lc_220_count*100) %>% 
  dplyr::select(!matches("_cou") & !matches("_sum")) %>% 
  #quito variables que no vaya a usar en el modelo 
  dplyr::select(-c(excel_file,excel_path, shp_file, shp_path, #species, #me cargo variables que no necesito para el modelo
                   grow_mean, #la variable grow tiene datos raros, as? que la voy a eliminar de nuestro set de datos
                   lc50, lc62, lc72, lc81, lc82,lc121, lc170,#, #Me cargo variables lc que en el summary son todo 0
                   #BR_todasr1
                   "individualCount...39", "individualCount...51"
  )) %>% #quito las bioregiones porque no las vamos a usar para el modelo
  filter(dens<50) %>%  #tres poligonos de >Spain con areas muy muy muy muy small
  # filter(dens>0 & dens<50) 
  select(locationID, locality, country, NUTS, BR_todasr1, area_km2, matches("harv"), dens, dens_r_trans, 
         matches("bio"), alt, snow, sun, Eucmean, hfp, 
         matches("lc"), 
         x_cen, y_cen)

sp$Eucmean<-ifelse(sp$Eucmean<0, 0.01, sp$Eucmean) #Como para wild boar transformamos la variable de Eucmean para eliminar los valores negativos por valores de0 muy peque?os
sp<-sp[!is.na(sp$dens),] #Remove NA values from RV: it does not remove any value

#Modeling phase
#1) data
colnames(sp)[15:33]<-gsub("_mean", "", colnames(sp)[15:33])
vars_std <- c(colnames(sp %>% as_tibble())[15:69])

#remove NA values from the dataset
sp<-sp[!is.na(sp$bio_01),] #Remove NA values from the spset
sp<-sp[!is.na(sp$snow),] #Remove NA values from the spset
sp<-sp[!is.na(sp$hfp),] #Remove NA values from the spset
sp<-sp[!is.na(sp$alt),] #Remove NA values from the spset
sp<-sp[!is.na(sp$lc10),] #Remove NA values from the spset

data<-sp  %>% mutate(across(all_of(vars_std),~ as.numeric(scale(.))))
summary(data)

#1.1 Creamos data set de calibrado y de validacion (80% / 20%)
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]; st_geometry(data_tra)<-NULL; st_geometry(data_val)<-NULL

#2.0 VIF <2
#2.1 seleccionamos variables para el vif, para eso miramos como es la distribucion de las variables
# hist(data_tra$Eucmean)
# hist(data_tra$alt) 
# hist(data_tra$snow)
# 
# hist(data_tra$lc10)
# hist(data_tra$lc11)
# hist(data_tra$lc12) #
# hist(data_tra$lc20) #
# hist(data_tra$lc30)
# hist(data_tra$lc40) 
# hist(data_tra$lc60) 
# hist(data_tra$lc70)
# hist(data_tra$lc90) 
# hist(data_tra$lc100)
# hist(data_tra$lc120) 
# hist(data_tra$lc130)
# 
# hist(data_tra$lc61) #
# hist(data_tra$lc71) #
# hist(data_tra$lc80) #
# hist(data_tra$lc110) #
# hist(data_tra$lc122) #
# hist(data_tra$lc140) #
# hist(data_tra$lc150) #
# hist(data_tra$lc152) #
# hist(data_tra$lc153) #
# hist(data_tra$lc160) #
# hist(data_tra$lc180) #
# hist(data_tra$lc190) #
# hist(data_tra$lc200) #
# hist(data_tra$lc201) #
# hist(data_tra$lc202) #
# hist(data_tra$lc210) #
# hist(data_tra$lc220) #

#2.1 hacemos el vif
variables<-data_tra %>% dplyr::select(matches("bio_"), alt, snow, Eucmean,lc10, lc11, lc12, lc30, lc40, lc60, lc70, lc90, 
                                      lc100, #lc120, 
                                      lc130, x_cen, y_cen)
resultado.vif0<-usdm::vifstep(variables, th=2); resultado.vif0

#2.2) Density Models: binomial negativa necesita que sea un numero entero, por eso se ha transformado la variable densidad *10000
#ZINB
count_vars <- c("bio_03","snow","lc10","lc30","lc40","lc60","lc70","lc100","lc130")
zi_vars <- count_vars
res <- stepwise_zeroinfl_alternate(
  data = data_tra,
  response = "dens_r_trans",
  count_fixed = c("x_cen","y_cen"),              # SIEMPRE dentro del count
  count_candidates = count_vars,
  zi_fixed = c("x_cen","y_cen"),                 # o c("x_cen","y_cen") si también los quieres fijos en ZI
  zi_candidates = zi_vars,
  start_in = "zi",                               # como tu gráfico: primero inflado, luego conteo
  dist = "negbin",
  criterion = "AIC",
  min_improve = 0                                # pon 2 si quieres “solo si mejora >2”
)

res$best_model
IbWildGoat_DensityModel<-res$best_model
# save(IbWildGoat_DensityModel, file="scripts/ModelResults/densityStepModel_IbWildGoat_Europe.Rdata")

#2.3) Calibration plots  validation data
data_val<-data_val %>%
  mutate(Pred_model = as.integer(predict(IbWildGoat_DensityModel, data_val,type = "response"))) %>% #predecir zeroinfl
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

# calibration plot
data_pos<-data_val %>% filter(dens>0)
cal_tot<-bind_rows(
  cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
)

ggplot(cal_tot, aes(x = pred_mean, y = obs_mean)) +
  geom_point(aes(size = n), shape = 21, fill = "grey70") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
  facet_wrap(~Type, scales = "free") +
  scale_size_continuous(name = "N per bin", range = c(2, 6)) +
  labs(
    x = "Predicted mean",
    y = "Observed mean"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

# s_class <- cut2(data_val$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean <- as.matrix(tapply(data_val$Pred_Red, s_class, mean)) # calculating mean values for predicted
# y_mean <- as.matrix(tapply(data_val$dens, s_class, mean)) # calculating mean values for observed
# # corbr1<-cor.test(data_val$Pred_Red, data_val$dens)
# data_val_w0<-data_val %>%filter(dens_r_trans>0)
# s_class_w0 <- cut2(data_val_w0$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean_w0 <- as.matrix(tapply(data_val_w0$Pred_Red, s_class_w0, mean)) # calculating mean values for predicted
# y_mean_w0 <- as.matrix(tapply(data_val_w0$dens, s_class_w0, mean)) # calculating mean values for observed

# cpE<-rbind(data.frame(s_mean=s_mean,   y_mean=y_mean, Type=as.factor("All data")),
#            data.frame(s_mean=s_mean_w0, y_mean=y_mean_w0, Type=as.factor("Data without 0")))
# ggplot(cpE, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=4)+ scale_colour_grey()+ scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="red", lty=2, cex=0.5)+
#   theme_bw()+
#   # xlim(0,max(s_mean_w0, max(na.omit(y_mean_w0))))+
#   # ylim(0,max(max(na.omit(y_mean_w0)), s_mean_w0))+
#   ylab("Observed Alpine ibex hunted bag density")+xlab("Predicted")+
#   theme(plot.title = element_text(hjust = 0.5))


#3.0) Spatial layer: 10x10km grid
grid_10km<-st_read(dsn="scripts/data/data_constrained_sp_20220902.gpkg", layer='grid_Goat') %>% dplyr::select(-Bioregion) %>% 
  rename(Bioregion=BR_T10, NUTS=NUT, Eucmean=Eucalyptus) %>% 
  mutate(NUTS=as.factor(NUTS), 
         Bioregion=as.factor(Bioregion),
         lc_50=as.numeric(lc_50),
         lc_62=as.numeric(lc_62),
         lc_72=as.numeric(lc_72), 
         lc_122=as.numeric(lc_122), 
         lc_153=as.numeric(lc_153), 
         Eucmean = replace_na(Eucmean, 0)) 
grid_10km<-grid_10km %>% 
  mutate(Eucmean = ifelse(grid_10km$Eucmean<0, 0, Eucmean)) %>% 
  mutate(x_cen=st_coordinates(st_centroid(grid_10km))[,1], y_cen=st_coordinates(st_centroid(grid_10km))[,2])
colnames(grid_10km)[32:65]<-gsub("lc_", "lc", colnames(grid_10km)[32:65])
colnames(grid_10km)[9:17]<-gsub("_", "_0", colnames(grid_10km)[9:17])


data %>% colnames()
means<-sp %>% as_tibble() %>% 
  select(15:69) %>%          # selecciona columnas
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
sdss<-sp %>% as_tibble() %>% 
  select(15:69) %>%          # selecciona columnas
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

library(modEvA) #MESS
summary(IbWildGoat_DensityModel) #Solo ponemos las variables que selecciona el modelo que son sobre las que vamos a hacer la prediccion
data_<-data %>% as.data.frame() %>% dplyr::select(c(x_cen, y_cen, bio_03, lc10, lc30, lc40, lc60, lc100, snow))

#10km
# model_vars<-c(x_cen, y_cen, bio_03, bio_10, bio_15, lc10, lc30, lc40, lc60, lc70, lc90, lc100, lc130, Eucmean)
grid_10km<-grid_10km %>% mutate(across(all_of(colnames(means)), ~ . - means[[cur_column()]])) %>% 
  mutate(across(all_of(colnames(means)), ~ . / sdss[[cur_column()]]))

grid_10km<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km)))
grid10km_df<-grid_10km %>% as.data.frame()  %>% dplyr::select("id_buena",
                                                              c(x_cen, y_cen, bio_03, lc10, lc30, lc40, lc60, lc100, snow))

mess_10km<-MESS(V=data_, P=grid10km_df, id.col=1)

grid_10km_<-grid_10km[c(-1:-73)]
mess_10km_<-mess_10km[c(-2:-10)]

grid10km_mess<-grid_10km %>%  full_join(mess_10km_, by="id_buena")

predRedDeer10km<-grid10km_mess %>% mutate(Pred_model = as.integer(predict(IbWildGoat_DensityModel, grid_10km, type = "response"))) %>% mutate(Pred_Red = Pred_model/10000) 
st_write(predRedDeer10km, 
         dsn="scripts/ModelResults/predictions/predicciones_20260203.gpkg",
         layer="ZINB_stepModel_grid10km_IbWildGoat", driver = "GPKG",append=TRUE)
# --------------------

# Moose #
# --------------------
# rm(list=ls())
data_sp<-st_read(dsn="scripts/data/predicciones_20220811.gpkg", 
                 layer='Moose_spatial_v5')

sp<-data_sp %>%
  mutate(lc10=lc_10_sum/lc_10_count*100, lc11=lc_11_sum/lc_11_count*100, lc12=lc_12_sum/lc_12_count*100, lc20=lc_20_sum/lc_20_count*100, lc30=lc_30_sum/lc_30_count*100, lc40=lc_40_sum/lc_40_count*100, lc50=lc_50_sum/lc_50_count*100, lc60=lc_60_sum/lc_60_count*100, lc61=lc_61_sum/lc_61_count*100, lc62=lc_62_sum/lc_62_count*100, lc70=lc_70_sum/lc_70_count*100, lc71=lc_71_sum/lc_71_count*100, lc72=lc_72_sum/lc_72_count*100, lc80=lc_80_sum/lc_80_count*100, lc81=lc_81_sum/lc_81_count*100, lc82=lc_82_sum/lc_82_count*100, lc90=lc_90_sum/lc_90_count*100, lc100=lc_100_sum/lc_100_count*100, lc110=lc_110_sum/lc_110_count*100,lc120=lc_120_sum/lc_120_count*100,lc121=lc_121_sum/lc_121_count*100,lc122=lc_122_sum/lc_122_count*100,lc130=lc_130_sum/lc_130_count*100,lc140=lc_140_sum/lc_140_count*100,lc150=lc_150_sum/lc_150_count*100,lc152=lc_152_sum/lc_152_count*100,lc153=lc_153_sum/lc_153_count*100,lc160=lc_160_sum/lc_160_count*100,lc170=lc_170_sum/lc_170_count*100,lc180=lc_180_sum/lc_180_count*100,lc190=lc_190_sum/lc_190_count*100,lc200=lc_200_sum/lc_200_count*100,lc201=lc_201_sum/lc_201_count*100,lc202=lc_202_sum/lc_202_count*100,lc210=lc_210_sum/lc_210_count*100,lc220=lc_220_sum/lc_220_count*100) %>% 
  dplyr::select(!matches("_cou") & !matches("_sum")) %>% 
  #quito variables que no vaya a usar en el modelo 
  dplyr::select(-c(excel_file,excel_path, shp_file, shp_path, #species, #me cargo variables que no necesito para el modelo
                   grow_mean, #la variable grow tiene datos raros, as? que la voy a eliminar de nuestro set de datos
                   lc50, lc62, lc72, lc81, lc82,lc121, lc170,#, #Me cargo variables lc que en el summary son todo 0
                   #BR_todasr1
                   "individualCount...39", "individualCount...51"
  )) %>% #quito las bioregiones porque no las vamos a usar para el modelo
  filter(dens<50) %>%  #tres poligonos de >Spain con areas muy muy muy muy small
  # filter(dens>0 & dens<50) 
  select(locationID, locality, country, NUTS, BR_todasr1, area_km2, matches("harv"), dens, dens_r_trans, 
         matches("bio"), alt, snow, sun, Eucmean, hfp, 
         matches("lc"), 
         x_cen, y_cen)

sp$Eucmean<-ifelse(sp$Eucmean<0, 0.01, sp$Eucmean) #Como para wild boar transformamos la variable de Eucmean para eliminar los valores negativos por valores de0 muy peque?os
sp<-sp[!is.na(sp$dens),] #Remove NA values from RV: it does not remove any value

#Modeling phase
#1) data
colnames(sp)[16:34]<-gsub("_mean", "", colnames(sp)[16:34])
vars_std <- c(colnames(sp %>% as_tibble())[16:70])

#remove NA values from the dataset
sp<-sp[!is.na(sp$bio_01),] #Remove NA values from the spset
sp<-sp[!is.na(sp$snow),] #Remove NA values from the spset
sp<-sp[!is.na(sp$hfp),] #Remove NA values from the spset
sp<-sp[!is.na(sp$alt),] #Remove NA values from the spset
sp<-sp[!is.na(sp$lc10),] #Remove NA values from the spset

data<-sp  %>% mutate(across(all_of(vars_std),~ as.numeric(scale(.))))
summary(data)

#1.1 Creamos data set de calibrado y de validacion (80% / 20%)
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]; st_geometry(data_tra)<-NULL; st_geometry(data_val)<-NULL

#2.0 VIF <2
#2.1 seleccionamos variables para el vif, para eso miramos como es la distribucion de las variables
# hist(data_tra$Eucmean)
# hist(data_tra$alt) 
# hist(data_tra$snow)
# 
# hist(data_tra$lc10)
# hist(data_tra$lc11)
# hist(data_tra$lc12) #
# hist(data_tra$lc20) #
# hist(data_tra$lc30)
# hist(data_tra$lc40) 
# hist(data_tra$lc60) 
# hist(data_tra$lc70)
# hist(data_tra$lc90) 
# hist(data_tra$lc100)
# hist(data_tra$lc120) 
# hist(data_tra$lc130)
# 
# hist(data_tra$lc61) #
# hist(data_tra$lc71) #
# hist(data_tra$lc80) #
# hist(data_tra$lc110) #
# hist(data_tra$lc122) #
# hist(data_tra$lc140) #
# hist(data_tra$lc150) #
# hist(data_tra$lc152) #
# hist(data_tra$lc153) #
# hist(data_tra$lc160) #
# hist(data_tra$lc180) #
# hist(data_tra$lc190) #
# hist(data_tra$lc200) #
# hist(data_tra$lc201) #
# hist(data_tra$lc202) #
# hist(data_tra$lc210) #
# hist(data_tra$lc220) #

#2.1 hacemos el vif
variables<-data_tra %>% dplyr::select(matches("bio_"), alt, snow, Eucmean, lc10, lc11, lc12, lc30, lc40, lc60, lc70, lc90, 
                                      lc100, #lc120, 
                                      lc130, x_cen, y_cen)
resultado.vif0<-usdm::vifstep(variables, th=2); resultado.vif0


#2.2) Density Models: binomial negativa necesita que sea un numero entero, por eso se ha transformado la variable densidad *10000
#ZINB
count_vars <- c("bio_15", "alt","lc10","lc30","lc40","lc60","lc70", "lc90","lc100","lc130")
zi_vars <- count_vars
res <- stepwise_zeroinfl_alternate(
  data = data_tra,
  response = "dens_r_trans",
  count_fixed = c("x_cen","y_cen"),              # SIEMPRE dentro del count
  count_candidates = count_vars,
  zi_fixed = c("x_cen","y_cen"),                 # o c("x_cen","y_cen") si también los quieres fijos en ZI
  zi_candidates = zi_vars,
  start_in = "zi",                               # como tu gráfico: primero inflado, luego conteo
  dist = "negbin",
  criterion = "AIC",
  min_improve = 0                                # pon 2 si quieres “solo si mejora >2”
)

res$best_model

Moose_DensityModel<-res$best_model
# save(Moose_DensityModel, file="scripts/ModelResults/densityStepModel_Moose_Europe.Rdata")

#2.3) Calibration plots  validation data
data_val<-data_val %>%
  mutate(Pred_model = as.integer(predict(Moose_DensityModel, data_val,type = "response"))) %>% #predecir zeroinfl
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

# calibration plot
data_pos<-data_val %>% filter(dens>0)
cal_tot<-bind_rows(
  cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
)

ggplot(cal_tot, aes(x = pred_mean, y = obs_mean)) +
  geom_point(aes(size = n), shape = 21, fill = "grey70") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
  facet_wrap(~Type, scales = "free") +
  scale_size_continuous(name = "N per bin", range = c(2, 6)) +
  labs(
    x = "Predicted mean",
    y = "Observed mean"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

# s_class <- cut2(data_val$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean <- as.matrix(tapply(data_val$Pred_Red, s_class, mean)) # calculating mean values for predicted
# y_mean <- as.matrix(tapply(data_val$dens, s_class, mean)) # calculating mean values for observed
# # corbr1<-cor.test(data_val$Pred_Red, data_val$dens)
# data_val_w0<-data_val %>%filter(dens_r_trans>0)
# s_class_w0 <- cut2(data_val_w0$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean_w0 <- as.matrix(tapply(data_val_w0$Pred_Red, s_class_w0, mean)) # calculating mean values for predicted
# y_mean_w0 <- as.matrix(tapply(data_val_w0$dens, s_class_w0, mean)) # calculating mean values for observed
# 
# cpE<-rbind(data.frame(s_mean=s_mean,   y_mean=y_mean, Type=as.factor("All data")),
#            data.frame(s_mean=s_mean_w0, y_mean=y_mean_w0, Type=as.factor("Data without 0")))
# ggplot(cpE, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=4)+ scale_colour_grey()+ scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="red", lty=2, cex=0.5)+
#   theme_bw()+
#   # xlim(0,max(s_mean_w0, max(na.omit(y_mean_w0))))+
#   # ylim(0,max(max(na.omit(y_mean_w0)), s_mean_w0))+
#   ylab("Observed Alpine ibex hunted bag density")+xlab("Predicted")+
#   theme(plot.title = element_text(hjust = 0.5))
# 


#3.0) Spatial layer: 10x10km grid
grid_10km<-st_read(dsn="scripts/data/data_constrained_sp_20220902.gpkg", layer='grid_Moose') %>% dplyr::select(-Bioregion) %>% 
  rename(Bioregion=BR_T10, NUTS=NUT, Eucmean=Eucalyptus) %>% 
  mutate(NUTS=as.factor(NUTS), 
         Bioregion=as.factor(Bioregion),
         lc_50=as.numeric(lc_50),
         lc_62=as.numeric(lc_62),
         lc_72=as.numeric(lc_72), 
         lc_122=as.numeric(lc_122), 
         lc_153=as.numeric(lc_153), 
         Eucmean = replace_na(Eucmean, 0)) 
grid_10km<-grid_10km %>% 
  mutate(Eucmean = ifelse(grid_10km$Eucmean<0, 0, Eucmean)) %>% 
  mutate(x_cen=st_coordinates(st_centroid(grid_10km))[,1], y_cen=st_coordinates(st_centroid(grid_10km))[,2])
colnames(grid_10km)[32:65]<-gsub("lc_", "lc", colnames(grid_10km)[32:65])
colnames(grid_10km)[9:17]<-gsub("_", "_0", colnames(grid_10km)[9:17])

data %>% colnames()
means<-sp %>% as_tibble() %>% 
  select(16:70) %>%          # selecciona columnas
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
sdss<-sp %>% as_tibble() %>% 
  select(16:70) %>%          # selecciona columnas
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

library(modEvA) #MESS
summary(Moose_DensityModel) #Solo ponemos las variables que selecciona el modelo que son sobre las que vamos a hacer la prediccion
data_<-data %>% as.data.frame() %>% dplyr::select(c(x_cen, y_cen, bio_15, lc10, lc30, lc40, lc60, lc70, lc90, lc130))

#10km
grid_10km<-grid_10km %>% mutate(across(all_of(colnames(means)), ~ . - means[[cur_column()]])) %>% 
  mutate(across(all_of(colnames(means)), ~ . / sdss[[cur_column()]]))

grid_10km<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km)))
grid10km_df<-grid_10km %>% as.data.frame()  %>% dplyr::select("id_buena",
                                                              c(x_cen, y_cen, bio_15, lc10, lc30, lc40, lc60, lc70, lc90, lc130))
mess_10km<-MESS(V=data_, P=grid10km_df, id.col=1)

grid_10km_<-grid_10km[c(-1:-73)]
mess_10km_<-mess_10km[c(-2:-11)]; names(mess_10km_)

grid10km_mess<-grid_10km %>%  full_join(mess_10km_, by="id_buena")

predRedDeer10km<-grid10km_mess %>% mutate(Pred_model = as.integer(predict(Moose_DensityModel, grid_10km, type = "response"))) %>% mutate(Pred_Red = Pred_model/10000) 
st_write(predRedDeer10km, 
         dsn="scripts/ModelResults/predictions/predicciones_20260203.gpkg",
         layer="ZINB_stepModel_grid10km_Moose", driver = "GPKG",append=TRUE)
# --------------------

# Mouflon #
# --------------------
# rm(list=ls())
data_sp<-st_read(dsn="scripts/data/predicciones_20220811.gpkg", 
                 layer='Mouflon_spatial_v5')
data_sp<-data_sp %>% mutate(x_cen=st_coordinates(st_centroid(data_sp))[,1], y_cen=st_coordinates(st_centroid(data_sp))[,2])
sp<-data_sp %>%
  mutate(lc10=lc_10_sum/lc_10_count*100, lc11=lc_11_sum/lc_11_count*100, lc12=lc_12_sum/lc_12_count*100, lc20=lc_20_sum/lc_20_count*100, lc30=lc_30_sum/lc_30_count*100, lc40=lc_40_sum/lc_40_count*100, lc50=lc_50_sum/lc_50_count*100, lc60=lc_60_sum/lc_60_count*100, lc61=lc_61_sum/lc_61_count*100, lc62=lc_62_sum/lc_62_count*100, lc70=lc_70_sum/lc_70_count*100, lc71=lc_71_sum/lc_71_count*100, lc72=lc_72_sum/lc_72_count*100, lc80=lc_80_sum/lc_80_count*100, lc81=lc_81_sum/lc_81_count*100, lc82=lc_82_sum/lc_82_count*100, lc90=lc_90_sum/lc_90_count*100, lc100=lc_100_sum/lc_100_count*100, lc110=lc_110_sum/lc_110_count*100,lc120=lc_120_sum/lc_120_count*100,lc121=lc_121_sum/lc_121_count*100,lc122=lc_122_sum/lc_122_count*100,lc130=lc_130_sum/lc_130_count*100,lc140=lc_140_sum/lc_140_count*100,lc150=lc_150_sum/lc_150_count*100,lc152=lc_152_sum/lc_152_count*100,lc153=lc_153_sum/lc_153_count*100,lc160=lc_160_sum/lc_160_count*100,lc170=lc_170_sum/lc_170_count*100,lc180=lc_180_sum/lc_180_count*100,lc190=lc_190_sum/lc_190_count*100,lc200=lc_200_sum/lc_200_count*100,lc201=lc_201_sum/lc_201_count*100,lc202=lc_202_sum/lc_202_count*100,lc210=lc_210_sum/lc_210_count*100,lc220=lc_220_sum/lc_220_count*100) %>% 
  dplyr::select(!matches("_cou") & !matches("_sum")) %>% 
  #quito variables que no vaya a usar en el modelo 
  dplyr::select(-c(excel_file,excel_path, shp_file, shp_path, #species, #me cargo variables que no necesito para el modelo
                   grow_mean, #la variable grow tiene datos raros, as? que la voy a eliminar de nuestro set de datos
                   lc50, lc62, lc72, lc81, lc82,lc121, lc170,#, #Me cargo variables lc que en el summary son todo 0
                   #BR_todasr1
                   "individualCount...39", "individualCount...51"
  )) %>% #quito las bioregiones porque no las vamos a usar para el modelo
  filter(dens<50) %>%  #tres poligonos de >Spain con areas muy muy muy muy small
  # filter(dens>0 & dens<50) 
  select(locationID, locality, country, NUTS, BR_todasr1, area_km2, matches("harv"), dens, dens_r_trans, 
         matches("bio"), alt, snow, sun, Eucmean, hfp, 
         matches("lc"), 
         x_cen, y_cen)

sp$Eucmean<-ifelse(sp$Eucmean<0, 0.01, sp$Eucmean) #Como para wild boar transformamos la variable de Eucmean para eliminar los valores negativos por valores de0 muy peque?os
sp<-sp[!is.na(sp$dens),] #Remove NA values from RV: it does not remove any value

#Modeling phase
#1) data
colnames(sp)[17:35]<-gsub("_mean", "", colnames(sp)[17:35])
vars_std <- c(colnames(sp %>% as_tibble())[17:71])

#remove NA values from the dataset
sp<-sp[!is.na(sp$bio_01),] #Remove NA values from the spset
sp<-sp[!is.na(sp$snow),] #Remove NA values from the spset
sp<-sp[!is.na(sp$hfp),] #Remove NA values from the spset
sp<-sp[!is.na(sp$alt),] #Remove NA values from the spset
sp<-sp[!is.na(sp$lc10),] #Remove NA values from the spset

data<-sp  %>% mutate(across(all_of(vars_std),~ as.numeric(scale(.))))
summary(data)


#1.1 Creamos data set de calibrado y de validacion (80% / 20%)
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]; st_geometry(data_tra)<-NULL; st_geometry(data_val)<-NULL

#2.0 VIF <2
#2.1 seleccionamos variables para el vif, para eso miramos como es la distribucion de las variables
# hist(data_tra$Eucmean)
# hist(data_tra$alt) 
# hist(data_tra$snow)
# 
# hist(data_tra$lc10)
# hist(data_tra$lc11)
# hist(data_tra$lc12) #
# hist(data_tra$lc20) #
# hist(data_tra$lc30)
# hist(data_tra$lc40) 
# hist(data_tra$lc60) 
# hist(data_tra$lc70)
# hist(data_tra$lc90) 
# hist(data_tra$lc100)
# hist(data_tra$lc120) 
# hist(data_tra$lc130)
# 
# hist(data_tra$lc61) #
# hist(data_tra$lc71) #
# hist(data_tra$lc80) #
# hist(data_tra$lc110) #
# hist(data_tra$lc122) #
# hist(data_tra$lc140) #
# hist(data_tra$lc150) #
# hist(data_tra$lc152) #
# hist(data_tra$lc153) #
# hist(data_tra$lc160) #
# hist(data_tra$lc180) #
# hist(data_tra$lc190) #
# hist(data_tra$lc200) #
# hist(data_tra$lc201) #
# hist(data_tra$lc202) #
# hist(data_tra$lc210) #
# hist(data_tra$lc220) #

#2.1 hacemos el vif
variables<-data_tra %>% dplyr::select(matches("bio_"), alt, snow, Eucmean,lc10, lc11, lc12, lc30, lc40, lc60, lc70, lc90, 
                                      lc100, #lc120, 
                                      lc130, x_cen, y_cen)
resultado.vif0<-usdm::vifstep(variables, th=2); resultado.vif0

#2.2) Density Models: binomial negativa necesita que sea un numero entero, por eso se ha transformado la variable densidad *10000
#ZINB
#modelo sin stepwise
# 1) ajustamos modelo con ZI solo intercepto  + Count variables
count_vars <- c("bio_07", "bio_14", "bio_19", "alt", "Eucmean", "lc10","lc30","lc40","lc60","lc70", "lc90","lc100","lc130")
zi_vars <- count_vars

res <- stepwise_zeroinfl_alternate(
  data = data_tra,
  response = "dens_r_trans",
  count_fixed = c("x_cen","y_cen"),              # SIEMPRE dentro del count
  count_candidates = count_vars,
  zi_fixed = c("x_cen","y_cen"),                 # o c("x_cen","y_cen") si también los quieres fijos en ZI
  zi_candidates = zi_vars,
  start_in = "zi",                               # como tu gráfico: primero inflado, luego conteo
  dist = "negbin",
  criterion = "AIC",
  min_improve = 0                                # pon 2 si quieres “solo si mejora >2”
)
res$best_model
Mouflon_DensityModel<-res$best_model
# save(Mouflon_DensityModel, file="scripts/ModelResults/densityStepModel_Mouflon_Europe.Rdata")

#2.3) Calibration plots  validation data
data_val<-data_val %>% 
  mutate(Pred_model = as.integer(predict(Mouflon_DensityModel, data_val,type = "response"))) %>% #predecir zeroinfl
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

data_pos<-data_val %>% filter(dens>0) 
cal_tot<-bind_rows(
  cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
)

ggplot(cal_tot, aes(x = pred_mean, y = obs_mean)) +
  geom_point(aes(size = n), shape = 21, fill = "grey70") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
  facet_wrap(~Type, scales = "free") +
  scale_size_continuous(name = "N per bin", range = c(2, 6)) +
  labs(
    x = "Predicted mean",
    y = "Observed mean"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

# s_class <- cut2(data_val$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean <- as.matrix(tapply(data_val$Pred_Red, s_class, mean)) # calculating mean values for predicted
# y_mean <- as.matrix(tapply(data_val$dens, s_class, mean)) # calculating mean values for observed
# # corbr1<-cor.test(data_val$Pred_Red, data_val$dens)
# data_val_w0<-data_val %>%filter(dens_r_trans>0)
# s_class_w0 <- cut2(data_val_w0$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean_w0 <- as.matrix(tapply(data_val_w0$Pred_Red, s_class_w0, mean)) # calculating mean values for predicted
# y_mean_w0 <- as.matrix(tapply(data_val_w0$dens, s_class_w0, mean)) # calculating mean values for observed
# 
# cpE<-rbind(data.frame(s_mean=s_mean,   y_mean=y_mean, Type=as.factor("All data")),
#            data.frame(s_mean=s_mean_w0, y_mean=y_mean_w0, Type=as.factor("Data without 0")))
# ggplot(cpE, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=4)+ scale_colour_grey()+ scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="red", lty=2, cex=0.5)+
#   theme_bw()+
#   # xlim(0,max(s_mean_w0, max(na.omit(y_mean_w0))))+
#   # ylim(0,max(max(na.omit(y_mean_w0)), s_mean_w0))+
#   ylab("Observed Alpine ibex hunted bag density")+xlab("Predicted")+
#   theme(plot.title = element_text(hjust = 0.5))


#3.0) Spatial layer: 10x10km grid
grid_10km<-st_read(dsn="scripts/data/data_constrained_sp_20220902.gpkg", layer='grid_Mouflon') %>% dplyr::select(-Bioregion) %>% 
  rename(Bioregion=BR_T10, NUTS=NUT, Eucmean=Eucalyptus) %>% 
  mutate(NUTS=as.factor(NUTS), 
         Bioregion=as.factor(Bioregion),
         lc_50=as.numeric(lc_50),
         lc_62=as.numeric(lc_62),
         lc_72=as.numeric(lc_72), 
         lc_122=as.numeric(lc_122), 
         lc_153=as.numeric(lc_153), 
         Eucmean = replace_na(Eucmean, 0)) 
grid_10km<-grid_10km %>% 
  mutate(Eucmean = ifelse(grid_10km$Eucmean<0, 0, Eucmean)) %>% 
  mutate(x_cen=st_coordinates(st_centroid(grid_10km))[,1], y_cen=st_coordinates(st_centroid(grid_10km))[,2])
colnames(grid_10km)[32:65]<-gsub("lc_", "lc", colnames(grid_10km)[32:65])
colnames(grid_10km)[9:17]<-gsub("_", "_0", colnames(grid_10km)[9:17])


data %>% colnames()
means<-sp %>% as_tibble() %>% 
  select(17:71) %>%          # selecciona columnas
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
sdss<-sp %>% as_tibble() %>% 
  select(17:71) %>%          # selecciona columnas
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

library(modEvA) #MESS
summary(Mouflon_DensityModel) #Solo ponemos las variables que selecciona el modelo que son sobre las que vamos a hacer la prediccion
data_<-data %>% as.data.frame() %>% dplyr::select(c(x_cen, y_cen, bio_07, bio_14, bio_19, alt, lc10, lc40, lc60, lc130, Eucmean))

#10km
grid_10km<-grid_10km %>% mutate(across(all_of(colnames(means)), ~ . - means[[cur_column()]])) %>% 
  mutate(across(all_of(colnames(means)), ~ . / sdss[[cur_column()]]))

grid_10km<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km)))
grid10km_df<-grid_10km %>% as.data.frame()  %>% dplyr::select("id_buena",
                                                              c(x_cen, y_cen, bio_07, bio_14, bio_19, alt, lc10, lc40, lc60, lc130, Eucmean))

mess_10km<-MESS(V=data_, P=grid10km_df, id.col=1)

grid_10km_<-grid_10km[c(-1:-73)]
mess_10km_<-mess_10km[c(-2:-12)]

grid10km_mess<-grid_10km %>%  full_join(mess_10km_, by="id_buena")

predRedDeer10km<-grid10km_mess %>% mutate(Pred_model = as.integer(predict(Mouflon_DensityModel, grid_10km, type = "response"))) %>% mutate(Pred_Red = Pred_model/10000) 
st_write(predRedDeer10km, 
         dsn="scripts/ModelResults/predictions/predicciones_20260203.gpkg",
         layer="ZINB_stepModel_grid10km_Mouflon", driver = "GPKG",append=TRUE)
# --------------------

# Northern chamois
# --------------------
data_sp<-st_read(dsn="scripts/data/predicciones_20220811.gpkg", 
                 layer='Chamois_spatial_v5')

sp<-data_sp %>%
  mutate(lc10=lc_10_sum/lc_10_count*100, lc11=lc_11_sum/lc_11_count*100, lc12=lc_12_sum/lc_12_count*100, lc20=lc_20_sum/lc_20_count*100, lc30=lc_30_sum/lc_30_count*100, lc40=lc_40_sum/lc_40_count*100, lc50=lc_50_sum/lc_50_count*100, lc60=lc_60_sum/lc_60_count*100, lc61=lc_61_sum/lc_61_count*100, lc62=lc_62_sum/lc_62_count*100, lc70=lc_70_sum/lc_70_count*100, lc71=lc_71_sum/lc_71_count*100, lc72=lc_72_sum/lc_72_count*100, lc80=lc_80_sum/lc_80_count*100, lc81=lc_81_sum/lc_81_count*100, lc82=lc_82_sum/lc_82_count*100, lc90=lc_90_sum/lc_90_count*100, lc100=lc_100_sum/lc_100_count*100, lc110=lc_110_sum/lc_110_count*100,lc120=lc_120_sum/lc_120_count*100,lc121=lc_121_sum/lc_121_count*100,lc122=lc_122_sum/lc_122_count*100,lc130=lc_130_sum/lc_130_count*100,lc140=lc_140_sum/lc_140_count*100,lc150=lc_150_sum/lc_150_count*100,lc152=lc_152_sum/lc_152_count*100,lc153=lc_153_sum/lc_153_count*100,lc160=lc_160_sum/lc_160_count*100,lc170=lc_170_sum/lc_170_count*100,lc180=lc_180_sum/lc_180_count*100,lc190=lc_190_sum/lc_190_count*100,lc200=lc_200_sum/lc_200_count*100,lc201=lc_201_sum/lc_201_count*100,lc202=lc_202_sum/lc_202_count*100,lc210=lc_210_sum/lc_210_count*100,lc220=lc_220_sum/lc_220_count*100) %>% 
  dplyr::select(!matches("_cou") & !matches("_sum")) %>% 
  #quito variables que no vaya a usar en el modelo 
  dplyr::select(-c(excel_file,excel_path, shp_file, shp_path, #species, #me cargo variables que no necesito para el modelo
                   grow_mean, #la variable grow tiene datos raros, as? que la voy a eliminar de nuestro set de datos
                   lc50, lc62, lc72, lc81, lc82,lc121, lc170,#, #Me cargo variables lc que en el summary son todo 0
                   #BR_todasr1
                   "individualCount...39", "individualCount...51"
  )) %>% #quito las bioregiones porque no las vamos a usar para el modelo
  filter(dens<50) %>%  #tres poligonos de >Spain con areas muy muy muy muy small
  select(locationID, locality, country, NUTS, BR_todasr1, area_km2, matches("harv"), dens, dens_r_trans, 
         matches("bio"), alt, snow, sun, Eucmean, hfp, 
         matches("lc"), 
         x_cen, y_cen)

sp$Eucmean<-ifelse(sp$Eucmean<0, 0.01, sp$Eucmean) #Como para wild boar transformamos la variable de Eucmean para eliminar los valores negativos por valores de0 muy peque?os
sp<-sp[!is.na(sp$dens),] #Remove NA values from RV: it does not remove any value

#Modeling phase
#1) data
colnames(sp)[16:34]<-gsub("_mean", "", colnames(sp)[16:34])
vars_std <- c(colnames(sp %>% as_tibble())[16:70])

#remove NA values from the dataset
sp<-sp[!is.na(sp$bio_01),] #Remove NA values from the spset
sp<-sp[!is.na(sp$snow),] #Remove NA values from the spset
sp<-sp[!is.na(sp$hfp),] #Remove NA values from the spset
sp<-sp[!is.na(sp$alt),] #Remove NA values from the spset
sp<-sp[!is.na(sp$lc10),] #Remove NA values from the spset

data<-sp  %>% mutate(across(all_of(vars_std),~ as.numeric(scale(.))))
summary(data)

#1.1 Creamos data set de calibrado y de validacion (80% / 20%)
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]; st_geometry(data_tra)<-NULL; st_geometry(data_val)<-NULL

#2.0 VIF <2
#2.1 seleccionamos variables para el vif, para eso miramos como es la distribucion de las variables
# hist(data_tra$Eucmean)
# hist(data_tra$alt) 
# hist(data_tra$snow)
# 
# hist(data_tra$lc10)
# hist(data_tra$lc11)
# hist(data_tra$lc12) #
# hist(data_tra$lc20) #
# hist(data_tra$lc30)
# hist(data_tra$lc40) 
# hist(data_tra$lc60) 
# hist(data_tra$lc70)
# hist(data_tra$lc90) 
# hist(data_tra$lc100)
# hist(data_tra$lc120) 
# hist(data_tra$lc130)
# 
# hist(data_tra$lc61) #
# hist(data_tra$lc71) #
# hist(data_tra$lc80) #
# hist(data_tra$lc110) #
# hist(data_tra$lc122) #
# hist(data_tra$lc140) #
# hist(data_tra$lc150) #
# hist(data_tra$lc152) #
# hist(data_tra$lc153) #
# hist(data_tra$lc160) #
# hist(data_tra$lc180) #
# hist(data_tra$lc190) #
# hist(data_tra$lc200) #
# hist(data_tra$lc201) #
# hist(data_tra$lc202) #
# hist(data_tra$lc210) #
# hist(data_tra$lc220) #

#2.1 hacemos el vif
variables<-data_tra %>% dplyr::select(matches("bio_"), alt, snow, Eucmean,lc10, lc11, lc12, lc30, lc40, lc60, lc70, lc90, 
                                      lc100, #lc120, 
                                      lc130, x_cen, y_cen)
resultado.vif0<-usdm::vifstep(variables, th=2); resultado.vif0


#2.2) Density Models: binomial negativa necesita que sea un numero entero, por eso se ha transformado la variable densidad *10000
#ZINB
count_vars <- c("bio_03", "bio_09", "bio_15", "bio_19", "lc10","lc30","lc40","lc60","lc70", "lc90","lc100","lc130")
zi_vars <- count_vars
res <- stepwise_zeroinfl_alternate(
  data = data_tra,
  response = "dens_r_trans",
  count_fixed = c("x_cen","y_cen"),              # SIEMPRE dentro del count
  count_candidates = count_vars,
  zi_fixed = c("x_cen","y_cen"),                 # o c("x_cen","y_cen") si también los quieres fijos en ZI
  zi_candidates = zi_vars,
  start_in = "zi",                               # como tu gráfico: primero inflado, luego conteo
  dist = "negbin",
  criterion = "AIC",
  min_improve = 0                                # pon 2 si quieres “solo si mejora >2”
)

res$best_model
NorthChamois_DensityModel<-res$best_model
# save(NorthChamois_DensityModel, file="scripts/ModelResults/densityStepModel_NorthChamois_Europe.Rdata")

#2.3) Calibration plots  validation data
data_val<-data_val %>% 
  mutate(Pred_model = as.integer(predict(NorthChamois_DensityModel, data_val,type = "response"))) %>% #predecir zeroinfl
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

data_pos<-data_val %>% filter(dens>0) 
cal_tot<-bind_rows(
  cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"), 
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
)

ggplot(cal_tot, aes(x = pred_mean, y = obs_mean)) +
  geom_point(aes(size = n), shape = 21, fill = "grey70") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
  facet_wrap(~Type, scales = "free") +
  scale_size_continuous(name = "N per bin", range = c(2, 6)) +
  labs(
    x = "Predicted mean",
    y = "Observed mean"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

# s_class <- cut2(data_val$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean <- as.matrix(tapply(data_val$Pred_Red, s_class, mean)) # calculating mean values for predicted
# y_mean <- as.matrix(tapply(data_val$dens, s_class, mean)) # calculating mean values for observed
# # corbr1<-cor.test(data_val$Pred_Red, data_val$dens)
# data_val_w0<-data_val %>%filter(dens_r_trans>0)
# s_class_w0 <- cut2(data_val_w0$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean_w0 <- as.matrix(tapply(data_val_w0$Pred_Red, s_class_w0, mean)) # calculating mean values for predicted
# y_mean_w0 <- as.matrix(tapply(data_val_w0$dens, s_class_w0, mean)) # calculating mean values for observed
# 
# cpE<-rbind(data.frame(s_mean=s_mean,   y_mean=y_mean, Type=as.factor("All data")),
#            data.frame(s_mean=s_mean_w0, y_mean=y_mean_w0, Type=as.factor("Data without 0")))
# ggplot(cpE, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=4)+ scale_colour_grey()+ scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="red", lty=2, cex=0.5)+
#   theme_bw()+
#   # xlim(0,max(s_mean_w0, max(na.omit(y_mean_w0))))+
#   # ylim(0,max(max(na.omit(y_mean_w0)), s_mean_w0))+
#   ylab("Observed Alpine ibex hunted bag density")+xlab("Predicted")+
#   theme(plot.title = element_text(hjust = 0.5))
# 

#3.0) Spatial layer: 10x10km grid
grid_10km<-st_read(dsn="scripts/data/data_constrained_sp_20220902.gpkg", layer='grid_Chamois') %>% dplyr::select(-Bioregion) %>% 
  rename(Bioregion=BR_T10, NUTS=NUT, Eucmean=Eucalyptus) %>% 
  mutate(NUTS=as.factor(NUTS), 
         Bioregion=as.factor(Bioregion),
         lc_50=as.numeric(lc_50),
         lc_62=as.numeric(lc_62),
         lc_72=as.numeric(lc_72), 
         lc_122=as.numeric(lc_122), 
         lc_153=as.numeric(lc_153), 
         Eucmean = replace_na(Eucmean, 0)) 
grid_10km<-grid_10km %>% 
  mutate(Eucmean = ifelse(grid_10km$Eucmean<0, 0, Eucmean)) %>% 
  mutate(x_cen=st_coordinates(st_centroid(grid_10km))[,1], y_cen=st_coordinates(st_centroid(grid_10km))[,2])
colnames(grid_10km)[32:65]<-gsub("lc_", "lc", colnames(grid_10km)[32:65])
colnames(grid_10km)[9:17]<-gsub("_", "_0", colnames(grid_10km)[9:17])


data %>% colnames()
means<-sp %>% as_tibble() %>% 
  select(16:70) %>%          # selecciona columnas
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
sdss<-sp %>% as_tibble() %>% 
  select(16:70) %>%          # selecciona columnas
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

library(modEvA) #MESS
summary(NorthChamois_DensityModel) #Solo ponemos las variables que selecciona el modelo que son sobre las que vamos a hacer la prediccion
data_<-data %>% as.data.frame() %>% dplyr::select(c(x_cen, y_cen, bio_03, bio_09, bio_15, bio_19))

#10km
grid_10km<-grid_10km %>% mutate(across(all_of(colnames(means)), ~ . - means[[cur_column()]])) %>% 
  mutate(across(all_of(colnames(means)), ~ . / sdss[[cur_column()]]))

grid_10km<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km)))
grid10km_df<-grid_10km %>% as.data.frame()  %>% dplyr::select("id_buena",
                                                              c(x_cen, y_cen, bio_03, bio_09, bio_15, bio_19))
mess_10km<-MESS(V=data_, P=grid10km_df, id.col=1)

grid_10km_<-grid_10km[c(-1:-73)]
mess_10km_<-mess_10km[c(-2:-7)]; names(mess_10km_)

grid10km_mess<-grid_10km %>%  full_join(mess_10km_, by="id_buena")

predRedDeer10km<-grid10km_mess %>% mutate(Pred_model = as.integer(predict(NorthChamois_DensityModel, grid_10km, type = "response"))) %>% mutate(Pred_Red = Pred_model/10000) 
st_write(predRedDeer10km, 
         dsn="scripts/ModelResults/predictions/predicciones_20260203.gpkg",
         layer="ZINB_stepModel_grid10km_NorthChamois", driver = "GPKG",append=TRUE)
# --------------------

# Pyrenean chamois #
# --------------------
data_sp<-st_read(dsn="scripts/data/predicciones_20220811.gpkg", 
                 layer='Isard_spatial_v5')

sp<-data_sp %>%
  mutate(lc10=lc_10_sum/lc_10_count*100, lc11=lc_11_sum/lc_11_count*100, lc12=lc_12_sum/lc_12_count*100, lc20=lc_20_sum/lc_20_count*100, lc30=lc_30_sum/lc_30_count*100, lc40=lc_40_sum/lc_40_count*100, lc50=lc_50_sum/lc_50_count*100, lc60=lc_60_sum/lc_60_count*100, lc61=lc_61_sum/lc_61_count*100, lc62=lc_62_sum/lc_62_count*100, lc70=lc_70_sum/lc_70_count*100, lc71=lc_71_sum/lc_71_count*100, lc72=lc_72_sum/lc_72_count*100, lc80=lc_80_sum/lc_80_count*100, lc81=lc_81_sum/lc_81_count*100, lc82=lc_82_sum/lc_82_count*100, lc90=lc_90_sum/lc_90_count*100, lc100=lc_100_sum/lc_100_count*100, lc110=lc_110_sum/lc_110_count*100,lc120=lc_120_sum/lc_120_count*100,lc121=lc_121_sum/lc_121_count*100,lc122=lc_122_sum/lc_122_count*100,lc130=lc_130_sum/lc_130_count*100,lc140=lc_140_sum/lc_140_count*100,lc150=lc_150_sum/lc_150_count*100,lc152=lc_152_sum/lc_152_count*100,lc153=lc_153_sum/lc_153_count*100,lc160=lc_160_sum/lc_160_count*100,lc170=lc_170_sum/lc_170_count*100,lc180=lc_180_sum/lc_180_count*100,lc190=lc_190_sum/lc_190_count*100,lc200=lc_200_sum/lc_200_count*100,lc201=lc_201_sum/lc_201_count*100,lc202=lc_202_sum/lc_202_count*100,lc210=lc_210_sum/lc_210_count*100,lc220=lc_220_sum/lc_220_count*100) %>% 
  dplyr::select(!matches("_cou") & !matches("_sum")) %>% 
  #quito variables que no vaya a usar en el modelo 
  dplyr::select(-c(excel_file,excel_path, shp_file, shp_path, #species, #me cargo variables que no necesito para el modelo
                   grow_mean, #la variable grow tiene datos raros, as? que la voy a eliminar de nuestro set de datos
                   lc50, lc62, lc72, lc81, lc82,lc121, lc170,#, #Me cargo variables lc que en el summary son todo 0
                   #BR_todasr1
                   "individualCount...39", "individualCount...51"
  )) %>% #quito las bioregiones porque no las vamos a usar para el modelo
  filter(dens<50) %>%  #tres poligonos de >Spain con areas muy muy muy muy small
  # filter(dens>0 & dens<50) 
  select(locationID, locality, country, NUTS, BR_todasr1, area_km2, matches("harv"), dens, dens_r_trans, 
         matches("bio"), alt, snow, sun, Eucmean, hfp, 
         matches("lc"), 
         x_cen, y_cen)

sp$Eucmean<-ifelse(sp$Eucmean<0, 0.01, sp$Eucmean) #Como para wild boar transformamos la variable de Eucmean para eliminar los valores negativos por valores de0 muy peque?os
sp<-sp[!is.na(sp$dens),] #Remove NA values from RV: it does not remove any value

#Modeling phase
#1) data
colnames(sp)[15:33]<-gsub("_mean", "", colnames(sp)[15:33])
vars_std <- c(colnames(sp %>% as_tibble())[15:69])

#remove NA values from the dataset
sp<-sp[!is.na(sp$bio_01),] #Remove NA values from the spset
sp<-sp[!is.na(sp$snow),] #Remove NA values from the spset
sp<-sp[!is.na(sp$hfp),] #Remove NA values from the spset
sp<-sp[!is.na(sp$alt),] #Remove NA values from the spset
sp<-sp[!is.na(sp$lc10),] #Remove NA values from the spset

data<-sp  %>% mutate(across(all_of(vars_std),~ as.numeric(scale(.))))
summary(data)

#1.1 Creamos data set de calibrado y de validacion (80% / 20%)
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]; st_geometry(data_tra)<-NULL; st_geometry(data_val)<-NULL

#2.0 VIF <2
#2.1 seleccionamos variables para el vif, para eso miramos como es la distribucion de las variables
# hist(data_tra$Eucmean)
# hist(data_tra$alt) 
# hist(data_tra$snow)
# 
# hist(data_tra$lc10)
# hist(data_tra$lc11)
# hist(data_tra$lc12) #
# hist(data_tra$lc20) #
# hist(data_tra$lc30)
# hist(data_tra$lc40) 
# hist(data_tra$lc60) 
# hist(data_tra$lc70)
# hist(data_tra$lc90) 
# hist(data_tra$lc100)
# hist(data_tra$lc120) 
# hist(data_tra$lc130)
# 
# hist(data_tra$lc61) #
# hist(data_tra$lc71) #
# hist(data_tra$lc80) #
# hist(data_tra$lc110) #
# hist(data_tra$lc122) #
# hist(data_tra$lc140) #
# hist(data_tra$lc150) #
# hist(data_tra$lc152) #
# hist(data_tra$lc153) #
# hist(data_tra$lc160) #
# hist(data_tra$lc180) #
# hist(data_tra$lc190) #
# hist(data_tra$lc200) #
# hist(data_tra$lc201) #
# hist(data_tra$lc202) #
# hist(data_tra$lc210) #
# hist(data_tra$lc220) #

#2.1 hacemos el vif
variables<-data_tra %>% dplyr::select(matches("bio_"), alt, snow, Eucmean,lc10, lc11, lc12, lc30, lc40, lc60, lc70, lc90, 
                                      lc100, #lc120, 
                                      lc130, x_cen, y_cen)
resultado.vif0<-usdm::vifstep(variables, th=2); resultado.vif0

#2.2) Density Models: binomial negativa necesita que sea un numero entero, por eso se ha transformado la variable densidad *10000
#ZINB
count_vars <- c("bio_03", "bio_09", "bio_15", "alt", "lc10","lc30","lc40","lc60", "lc70", "lc90","lc100","lc130")
zi_vars <- count_vars

res <- stepwise_zeroinfl_alternate(
  data = data_tra,
  response = "dens_r_trans",
  count_fixed = c("x_cen","y_cen"),              # SIEMPRE dentro del count
  count_candidates = count_vars,
  zi_fixed = c("x_cen","y_cen"),                 # o c("x_cen","y_cen") si también los quieres fijos en ZI
  zi_candidates = zi_vars,
  start_in = "zi",                               # como tu gráfico: primero inflado, luego conteo
  dist = "negbin",
  criterion = "AIC",
  min_improve = 0                                # pon 2 si quieres “solo si mejora >2”
)

res$best_model

PyrChamois_DensityModel<-res$best_model
# save(PyrChamois_DensityModel, file="scripts/ModelResults/densityStepModel_PyrChamois_Europe.Rdata")

#2.3) Calibration plots  validation data
data_val<-data_val %>% 
  mutate(Pred_model = as.integer(predict(PyrChamois_DensityModel, data_val,type = "response"))) %>% #predecir zeroinfl
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

data_pos<-data_val %>% filter(dens>0) 
cal_tot<-bind_rows(
  cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
)

ggplot(cal_tot, aes(x = pred_mean, y = obs_mean)) +
  geom_point(aes(size = n), shape = 21, fill = "grey70") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
  facet_wrap(~Type, scales = "free") +
  scale_size_continuous(name = "N per bin", range = c(2, 6)) +
  labs(
    x = "Predicted mean",
    y = "Observed mean"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

# s_class <- cut2(data_val$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean <- as.matrix(tapply(data_val$Pred_Red, s_class, mean)) # calculating mean values for predicted
# y_mean <- as.matrix(tapply(data_val$dens, s_class, mean)) # calculating mean values for observed
# # corbr1<-cor.test(data_val$Pred_Red, data_val$dens)
# data_val_w0<-data_val %>%filter(dens_r_trans>0)
# s_class_w0 <- cut2(data_val_w0$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean_w0 <- as.matrix(tapply(data_val_w0$Pred_Red, s_class_w0, mean)) # calculating mean values for predicted
# y_mean_w0 <- as.matrix(tapply(data_val_w0$dens, s_class_w0, mean)) # calculating mean values for observed
# 
# cpE<-rbind(data.frame(s_mean=s_mean,   y_mean=y_mean, Type=as.factor("All data")),
#            data.frame(s_mean=s_mean_w0, y_mean=y_mean_w0, Type=as.factor("Data without 0")))
# ggplot(cpE, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=4)+ scale_colour_grey()+ scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="red", lty=2, cex=0.5)+
#   theme_bw()+
#   # xlim(0,max(s_mean_w0, max(na.omit(y_mean_w0))))+
#   # ylim(0,max(max(na.omit(y_mean_w0)), s_mean_w0))+
#   ylab("Observed Alpine ibex hunted bag density")+xlab("Predicted")+
#   theme(plot.title = element_text(hjust = 0.5))


#3.0) Spatial layer: 10x10km grid
grid_10km<-st_read(dsn="scripts/data/data_constrained_sp_20220902.gpkg", 
                   layer='grid_Isard') %>% dplyr::select(-Bioregion) %>% 
  rename(Bioregion=BR_T10, NUTS=NUT, Eucmean=Eucalyptus) %>% 
  mutate(NUTS=as.factor(NUTS), 
         Bioregion=as.factor(Bioregion),
         lc_50=as.numeric(lc_50),
         lc_62=as.numeric(lc_62),
         lc_72=as.numeric(lc_72), 
         lc_122=as.numeric(lc_122), 
         lc_153=as.numeric(lc_153), 
         Eucmean = replace_na(Eucmean, 0)) 
grid_10km<-grid_10km %>% 
  mutate(Eucmean = ifelse(grid_10km$Eucmean<0, 0, Eucmean)) %>% 
  mutate(x_cen=st_coordinates(st_centroid(grid_10km))[,1], y_cen=st_coordinates(st_centroid(grid_10km))[,2])
colnames(grid_10km)[32:65]<-gsub("lc_", "lc", colnames(grid_10km)[32:65])
colnames(grid_10km)[9:17]<-gsub("_", "_0", colnames(grid_10km)[9:17])

data %>% colnames()
means<-sp %>% as_tibble() %>% 
  select(15:69) %>%          # selecciona columnas
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
sdss<-sp %>% as_tibble() %>% 
  select(15:69) %>%          # selecciona columnas
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

library(modEvA) #MESS
summary(PyrChamois_DensityModel) #Solo ponemos las variables que selecciona el modelo que son sobre las que vamos a hacer la prediccion
data_<-data %>% as.data.frame() %>% dplyr::select(c(x_cen, y_cen, bio_09, alt, lc60, lc70))

#10km
grid_10km<-grid_10km %>% mutate(across(all_of(colnames(means)), ~ . - means[[cur_column()]])) %>% 
  mutate(across(all_of(colnames(means)), ~ . / sdss[[cur_column()]]))

grid_10km<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km)))
grid10km_df<-grid_10km %>% as.data.frame()  %>% dplyr::select("id_buena",
                                                              c(x_cen, y_cen, bio_09, alt, lc60, lc70))

mess_10km<-MESS(V=data_, P=grid10km_df, id.col=1)

grid_10km_<-grid_10km[c(-1:-73)]
mess_10km_<-mess_10km[c(-2:-7)]; names(mess_10km_)

grid10km_mess<-grid_10km %>%  full_join(mess_10km_, by="id_buena")

predRedDeer10km<-grid10km_mess %>% mutate(Pred_model = as.integer(predict(PyrChamois_DensityModel, grid_10km, type = "response"))) %>% mutate(Pred_Red = Pred_model/10000) 
st_write(predRedDeer10km, 
         dsn="scripts/ModelResults/predictions/predicciones_20260203.gpkg",
         layer="ZINB_stepModel_grid10km_PyrChamois", driver = "GPKG",append=TRUE)
# --------------------

# Reindeer #
# --------------------
data_sp<-st_read(dsn="scripts/data/predicciones_20220811.gpkg",
                 layer='Reindeer_spatial_v5')

sp<-data_sp %>%
  mutate(lc10=lc_10_sum/lc_10_count*100, lc11=lc_11_sum/lc_11_count*100, lc12=lc_12_sum/lc_12_count*100, lc20=lc_20_sum/lc_20_count*100, lc30=lc_30_sum/lc_30_count*100, lc40=lc_40_sum/lc_40_count*100, lc50=lc_50_sum/lc_50_count*100, lc60=lc_60_sum/lc_60_count*100, lc61=lc_61_sum/lc_61_count*100, lc62=lc_62_sum/lc_62_count*100, lc70=lc_70_sum/lc_70_count*100, lc71=lc_71_sum/lc_71_count*100, lc72=lc_72_sum/lc_72_count*100, lc80=lc_80_sum/lc_80_count*100, lc81=lc_81_sum/lc_81_count*100, lc82=lc_82_sum/lc_82_count*100, lc90=lc_90_sum/lc_90_count*100, lc100=lc_100_sum/lc_100_count*100, lc110=lc_110_sum/lc_110_count*100,lc120=lc_120_sum/lc_120_count*100,lc121=lc_121_sum/lc_121_count*100,lc122=lc_122_sum/lc_122_count*100,lc130=lc_130_sum/lc_130_count*100,lc140=lc_140_sum/lc_140_count*100,lc150=lc_150_sum/lc_150_count*100,lc152=lc_152_sum/lc_152_count*100,lc153=lc_153_sum/lc_153_count*100,lc160=lc_160_sum/lc_160_count*100,lc170=lc_170_sum/lc_170_count*100,lc180=lc_180_sum/lc_180_count*100,lc190=lc_190_sum/lc_190_count*100,lc200=lc_200_sum/lc_200_count*100,lc201=lc_201_sum/lc_201_count*100,lc202=lc_202_sum/lc_202_count*100,lc210=lc_210_sum/lc_210_count*100,lc220=lc_220_sum/lc_220_count*100) %>%
  dplyr::select(!matches("_cou") & !matches("_sum")) %>%
  #quito variables que no vaya a usar en el modelo
  dplyr::select(-c(excel_file,excel_path, shp_file, shp_path, #species, #me cargo variables que no necesito para el modelo
                   grow_mean, #la variable grow tiene datos raros, as? que la voy a eliminar de nuestro set de datos
                   lc50, lc62, lc72, lc81, lc82,lc121, lc170,#, #Me cargo variables lc que en el summary son todo 0
                   #BR_todasr1
                   "individualCount...39", "individualCount...51"
  )) %>% #quito las bioregiones porque no las vamos a usar para el modelo
  filter(dens<50) %>%  #tres poligonos de >Spain con areas muy muy muy muy small
  # filter(dens>0 & dens<50)
  select(locationID, locality, country, NUTS, BR_todasr1, area_km2, matches("harv"), dens, dens_r_trans,
         matches("bio"), alt, snow, sun, Eucmean, hfp,
         matches("lc"),
         x_cen, y_cen)

sp$Eucmean<-ifelse(sp$Eucmean<0, 0.01, sp$Eucmean) #Como para wild boar transformamos la variable de Eucmean para eliminar los valores negativos por valores de0 muy peque?os
sp<-sp[!is.na(sp$dens),] #Remove NA values from RV: it does not remove any value

#Modeling phase
#1) data
colnames(sp)[16:34]<-gsub("_mean", "", colnames(sp)[16:34])
vars_std <- c(colnames(sp %>% as_tibble())[16:70])

#remove NA values from the dataset
sp<-sp[!is.na(sp$bio_01),] #Remove NA values from the spset
sp<-sp[!is.na(sp$snow),] #Remove NA values from the spset
sp<-sp[!is.na(sp$hfp),] #Remove NA values from the spset
sp<-sp[!is.na(sp$alt),] #Remove NA values from the spset
sp<-sp[!is.na(sp$lc10),] #Remove NA values from the spset

data<-sp  %>% mutate(across(all_of(vars_std),~ as.numeric(scale(.))))
summary(data)

#1.1 Creamos data set de calibrado y de validacion (80% / 20%)
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]; st_geometry(data_tra)<-NULL; st_geometry(data_val)<-NULL

#2.0 VIF <2
#2.1 seleccionamos variables para el vif, para eso miramos como es la distribucion de las variables
# hist(data_tra$Eucmean)
# hist(data_tra$alt)
# hist(data_tra$snow)
#
# hist(data_tra$lc10)
# hist(data_tra$lc11)
# hist(data_tra$lc12) #
# hist(data_tra$lc20) #
# hist(data_tra$lc30)
# hist(data_tra$lc40)
# hist(data_tra$lc60)
# hist(data_tra$lc70)
# hist(data_tra$lc90)
# hist(data_tra$lc100)
# hist(data_tra$lc120)
# hist(data_tra$lc130)
#
# hist(data_tra$lc61) #
# hist(data_tra$lc71) #
# hist(data_tra$lc80) #
# hist(data_tra$lc110) #
# hist(data_tra$lc122) #
# hist(data_tra$lc140) #
# hist(data_tra$lc150) #
# hist(data_tra$lc152) #
# hist(data_tra$lc153) #
# hist(data_tra$lc160) #
# hist(data_tra$lc180) #
# hist(data_tra$lc190) #
# hist(data_tra$lc200) #
# hist(data_tra$lc201) #
# hist(data_tra$lc202) #
# hist(data_tra$lc210) #
# hist(data_tra$lc220) #

#2.1 hacemos el vif
variables<-data_tra %>% dplyr::select(matches("bio_"), alt, snow, Eucmean,lc10, lc11, #lc12,
                                      lc30, lc40, lc60, lc70, lc90,
                                      lc100, #lc120,
                                      lc130, x_cen, y_cen)
resultado.vif0<-usdm::vifstep(variables, th=2); resultado.vif0


#2.2) Density Models: binomial negativa necesita que sea un numero entero, por eso se ha transformado la variable densidad *10000
#ZINB
#modelo sin stepwise
count_vars <- c("lc30","lc40", "lc60", "lc130")
zi_vars <- count_vars
res <- stepwise_zeroinfl_alternate(
  data = data_tra,
  response = "dens_r_trans",
  count_fixed = c("x_cen","y_cen"),              # SIEMPRE dentro del count
  count_candidates = count_vars,
  zi_fixed = c("x_cen","y_cen"),                 # o c("x_cen","y_cen") si también los quieres fijos en ZI
  zi_candidates = zi_vars,
  start_in = "zi",                               # como tu gráfico: primero inflado, luego conteo
  dist = "negbin",
  criterion = "AIC",
  min_improve = 0                                # pon 2 si quieres “solo si mejora >2”
)

res$best_model

Reindeer_DensityModel<-res$best_model
# save(Reindeer_DensityModel, file="scripts/ModelResults/densityStepModel_Reindeer_poisson_Europe.Rdata")

#2.3) Calibration plots  validation data
data_val<-data_val %>%
mutate(Pred_model = as.integer(predict(Reindeer_DensityModel, data_val,type = "response"))) %>% #predecir zeroinfl
mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

data_pos<-data_val %>% filter(dens>0)
cal_tot<-bind_rows(
  cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0")
)

ggplot(cal_tot, aes(x = pred_mean, y = obs_mean)) +
  geom_point(aes(size = n), shape = 21, fill = "grey70") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
  facet_wrap(~Type, scales = "free") +
  scale_size_continuous(name = "N per bin", range = c(2, 6)) +
  labs(
    x = "Predicted mean",
    y = "Observed mean"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

# s_class <- cut2(data_val$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean <- as.matrix(tapply(data_val$Pred_Red, s_class, mean)) # calculating mean values for predicted
# y_mean <- as.matrix(tapply(data_val$dens, s_class, mean)) # calculating mean values for observed
# # corbr1<-cor.test(data_val$Pred_Red, data_val$dens)
# data_val_w0<-data_val %>%filter(dens_r_trans>0)
# s_class_w0 <- cut2(data_val_w0$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean_w0 <- as.matrix(tapply(data_val_w0$Pred_Red, s_class_w0, mean)) # calculating mean values for predicted
# y_mean_w0 <- as.matrix(tapply(data_val_w0$dens, s_class_w0, mean)) # calculating mean values for observed
#
# 
# cpE<-rbind(data.frame(s_mean=s_mean,   y_mean=y_mean, Type=as.factor("All data")),
#            data.frame(s_mean=s_mean_w0, y_mean=y_mean_w0, Type=as.factor("Data without 0")))
# ggplot(cpE, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=4)+ scale_colour_grey()+ scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="red", lty=2, cex=0.5)+
#   theme_bw()+
#   # xlim(0,max(s_mean_w0, max(na.omit(y_mean_w0))))+
#   # ylim(0,max(max(na.omit(y_mean_w0)), s_mean_w0))+
#   ylab("Observed Alpine ibex hunted bag density")+xlab("Predicted")+
#   theme(plot.title = element_text(hjust = 0.5))

#3.0) Spatial layer: 10x10km grid
grid_10km<-st_read(dsn="scripts/data/data_constrained_sp_20220902.gpkg", layer='grid_Reindeer') %>% dplyr::select(-Bioregion) %>%
  rename(Bioregion=BR_T10, NUTS=NUT, Eucmean=Eucalyptus) %>%
  mutate(NUTS=as.factor(NUTS),
         Bioregion=as.factor(Bioregion),
         lc_50=as.numeric(lc_50),
         lc_62=as.numeric(lc_62),
         lc_72=as.numeric(lc_72),
         lc_122=as.numeric(lc_122),
         lc_153=as.numeric(lc_153),
         Eucmean = replace_na(Eucmean, 0))
grid_10km<-grid_10km %>%
  mutate(Eucmean = ifelse(grid_10km$Eucmean<0, 0, Eucmean)) %>%
  mutate(x_cen=st_coordinates(st_centroid(grid_10km))[,1], y_cen=st_coordinates(st_centroid(grid_10km))[,2])
colnames(grid_10km)[32:65]<-gsub("lc_", "lc", colnames(grid_10km)[32:65])
colnames(grid_10km)[9:17]<-gsub("_", "_0", colnames(grid_10km)[9:17])


data %>% colnames()
means<-sp %>% as_tibble() %>%
  select(16:70) %>%          # selecciona columnas
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
sdss<-sp %>% as_tibble() %>%
  select(16:70) %>%          # selecciona columnas
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

library(modEvA) #MESS
summary(Reindeer_DensityModel) #Solo ponemos las variables que selecciona el modelo que son sobre las que vamos a hacer la prediccion
data_<-data %>% as.data.frame() %>% dplyr::select(c(x_cen, y_cen))

#10km
grid_10km<-grid_10km %>% mutate(across(all_of(colnames(means)), ~ . - means[[cur_column()]])) %>%
  mutate(across(all_of(colnames(means)), ~ . / sdss[[cur_column()]]))

grid_10km<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km)))
grid10km_df<-grid_10km %>% as.data.frame()  %>% dplyr::select("id_buena",
                                                              c(x_cen, y_cen))

mess_10km<-MESS(V=data_, P=grid10km_df, id.col=1)

grid_10km_<-grid_10km[c(-1:-73)]
mess_10km_<-mess_10km[c(-2:-3)]

grid10km_mess<-grid_10km %>%  full_join(mess_10km_, by="id_buena")

predRedDeer10km<-grid10km_mess %>% mutate(Pred_model = as.integer(predict(Reindeer_DensityModel, grid_10km, type = "response"))) %>% mutate(Pred_Red = Pred_model/10000)
st_write(predRedDeer10km,
         dsn="scripts/ModelResults/predictions/predicciones_20260203.gpkg",
         layer="ZINB_stepModel_grid10km_Reindeer", driver = "GPKG",append=TRUE)
# --------------------

# Sika deer #
# --------------------
data_sp<-st_read(dsn="scripts/data/predicciones_20220811.gpkg", 
                 layer='Sd_spatial_v5')

sp<-data_sp %>%
  mutate(lc10=lc_10_sum/lc_10_count*100, lc11=lc_11_sum/lc_11_count*100, lc12=lc_12_sum/lc_12_count*100, lc20=lc_20_sum/lc_20_count*100, lc30=lc_30_sum/lc_30_count*100, lc40=lc_40_sum/lc_40_count*100, lc50=lc_50_sum/lc_50_count*100, lc60=lc_60_sum/lc_60_count*100, lc61=lc_61_sum/lc_61_count*100, lc62=lc_62_sum/lc_62_count*100, lc70=lc_70_sum/lc_70_count*100, lc71=lc_71_sum/lc_71_count*100, lc72=lc_72_sum/lc_72_count*100, lc80=lc_80_sum/lc_80_count*100, lc81=lc_81_sum/lc_81_count*100, lc82=lc_82_sum/lc_82_count*100, lc90=lc_90_sum/lc_90_count*100, lc100=lc_100_sum/lc_100_count*100, lc110=lc_110_sum/lc_110_count*100,lc120=lc_120_sum/lc_120_count*100,lc121=lc_121_sum/lc_121_count*100,lc122=lc_122_sum/lc_122_count*100,lc130=lc_130_sum/lc_130_count*100,lc140=lc_140_sum/lc_140_count*100,lc150=lc_150_sum/lc_150_count*100,lc152=lc_152_sum/lc_152_count*100,lc153=lc_153_sum/lc_153_count*100,lc160=lc_160_sum/lc_160_count*100,lc170=lc_170_sum/lc_170_count*100,lc180=lc_180_sum/lc_180_count*100,lc190=lc_190_sum/lc_190_count*100,lc200=lc_200_sum/lc_200_count*100,lc201=lc_201_sum/lc_201_count*100,lc202=lc_202_sum/lc_202_count*100,lc210=lc_210_sum/lc_210_count*100,lc220=lc_220_sum/lc_220_count*100) %>% 
  dplyr::select(!matches("_cou") & !matches("_sum")) %>% 
  #quito variables que no vaya a usar en el modelo 
  dplyr::select(-c(excel_file,excel_path, shp_file, shp_path, #species, #me cargo variables que no necesito para el modelo
                   grow_mean, #la variable grow tiene datos raros, as? que la voy a eliminar de nuestro set de datos
                   lc50, lc62, lc72, lc81, lc82,lc121, lc170,#, #Me cargo variables lc que en el summary son todo 0
                   #BR_todasr1
                   "individualCount...39", "individualCount...51"
  )) %>% #quito las bioregiones porque no las vamos a usar para el modelo
  filter(dens<50) %>%  #tres poligonos de >Spain con areas muy muy muy muy small
  # filter(dens>0 & dens<50) 
  select(locationID, locality, country, NUTS, BR_todasr1, area_km2, matches("harv"), dens, dens_r_trans, 
         matches("bio"), alt, snow, sun, Eucmean, hfp, 
         matches("lc"), 
         x_cen, y_cen)

sp$Eucmean<-ifelse(sp$Eucmean<0, 0.01, sp$Eucmean) #Como para wild boar transformamos la variable de Eucmean para eliminar los valores negativos por valores de0 muy peque?os
sp<-sp[!is.na(sp$dens),] #Remove NA values from RV: it does not remove any value

#Modeling phase
#1) data
colnames(sp)[15:33]<-gsub("_mean", "", colnames(sp)[15:33])
vars_std <- c(colnames(sp %>% as_tibble())[15:69])

#remove NA values from the dataset
sp<-sp[!is.na(sp$bio_01),] #Remove NA values from the spset
sp<-sp[!is.na(sp$snow),] #Remove NA values from the spset
sp<-sp[!is.na(sp$hfp),] #Remove NA values from the spset
sp<-sp[!is.na(sp$alt),] #Remove NA values from the spset
sp<-sp[!is.na(sp$lc10),] #Remove NA values from the spset

data<-sp  %>% mutate(across(all_of(vars_std),~ as.numeric(scale(.))))
summary(data)

#1.1 Creamos data set de calibrado y de validacion (80% / 20%)
set.seed(500)
s_model <- sample(nrow(data), nrow(data)*0.8) ; data_tra <- data[s_model,]; data_val <- data[-s_model,]; st_geometry(data_tra)<-NULL; st_geometry(data_val)<-NULL

#2.0 VIF <2
#2.1 seleccionamos variables para el vif, para eso miramos como es la distribucion de las variables
# hist(data_tra$Eucmean)
# hist(data_tra$alt) 
# hist(data_tra$snow)
# 
# hist(data_tra$lc10)
# hist(data_tra$lc11)
# hist(data_tra$lc12) #
# hist(data_tra$lc20) #
# hist(data_tra$lc30)
# hist(data_tra$lc40) 
# hist(data_tra$lc60) 
# hist(data_tra$lc70)
# hist(data_tra$lc90) 
# hist(data_tra$lc100)
# hist(data_tra$lc120) 
# hist(data_tra$lc130)
# 
# hist(data_tra$lc61) #
# hist(data_tra$lc71) #
# hist(data_tra$lc80) #
# hist(data_tra$lc110) #
# hist(data_tra$lc122) #
# hist(data_tra$lc140) #
# hist(data_tra$lc150) #
# hist(data_tra$lc152) #
# hist(data_tra$lc153) #
# hist(data_tra$lc160) #
# hist(data_tra$lc180) #
# hist(data_tra$lc190) #
# hist(data_tra$lc200) #
# hist(data_tra$lc201) #
# hist(data_tra$lc202) #
# hist(data_tra$lc210) #
# hist(data_tra$lc220) #

#2.1 hacemos el vif
variables<-data_tra %>% dplyr::select(matches("bio_"), alt, snow, Eucmean,lc10, lc11, lc12, lc30, lc40, lc60, lc70, lc90, 
                                      lc100, #lc120, 
                                      lc130, x_cen, y_cen)
resultado.vif0<-usdm::vifstep(variables, th=2); resultado.vif0


#2.2) Density Models: binomial negativa necesita que sea un numero entero, por eso se ha transformado la variable densidad *10000
#ZINB
count_vars <- c("bio_03", "bio_10", "bio_15", "lc10","lc30","lc40","lc60","lc70", "lc90","lc100","lc130")
zi_vars <- count_vars
res <- stepwise_zeroinfl_alternate(
  data = data_tra,
  response = "dens_r_trans",
  count_fixed = c("x_cen","y_cen"),              # SIEMPRE dentro del count
  count_candidates = count_vars,
  zi_fixed = c("x_cen","y_cen"),                 # o c("x_cen","y_cen") si también los quieres fijos en ZI
  zi_candidates = zi_vars,
  start_in = "zi",                               # como tu gráfico: primero inflado, luego conteo
  dist = "negbin",
  criterion = "AIC",
  min_improve = 0                                # pon 2 si quieres “solo si mejora >2”
)

res$best_model

SikaDeer_DensityModel<-res$best_model
# save(SikaDeer_DensityModel, file="scripts/ModelResults/densityStepModel_SikaDeer_Europe.Rdata")

#2.3) Calibration plots  validation data
data_val<-data_val %>% 
  mutate(Pred_model = as.integer(predict(SikaDeer_DensityModel, data_val,type = "response"))) %>% #predecir zeroinfl
  mutate(Pred_Red = Pred_model/10000) ; range(unique(data_val$Pred_Red))

data_pos<-data_val %>% filter(dens>0) 
cal_tot<-bind_rows(
  cal_data(pred = data_val$Pred_Red, obs=data_val$dens) %>% mutate(Type = "All data"),
  cal_data(pred = data_pos$Pred_Red, obs=data_pos$dens) %>% mutate(Type = "Data without 0 ")
)
ggplot(cal_tot, aes(x = pred_mean, y = obs_mean)) +
  geom_point(aes(size = n), shape = 21, fill = "grey70") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
  facet_wrap(~Type, scales = "free") +
  scale_size_continuous(name = "N per bin", range = c(2, 6)) +
  labs(
    x = "Predicted mean",
    y = "Observed mean"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )


# s_class <- cut2(data_val$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean <- as.matrix(tapply(data_val$Pred_Red, s_class, mean)) # calculating mean values for predicted
# y_mean <- as.matrix(tapply(data_val$dens, s_class, mean)) # calculating mean values for observed
# # corbr1<-cor.test(data_val$Pred_Red, data_val$dens)
# data_val_w0<-data_val %>%filter(dens_r_trans>0)
# s_class_w0 <- cut2(data_val_w0$Pred_Red, g=9) # defining bins (percentiles) on the predicted
# s_mean_w0 <- as.matrix(tapply(data_val_w0$Pred_Red, s_class_w0, mean)) # calculating mean values for predicted
# y_mean_w0 <- as.matrix(tapply(data_val_w0$dens, s_class_w0, mean)) # calculating mean values for observed
# 
# cpE<-rbind(data.frame(s_mean=s_mean,   y_mean=y_mean, Type=as.factor("All data")),
#            data.frame(s_mean=s_mean_w0, y_mean=y_mean_w0, Type=as.factor("Data without 0")))
# ggplot(cpE, aes(x=s_mean, y=y_mean, col=Type, shape=Type))+
#   geom_point(size=4)+ scale_colour_grey()+ scale_shape_manual(values=c(16, 17))+
#   geom_abline(col="red", lty=2, cex=0.5)+
#   theme_bw()+
#   # xlim(0,max(s_mean_w0, max(na.omit(y_mean_w0))))+
#   # ylim(0,max(max(na.omit(y_mean_w0)), s_mean_w0))+
#   ylab("Observed Alpine ibex hunted bag density")+xlab("Predicted")+
#   theme(plot.title = element_text(hjust = 0.5))


#3.0) Spatial layer: 10x10km grid
grid_10km<-st_read(dsn="scripts/data/data_constrained_sp_20220902.gpkg", layer='grid_Sd') %>% dplyr::select(-Bioregion) %>% 
  rename(Bioregion=BR_T10, NUTS=NUT, Eucmean=Eucalyptus) %>% 
  mutate(NUTS=as.factor(NUTS), 
         Bioregion=as.factor(Bioregion),
         lc_50=as.numeric(lc_50),
         lc_62=as.numeric(lc_62),
         lc_72=as.numeric(lc_72), 
         lc_122=as.numeric(lc_122), 
         lc_153=as.numeric(lc_153), 
         Eucmean = replace_na(Eucmean, 0), 
         Eucmean = ifelse(Eucmean<0, 0, Eucmean)) 
grid_10km<-grid_10km %>% 
  mutate(x_cen=st_coordinates(st_centroid(grid_10km))[,1], y_cen=st_coordinates(st_centroid(grid_10km))[,2])
colnames(grid_10km)[32:65]<-gsub("lc_", "lc", colnames(grid_10km)[32:65])
colnames(grid_10km)[9:17]<-gsub("_", "_0", colnames(grid_10km)[9:17])

data %>% colnames()
means<-sp %>% as_tibble() %>% 
  select(15:69) %>%          # selecciona columnas
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
sdss<-sp %>% as_tibble() %>% 
  select(15:69) %>%          # selecciona columnas
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

library(modEvA) #MESS
summary(SikaDeer_DensityModel) #Solo ponemos las variables que selecciona el modelo que son sobre las que vamos a hacer la prediccion
data_<-data %>% as.data.frame() %>% dplyr::select(c(x_cen, y_cen, bio_03, bio_10, bio_15, lc60))

#10km
grid_10km<-grid_10km %>% mutate(across(all_of(colnames(means)), ~ . - means[[cur_column()]])) %>% 
  mutate(across(all_of(colnames(means)), ~ . / sdss[[cur_column()]]))

grid_10km<-grid_10km %>% mutate(id_buena=seq(1:nrow(grid_10km)))
grid10km_df<-grid_10km %>% as.data.frame()  %>% dplyr::select("id_buena",
                                                              c(x_cen, y_cen, bio_03, bio_10, bio_15, lc60))

mess_10km<-MESS(V=data_, P=grid10km_df, id.col=1)

grid_10km_<-grid_10km[c(-1:-73)]
mess_10km_<-mess_10km[c(-2:-7)];names(mess_10km_)

grid10km_mess<-grid_10km %>%  full_join(mess_10km_, by="id_buena")

predRedDeer10km<-grid10km_mess %>% mutate(Pred_model = as.integer(predict(SikaDeer_DensityModel, grid_10km, type = "response"))) %>% mutate(Pred_Red = Pred_model/10000) 
st_write(predRedDeer10km, 
         dsn="scripts/ModelResults/predictions/predicciones_20260203.gpkg",
         layer="ZINB_stepModel_grid10km_SikaDeer", driver = "GPKG",append=TRUE)
# --------------------