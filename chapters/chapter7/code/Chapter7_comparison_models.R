# -----------------------------------------
# From the work: Illanas, et al (in prep). What a camera trap network reveals 
# about wild ungulate communities in Spain: large-scale monitoring using 
# hierarchical N-mixture models
# -----------------------------------------
# script: comparison_models.R
# This script was used for making the multispecies models
# author: Sonia Illanas
# v1.0.0
# date last modification: 2026/03/17
# -----------------------------------------

# 0. Setup. Load required packages and data
# ------------------------------------------
# load packages
library(tidyverse)   # v2.0.0
library(spAbundance) # v0.2.2
library(brms)        # v2.23.0
library(bayesplot)   # v1.13.0
library(tidybayes)   # v3.0.7
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)    # v0.8.0
install_cmdstan(cores = 4) # select number of cores for a faster compilation
check_cmdstan_toolchain() # check that it's properly installed. If okay, brm (backend = "cmdstanr") can be used

rm(list=ls())
set.seed(1993)

# load data (needed for standardization validation covariates values)
load(file = "~/Documentos/PhD/5_Estancia/2_WIP/0_DataPrepared/modelingData/data.model.Rdata")
# str(data.model)

scale_covs <- function(df) {
  df %>%
    as_tibble() %>%
    mutate(
      # use the covariates used for modeling 
      precip  = (precip - mean(data.model$abund.covs$precip)) / sd(data.model$abund.covs$precip),
      precip2 = (precip2 - mean(data.model$abund.covs$precip2)) / sd(data.model$abund.covs$precip2),
      elev    = (elev - mean(data.model$abund.covs$elev)) / sd(data.model$abund.covs$elev),
      elev2   = (elev2 - mean(data.model$abund.covs$elev2, na.rm=TRUE)) / sd(data.model$abund.covs$elev2, na.rm=TRUE),
      slope   = (slope - mean(data.model$abund.covs$slope)) / sd(data.model$abund.covs$slope),
      forest  = (forest * 100 - mean(data.model$abund.covs$forest)) / sd(data.model$abund.covs$forest),
      shrb    = (shrb * 100 - mean(data.model$abund.covs$shrb)) / sd(data.model$abund.covs$shrb)
    )
}

# prepare comparison data (long format)
data_prep <- function(data, type = c("hunt", "gbif")) {
  data %>%
    as_tibble() %>%
    select(-any_of("geom")) %>%
    pivot_longer(
      cols = contains(c(".med", ".se")),
      names_to = c("species", ".value"),
      names_sep = "\\."
    ) %>%
    rename(lambda_pred = med, se_lambda = se) %>%
    mutate(
      value_obs = case_when(
        type == "hunt" ~ case_when(
          species == "cp" ~ Cp.max, species == "cc" ~ Cc.max, species == "ce" ~ Ce.max,
          species == "dd" ~ Dd.max, species == "oa" ~ Oa.max, species == "rp" ~ Rp.max,
          species == "wb" ~ Ss.max, TRUE ~ 0
        ) / area_km2,
        type == "gbif" ~ case_when(
          TotP > 0 & species == "cp" ~ cp / TotP,
          TotP > 0 & species == "cc" ~ cc / TotP,
          TotP > 0 & species == "ce" ~ ce / TotP,
          TotP > 0 & species == "dd" ~ dd / TotP,
          TotP > 0 & species == "oa" ~ oa / TotP,
          TotP > 0 & species == "rp" ~ rp / TotP,
          TotP > 0 & species == "wb" ~ ss / TotP,
          TRUE ~ 0 # if TotP == 0 or does not match species = 0
        ),
        TRUE ~ 0
      ),
      log_lambda = log(lambda_pred + 0.001),
      se_log = se_lambda / (lambda_pred + 0.001)#,
      # species = factor(species, levels = order_species)
    ) %>%
    select(!matches("low|high|mean|max"))
}


# load model results for making predictions
load(file = "~/Documentos/PhD/5_Estancia/2_WIP/2_Results/MSNM_20260312_spAbundance/msNmix_ms07_2_largerSamples.Rdata")


# internal comparisons
# -----------------------
load (file = "~/Documentos/PhD/5_Estancia/2_WIP/0_DataPrepared/modelingData/data.internal.Rdata")

hb.pred <- scale_covs(hb)

hb.pred<-hb.pred[hb.pred$precip<3 & hb.pred$precip>-3,]
hb.pred<-hb.pred[hb.pred$precip2<3 & hb.pred$precip2>-3,]
hb.pred<-hb.pred[hb.pred$elev<3 & hb.pred$elev>-3,]
hb.pred<-hb.pred[hb.pred$elev2<3 & hb.pred$elev2>-3,]
hb.pred<-hb.pred[hb.pred$slope<3 & hb.pred$slope>-3,]
hb.pred<-hb.pred[hb.pred$forest<3 & hb.pred$forest>-3,]
hb.pred<-hb.pred[hb.pred$shrb<3 & hb.pred$shrb>-3,]

summary(hb.pred)
colnames(hb.pred)

X.0 <- model.matrix(~ elev + elev2 + precip + precip2 + slope + forest + shrb,
                    data = hb.pred)

out.pred<-predict(out.ms07_2, X.0)
out.pred$mu.0.samples %>% dim()


mu.0.quants <- apply(out.pred$mu.0.samples, c(2,3), quantile, prob = c(0.025, 0.5, 0.975))
mu.0.quants %>% dim()

hb.pred<-cbind(hb.pred, t(as.data.frame(mu.0.quants[2, , ])))
rownames(hb.pred) <- NULL
sp.codes<-dimnames(data.model$y)[[1]]
colnames(hb.pred)[26:32]<-c("cp.med", "cc.med", "ce.med", "dd.med", "oa.med", "rp.med", "wb.med")

se<-apply(out.pred$mu.0.samples, c(2,3), sd)
hb.pred<-cbind(hb.pred, t(as.data.frame(se)))
rownames(hb.pred) <- NULL
colnames(hb.pred)[33:39]<-c("cp.se", "cc.se", "ce.se", "dd.se", "oa.se", "rp.se", "wb.se")

# hunting data
df_hunt <- hb.pred %>% data_prep(type = "hunt") 
df_hunt<-df_hunt %>% filter(value_obs >0)
df_hunt$log_obs <- log(df_hunt$value_obs + 0.001)
df_hunt$obs_z <- (df_hunt$log_obs - mean(df_hunt$log_obs)) / sd(df_hunt$log_obs)

m_hunt <- brm(
  log(lambda_pred) | se(se_log, sigma = TRUE) ~ obs_z + (1 + obs_z || species),
  data = df_hunt, family = gaussian(),
  prior = c(prior(normal(0, 2), class = "b"), 
            prior(exponential(1), class = "sd")),
  control = list(adapt_delta = 0.99, 
                 max_treedepth = 15), 
  cores = 12, 
  threads = threading(4),
  chains = 3,
  thin = 1,
  iter = 5000, 
  warmup = 1000, 
  backend = "cmdstanr"
)
# save(m_hunt, file = "m_hunt_WITTHOUT0.Rdata")

# gbif data
df_gbif <- hb.pred %>% data_prep(type = "gbif") 

df_gbif$obs_z <- (df_gbif$value_obs - mean(df_gbif$value_obs, na.rm=TRUE)) / sd(df_gbif$value_obs, na.rm=TRUE)

m_gbif <- brm(
  log(lambda_pred) | se(se_log, sigma = TRUE) ~ obs_z + (1 + obs_z || species),
  data = df_gbif, family = gaussian(),
  prior = c(prior(normal(0, 2), class = "b"), 
            prior(exponential(1), class = "sd")),
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15), 
  cores = 12, 
  threads = threading(4),
  chains = 3,
  thin = 1,
  iter = 5000, 
  warmup = 1000, 
  backend = "cmdstanr"
)
# save(m_gbif, file = "m_gbif.Rdata")

performance::r2_bayes(m_hunt)
performance::r2_bayes(m_gbif)

# -------------------------------

# Prepare data - external validation 20 km
# 20km
# -------------------------------
load(file = "~/Documentos/PhD/5_Estancia/2_WIP/0_DataPrepared/modelingData/data.external20.Rdata")

hb.pred <- scale_covs(hb20)

hb.pred<-hb.pred[hb.pred$precip<3 & hb.pred$precip>-3,]
hb.pred<-hb.pred[hb.pred$precip2<3 & hb.pred$precip2>-3,]
hb.pred<-hb.pred[hb.pred$elev<3 & hb.pred$elev>-3,]
hb.pred<-hb.pred[hb.pred$elev2<3 & hb.pred$elev2>-3,]
hb.pred<-hb.pred[hb.pred$slope<3 & hb.pred$slope>-3,]
hb.pred<-hb.pred[hb.pred$forest<3 & hb.pred$forest>-3,]
hb.pred<-hb.pred[hb.pred$shrb<3 & hb.pred$shrb>-3,]


summary(hb.pred)
colnames(hb.pred)

rand<-sample(nrow(hb.pred), nrow(hb.pred)*0.3) 
hb.pred.s <- hb.pred[rand,]

X.0 <- model.matrix(~ elev + elev2 + precip + precip2 + slope + forest + shrb,
                    data = hb.pred.s)

out.pred<-predict(out.ms07_2, X.0)
out.pred$mu.0.samples %>% dim()

mu.0.quants <- apply(out.pred$mu.0.samples, c(2,3), quantile, prob = c(0.025, 0.5, 0.975))
mu.0.quants %>% dim()

hb.pred.s<-cbind(hb.pred.s, t(as.data.frame(mu.0.quants[2, , ])))
rownames(hb.pred.s) <- NULL
sp.codes
colnames(hb.pred.s)[26:32]<-c("cp.med", "cc.med", "ce.med", "dd.med", "oa.med", "rp.med", "wb.med")

se<-apply(out.pred$mu.0.samples, c(2,3), sd)
hb.pred.s<-cbind(hb.pred.s, t(as.data.frame(se)))
rownames(hb.pred.s) <- NULL
colnames(hb.pred.s)[33:39]<-c("cp.se", "cc.se", "ce.se", "dd.se", "oa.se", "rp.se", "wb.se")

# hy data
df_hunt20 <- hb.pred.s %>% data_prep(type = "hunt")
df_hunt20<-df_hunt20 %>% filter(value_obs >0)
df_hunt20$log_obs <- log(df_hunt20$value_obs + 0.001)
df_hunt20$obs_z <- (df_hunt20$log_obs - mean(df_hunt20$log_obs)) / sd(df_hunt20$log_obs)

m_hunt20 <- brm(
  log(lambda_pred) | se(se_log, sigma = TRUE) ~ obs_z + (1 + obs_z || species),
  data = df_hunt20, 
  family = gaussian(),
  prior = c(prior(normal(0, 2), class = "b"), 
            prior(exponential(1), class = "sd")),
  control = list(adapt_delta = 0.99, 
                 max_treedepth = 15), 
  cores = 12, 
  threads = threading(6),
  chains = 3,
  thin = 1,
  iter = 5000, 
  warmup = 1000, 
  backend = "cmdstanr"
)
# plot(m_hunt20)
# save(m_hunt20, file = "m_hunt20WITHOUT0.Rdata")

# gbif data
df_gbif20 <- hb.pred.s %>% data_prep(type = "gbif") 
df_gbif20$obs_z <- (df_gbif20$value_obs - mean(df_gbif20$value_obs, na.rm=TRUE)) / sd(df_gbif20$value_obs, na.rm=TRUE)

m_gbif20 <- brm(
  log(lambda_pred) | se(se_log, sigma = TRUE) ~ obs_z + (1 + obs_z || species),
  data = df_gbif20, family = gaussian(),
  prior = c(prior(normal(0, 2), class = "b"), 
            prior(exponential(1), class = "sd")),
  control = list(adapt_delta = 0.99, 
                 max_treedepth = 15),
  cores = 12, 
  threads = threading(4),
  chains = 3,
  thin = 1,
  iter = 5000, 
  warmup = 1000, 
  backend = "cmdstanr"
)
# save(m_gbif20, file = "m_gbif20.Rdata")


performance::r2_bayes(m_hunt20)
performance::r2_bayes(m_gbif20)

# -------------------------------


# Prepare data - external validation 20 km
# 50km
# -------------------------------
load(file = "~/Documentos/PhD/5_Estancia/2_WIP/0_DataPrepared/modelingData/data.external50.Rdata")

hb.pred <- scale_covs(hb50)

hb.pred<-hb.pred[hb.pred$precip<3 & hb.pred$precip>-3,]
hb.pred<-hb.pred[hb.pred$precip2<3 & hb.pred$precip2>-3,]
hb.pred<-hb.pred[hb.pred$elev<3 & hb.pred$elev>-3,]
hb.pred<-hb.pred[hb.pred$elev2<3 & hb.pred$elev2>-3,]
hb.pred<-hb.pred[hb.pred$slope<3 & hb.pred$slope>-3,]
hb.pred<-hb.pred[hb.pred$forest<3 & hb.pred$forest>-3,]
hb.pred<-hb.pred[hb.pred$shrb<3 & hb.pred$shrb>-3,]

summary(hb.pred)
colnames(hb.pred)

rand<-sample(nrow(hb.pred), nrow(hb.pred)*0.3) 
hb.pred.s <- hb.pred[rand,]

X.0 <- model.matrix(~ elev + elev2 + precip + precip2 + slope + forest + shrb,
                    data = hb.pred.s)


# out.ms07_2$beta.samples %>% colnames()
out.pred<-predict(out.ms07_2, X.0)
out.pred$mu.0.samples %>% dim()

mu.0.quants <- apply(out.pred$mu.0.samples, c(2,3), quantile, prob = c(0.025, 0.5, 0.975))
mu.0.quants %>% dim()

hb.pred.s<-cbind(hb.pred.s, t(as.data.frame(mu.0.quants[2, , ])))
rownames(hb.pred.s) <- NULL
sp.codes
colnames(hb.pred.s)[26:32]<-c("cp.med", "cc.med", "ce.med", "dd.med", "oa.med", "rp.med", "wb.med")

se<-apply(out.pred$mu.0.samples, c(2,3), sd)
hb.pred.s<-cbind(hb.pred.s, t(as.data.frame(se)))
rownames(hb.pred.s) <- NULL
colnames(hb.pred.s)[33:39]<-c("cp.se", "cc.se", "ce.se", "dd.se", "oa.se", "rp.se", "wb.se")

# hy data
df_hunt50 <- hb.pred.s %>% data_prep(type = "hunt")
df_hunt50<-df_hunt50 %>% filter(value_obs >0)
df_hunt50$log_obs <- log(df_hunt50$value_obs + 0.001)
df_hunt50$obs_z <- (df_hunt50$log_obs - mean(df_hunt50$log_obs, na.rm = TRUE)) / sd(df_hunt50$log_obs, na.rm = TRUE)

m_hunt50 <- brm(
  log(lambda_pred) | se(se_log, sigma = TRUE) ~ obs_z + (1 + obs_z || species),
  data = df_hunt50, 
  family = gaussian(),
  prior = c(prior(normal(0, 2), class = "b"), 
            prior(exponential(1), class = "sd")),
  control = list(adapt_delta = 0.99, 
                 max_treedepth = 15),
  cores = 12, 
  threads = threading(6),
  chains = 3,
  thin = 1,
  iter = 5000, 
  warmup = 1000, 
  backend = "cmdstanr", 
  save_pars = save_pars(all = TRUE)
)
# save(m_hunt50, file = "m_hunt50.Rdata")


# gbif data
df_gbif50 <- hb.pred.s %>% data_prep(type = "gbif") 
df_gbif50$obs_z <- (df_gbif50$value_obs - mean(df_gbif50$value_obs, na.rm=TRUE)) / sd(df_gbif50$value_obs, na.rm=TRUE)

m_gbif50 <- brm(
  log(lambda_pred) | se(se_log, sigma = TRUE) ~ obs_z + (1 + obs_z || species),
  data = df_gbif50, family = gaussian(),
  prior = c(prior(normal(0, 2), class = "b"), 
            prior(exponential(1), class = "sd")),
  control = list(adapt_delta = 0.99, 
                 max_treedepth = 15),
  cores = 12, 
  threads = threading(6),
  chains = 3,
  thin = 1,
  iter = 5000, 
  warmup = 1000,
  backend = "cmdstanr", 
  save_pars = save_pars(all = TRUE)
)
# save(m_gbif50, file = "m_gbif50.Rdata")

performance::r2_bayes(m_hunt50)
performance::r2_bayes(m_gbif50)
# -------------------------------
