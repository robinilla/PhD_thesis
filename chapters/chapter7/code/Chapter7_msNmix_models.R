# -----------------------------------------
# From the work: Illanas, et al (in prep). 
# Large-scale modeling of wild ungulate relative abundance 
# from a national camera trap network
# -----------------------------------------
# script: msNmix_models.R
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
library(flextable)   # v0.9.11
rm(list=ls())

# load data
load(file = "~/modelingData/data.model.Rdata")
str(data.model)

# prepare data for modeling 
# standardize covariates of latent and observational process
data.sp <- list (y = data.model$y, 
                 abund.covs = scale(data.model$abund.covs)[,1:8], 
                 det.covs = list(
                    height.CT = ((data.model$det.covs$height.CT) - mean(data.model$det.covs$height.CT[,1])) /sd(data.model$det.covs$height.CT[,1]), 
                    slope25m = ((data.model$det.covs$slope25m) - mean(data.model$det.covs$slope25m[,1])) /sd(data.model$det.covs$slope25m[,1]),
                    week_s = ((data.model$det.covs$week_s) - mean(data.model$det.covs$week_s[,1])) /sd(data.model$det.covs$week_s[,1]),
                    week2 = ((data.model$det.covs$week2) - mean(data.model$det.covs$week2[,1])) /sd(data.model$det.covs$week2[,1]),
                    effort = ((data.model$det.covs$effort) - mean(data.model$det.covs$effort[,1])) /sd(data.model$det.covs$effort[,1])
                    )
                  )

inits <- list(alpha.comm = 0,
              beta.comm = 0,
              beta = 0, 
              alpha = 0, 
              tau.sq.alpha = 0.5, 
              tau.sq.beta = 1,
              kappa = 1,
              sigma.sq.mu = 0.5,
              N = apply(data.one.sp$y, c(1, 2), max, na.rm = TRUE))

priors <- list(alpha.comm.normal = list(mean = -1.5, var = 1.5), 
               beta.comm.normal = list(mean = 0, var = 2.72),
               tau.sq.beta.ig = list(a = 0.1, b = 0.1),
               tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
               kappa.unif = c(0, 50), 
               sigma.sq.mu.ig = list(3, 2)                
)

n.batch<- 5000
batch.length<-50 
n.burn <- n.batch*batch.length*0.4
n.thin <-  10 
n.chains <- 3
(n.batch*batch.length-n.burn)/n.thin * n.chains
tuning <- list(beta = 0.5, alpha = 0.5, kappa = 0.5, beta.star = 0.5)
# accept.rate = 0.43 by default, so we do not specify it.

# Run msNmix models
# ----------------------------
out.micro<-
  msNMix(abund.formula = ~ forest + shrb,
         det.formula = ~ height.CT + slope25m + week_s + week2 +  effort,
         data = data.one.sp,
         inits = inits,
         priors = priors,
         n.batch = n.batch,
         batch.length = batch.length,
         tuning = tuning,
         n.omp.threads = 8,
         n.report = 1000,
         family = 'Poisson',
         verbose = TRUE,
         n.burn = n.burn,
         n.thin = n.thin,
         n.chains = n.chains)
summary(out.micro)
# save(out.micro, file = "msNmix_micro.Rdata")

out.macro<-
  msNMix(abund.formula = ~ elev + elev2 + precip + precip2 + slope,
         det.formula = ~ height.CT + slope25m + week_s + week2 +  effort,
         data = data.one.sp,
         inits = inits,
         priors = priors,
         n.batch = n.batch,
         batch.length = batch.length,
         tuning = tuning,
         n.omp.threads = 14,
         n.report = 1000,
         family = 'Poisson',
         verbose = TRUE,
         n.burn = n.burn,
         n.thin = n.thin,
         n.chains = n.chains)
summary(out.macro)
# save(out.macro, file = "msNmix_macro.Rdata")

out.comb<-
  msNMix(abund.formula = ~ elev + elev2 + precip + precip2 + slope + forest + shrb,
         det.formula = ~ height.CT + slope25m + week_s + week2 +  effort,
         data = data.sp,
         inits = inits,
         priors = priors,
         n.batch = n.batch,
         batch.length = batch.length,
         tuning = tuning,
         n.omp.threads = 8,
         n.report = 1000,
         family = 'Poisson',
         verbose = TRUE,
         n.burn = n.burn,
         n.thin = n.thin,
         n.chains = n.chains)
summary(out.comb)
# save(out.comb, file = "msNmix_out.comb.Rdata")
