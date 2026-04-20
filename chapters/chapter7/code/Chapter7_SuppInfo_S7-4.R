# -----------------------------------------
# From the work: Illanas, et al (in prep). 
# Large-scale modeling of wild ungulate relative abundance 
# from a national camera trap network
# -----------------------------------------
# script: SuppInfo_S4.R
# This script takes multispecies model output 
# and replicate supplementary information S4
# author: Sonia Illanas
# v1.0.0
# date last modification: 2026/03/17
# -----------------------------------------

# 0. Setup. Load required packages and data
# ------------------------------------------
# load packages
library(tidyverse)   # v2.0.0
library(spAbundance) # v0.2.2
library(coda)        # v0.19-4.1
library(scales)      # v1.4.0
library(ggh4x)       # v0.3.1

load(file = "~/Documentos/PhD/5_Estancia/2_WIP/2_Results/MSNM_20260312_spAbundance/msNmix_ms07_2_largerSamples.Rdata")
out <- out.ms07_2

levels.sp <- c("C. pyrenaica",  "C. capreolus",  "C. elaphus", "D. dama", "O. aries", "R. pyrenaica", "S. scrofa")

# S4 - TRACEPLOTS, RHAT, ESS
# -------------------------
# ABUND - BETA
n_chains <- 3
n_samples <- nrow(out$beta.samples) / n_chains

samples_df <- as.data.frame(out$beta.samples) %>%
  mutate(iteration = rep(1:n_samples, n_chains),
         chain = factor(rep(1:n_chains, each = n_samples))) %>%
  pivot_longer(-c(iteration, chain), names_to = "parameter", values_to = "value") %>% 
  tidyr::extract(parameter, 
                 into = c("parameter", "species"), 
                 regex = "(.*)-(.*)") 

samples_df <- samples_df %>% 
  mutate(species = paste(substr(samples_df$species, 1, 1), 
                         ". ", 
                         str_extract(samples_df$species, "(?<=\\s).*"), 
                         sep = "")) %>% 
  mutate(species = factor(samples_df$species, levels = levels.sp))

# 2. Extract Rhat & ESS 
# transformt to mcmc.list to allow coda understand the chains
stats <- data.frame(
  parameter = colnames(out$beta.samples),
  rhat = out$rhat$beta,
  ess = out$ESS$beta ) %>%
  mutate(label = paste0("Rhat: ", round(rhat, 4), "\nESS: ", round(ess, 0))) %>%
  tidyr::extract(parameter, 
                 into = c("parameter", "species"), 
                 regex = "(.*)-(.*)") %>%
  mutate(species = paste(substr(stats$species, 1, 1), ". ", 
                         str_extract(stats$species, "(?<=\\s).*"), 
                         sep="")) %>% 
  mutate(species = factor(samples_df$species, levels = levels.sp))


trace_plot<-function(data_samples){
  ggplot(data_samples, aes(x = iteration, y = value, color = chain)) +
  geom_line(alpha = 0.6, linewidth = 0.3) +
  geom_text(data = stats, 
            aes(x = 500, y = Inf, label = label),
            inherit.aes = FALSE, 
            hjust = 0, vjust = 2.9,
            size = 2.2, 
            fontface = "bold.italic", 
            color = "white") +
  theme_minimal() +
  facet_grid2(parameter ~ species, 
              scales = "free_y", 
              independent = "y", 
              switch = NULL) +
  theme(
    axis.title = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold.italic"),
    strip.text.y = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 6)
  ) +
  labs(x = "Samples per chain", 
       y = "Parameter value", 
       color = "Chain")
}

traceAb <- trace_plot(samples_df)
# trace
# ggsave(file = "/FigureS4.1_TraceRhatESS_Abun.png",
#        traceAb,
#        width = 25.2,
#        height = 17.6,
#        units = "cm",
#        dpi = 300)


# DET - alpha
# convert the output to mcmc.list
n_samples <- nrow(out$alpha.samples) / n_chains
samples_df <- as.data.frame(out$alpha.samples) %>%
  mutate(iteration = rep(1:n_samples, n_chains),
         chain = factor(rep(1:n_chains, each = n_samples))
         ) %>%
  pivot_longer(-c(iteration, chain), names_to = "parameter", values_to = "value") %>% 
  tidyr::extract(parameter, 
                 into = c("parameter", "species"), 
                 regex = "(.*)-(.*)")

samples_df <- samples_df %>%
  mutate(species = paste(substr(samples_df$species, 1, 1), 
                         ". ", 
                         str_extract(samples_df$species, "(?<=\\s).*"), 
                         sep=""))  %>% 
  mutate(species = factor(samples_df$species, levels = levels.sp))


samples_df$parameter<-gsub("height.CT", "canopy.h", samples_df$parameter)
samples_df$parameter<-gsub("week_s", "week", samples_df$parameter)
samples_df$parameter<-factor(samples_df$parameter, 
                             levels = c("(Intercept)", "slope25m", "canopy.h", "week", "week2", "effort"    
                             ))

# 2. Extract Rhat & ESS 
stats <- data.frame(
  parameter = colnames(out$alpha.samples),
  rhat = out$rhat$alpha,
  ess = out$ESS$alpha
) %>%
  mutate(label = paste0("Rhat: ", round(rhat, 4), "\nESS: ", round(ess, 0))) %>% 
  tidyr::extract(parameter, 
                 into = c("parameter", "species"), 
                 regex = "(.*)-(.*)")

stats <- stats %>% 
  mutate(species = paste(substr(stats$species, 1, 1), ". ", 
                         str_extract(stats$species, "(?<=\\s).*"), 
                         sep=""))  %>% 
  mutate(species = factor(samples_df$species, levels = levels.sp))


stats$parameter<-gsub("height.CT", "canopy.h", stats$parameter)
stats$parameter<-gsub("week_s", "week", stats$parameter)
stats$parameter<-factor(stats$parameter, 
                        levels = c("(Intercept)", "slope25m", "canopy.h", "week", "week2", "effort"    
                        ))

traceDet <- trace_plot(samples_df)
# trace
# ggsave(file = "FigureS4.2_TraceRhatESS_Det.png",
#        traceDet,
#        width = 25.2,
#        height = 17.6,
#        units = "cm",
#        dpi = 300)
# -----------------------------