# -----------------------------------------
# From the work: Illanas, et al (in prep). 
# Large-scale modeling of wild ungulate relative abundance 
# from a national camera trap network
# -----------------------------------------
# script: Figure3_&_4.R
# This script takes model output coefficient of the msNmix model with better WAIC 
# and replicate figure 3 and 4
# author: Sonia Illanas
# v1.0.0
# date last modification: 2026/03/17
# -----------------------------------------

# 0. Setup. Load required packages and data
# ------------------------------------------
# load packages
library(tidyverse)   # v2.0.0
library(spAbundance) # v0.2.2
library(scales)      # v1.4.0 - needed for label_number
library(patchwork)   # v.1.3.0

# load model results
load(file = "msNmix_modelResults.Rdata")
out<-out.comb; rm(out.comb)

# 1. Prepare model coefficients values
# --------------------------------------------
a <- rbind(
      tibble(
        median = apply(out$beta.samples, 2, median) ,
        low = apply(out$beta.samples, 2, quantile, 0.025),
        high = apply(out$beta.samples, 2, quantile, 0.975),
        low2 = apply(out$beta.samples, 2, quantile, 0.125),
        high2 = apply(out$beta.samples, 2, quantile, 0.875),
        vars = out$beta.samples %>% as_tibble( ) %>% names()
      ) %>% separate_wider_delim(vars, delim = "-", names = c("vars", "species")) %>% 
        mutate( level = "species", class = "Abundance"),
      
      tibble(
        median = apply(out$alpha.samples, 2, median) ,
        low = apply(out$alpha.samples, 2, quantile, 0.025),
        high = apply(out$alpha.samples, 2, quantile, 0.975),
        low2 = apply(out$alpha.samples, 2, quantile, 0.125),
        high2 = apply(out$alpha.samples, 2, quantile, 0.875),
        vars = out$alpha.samples %>% as_tibble( ) %>% names()
      ) %>% separate_wider_delim(vars, delim = "-", names = c("vars", "species")) %>% 
        mutate( level = "species", class = "Detectability"),
      
      tibble(
        median = apply(out$beta.comm.samples, 2, median) ,
        low = apply(out$beta.comm.samples, 2, quantile, 0.025),
        high = apply(out$beta.comm.samples, 2, quantile, 0.975),
        low2 = apply(out$beta.comm.samples, 2, quantile, 0.125),
        high2 = apply(out$beta.comm.samples, 2, quantile, 0.875),
        vars = out$beta.comm.samples %>% as_tibble( ) %>% names()
      ) %>% mutate(species = "Community", level = "community", class = "Abundance"), 
      
      tibble(
        median = apply(out$alpha.comm.samples, 2, median) ,
        low = apply(out$alpha.comm.samples, 2, quantile, 0.025),
        high = apply(out$alpha.comm.samples, 2, quantile, 0.975),
        low2 = apply(out$alpha.comm.samples, 2, quantile, 0.125),
        high2 = apply(out$alpha.comm.samples, 2, quantile, 0.875),
        vars = out$alpha.comm.samples %>% as_tibble( ) %>% names()
      ) %>% mutate(species = "Community", level = "community", class = "Detectability")
    )

# rename variables and factors for plotting
a <- a %>% mutate(species = factor(species, levels = c("Capra pyrenaica", "Capreolus capreolus", "Cervus elaphus", "Dama dama", "Ovis aries", "Rupicapra pyrenaica", "Sus scrofa", "Community")), 
                   vars = case_match(vars, "shrb"      ~ "shrub",
                                           "precip"    ~ "prec",
                                           "height.CT" ~ "canopy.h",
                                           "week_s"    ~ "week",
                                           .default    = vars),
                   vars = factor(vars, levels = c("(Intercept)", "elev", "elev2", "prec", "prec2", 
                                                  "slope", "forest", "shrub", "slope25m", 
                                                  "canopy.h", "week", "week2", "effort")),
                  overlap = case_when(low2 < 0 & high2 > 0 ~ "CI 75",
                                      low < 0 & high > 0   ~ "CI 95",
                                      TRUE ~ "No")
                  )


# Figure 3
# -------------------------
DetAb<-
  ggplot(a %>% filter(vars=="(Intercept)") %>% filter(level=="species") %>% 
           mutate(median = ifelse(class=="Abundance", exp(median), plogis(median)),
                  low = ifelse(class=="Abundance", exp(low), plogis(low)),
                  high = ifelse(class=="Abundance", exp(high), plogis(high))) %>% 
           mutate(species = reorder(species, median)),
         aes(y = species,
             color = class))+
  geom_errorbar(aes(xmin = low, xmax = high), width = 0.0, size = 0.5) +
  geom_point(aes(x = median), size= 1.2) +
  theme_classic() +
  theme(axis.title = element_blank(), 
        axis.text.y = element_text(face = "italic", hjust = 1, vjust = 0.5), 
        legend.position = "none", 
        strip.background = element_blank(), 
        strip.text = element_text(face = "bold", size = 10))+
  facet_wrap(~class, scale = "free_x", 
             labeller = as_labeller(c("Abundance" = "Relative abundance", 
                                      "Detectability" = "Detection probability")))+
  scale_color_manual(values = c("Abundance" = "#94CA94",
                                "Detectability" = "#7C4D87"))

# see the plot in the plot pane
# DetAb

# save the plot in a file
ggsave(file = "Figure3.png",
       DetAb,
       width = 14,
       height = 8,
       units = "cm",
       dpi = 300)
# -------------------------


# Figure 4 
# -------------------------
abCoef<-ggplot(a %>% filter(class =="Abundance") %>% 
               filter(vars != "(Intercept)"), aes(y = species, col = overlap, alpha = overlap))+
  geom_vline(xintercept = 0, col="grey50", lty = 3, linewidth= 0.8) +
  geom_errorbar(aes(xmin = low,
                    xmax = high),
                width = 0.0,
                size = 0.8) +
  geom_point(aes(x = median), size= 3) +
  theme_bw() +
  theme(axis.title = element_blank(), 
        axis.text.y = element_text(h = 1, v = 0.5, face = "italic"), 
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"), 
        strip.text.y = element_blank(),
        strip.placement = "outside") +
  facet_grid ( level ~ vars, scales = "free", space = "free_y", switch = "x") +
  scale_color_manual(values = c("#7C4D87", "#bc80bd", "#c4cf7c")) +
  scale_alpha_manual(values = c(1, 1, 1, 1, 1, 1, 1, 1)) +
  scale_x_continuous(labels = label_number(accuracy = 0.1) # Fuerza siempre 1 decimal
  )

detCoef<-ggplot(a %>% filter(class =="Detectability") %>% 
                filter(vars != "(Intercept)"), aes(y = species, col = overlap, alpha = overlap))+
  geom_vline(xintercept = 0, col="grey50", lty = 3, linewidth= 0.8) +
  geom_errorbar(aes(xmin = low,
                    xmax = high),
                position = position_dodge(width = 0.4),
                width = 0.0,
                size = 0.8) +
  geom_point(aes(x = median), size= 3) +
  theme_bw() +
  theme(axis.title = element_blank(), 
        axis.text.y = element_text(h = 1, v = 0.5, face = "italic"), 
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"), 
        strip.text.y = element_blank(),
        strip.placement = "outside") +
  facet_grid ( level ~ vars, scales = "free", space = "free_y", switch = "x")+
  scale_color_manual(values = c("#7C4D87", "#bc80bd", "#c4cf7c")) +
  scale_alpha_manual(values = c(1, 1, 1, 1, 1, 1, 1, 1)) +
  scale_x_continuous(labels = label_number(accuracy = 0.1) # Fuerza siempre 1 decimal
  )

# Plot together abundance and detectability coefficients plots
bothCoef <- wrap_plots(
  abCoef + theme(plot.tag = element_text(face = "bold", size = 15), 
               strip.text.x = element_text(size = 14), 
               axis.text.y =  element_text(size = 14)),
  detCoef + theme(plot.tag = element_text(face = "bold", size = 15), 
                strip.text.x = element_text(size = 14), 
                axis.text.y = element_text(size = 14)),
  ncol = 1) +
  plot_annotation( tag_levels = "a", tag_suffix = ")" )

# save the plot in a file
ggsave(file = "Figure4.png",
       bothCoef,
       width = 30,
       height = 20,
       units = "cm",
       dpi = 300)
