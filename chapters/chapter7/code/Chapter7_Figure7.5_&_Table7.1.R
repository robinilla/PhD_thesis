# -----------------------------------------
# From the work: Illanas, et al (in prep). 
# Large-scale modeling of wild ungulate relative abundance 
# from a national camera trap network
# -----------------------------------------
# script: Figure5_&_Table2.R
# This script takes model output of XXXXX 
# and replicate figure 5 and Table 2
# author: Sonia Illanas
# v1.0.0
# date last modification: 2026/03/17
# -----------------------------------------

# 0. Setup. Load required packages and data
# ------------------------------------------
# load packages
library(tidyverse)   # v2.0.0
library(brms)        # v2.23.0
library(performance) # v0.14.0
library(flextable)   # v0.9.11

# load model results
setwd("~/Documentos/PhD/5_Estancia/2_WIP/2_Results/MSNM_20260312_spAbundance/m_HY_val/")

load(file = "m_huntWITTHOUT0.Rdata")
load(file = "m_hunt20WITHOUT0.Rdata")
load(file = "m_hunt50WITHOUT0.Rdata")

load(file = "m_gbif.Rdata")
load(file = "m_gbif20.Rdata")
load(file = "m_gbif50.Rdata")


# Prepare model coefficents to export Table 2
# -------------------------------------------
extract_model_stats<-function(data){
  draws <- as_draws_df(data)
  r2 <- r2_bayes(data)
  tibble( beta = round(mean(draws[[2]]), 3),
          beta_95CI = paste(round(quantile(draws[[2]], 0.025), 3), 
                            "-",
                            round(quantile(draws[[2]], 0.975), 3)), 
          SD_intercept_sp = round(mean(draws$sd_especie__Intercept), 5),
          SD_slope_sp     = round(mean(draws$sd_especie__obs_z), 3),
          sigma           = round(mean(draws$sigma), 3),
          R2_conditional  = round(r2$R2_Bayes, 3),
          R2_marginal     = round(r2$R2_Bayes_marginal, 3)
          )
  }


table2<-rbind(extract_model_stats(m_hunt), 
              extract_model_stats(m_hunt20), 
              extract_model_stats(m_hunt50), 
              extract_model_stats(m_gbif), 
              extract_model_stats(m_gbif20), 
              extract_model_stats(m_gbif50)
              ) %>% 
            mutate(model = c("H1", "H2", "H3", "G1", "G2", "G3"),
                   spatial_extent = rep(c("internal", "external (20 km)" , "external (50 km)"), 2)
                   ) 
table2 <- table2 %>% select(model, spatial_extent, colnames(table2))

ft <- flextable(table2) %>% autofit() %>% theme_booktabs()

# Save Table 2 to word 
# save_as_docx(ft, path = "Table2.docx")


# Make Figure 5
# --------------------------------------
sp_names <- c( "rp" = "Rupicapra pyrenaica", "wb" = "Sus scrofa", "ce" = "Cervus elaphus", "cp" = "Capra pyrenaica", "cc" = "Capreolus capreolus", "oa" = "Ovis aries", "dd" = "Dama dama")
sort_pred<-c("Internal validation - site", "External validation - 20km", "External validation - 50km")

df.post <- 
  rbind(
  rbind(as_tibble(coef(m_hunt, probs = c(0.025, 0.125, 0.875, 0.975))$especie[,,2],   rownames = "species") %>% mutate(pred = "Internal validation - site"),
        as_tibble(coef(m_hunt20, probs = c(0.025, 0.125, 0.875, 0.975))$especie[,,2], rownames = "species") %>% mutate(pred = "External validation - 20km"),
        as_tibble(coef(m_hunt50, probs = c(0.025, 0.125, 0.875, 0.975))$especie[,,2], rownames = "species") %>% mutate(pred = "External validation - 50km")
        ) %>% 
        mutate(source = "Hunting yield"),
  rbind(as_tibble(coef(m_gbif, probs = c(0.025, 0.125, 0.875, 0.975))$especie[,,2],   rownames = "species") %>% mutate(pred = "Internal validation - site"),
        as_tibble(coef(m_gbif20, probs = c(0.025, 0.125, 0.875, 0.975))$especie[,,2], rownames = "species") %>% mutate(pred = "External validation - 20km"),
        as_tibble(coef(m_gbif50, probs = c(0.025, 0.125, 0.875, 0.975))$especie[,,2], rownames = "species") %>% mutate(pred = "External validation - 50km")) %>% 
        mutate(source = "GBIF occurrence frequency")) %>%  
  mutate(species = sp_names[species]) %>% 
  mutate( pred = factor(pred, levels = sort_pred)) %>%
  mutate(species = fct_reorder(species, Estimate, .fun = mean))

post<-
  ggplot(data = df.post, aes(y = species, x = Estimate, shape = pred))+
  geom_vline(xintercept = 0, col="grey50", lty = 3, linewidth= 0.5)+
  
  geom_errorbar(aes(xmin = Q2.5, xmax = Q97.5),
                width = 0,
                size = 0.3,
                col = "grey70", 
                position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(xmin = Q12.5,xmax = Q87.5),
                width = 0,
                size = 0.5, 
                col = "grey40", 
                position = position_dodge(width = 0.6)) +  
  
  geom_point(aes(x = Estimate), size = 1.2, position = position_dodge(width = 0.6)) +
  xlab("Coefficient slope")+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(face = "italic", size = 7),
        axis.text.x = element_text(size = 8), 
        axis.title.x = element_text(size = 8, 
                                    margin = margin(t = 10)),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7), 
        legend.position = "bottom",
        plot.title = element_blank(), 
        legend.title = element_blank()
  )+
  facet_grid( ~ source, scale = "free_x")

# ggsave(file = "Figure5.png",
#        post,
#        width = 14,
#        height = 10,
#        units = "cm",
#        dpi = 300)