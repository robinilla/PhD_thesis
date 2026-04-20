# -----------------------------------------
# From the work: Illanas, et al (in prep). 
# Large-scale modeling of wild ungulate relative abundance 
# from a national camera trap network
# -----------------------------------------
# script: SuppInfo_S5.R
# This script takes model output of XXXXX 
# and replicate supplementary information S5
# author: Sonia Illanas
# v1.0.0
# date last modification: 2026/03/17
# -----------------------------------------

# 0. Setup. Load required packages and data
# ------------------------------------------
# load packages
library(tidyverse)   # v2.0.0
library(brms)        # v2.23.0
# library(performance) # v0.14.0
library(tidybayes)   # v3.0.7
library(ggh4x)       #v0.3.1


# load model results
load(file = "m_huntWITTHOUT0.Rdata")
load(file = "m_hunt20WITHOUT0.Rdata")
load(file = "m_hunt50WITHOUT0.Rdata")

sp_names <- c( "cc" = "Capreolus capreolus", "ce" = "Cervus elaphus", 
               "cp" = "Capra pyrenaica", "dd" = "Dama dama",   "oa" = "Ovis aries",
               "rp" = "Rupicapra pyrenaica", "wb" = "Sus scrofa")


# make a function to extract and plot all the data from model output
traceS5<-
function(model_output, validation = "External validation - 50 km", source = "HYs"){
  
  # extract fixed effect coefficients 
  samples_df <- model_output %>%
    gather_draws(`b_.*`, regex = TRUE) %>% 
    mutate( chain = factor(.chain), 
            iteration = .iteration,
            parameter = gsub("b_", "", .variable) # Clean predictor names
    )
  
  # extract ranef effect coefficients 
  samples_ranef <- model_output %>%
    # check that name matches with your parameter names r_especie[nivel, parámetro]
    spread_draws(r_especie[especie, parameter]) %>%
    mutate(
      chain = factor(.chain),
      iteration = .iteration,
      value = r_especie,
      #  change parameter names 
      parameter = gsub("Intercept", "(Intercept)", parameter),
      parameter = gsub("obs_z", "slope_obs", parameter)
    )

  # extract Rhat & ESS for fixed effects
  stats <- rhat(model_output) %>%
    enframe(name = "variable", value = "rhat") %>%
    filter(grepl("b_", variable)) %>%
    mutate(
      parameter = gsub("b_", "", variable),
      ess = neff_ratio(model_output)[variable] * ndraws(model_output),
      label = paste0("Rhat: ", round(rhat, 3), "\nESS: ", round(ess, 0))
    )
    
  # extract Rhat & ESS for random effects 
  # species level in our case
  stats_ranef <- model_output %>%
    as_draws_df() %>%
    summarise_draws() %>%
    filter(grepl("^r_especie", variable)) %>%
    mutate(
      # extract data from the variables
      content = str_extract(variable, "(?<=\\[).*(?=\\])"),
      especie = str_remove(content, ",.*"),
      parameter = str_remove(content, ".*,"),
      #  change parameter names 
      parameter = gsub("Intercept", "(Intercept)", parameter),
      parameter = gsub("obs_z", "slope_obs", parameter),
      label = paste0("Rhat: ", round(rhat, 3), "\nESS: ", round(ess_bulk, 0))
    )
  
  # fix effects trace plot
  fix_ef<-
    ggplot(samples_df, aes(x = iteration, y = .value, color = chain)) +
    geom_line(alpha = 0.6, linewidth = 0.3) +
    # add Rhat/ESS
    geom_text(data = stats, 
              aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, 
              hjust = -0.1, vjust = 1.5,
              size = 2.2, fontface = "bold.italic") +
    facet_wrap( ~ parameter, scales = "free_y") +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      axis.title = element_blank()
    ) + 
    labs(title = paste("Traceplots for fixed effects -", validation, source))
  
  ranef_trace<-
    ggplot(samples_ranef, aes(x = iteration, y = value, color = chain)) +
    geom_line(alpha = 0.5, linewidth = 0.2) +
    # Etiquetas de diagnóstico en cada panel
    geom_text(data = stats_ranef, 
              aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, 
              hjust = -0.1, vjust = 1.5,
              size = 2.2, fontface = "bold.italic") +
    facet_grid2(parameter ~ especie, 
                independent = "y", scales = "free_y",
                labeller = labeller(especie = as_labeller(sp_names), 
                                    parameter = as_labeller(c("(Intercept)"="(Intercept)", 
                                                              "slope_obs"="slope_obs")))) +
    theme_minimal() +
    theme(
      legend.position = "none",
      strip.text.y = element_text(face = "bold", size = 7),
      strip.text.x = element_text(face = "bold.italic", size = 7),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 6),
      axis.title = element_blank()
    ) +
    labs(title = paste("Traceplots for species-level random effects -", validation, source))
  
  return(list(
    fixed_trace = fix_ef,
    ranef_trace = ranef_trace
  ))
  
}

# Extract and save traceplots from model output to a global environment object
# plot.trace<-traceS5(m_hunt,   validation = "Within study area", source = "(HYs)")
# plot.trace<-traceS5(m_hunt20, validation = "Nearby study area - 20 km", source = "(HYs)")
# plot.trace<-traceS5(m_hunt50, validation = "Nearby study area - 50 km", source = "(HYs)")

# save plots to a file
ggsave(file = "m_hunt_TRACE_fix.png",
       plot.trace$fixed_trace,
       width = 10,
       height = 6.5,
       units = "cm",
       dpi = 300)

ggsave(file = "m_hunt_TRACE.png",
       plot.trace$ranef_trace,
       width = 32,
       height = 12,
       units = "cm",
       dpi = 300)
