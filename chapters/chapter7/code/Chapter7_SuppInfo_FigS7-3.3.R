# -----------------------------------------
# From the work: Illanas, et al (in prep). 
# Large-scale modeling of wild ungulate relative abundance 
# from a national camera trap network
# -----------------------------------------
# script: SuppInfo_FigS2.3.R
# This script was used for making Figure S2.3 
# author: Sonia Illanas
# v1.0.0
# date last modification: 2026/03/17
# -----------------------------------------

# 0. Setup. Load required packages and data
# ------------------------------------------

# load packages
library(tidyverse)   # v2.0.0
library(ggpubr)      # v0.6.0

# load data covariates information
load(file = "data.model.Rdata")
vars<-tibble(deployment = dimnames(data.model$y)[[2]], 
             data.model$abund.covs) %>% 
             separate(deployment, into = c("loc", "id"), sep = "\\.") %>% 
            mutate(location = sprintf("L%02d", as.numeric(location)))


plot.vars<-function(data, covariate, label.y, color){
ggplot(data, aes(x = location, y = {{covariate}})) +
  geom_boxplot(alpha = 0.7, 
               fill = color, 
               outlier.color = "darkred", 
               outlier.shape = 16) +
  scale_fill_brewer(palette = "Set3") + 
  labs(y = label.y) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size = 8), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90,
                                   vjust = 1,
                                   hjust = 1,
                                   margin = margin(t = 10))
        )
  }

var.var<-
ggarrange(
  plot.vars(data = vars, covariate = precip,label.y = "Precipitation (mm)", color = "#728EAF"),
  plot.vars(data = vars, covariate = elev,  label.y = "Elevation (m)",  color = "#62B9AD"),
  plot.vars(data = vars, covariate = slope,  label.y = "Slope (%)", color = "#fec46c"),
  plot.vars(data = vars, covariate = forest,  label.y = "Forest cover (%)", color = "#3da75a"),
  plot.vars(data = vars, covariate = shrb,  label.y = "Shrub cover (%)", color = "#8cab52"),
nrow = 5, ncol= 1)
ggsave(file = "~/FigureS7-3.3.png",
       plot = var.var,
       width = 8, 
       height = 24,
       units = "cm",
       dpi = 300)
