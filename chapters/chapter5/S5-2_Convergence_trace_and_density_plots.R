# ---------------------------------------------
#  PhD Thesis: Large-scale monitoring of wild mammal abundance: 
#  modeling frameworks for structured and unstructured data
#  Chapter 5. Accounting for spatial dependence to smooth regional patterns of wild boar abundance across Europe
#  Author: Sonia Illanas
#  Institution: Institute of Game and Wildlife Research
#  Date of last modification: 12/12/2026
# ---------------------------------------------
## R version 4.5.2
## tidyverse version: 2.0.0
## spAbundance version: 0.2.2

library(spAbundance)
library(tidyverse)

load(file="out_20251125_MCMCparamTest12_out4b_EU_largeSAMPLES.Rdata")

mcmc.model.params4<-function(model){
  collapsed <- mcmc(cbind(model$beta.samples, 
                          model$theta.samples, 
                          model$kappa.samples))# y es class == "mcmc"
  
  # number of iterations and parameters
  n_total <- niter(collapsed)  # total nrow
  n_chains <- 3
  n_iter_per_chain <- n_total / n_chains
  
  # 
  if (n_total %% n_chains != 0) stop("No se puede dividir equitativamente en 3 cadenas.")
  
  # Number of chains
  chain_list <- lapply(1:n_chains, function(i) {
    start <- ((i - 1) * n_iter_per_chain + 1)
    end <- i * n_iter_per_chain
    mcmc(as.matrix(collapsed)[start:end, ])
  })
  
  # Rebuilt mcmc.list
  samples_fixed <- mcmc.list(chain_list)
  
  # Convert to an array compatible with bayesplot
  mcmc_array <- as.array(samples_fixed)
  mcmc_array <- aperm(mcmc_array, perm = c(1, 3, 2))
  
  # Verify:
  # dim(mcmc_array)
  
  # parameters name available
  dimnames(mcmc_array)[[3]]<-gsub("\\(Intercept)", "Intercept", 
                                  gsub("cforest", "forest", 
                                       gsub("bio_12", "prec", dimnames(mcmc_array)[[3]])))
  params <- dimnames(mcmc_array)[[3]]
  
  summary_stats<-tibble(params=params, 
                        Rhat=c(round(modelo$rhat$beta, 2), 
                               round(modelo$rhat$theta, 2), 
                               round(modelo$rhat$kappa, 2)), 
                        ESS=c(round(modelo$ESS$beta,0), 
                              round(modelo$ESS$theta, 0), 
                              round(modelo$ESS$kappa, 0)))
  
  # make traceplots for each parameter
  traceplots <- lapply(params, function(p) {
    bayesplot::mcmc_trace(mcmc_array, pars = p) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(breaks = function(x) range(x), 
                         labels = scales::number_format(accuracy = 0.01))+
      theme_classic()+
      theme(legend.position = "none", 
            axis.text = element_text(family = "sans", size=7), 
            axis.title.y = element_text(family = "sans", size=9))
  })
  
  # make density plots  for each parameter
  densityplots <- lapply(params, function(p) {
    bayesplot::mcmc_dens(mcmc_array, pars = p) + xlab(NULL) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 2))+
      theme_classic()+
      theme(legend.position = "none")
  })
  
  for (i in 1:length(densityplots)){
    gb <- ggplot_build(densityplots[[i]])
    x_range <- gb$layout$panel_params[[1]]$x.range
    y_range <- gb$layout$panel_params[[1]]$y.range
    
    x_pos <- x_range[1] + 0.88 * (x_range[2] - x_range[1])
    y_pos <- y_range[1] + 0.54 * (y_range[2] - y_range[1])
    
    densityplots[[i]]<-densityplots[[i]]+
      annotate("text", 
               label = paste("Rhat : ",  summary_stats[i,2:3]$Rhat, "\n", "ESS : ", summary_stats[i,2:3]$ESS),
               color="#005b96", 
               size=2.5,
               x=x_pos, 
               y=y_pos)
  }
  
  # vertical visual placement
  trace_col <- patchwork::wrap_plots(traceplots, ncol = 1)
  dens_col  <- patchwork::wrap_plots(densityplots, ncol = 1)
  
  # Combine columns: traceplots | density plots
  final_plot <- trace_col | dens_col
  final_plot
}

b<-mcmc.model.params4(out4b_EU)

ggsave(filename = "FigureS5-2.1_Traceplots.png",
       plot = b,
       device = NULL,
       # path = "",
       scale = 1,
       width=18,
       height=26,
       units="cm",
       dpi = 300)
