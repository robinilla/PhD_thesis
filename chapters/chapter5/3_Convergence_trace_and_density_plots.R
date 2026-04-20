library(spAbundance)
library(tidyverse)

load(file="J:/IREC_Sonia/2_WIP/Report202512/2_modelResults/out_20251125_MCMCparamTest12_out4b_EU_largeSAMPLES.Rdata")

mcmc.model.params4<-function(modelo){# Supón que tienes un objeto tipo mcmc colapsado
  collapsed <- mcmc(cbind(modelo$beta.samples, 
                          modelo$theta.samples, 
                          modelo$kappa.samples))# y es class == "mcmc"
  
  # Ver el número de iteraciones y parámetros
  n_total <- niter(collapsed)  # Total de filas
  n_chains <- 3
  n_iter_per_chain <- n_total / n_chains
  
  # Verifica que la división es exacta
  if (n_total %% n_chains != 0) stop("No se puede dividir equitativamente en 3 cadenas.")
  
  # Divide en 3 cadenas
  chain_list <- lapply(1:n_chains, function(i) {
    start <- ((i - 1) * n_iter_per_chain + 1)
    end <- i * n_iter_per_chain
    mcmc(as.matrix(collapsed)[start:end, ])
  })
  
  # Reconstruye el mcmc.list
  samples_fixed <- mcmc.list(chain_list)
  
  # Convertir a array compatible con bayesplot
  mcmc_array <- as.array(samples_fixed)
  mcmc_array <- aperm(mcmc_array, perm = c(1, 3, 2))
  
  # Verifica:
  # dim(mcmc_array)
  
  # Ver los nombres de los parámetros disponibles
  dimnames(mcmc_array)[[3]]<-gsub("\\(Intercept)", "Intercept", 
                                  gsub("cforest", "forest", 
                                       gsub("bio_12", "prec", dimnames(mcmc_array)[[3]])))
  params <- dimnames(mcmc_array)[[3]]
  
  # O selecciona solo algunos:
  # params <- c("beta[1]", "beta[2]", ...)
  
  summary_stats<-tibble(params=params, 
                        Rhat=c(round(modelo$rhat$beta, 2), 
                               round(modelo$rhat$theta, 2), 
                               round(modelo$rhat$kappa, 2)), 
                        ESS=c(round(modelo$ESS$beta,0), 
                              round(modelo$ESS$theta, 0), 
                              round(modelo$ESS$kappa, 0)))
  
  # Generar traceplots por parámetro
  traceplots <- lapply(params, function(p) {
    bayesplot::mcmc_trace(mcmc_array, pars = p) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(breaks = function(x) range(x), 
                         labels = scales::number_format(accuracy = 0.01))+
      theme_classic()+
      theme(legend.position = "none", 
            axis.text = element_text(family = "sans", size=7), 
            axis.title.y = element_text(family = "sans", size=9)) #+ ggtitle(p)
  })
  
  # Generar density plots por parámetro
  densityplots <- lapply(params, function(p) {
    bayesplot::mcmc_dens(mcmc_array, pars = p) + xlab(NULL) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 2))+
      # scale_y_continuous(breaks = function(x) range(x), 
      #                    labels = scales::number_format(accuracy = 1))+
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
  
  # Apilar verticalmente
  trace_col <- patchwork::wrap_plots(traceplots, ncol = 1)
  dens_col  <- patchwork::wrap_plots(densityplots, ncol = 1)
  
  # Combinar columnas: traceplots | density plots
  final_plot <- trace_col | dens_col
  final_plot
}

b<-mcmc.model.params4(out4b_EU)

ggsave(filename = "J:/IREC_Sonia/2_WIP/Report202512/4_plots/Figure_Traceplots1.2.png",
       plot = b,
       device = NULL,
       # path = "",
       scale = 1,
       width=18,
       height=26,
       units="cm",
       dpi = 300)
