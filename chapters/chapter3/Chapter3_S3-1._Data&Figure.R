library(tidyverse)
library(ggbreak)
library(patchwork)
rm(list=ls())
# search date = 29/01/2026
# https://www.inaturalist.org/observations?captive=false&iconic_taxa=Insecta&view=species


# world data
# ----------------
df<-tibble(grupo         = c("Birds",  "Amphibians", "Reptiles", "Mammals", 
                             "Ray-Finned Fishes", 
                             "Mollusks", "Arachnids", "Insects"), 
           N_especies    = c(11193,   6130, 9407, 4702, 17016, 20724, 17715, 225878), 
           Observaciones = c(40435862,3816938, 6078625, 117908, 3454080, 4768580, 8714362, 78370214),
           Observadores =  c(1161168, 531536, 737041, 653064, 218590, 484992, 888513, 1985595)
       )

df<-df %>% mutate(R_Obs_Obv=Observaciones/Observadores*100) %>% 
           mutate(R_Obs_spp = Observaciones/N_especies*100) %>% 
           mutate (R_spp = N_especies / sum(df$N_especies)*100) %>% 
           mutate (R_Obs = Observaciones / sum(df$Observaciones)*100) %>% 
           mutate (R_Obv = Observadores / sum(df$Observadores)*100)


df2 <- df %>% select(grupo, R_Obs_Obv:R_Obv) %>% 
              pivot_longer(cols = R_Obs_Obv:R_Obv)

df2$name<-factor(df2$name, 
                 levels = c("R_Obs", "R_Obs_Obv", "R_Obs_spp", "R_Obv", "R_spp"), 
                 labels = c("Ratio observations", 
                            "Observation/Observer", "Observation/Spp", "Ratio observers", "Ratio species"))
# "#c799b0", "#b97c9b",

p1<-
ggplot(df2 %>% filter(name!="Observation/Spp") %>% filter(name!="Observation/Observer") , 
       aes(y = grupo, x = value, colour = grupo, shape = name)) +
    geom_point(size = 3.5) +
    scale_colour_manual(values = c("#62B9AD", "#D8D2C4","#FAEE62", "#B05C6E", "#7C4D87", "#728EAF", "#3F4A59", "#AAE280")) +
  theme_classic()+
  scale_x_continuous(limits = c(0, 80),
                     # labels = scales::label_number(scale = 1e-3, suffix = "k"), 
                     breaks = seq(0, 80, 5))+
  scale_x_break(c(20, 23), space = 0.6, scales=0.1, ticklabels = c(30, 40))+
  scale_x_break(c(32, 50), space = 0.6, scales=0.2, ticklabels = c(50, 75))+
  theme(
    axis.title = element_blank(),
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(6, "pt"),
    axis.text.x = element_text(hjust=1, size=10), 
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x.top = element_blank(), #quitar el eje superior
    axis.ticks.x.top = element_blank(),
    axis.line.x.top  = element_blank()
    )+ guides(color="none")
p1

p2<-ggplot(df2 %>% filter(name=="Observation/Spp" | name=="Observation/Observer") , 
       aes(y = grupo, x = value, colour = grupo, shape = name)) +
  geom_point(size = 3.5) +
  scale_colour_manual(values = c("#62B9AD", "#D8D2C4","#FAEE62", "#B05C6E", "#7C4D87", "#728EAF", "#3F4A59", "#AAE280")) +
  # scale_shape_manual(values=c(3, 21))+
  scale_shape_manual(values=c(8, 18))+
  theme_classic()+
  scale_x_continuous(limits = c(0, 375000),
                     labels = scales::label_number(scale = 1e-3, suffix = "k"),
                     breaks = seq(0, 375000, 15000))+
  # scale_x_break(c(73e+04, 3.7e+05)) +
  scale_x_break(c(4000, 20000), space = 0.6, scales=1.1)+
  scale_x_break(c(75000, 360000), space = 0.8, scales=1.2)+
  theme(
    axis.title = element_blank(),
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(6, "pt"),
    axis.text.x = element_text(hjust=1, size=10), 
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x.top = element_blank(), #quitar el eje superior
    axis.ticks.x.top = element_blank(),
    axis.line.x.top  = element_blank(),
    plot.margin = margin(t=20, r=10, b=20, l=40)
  )+ guides(color="none")
p2

leg1 <- cowplot::get_legend(p1)
leg2 <- cowplot::get_legend(p2)

leg_panel1 <- ggplot() + 
  ggimage::geom_subview(x = 0.6, y = 0.5, subview = leg1) +
  #   ggimage::geom_subview(x = 0.77, y = 0.5, subview = leg2) +
  theme_void()

leg_panel2 <- ggplot() + 
  # ggimage::geom_subview(x = 0.33, y = 0.5, subview = leg1) +
  ggimage::geom_subview(x = 0.6, y = 0.5, subview = leg2) +
  theme_void()

## redraw the figure
pf<-(ggplotify::as.ggplot(print(p1 + theme(legend.position = "none"))) /
     ggplotify::as.ggplot(print(p2 + theme(legend.position = "none"))) / 
     leg_panel1 / leg_panel2)+ 
     # patchwork::plot_spacer()) + 
    plot_layout(heights = c(1, 1, 0.1, 0.1)) 
# pf

ggsave(filename = "F:/IREC/PhD/D/Docs/Sonia/0_PhD/7_DocTesis/Chapter_ENETWILD/Figure_iNat.png", # filename of the plot
       plot = pf,           # plot to save
       width = 13, height = 17,   # dimensions
       units = "cm",              # dimension units
       dpi = 300)
# ----------------

# Europe data
# --------------------
df<-tibble(grupo         = c("Birds",  "Amphibians", "Reptiles", "Mammals", 
                             "Ray-Finned Fishes", 
                             "Mollusks", "Arachnids", "Insects"), 
           
           N_especies    = c(1133,   104, 216, 331, 
                             984, 
                             2794, 2866, 34707), 
           Observaciones = c(7751531, 515729, 540757, 793651, 
                             345952, 
                             11418889, 1856052, 18773864),
           Observadores =  c(227605, 95300, 96297, 105838, 
                             34445, 
                             125808, 173987, 449692)
)

df<-df %>% mutate(R_Obs_Obv=Observaciones/Observadores*100) %>% 
  mutate(R_Obs_spp = Observaciones/N_especies*100) %>% 
  mutate (R_spp = N_especies / sum(df$N_especies)*100) %>% 
  mutate (R_Obs = Observaciones / sum(df$Observaciones)*100) %>% 
  mutate (R_Obv = Observadores / sum(df$Observadores)*100)


df2 <- df %>% select(grupo, R_Obs_Obv:R_Obv) %>% 
  pivot_longer(cols = R_Obs_Obv:R_Obv)

df2$name<-factor(df2$name, 
                 levels = c("R_Obs", "R_Obs_Obv", "R_Obs_spp", "R_Obv", "R_spp"), 
                 labels = c("Ratio observations", 
                            "Observation/Observer", "Observation/Spp", "Ratio observers", "Ratio species"))
# "#c799b0", "#b97c9b",
p11<-
  ggplot(df2 %>% filter(name!="Observation/Spp") %>% filter(name!="Observation/Observer") , 
         aes(y = grupo, x = value, colour = grupo, shape = name)) +
  geom_point(size = 3.5) +
  scale_colour_manual(values = c("#62B9AD", "#D8D2C4","#FAEE62", "#B05C6E", "#7C4D87", "#728EAF", "#3F4A59", "#AAE280")) +
  theme_classic()+
  scale_x_continuous(limits = c(0, 81),
                     # labels = scales::label_number(scale = 1e-3, suffix = "k"), 
                     breaks = seq(0, 81, 15))+
  theme(
    axis.title = element_blank(),
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(6, "pt"),
    axis.text.x = element_text(hjust=1, size=10), 
    legend.title = element_blank(),
    legend.position = "bottom"
  )+ guides(color="none")
# p11

p22<-ggplot(df2 %>% filter(name=="Observation/Spp" | name=="Observation/Observer") , 
           aes(y = grupo, x = value, colour = grupo, shape = name)) +
  geom_point(size = 3.5) +
  scale_colour_manual(values = c("#62B9AD", "#D8D2C4","#FAEE62", "#B05C6E", "#7C4D87", "#728EAF", "#3F4A59", "#AAE280")) +
  # scale_shape_manual(values=c(3, 21))+
  scale_shape_manual(values=c(8, 18))+
  theme_classic()+
  scale_x_continuous(limits = c(0, 255000),
                     labels = scales::label_number(scale = 1e-3, suffix = "k"),
                     breaks = seq(0, 255000, 15000))+
  # scale_x_break(c(73e+04, 3.7e+05)) +
  scale_x_break(c(5000, 30000), space = 0.6, scales=1.1)+
  scale_x_break(c(75000, 230000), space = 0.8, scales=0.5)+
  theme(
    axis.title = element_blank(),
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(6, "pt"),
    axis.text.x = element_text(hjust=1, size=10), 
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x.top = element_blank(), #quitar el eje superior
    axis.ticks.x.top = element_blank(),
    axis.line.x.top  = element_blank(),
    plot.margin = margin(t=20, r=10, b=20, l=40)
  )+ guides(color="none")
# p22

leg11 <- cowplot::get_legend(p11)
leg22 <- cowplot::get_legend(p22)

leg_panel11 <- ggplot() + 
  ggimage::geom_subview(x = 0.6, y = 0.5, subview = leg11) +
  #   ggimage::geom_subview(x = 0.77, y = 0.5, subview = leg2) +
  theme_void()

leg_panel22 <- ggplot() + 
  # ggimage::geom_subview(x = 0.33, y = 0.5, subview = leg1) +
  ggimage::geom_subview(x = 0.6, y = 0.5, subview = leg22) +
  theme_void()

## redraw the figure
pf2<-(ggplotify::as.ggplot(print(p11 + theme(legend.position = "none"))) /
       ggplotify::as.ggplot(print(p22 + theme(legend.position = "none"))) / 
       leg_panel11 / leg_panel22)+ 
  # patchwork::plot_spacer()) + 
  plot_layout(heights = c(1, 1, 0.1, 0.1)) 
pf2

pf+pf2 +
  plot_layout(heights = c(1, 1, 0.1, 0.1), widths = c(1, 1)) 
# ----------------

# Asia data
# --------------------
df<-tibble(grupo         = c("Birds",  "Amphibians", "Reptiles", "Mammals", 
                             "Ray-Finned Fishes", 
                             "Mollusks", "Arachnids", "Insects"), 
           
           N_especies    = c(3347,   1632, 2721, 1464, 
                             6652, 
                             7006, 5172, 71651), 
           Observaciones = c(3801579,394819, 609171, 408622, 
                             552687, 
                             503642, 943137, 8631570),
           Observadores =  c(101469, 42034, 70777, 53622, 
                             27708, 
                             49074, 71225, 176410)
)

df<-df %>% mutate(R_Obs_Obv=Observaciones/Observadores*100) %>% 
  mutate(R_Obs_spp = Observaciones/N_especies*100) %>% 
  mutate (R_spp = N_especies / sum(df$N_especies)*100) %>% 
  mutate (R_Obs = Observaciones / sum(df$Observaciones)*100) %>% 
  mutate (R_Obv = Observadores / sum(df$Observadores)*100)


df2 <- df %>% select(grupo, R_Obs_Obv:R_Obv) %>% 
  pivot_longer(cols = R_Obs_Obv:R_Obv)

df2$name<-factor(df2$name, 
                 levels = c("R_Obs", "R_Obs_Obv", "R_Obs_spp", "R_Obv", "R_spp"), 
                 labels = c("Ratio observations", 
                            "Observation/Observer", "Observation/Spp", "Ratio observers", "Ratio species"))

# --------------------

# Africa data
# --------------------
df<-tibble(grupo         = c("Birds",  "Amphibians", "Reptiles", "Mammals", 
                             "Ray-Finned Fishes", 
                             "Mollusks", "Arachnids", "Insects"), 
           
           N_especies    = c(2568,   962, 1808, 1102, 
                             3506, 
                             3734, 2876, 28883), 
           Observaciones = c(1551672,118743, 303281, 473758, 
                             191987, 
                             131392, 234776, 1748254),
           Observadores =  c(33761, 14007, 29600, 24351, 
                             8048, 
                             13842, 22772, 46462)
)

df<-df %>% mutate(R_Obs_Obv=Observaciones/Observadores*100) %>% 
  mutate(R_Obs_spp = Observaciones/N_especies*100) %>% 
  mutate (R_spp = N_especies / sum(df$N_especies)*100) %>% 
  mutate (R_Obs = Observaciones / sum(df$Observaciones)*100) %>% 
  mutate (R_Obv = Observadores / sum(df$Observadores)*100)


df2 <- df %>% select(grupo, R_Obs_Obv:R_Obv) %>% 
  pivot_longer(cols = R_Obs_Obv:R_Obv)

df2$name<-factor(df2$name, 
                 levels = c("R_Obs", "R_Obs_Obv", "R_Obs_spp", "R_Obv", "R_spp"), 
                 labels = c("Ratio observations", 
                            "Observation/Observer", "Observation/Spp", "Ratio observers", "Ratio species"))

# --------------------

# North America data
# --------------------
df<-tibble(grupo         = c("Birds",  "Amphibians", "Reptiles", "Mammals", 
                             "Ray-Finned Fishes", 
                             "Mollusks", "Arachnids", "Insects"), 
           
           N_especies    = c(2417,   1146, 1939, 1007, 
                             3648, 
                             5537, 4959, 67506), 
           Observaciones = c(21720557,2335287, 3797285, 3766943, 
                             1320167, 
                             2169075, 4259561, 39630272),
           Observadores =  c(721650, 340725, 494798, 445438, 
                             134059, 
                             265698, 534294, 1159504)
)

df<-df %>% mutate(R_Obs_Obv=Observaciones/Observadores*100) %>% 
  mutate(R_Obs_spp = Observaciones/N_especies*100) %>% 
  mutate (R_spp = N_especies / sum(df$N_especies)*100) %>% 
  mutate (R_Obs = Observaciones / sum(df$Observaciones)*100) %>% 
  mutate (R_Obv = Observadores / sum(df$Observadores)*100)


df2 <- df %>% select(grupo, R_Obs_Obv:R_Obv) %>% 
  pivot_longer(cols = R_Obs_Obv:R_Obv)

df2$name<-factor(df2$name, 
                 levels = c("R_Obs", "R_Obs_Obv", "R_Obs_spp", "R_Obv", "R_spp"), 
                 labels = c("Ratio observations", 
                            "Observation/Observer", "Observation/Spp", "Ratio observers", "Ratio species"))

# --------------------

# South America data
# --------------------
df<-tibble(grupo         = c("Birds",  "Amphibians", "Reptiles", "Mammals", 
                             "Ray-Finned Fishes", 
                             "Mollusks", "Arachnids", "Insects"), 
           
           N_especies    = c(3595,   2137, 1864, 1175, 
                             3415, 
                             2827, 3449, 44129), 
           Observaciones = c(2858849,265037, 344415, 285115, 
                             126687, 
                             141117, 633553, 4685417),
           Observadores =  c(92716, 38285, 53682, 42165, 
                             12686, 
                             29175, 69312, 143725)
)

df<-df %>% mutate(R_Obs_Obv=Observaciones/Observadores*100) %>% 
  mutate(R_Obs_spp = Observaciones/N_especies*100) %>% 
  mutate (R_spp = N_especies / sum(df$N_especies)*100) %>% 
  mutate (R_Obs = Observaciones / sum(df$Observaciones)*100) %>% 
  mutate (R_Obv = Observadores / sum(df$Observadores)*100)


df2 <- df %>% select(grupo, R_Obs_Obv:R_Obv) %>% 
  pivot_longer(cols = R_Obs_Obv:R_Obv)

df2$name<-factor(df2$name, 
                 levels = c("R_Obs", "R_Obs_Obv", "R_Obs_spp", "R_Obv", "R_spp"), 
                 labels = c("Ratio observations", 
                            "Observation/Observer", "Observation/Spp", "Ratio observers", "Ratio species"))

# --------------------

# Oceania data
# --------------------
df<-tibble(grupo         = c("Birds",  "Amphibians", "Reptiles", "Mammals", 
                             "Ray-Finned Fishes", 
                             "Mollusks", "Arachnids", "Insects"), 
           
           N_especies    = c(2034,   378, 1428, 559, 
                             3926, 
                             6283, 2814, 32019), 
           Observaciones = c(2509263,149296, 407288, 350818, 
                             755779, 
                             641238, 737431, 4608300),
           Observadores =  c(83593, 25411, 51555, 38530, 
                             20749, 
                             31642, 53200, 102895)
)

df<-df %>% mutate(R_Obs_Obv=Observaciones/Observadores*100) %>% 
  mutate(R_Obs_spp = Observaciones/N_especies*100) %>% 
  mutate (R_spp = N_especies / sum(df$N_especies)*100) %>% 
  mutate (R_Obs = Observaciones / sum(df$Observaciones)*100) %>% 
  mutate (R_Obv = Observadores / sum(df$Observadores)*100)


df2 <- df %>% select(grupo, R_Obs_Obv:R_Obv) %>% 
  pivot_longer(cols = R_Obs_Obv:R_Obv)

df2$name<-factor(df2$name, 
                 levels = c("R_Obs", "R_Obs_Obv", "R_Obs_spp", "R_Obv", "R_spp"), 
                 labels = c("Ratio observations", 
                            "Observation/Observer", "Observation/Spp", "Ratio observers", "Ratio species"))

# --------------------

