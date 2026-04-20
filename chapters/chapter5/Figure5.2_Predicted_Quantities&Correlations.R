# ---------------------------------------------
#  PhD Thesis: Large-scale monitoring of wild mammal abundance: 
#  modeling frameworks for structured and unstructured data
#  Chapter 5. Accounting for Spatial Dependence to Smooth
#  Regional Patterns of Wild Boar Abundance across Europe
#  Author: Sonia Illanas
#  Institution: Institute of Game and Wildlife Research
#  Date of last modification: 12/12/2025
# ---------------------------------------------
## R version 4.5.2
## tidyverse version: 2.0.0
## sf version: 1.0-21
## ggbreak version: 0.1.6
## tmap version: 4.2
## ggh4x version: 0.3.1
## corrplot version: 0.95
## grDevices version: 4.5.2
## Hmisc version: 5.2-4

library(sf)
library(tidyverse)
library(ggbreak) 
library(tmap)
tmap::tmap_mode("view")

# load 10km x 10km spatial predictions
grid10km_env<-
  st_read(dsn = "outSpeciesPredictions_2025December.gpkg", 
          layer="out_20251202_grid10km_environmental2")
grid10km_spt<-
  st_read(dsn = "outSpeciesPredictions_2025December.gpkg", 
          layer="out_20251202_grid10km")

wb_model<-
  st_read(dsn="predicciones_20220811.gpkg", 
          layer="DensityModel_grid10km_WildBoar_BR_final4_MESS_IC_v1.1")

wb_model<-st_transform(wb_model, 
                       st_crs("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=km +no_defs"))
wb_model <- wb_model %>% mutate(Var1=round(x_cen/1000, digits=3),
                                Var2=round(y_cen/1000, digits=3))

crop<-st_read("wb_outside_range_clean3.shp") %>% select(geometry) %>% mutate(ID=1)
tm_shape(crop) + tm_fill()

countries<-st_read("NUTS_RG_01M_2021_3035_LEVL_0.shp")
countries<-countries %>% select(NUTS_ID:NUTS_NAME)
grid10km_env_pt<-st_centroid(st_transform(grid10km_env, 3035))
grid10km_spt_pt<-st_centroid(st_transform(grid10km_spt, 3035))
wb_model_pt<-st_centroid(st_transform(wb_model, 3035))


join_env<-st_join(grid10km_env_pt, countries, join = st_nearest_feature)
join_spt<-st_join(grid10km_spt_pt, countries, join = st_nearest_feature)
join_wb<-st_join(wb_model_pt, countries, join = st_nearest_feature)

idx_i<-lengths(st_intersects(join_env, crop))>0
join_env2<-join_env[!idx_i,]
idx_i<-lengths(st_intersects(join_spt, crop))>0
join_spt2<-join_spt[!idx_i,]
idx_i<-lengths(st_intersects(join_wb, crop))>0
join_wb2<-join_wb[!idx_i,]

a <- tibble(
  median  = c(sum(join_env2$spA.m4b),
              sum(join_spt2$spA.m4b), 
              sum(join_wb2$Pred_Red*100)),
  lowCI   = c(sum(join_env2$spA.m4b.lowCI),
              sum(join_spt2$spA.m4b.lowCI), 
              sum(join_wb2$pred_lower*100)),
  highCI  = c(sum(join_env2$spA.m4b.highCI),
              sum(join_spt2$spA.m4b.highCI), 
              sum(join_wb2$pred_upper*100)),
  model   = c("environmental", "spatial", "lastReport")
)

a$model<-factor(a$model, levels=c("environmental", "spatial", "lastReport"), labels = c("ENVIRONMENTAL", "SPATIAL", "ENETWILD'22"))

# Define break limits
low_end   <- 60000000
high_start <- 150000000

# maximum value 
y_max <- max(a$highCI, a$median, na.rm = TRUE)

plot.a<-
ggplot(a, aes(x = model, y = median, colour = model)) +
  geom_errorbar(aes(ymin = lowCI, ymax = highCI),
                width = 0.2, alpha = 0.8, linewidth = 0.6) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = c("#b3cde0", "#005b96", "#011f4b")) +
  # Define manually breaks
  scale_y_continuous(
    breaks = c(
      seq(0, low_end, by = 20000000),               # 0, 20M, 40M, 60M
      seq(high_start, y_max, by = 10000000)         # 154M, 164M, 174M, ...
    ),
    labels = scales::label_number(scale = 1e-6, suffix = "M"),
    minor_breaks = NULL                              
  ) +
  scale_y_break(c(low_end, high_start), scales = 0.1) +
  ylab("Total predicted \n harvested wild boar") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.ticks = element_line(color = "#011f4b", linewidth = 1),
    axis.ticks.length = unit(6, "pt"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text.x = element_text(angle = 90, hjust = 1, size=8), 
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right  = element_blank()
    
  )
ggsave(filename = "Figure_SumQuantityPredicted2.1.png",
       plot = plot.a,
       device = NULL,
       scale = 1,
       width=8,
       height=6,
       units="cm",
       dpi = 300)

b<-rbind(
rbind(
  join_env2 %>% as_tibble() %>% filter(CNTR_CODE=="ES") %>% select(spA.m4b:NUTS_NAME) %>% rename(median=spA.m4b ) %>%  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% mutate(country="Spain", model="environmental"),
  join_spt2 %>% as_tibble() %>% filter(CNTR_CODE=="ES") %>% select(spA.m4b:NUTS_NAME) %>% rename(median=spA.m4b ) %>%  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% mutate(country="Spain", model="spatial"),
  join_wb2  %>% as_tibble() %>% filter(CNTR_CODE=="ES") %>% select(Pred_Red:pred_upper, NUTS_ID:NUTS_NAME) %>% rename(median=Pred_Red, spA.m4b.lowCI=pred_lower, spA.m4b.highCI=pred_upper) %>%  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% mutate(median=median*100, spA.m4b.lowCI=spA.m4b.lowCI*100, spA.m4b.highCI=spA.m4b.highCI*100, country="Spain", model="lastReport")),

rbind(
  join_env2 %>% as_tibble() %>% filter(CNTR_CODE=="IT") %>% select(spA.m4b:NUTS_NAME) %>% rename(median=spA.m4b ) %>%  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% mutate(country="Italy", model="environmental"),
  join_spt2 %>% as_tibble() %>% filter(CNTR_CODE=="IT") %>% select(spA.m4b:NUTS_NAME) %>% rename(median=spA.m4b ) %>%  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% mutate(country="Italy", model="spatial"),
  join_wb2  %>% as_tibble() %>% filter(CNTR_CODE=="IT") %>% select(Pred_Red:pred_upper, NUTS_ID:NUTS_NAME) %>% rename(median=Pred_Red, spA.m4b.lowCI=pred_lower, spA.m4b.highCI=pred_upper) %>%  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% mutate(median=median*100, spA.m4b.lowCI=spA.m4b.lowCI*100, spA.m4b.highCI=spA.m4b.highCI*100, country="Italy", model="lastReport")),

rbind(
  join_env2 %>% as_tibble() %>% filter(CNTR_CODE=="FR") %>% select(spA.m4b:NUTS_NAME) %>% rename(median=spA.m4b ) %>%  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% mutate(country="France", model="environmental"),
  join_spt2 %>% as_tibble() %>% filter(CNTR_CODE=="FR") %>% select(spA.m4b:NUTS_NAME) %>% rename(median=spA.m4b ) %>%  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% mutate(country="France", model="spatial"),
  join_wb2  %>% as_tibble() %>% filter(CNTR_CODE=="FR") %>% select(Pred_Red:pred_upper, NUTS_ID:NUTS_NAME) %>% rename(median=Pred_Red, spA.m4b.lowCI=pred_lower, spA.m4b.highCI=pred_upper) %>%  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% mutate(median=median*100, spA.m4b.lowCI=spA.m4b.lowCI*100, spA.m4b.highCI=spA.m4b.highCI*100, country="France", model="lastReport")),

rbind(
  join_env2 %>% as_tibble() %>% filter(CNTR_CODE=="PL") %>% select(spA.m4b:NUTS_NAME) %>% rename(median=spA.m4b ) %>%  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% mutate(country="Poland", model="environmental"),
  join_spt2 %>% as_tibble() %>% filter(CNTR_CODE=="PL") %>% select(spA.m4b:NUTS_NAME) %>% rename(median=spA.m4b ) %>%  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% mutate(country="Poland", model="spatial"),
  join_wb2  %>% as_tibble() %>% filter(CNTR_CODE=="PL") %>% select(Pred_Red:pred_upper, NUTS_ID:NUTS_NAME) %>% rename(median=Pred_Red, spA.m4b.lowCI=pred_lower, spA.m4b.highCI=pred_upper) %>%  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% mutate(median=median*100, spA.m4b.lowCI=spA.m4b.lowCI*100, spA.m4b.highCI=spA.m4b.highCI*100, country="Poland", model="lastReport")
  )
) %>% rename(lowCI=spA.m4b.lowCI, highCI=spA.m4b.highCI)

b$model<-factor(b$model, levels=c("environmental", "spatial", "lastReport"), labels = c("ENVIRONMENTAL", "SPATIAL", "ENETWILD'22"))

library(ggh4x)
plot.b<-
  ggplot(b, aes(x = model, y = median, colour = model)) +
  geom_errorbar(aes(ymin = lowCI, ymax = highCI),
                width = 0.2, 
                alpha = 0.8,
                linewidth=0.6) +
  geom_point(size = 1.5) +
  scale_colour_manual(values=c("#b3cde0", "#005b96", "#011f4b")) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M"), 
                     breaks = seq(0, 6500000, 1000000)
  )+
  ylab("Total predicted \n harvested wild boar")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.ticks = element_line(color = "#011f4b",
                                  linewidth = 1), 
        axis.title.y = element_text(margin = margin(r = 10)), 
        axis.text.x = element_text(angle = 90, hjust = 1, size=8), 
        axis.text.y.right  = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right  = element_blank())+
  facet_wrap(~country, scales="free_y")

plot.b
ggsave(filename = "Figure5.2A_SumQuantityPredictedb_1.1.png",
       plot = plot.b,
       device = NULL,
       scale = 1,
       width=8,
       height=8,
       units="cm",
       dpi = 300)

library(patchwork)
plots<-plot.a | plot.b
ggsave(filename = "Figure5.2A_SumQuantityPredicted_3.1.png",
       plot = plots,
       device = NULL,
       scale = 1,
       width=20,
       height=10,
       units="cm",
       dpi = 300)


env2<-join_env2 %>% as_tibble() %>% select(spA.m4b:spA.m4b.highCI, Var1, Var2, CNTR_CODE)
colnames(env2)<-gsub("spA.m4b", "env", colnames(env2))
spt2<-join_spt2 %>% as_tibble() %>% select(spA.m4b:spA.m4b.highCI, Var1, Var2, CNTR_CODE)
colnames(spt2)<-gsub("spA.m4b", "spt", colnames(spt2))
last<-join_wb2 %>%  as_tibble() %>% select(Pred_model, Pred_Red, Var1, Var2, CNTR_CODE)

df<-env2 %>% mutate(Var1=round(Var1, digits=3), Var2=round(Var2, digits=3)) %>% 
      left_join(spt2  %>% mutate(Var1=round(Var1, digits=3), Var2=round(Var2, digits=3)), 
                by=c("Var1", "Var2", "CNTR_CODE")) %>% 
      left_join(last, by=c("Var1", "Var2", "CNTR_CODE")) %>% rename(wb_last=Pred_Red)

library(corrplot)
library(grDevices)
library(Hmisc)
df.g<-df %>% select(env, spt, wb_last) %>% mutate(wb_last=wb_last*100) %>% rename(ENVIRONMENTAL=env, 
                                                                                SPATIAL=spt, `ENETWILD'22`=wb_last)
M<-cor(df.g, method = "spearman")
M.p<-cor.mtest(df.g, conf.level=0.95)
png(filename = "Figure5.2B_Correlation.png",
       width=8,
       height=6,
       units="cm",
       res = 300)
corrplot::corrplot(M, 
                   method= 'shade',
                   type="lower", 
                   diag=FALSE,
                   sig.level = -1,
                   # insig = "p-value",
                   tl.col="#011f4b",
                   tl.cex = 0.5,
                   tl.srt = 0,
                   tl.offset = 0.6,
                   addCoef.col = "white",
                   col=colorRampPalette(c("white", "#b3cde0", "#005b96", "#011f4b"))(5), 
                   cl.cex=0.7, 
                   cl.offset = 0.2,
                   cl.ratio=0.5)
dev.off()


df.es<-df %>% filter(CNTR_CODE=="ES") %>% select(env, spt, wb_last) %>% mutate(wb_last=wb_last*100) %>% rename(environ=env, spatial=spt, lastReport=wb_last)
df.fr<-df %>% filter(CNTR_CODE=="FR") %>% select(env, spt, wb_last) %>% mutate(wb_last=wb_last*100) %>% rename(environ=env, spatial=spt, lastReport=wb_last)
df.it<-df %>% filter(CNTR_CODE=="IT") %>% select(env, spt, wb_last) %>% mutate(wb_last=wb_last*100) %>% rename(environ=env, spatial=spt, lastReport=wb_last)
df.pl<-df %>% filter(CNTR_CODE=="PL") %>% select(env, spt, wb_last) %>% mutate(wb_last=wb_last*100) %>% rename(environ=env, spatial=spt, lastReport=wb_last)

M.es<-cor(df.es, method = "spearman")
M.fr<-cor(df.fr, method = "spearman")
M.it<-cor(df.it, method = "spearman")
M.pl<-cor(df.pl, method = "spearman")
p.es<-cor.mtest(df.es, conf.level=0.95)
p.fr<-cor.mtest(df.fr, conf.level=0.95)
p.it<-cor.mtest(df.it, conf.level=0.95)
p.pl<-cor.mtest(df.pl, conf.level=0.95)

png(filename = "Figure5.2B_Correlation.plot.es.png",
    width=8,
    height=6,
    units="cm",
    res = 300)
corrplot::corrplot(M.es, 
                   method= 'shade',
                   type="lower", 
                   diag=FALSE,
                   sig.level = -1,
                   # insig = "p-value",
                   tl.col="#011f4b",
                   tl.cex = 0.5,
                   tl.srt = 0,
                   tl.offset = 0.6,
                   addCoef.col = "white",
                   col=colorRampPalette(c("white", "#b3cde0", "#005b96", "#011f4b"))(5), 
                   cl.cex=0.7, 
                   cl.offset = 0.2,
                   cl.ratio=0.5)
dev.off()

png(filename = "Figure5.2B_Correlation.plot.pl.png",
    width=8,
    height=6,
    units="cm",
    res = 300)
corrplot::corrplot(M.pl, 
                   method= 'shade',
                   type="lower", 
                   diag=FALSE,
                   sig.level = -1,
                   # insig = "p-value",
                   tl.col="#011f4b",
                   tl.cex = 0.5,
                   tl.srt = 0,
                   tl.offset = 0.6,
                   addCoef.col = "white",
                   col=colorRampPalette(c("white", "#b3cde0", "#005b96", "#011f4b"))(5), 
                   cl.cex=0.7, 
                   cl.offset = 0.2,
                   cl.ratio=0.5)
dev.off()

png(filename = "Figure5.2B_Correlation.plot.fr.png",
    width=8,
    height=6,
    units="cm",
    res = 300)
corrplot::corrplot(M.fr, 
                   method= 'shade',
                   type="lower", 
                   diag=FALSE,
                   sig.level = -1,
                   # insig = "p-value",
                   tl.col="#011f4b",
                   tl.cex = 0.5,
                   tl.srt = 0,
                   tl.offset = 0.6,
                   addCoef.col = "white",
                   col=colorRampPalette(c("white", "#b3cde0", "#005b96", "#011f4b"))(5), 
                   cl.cex=0.7, 
                   cl.offset = 0.2,
                   cl.ratio=0.5)
dev.off()


png(filename = "Figure5.2B_Correlation.plot.it2.png",
    width=8,
    height=6,
    units="cm",
    res = 300)
corrplot::corrplot(M.it, 
                   method= 'shade',
                   type="lower", 
                   diag=FALSE,
                   sig.level = -1,
                   insig = "p-value",
                   tl.col="#011f4b",
                   tl.cex = 0.5,
                   tl.srt = 0,
                   tl.offset = 0.6,
                   addCoef.col = "white",
                   col=colorRampPalette(c("white", "#b3cde0", "#005b96", "#011f4b"))(5), 
                   cl.cex=0.7, 
                   cl.offset = 0.2,
                   cl.ratio=0.5)
dev.off()
