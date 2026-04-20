# Run the script 4_Predicted_Quantities_Plot til 81 line
a<-join_env2 %>% as_tibble() %>% 
  select(spA.m4b:NUTS_NAME, -LEVL_CODE, -fasterize) %>% 
  rename(environmental.median=spA.m4b, 
         environmental.lowCI=spA.m4b.lowCI,
         environmental.highCI=spA.m4b.highCI) %>% 
  group_by(CNTR_CODE, NAME_LATN) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
left_join(
join_spt2 %>% as_tibble() %>% 
  select(spA.m4b:NUTS_NAME, -LEVL_CODE, -fasterize) %>% 
  rename(spatial.median=spA.m4b, 
         spatial.lowCI=spA.m4b.lowCI,
         spatial.highCI=spA.m4b.highCI) %>% 
  group_by(CNTR_CODE, NAME_LATN) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)), 
by=c("CNTR_CODE", "NAME_LATN")) %>% 

left_join(
join_wb2  %>% 
  select(Pred_Red:pred_upper, NUTS_ID:NUTS_NAME, -LEVL_CODE) %>% 
  rename(ENETW22.median=Pred_Red, ENETW22.lowCI=pred_lower, ENETW22.highCI=pred_upper) %>%  
  group_by(CNTR_CODE, NAME_LATN) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  mutate(ENETW22.median=ENETW22.median*100, 
         ENETW22.lowCI=ENETW22.lowCI*100, 
         ENETW22.highCI=ENETW22.highCI*100), 
by=c("CNTR_CODE", "NAME_LATN")
)



world<-st_read("J:/IREC/politico/TM_WORLD_BORDERS_0.3_crs3035.shp") %>% rename(CNTR_CODE=ISO2) %>% select(CNTR_CODE, ISO3, NAME)


a<-a %>% left_join(world, "CNTR_CODE") %>% 
  select(CNTR_CODE, NAME_LATN, ISO3, NAME,  environmental.median:ENETW22.highCI)

xlsx::write.xlsx(data.frame(a), 
                 file="J:/IREC_Sonia/2_WIP/Report202512/4_Outputs_HYResults/PredictedQuantitiesCountrySummary.xlsx", 
                 sheetName = "CountriesSummary", 
                 row.names = FALSE)
