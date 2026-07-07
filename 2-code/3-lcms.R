
lcms_report_positive = read.csv("1-data/lcms/JK/POS/Kaizad_LCMS_POS_long.csv")
lcms_report_negative = read.csv("1-data/lcms/JK/NEG/Kaizad_LCMS_NEG_long.csv")


lcms_combined_long = 
  lcms_report_negative %>% mutate(mode = "negative") %>% 
  bind_rows(lcms_report_positive %>% mutate(mode = "positive")) %>% 
  mutate(formula = str_remove_all(formula, " "))

lcms_data_long = 
  lcms_combined_long %>% 
  mutate(region = case_when(site == "MSM" ~ "CB",
                            site == "OWC" ~ "LE")) %>% 
  dplyr::select(formula, intensity, sample_name, sample_type, region, site, transect_location, horizon, mode)

lcms_metadata = 
  lcms_combined_long %>% 
  dplyr::select(-c(intensity, sample_name, sample_type, site, transect_location, horizon, orbitrap_id, X)) %>% 
  dplyr::select(-feature) %>% 
  distinct() %>% 
  group_by(formula) %>% 
  mutate(n = n()) %>% 
  mutate(mode = case_when(n == 2 ~ "negative, positive", TRUE ~ mode)) %>% 
  distinct()
#  separate(feature, sep = "_", into = c("x", "mode2", "y", "mz")) %>% 
#  dplyr::select(-c(x, mode2, y)) %>% 
  


lcms_data_trt = 
  lcms_data_long %>% 
  distinct(formula, region, site, transect_location, horizon, sample_name, mode)



compute_lcms_relabund = function(lcms_data_long, lcms_metadata){
  
  lcms_relabundance = 
    lcms_data_long %>% 
    #  filter(horizon != "B") %>% 
    # add the Class column to the data
    left_join(dplyr::select(lcms_metadata, formula, Class), by = "formula") %>% 
    # calculate abundance of each Class as the sum of all counts
    group_by(sample_name, mode, Class, region, transect_location, horizon) %>%
    dplyr::summarise(abund = sum(presence)) %>%
    ungroup %>% 
    # create a new column for total counts per core assignment
    # and then calculate relative abundance  
    group_by(mode, sample_name) %>% 
    dplyr::mutate(total = sum(abund),
                  relabund  = ((abund/total)*100))
  
  icr_relabundance_wide = 
    icr_relabundance %>% 
    ungroup() %>% 
    mutate(Class = factor(Class, 
                          levels = c("aliphatic", "unsaturated/lignin", 
                                     "aromatic", "condensed aromatic"))) %>% 
    dplyr::select(-c(abund, total)) %>% 
    pivot_wider(names_from = "Class", values_from = "relabund") %>% 
    drop_na() %>% 
    force()
  
  list(icr_relabundance = icr_relabundance,
       icr_relabundance_wide = icr_relabundance_wide)
  
}
