library(tidyverse)

lcms_long_pos = read.csv("1-data/lcms/LCMS_POS_long.csv")
lcms_long_neg = read.csv("1-data/lcms/LCMS_NEG_long.csv")

lcms_combine_long = function(lcms_long_pos, lcms_long_neg){
  
  combined = 
    lcms_long_neg %>% 
    mutate(mode = "negative") %>% 
    bind_rows(lcms_long_pos %>% mutate(mode = "positive")) %>% 
    dplyr::select(feature, formula, intensity, sample_name, orbitrap_id, sample_type, site, transect_location, horizon)
  
  combined_trt = 
    combined %>% 
    dplyr::select(-feature, -intensity) %>% 
    distinct()
  
}
lcms_make_metadata = function(lcms_long_pos, lcms_long_neg){
  
  lcms_metadata = 
    lcms_long_neg %>% 
    mutate(mode = "negative") %>% 
    bind_rows(lcms_long_pos %>% mutate(mode = "positive")) %>% 
    dplyr::select(formula, C, H, O, N, S, P, HC, OC, AImod, NOSC, contains("Class"), mode) %>% 
    distinct() %>% 
    group_by(formula) %>% 
    dplyr::mutate(mode = paste(mode, collapse = ", ")) %>% 
    distinct()
  
}

lcms_compute_relabund = function(lcms_long_combined, lcms_metadata){
  
  ra_id <- 
    lcms_long_combined %>%
    group_by(sample_type, site, transect_location, horizon, Class_detailed) %>%
    dplyr::summarize(inten = sum(intensity, na.rm = TRUE), .groups = "drop") %>%
    group_by(sample_name, orbitrap_id, sample_type, site) %>%
    dplyr::mutate(total = sum(inten)) %>%
    mutate(rel_abund = (inten/total)*100) %>%
    ungroup() %>%
    mutate(Class_chem = factor(Class_chem, levels = c("lipid", "amino sugar", "protein", "carbohydrate", "unsatHC", "lignin", "tannin", "condensed HC", "other"))) %>%
    mutate(short_id = str_extract_all(orbitrap_id, "C\\d+")) %>%
    mutate(short_id = str_extract_all(short_id, "\\d+")) %>%
    mutate(short_id = ifelse(short_id %in% c(1:9), gsub("C", "C0", orbitrap_id), orbitrap_id)) %>%
    rename(plot_id = short_id)
  
  
}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

