library(tidyverse)

lcms_long_pos = read.csv("1-data/lcms/LCMS_POS_long.csv")
lcms_long_neg = read.csv("1-data/lcms/LCMS_NEG_long.csv")

lcms_combine_long = function(lcms_long_pos, lcms_long_neg){
  
  combined = 
    lcms_long_neg %>% 
    mutate(mode = "negative") %>% 
    bind_rows(lcms_long_pos %>% mutate(mode = "positive")) %>% 
    dplyr::select(mode, feature, formula, intensity, sample_name, orbitrap_id, sample_type, site, transect_location, horizon)
  
  combined_trt = 
    combined %>% 
    #dplyr::select(-feature, -intensity) %>% 
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
  ## two different ways to calculate relative abundance
  
  ## 1. based on intensities - sum all intensities within each molecular class
  ## do this for each core, then summarize by treatment
  relabund_intensities_ALL_CORES <- 
    lcms_long_combined %>%
    left_join(lcms_metadata) %>% 
    group_by(mode, orbitrap_id, sample_type, site, transect_location, horizon, Class_detailed) %>%
    dplyr::summarize(inten = sum(intensity, na.rm = TRUE), .groups = "drop") %>%
    group_by(mode, orbitrap_id) %>%
    dplyr::mutate(total = sum(inten)) %>%
    mutate(rel_abund = (inten/total)*100) %>%
    ungroup()
  
  relabund_intensities_TREATMENT = 
    relabund_intensities_ALL_CORES %>% 
    group_by(mode, sample_type, site, transect_location, horizon, Class_detailed) %>% 
    dplyr::summarize(mean = mean(rel_abund, na.rm = TRUE), .groups = "drop") %>% 
    ungroup()
  
  ## 2. based on presence/absence of molecular assignments
  ## do this for each core, then summarize by treatment
  relabund_presence_ALL_CORES <- 
    lcms_long_combined %>%
    mutate(presence = 1) %>% 
    distinct(mode, formula, orbitrap_id, sample_type, site, transect_location, horizon, presence) %>% 
    left_join(lcms_metadata) %>% 
    group_by(mode, orbitrap_id, sample_type, site, transect_location, horizon, Class_detailed) %>%
    dplyr::summarize(inten = sum(presence, na.rm = TRUE), .groups = "drop") %>%
    group_by(mode, orbitrap_id) %>%
    dplyr::mutate(total = sum(inten)) %>%
    mutate(rel_abund = (inten/total)*100) %>%
    ungroup()
  
  relabund_presence_TREATMENT = 
    relabund_presence_ALL_CORES %>% 
    group_by(mode, sample_type, site, transect_location, horizon, Class_detailed) %>% 
    dplyr::summarize(mean = mean(rel_abund, na.rm = TRUE), .groups = "drop") %>% 
    ungroup()
  
  list(relabund_intensities_ALL_CORES = relabund_intensities_ALL_CORES,
       relabund_intensities_TREATMENT = relabund_intensities_TREATMENT,
       relabund_presence_ALL_CORES = relabund_presence_ALL_CORES,
       relabund_presence_TREATMENT = relabund_presence_TREATMENT)
}

lcms_wide = function(lcms_long_combined){
  
  wide_pos_features = 
    lcms_long_combined %>% 
    filter(mode == "positive") %>% 
    dplyr::select(-formula) %>% 
    pivot_wider(names_from = "feature", values_from = "intensity") %>% 
    replace(is.na(.), 0)
  
  wide_pos_formula = 
    lcms_long_combined %>% 
    filter(mode == "positive") %>%
    distinct(formula, orbitrap_id, sample_type, site, transect_location, horizon) %>% 
    mutate(presence = 1) %>%  
    pivot_wider(names_from = "formula", values_from = "presence") %>% 
    replace(is.na(.), 0)
  
  wide_neg_features = 
    lcms_long_combined %>% 
    filter(mode == "negative") %>% 
    dplyr::select(-formula) %>% 
    pivot_wider(names_from = "feature", values_from = "intensity") %>% 
    replace(is.na(.), 0)
  
  wide_neg_formula = 
    lcms_long_combined %>% 
    filter(mode == "negative") %>%
    distinct(formula, orbitrap_id, sample_type, site, transect_location, horizon) %>% 
    mutate(presence = 1) %>%  
    pivot_wider(names_from = "formula", values_from = "presence") %>% 
    replace(is.na(.), 0)
  
  list(wide_pos_features = wide_pos_features,
       wide_pos_formula = wide_pos_formula,
       wide_neg_features = wide_neg_features,
       wide_neg_formula = wide_neg_formula)
}



