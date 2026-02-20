library(tidyverse)
theme_set(theme_bw(base_size = 14))

#############################
#############################

# 1. PROCESSING FUNCTIONS -------------------------------------------------
## LEVEL I FUNCTIONS -------------------------------------------------------
## for metadata file
apply_filter_report = function(report){
  report %>% 
    # filter appropriate mass range
    filter(Mass>200 & Mass<800) %>% 
    # remove isotopes
    filter(C13==0) %>% 
    # remove peaks without C assignment
    filter(C>0)
}
compute_indices = function(dat){
  dat %>% 
    dplyr::select(Mass, C:P) %>% 
    dplyr::mutate(AImod = round((1+C-(0.5*O)-S-(0.5*(N+P+H)))/(C-(0.5*O)-S-N-P),4),
                  NOSC =  round(4-(((4*C)+H-(3*N)-(2*O)-(2*S))/C),4),
                  HC = round(H/C,2),
                  OC = round(O/C,2),
                  DBE_AI = 1+C-O-S-0.5*(N+P+H),
                  DBE =  1 + ((2*C-H + N + P))/2,
                  DBE_C = round(DBE_AI/C,4)) %>% 
    dplyr::select(-c(C:P))
}
compute_mol_formula = function(dat){
  dat %>% 
    dplyr::select(Mass, C:P) %>% 
    dplyr::mutate(formula_c = if_else(C>0,paste0("C",C),as.character(NA)),
                  formula_h = if_else(H>0,paste0("H",H),as.character(NA)),
                  formula_o = if_else(O>0,paste0("O",O),as.character(NA)),
                  formula_n = if_else(N>0,paste0("N",N),as.character(NA)),
                  formula_s = if_else(S>0,paste0("S",S),as.character(NA)),
                  formula_p = if_else(P>0,paste0("P",P),as.character(NA)),
                  formula = paste0(formula_c,formula_h, formula_o, formula_n, formula_s, formula_p),
                  formula = str_replace_all(formula,"NA","")) %>% 
    dplyr::select(Mass, formula, C:P)
}
assign_class_seidel = function(meta_clean, meta_indices){
  meta_clean %>%
    left_join(meta_indices, by = "Mass") %>% 
    mutate(Class = case_when(AImod>0.66 ~ "condensed aromatic",
                             AImod<=0.66 & AImod > 0.50 ~ "aromatic",
                             AImod <= 0.50 & HC < 1.5 ~ "unsaturated/lignin",
                             HC >= 1.5 ~ "aliphatic"),
           Class = replace_na(Class, "other"),
           Class_detailed = case_when(AImod>0.66 ~ "condensed aromatic",
                                      AImod<=0.66 & AImod > 0.50 ~ "aromatic",
                                      AImod <= 0.50 & HC < 1.5 ~ "unsaturated/lignin",
                                      HC >= 2.0 & OC >= 0.9 ~ "carbohydrate",
                                      HC >= 2.0 & OC < 0.9 ~ "lipid",
                                      HC < 2.0 & HC >= 1.5 & N==0 ~ "aliphatic",
                                      HC < 2.0 & HC >= 1.5 & N > 0 ~ "aliphatic+N")) %>% 
    dplyr::select(Mass, Class, Class_detailed)
}

## for data file
compute_presence = function(dat){
  dat %>% 
    pivot_longer(-("Mass"), values_to = "presence", names_to = "analysis_ID") %>% 
    # convert intensities to presence==1/absence==0  
    dplyr::mutate(presence = if_else(presence>0,1,0)) %>% 
    # keep only peaks present
    filter(presence>0)
}
apply_replication_filter = function(data_long_key, ...){
  max_replicates = 
    data_long_key %>% 
    ungroup() %>% 
    group_by(...) %>% 
    distinct(sample_name) %>% 
    dplyr::summarise(reps = n())
  
  # second, join the `reps` file to the long_key file
  # and then use the replication filter  
  data_long_key %>% 
    group_by(formula, ...) %>% 
    dplyr::mutate(n = n()) %>% 
    left_join(max_replicates) %>% 
    ungroup() %>% 
    mutate(keep = n >= (2/3)*reps) %>% 
    filter(keep) %>% 
    dplyr::select(-keep, -reps)
  
}


## LEVEL II FUNCTIONS ------------------------------------------------------

make_fticr_meta = function(icr_report){
  # filter peaks
  fticr_report = (apply_filter_report(icr_report))
  
  meta_clean = 
    fticr_report %>% 
    # select only the relevant columns for the formula assignments
    dplyr::select(Mass, C, H, O, N, S, P, El_comp)
  
  # calculate indices and assign classes - Seidel 2014 and Seidel 2017
  meta_indices = compute_indices(meta_clean)
  meta_formula = compute_mol_formula(meta_clean)
  meta_class = assign_class_seidel(meta_clean, meta_indices)
  
  # output
  meta2 = meta_formula %>% 
    left_join(meta_class, by = "Mass") %>% 
    left_join(meta_indices, by = "Mass") %>% dplyr::select(-Mass) %>% distinct(.)
  
  list(meta2 = meta2,
       meta_formula = meta_formula)
}
make_fticr_data = function(icr_report, sample_key){
  
  # a. convert intensities to compute presence/absence
  # b. pull samples vs. blanks vs. instrument blanks
  # c. do a blank correction 
  # d. apply replication filter (keep only peaks seen in > 2/3 of all replicates)
  # e. outputs: 
  # -- (1) peaks per sample -- use for relative abundance and stats, 
  # -- (2) peaks per treatment -- use for Van Krevelen plots and unique peak analysis
  
  
  fticr_report = (apply_filter_report(icr_report))
  mass_to_formula = make_fticr_meta(icr_report)$meta_formula
  
  data_columns = fticr_report %>% dplyr::select(Mass, starts_with("DOC"))
  
  data_presence = 
    compute_presence(data_columns) %>% 
    left_join(mass_to_formula, by = "Mass") %>% 
    dplyr::select(formula, analysis_ID, presence) 
  
  data_long_key = 
    data_presence %>% 
    mutate(analysis_ID = str_extract(analysis_ID, "DOC_CMPS_KFP_[0-9]{4}")) %>% 
    left_join(sample_key)
  
  # apply replication filter (to the samples only)
  data_long_key_repfiltered = 
    data_long_key %>% 
    filter(!is.na(region)) %>% 
    apply_replication_filter(., region, site, transect_location, horizon)
  
  # blanks - filter and solution
  blanks = 
    data_long_key %>% 
    filter(grepl("blank", sample_name))
  
  # instrument blanks - raw blanks, etc.
  instrument_blanks = 
    data_long_key %>% 
    filter(grepl("blank", sample_type))
  
  # do a blank correction
  blank_peaks = 
    blanks %>% 
    dplyr::distinct(formula) %>% 
    mutate(blank = TRUE)
  
  data_long_blank_corrected = 
    data_long_key_repfiltered %>% 
    left_join(blank_peaks) %>% 
    filter(is.na(blank)) %>% 
    dplyr::select(-blank)
  
  # get summary of peaks by treatment
  data_long_trt = 
    data_long_blank_corrected %>% 
    distinct(formula, region, site, transect_location, horizon)
  
  list(data_long_trt = data_long_trt,
       data_long_blank_corrected = data_long_blank_corrected)
  
}

#


#############################
#############################

#
## DATA FUNCTIONS ----------------------------------------------------------


icr_process_metadata = function(icr_report_negative, icr_report_positive){
  
  meta_negative = make_fticr_meta(icr_report = icr_report_negative)$meta2
  meta_positive = make_fticr_meta(icr_report = icr_report_positive)$meta2
  
  meta_combined = 
    meta_negative %>% mutate(mode = "negative") %>% 
    bind_rows(meta_positive %>% mutate(mode = "positive")) %>% 
    group_by(formula) %>% 
    mutate(n = n()) %>% 
    mutate(mode = case_when(n == 2 ~ "negative, positive", TRUE ~ mode)) %>% 
    distinct()
}

icr_process_data_trt = function(icr_report_negative, icr_report_positive, sample_key){
  
  data_negative = make_fticr_data(icr_report_negative, sample_key)$data_long_trt
  data_positive = make_fticr_data(icr_report_positive, sample_key)$data_long_trt
  
  icr_data_trt = 
    data_negative %>% mutate(mode = "negative") %>% 
    bind_rows(data_positive %>% mutate(mode = "positive"))
  
}

icr_process_data_cores = function(icr_report_negative, icr_report_positive, sample_key){
  
  data_negative_longform = make_fticr_data(icr_report_negative, sample_key)$data_long_blank_corrected
  data_positive_longform = make_fticr_data(icr_report_positive, sample_key)$data_long_blank_corrected
  
  icr_data_cores = 
    data_negative_longform %>% mutate(mode = "negative") %>% 
    bind_rows(data_positive_longform %>% mutate(mode = "positive")) %>% 
    dplyr::select(-analysis_ID, -n)
  
}


compute_icr_relabund = function(icr_data_cores, icr_metadata){
  
  icr_relabundance = 
    icr_data_cores %>% 
    #  filter(horizon != "B") %>% 
    # add the Class column to the data
    left_join(dplyr::select(icr_metadata, formula, Class), by = "formula") %>% 
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
