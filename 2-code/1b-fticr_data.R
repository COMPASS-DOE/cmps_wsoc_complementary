


# icr_report_positive = read.csv("1-data/fticr/DOC 367-420 Positive Report.csv")
# icr_report_negative = read.csv("1-data/fticr/DOC 367-420 Negative Report.csv")
# icr_sample_key = read.csv("1-data/sample_key.csv", na = "") %>% dplyr::select(-notes)

icr_process_metadata = function(icr_report_negative, icr_report_positive){
  
  meta_negative = make_fticr_meta(icr_report = icr_report_negative)$meta2
  meta_positive = make_fticr_meta(icr_report = icr_report_positive)$meta2
  
  meta_combined = 
    meta_negative %>% mutate(mode = "negative") %>% 
    bind_rows(meta_positive %>% mutate(mode = "positive")) %>% 
    group_by(formula) %>% 
    mutate(n = n()) %>% 
    mutate(mode = case_when(n == 2 ~ "both", TRUE ~ mode)) %>% 
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




compute_icr_pca = function(icr_relabundance, sample_key){
  
  fit_pca_function = function(dat){
    relabund_pca=
      dat %>% 
      ungroup %>% 
      #  dplyr::select(-c(abund, total)) %>% 
      #  spread(Class, relabund) %>% 
      #  replace(.,is.na(.),0)  %>% 
      dplyr::select(-1)
    
    num = 
      relabund_pca %>% 
      dplyr::select(where(is.numeric))
    
    grp = 
      relabund_pca %>% 
      dplyr::select(-where(is.numeric)) %>% 
      dplyr::mutate(row = row_number())
    
    pca_int = prcomp(num, scale. = T)
    
    list(num = num,
         grp = grp,
         pca_int = pca_int)
  }
  
  relabund_wide = 
    pivot_wider(icr_relabundance %>% dplyr::select(-abund, -total), names_from = "Class", values_from = "relabund") %>% 
    left_join(sample_key %>% dplyr::select(-analysis_ID)) %>% 
    filter(horizon != "B")
  
  pca_negative = fit_pca_function(relabund_wide %>% filter(mode == "negative"))
  
  ## plot PCA
  ggbiplot(pca_negative$pca_int, obs.scale = 1, var.scale = 1,
           groups = as.character(pca_negative$grp$transect_location), 
           ellipse = FALSE, circle = FALSE, 
           var.axes = TRUE, alpha = 0) +
    geom_convexhull(aes(group = groups,
                        fill = groups), 
                    alpha = 0.2)+
    geom_point(size = 5,stroke = 1, #alpha = 0.5,
               aes(shape = as.character(pca_negative$grp$site),
                   fill = groups))+
    labs(shape = "Site",
         title = "Surface horizons only")+
    scale_shape_manual(values = c(21, 23))+
    scale_fill_manual(values = pnw_palette("Bay", 3))+
    xlim(-3, 3)+
    ylim(-2,2)+
    guides(fill = guide_legend(override.aes = list(shape = 21)))
  
  
  pca_positive = fit_pca_function(relabund_wide %>% filter(mode == "positive"))
  
  ## plot PCA
  ggbiplot(pca_positive$pca_int, obs.scale = 1, var.scale = 1,
           groups = as.character(pca_positive$grp$transect_location), 
           ellipse = FALSE, circle = FALSE, 
           var.axes = TRUE, alpha = 0) +
    geom_convexhull(aes(group = groups,
                        fill = groups), 
                    alpha = 0.2)+
    geom_point(size = 5,stroke = 1, #alpha = 0.5,
               aes(shape = as.character(pca_positive$grp$site),
                   fill = groups))+
    labs(shape = "Site",
         title = "Surface horizons only")+
    scale_shape_manual(values = c(21, 23))+
    scale_fill_manual(values = pnw_palette("Bay", 3))+
    xlim(-3, 3)+
    ylim(-2,2)+
    guides(fill = guide_legend(override.aes = list(shape = 21)))
  
  
}


  



###############
###############

x = function(){
gg_vk_domains = 
  gg_vankrev(meta_negative, aes(x = OC, y = HC, color = Class))+
  scale_color_manual(values = PNWColors::pnw_palette("Sunset2"))

gg_vankrev(meta_positive, aes(x = OC, y = HC, color = Class))+
  scale_color_manual(values = PNWColors::pnw_palette("Sunset2"))


gg_vankrev(meta_combined, aes(x = OC, y = HC, color = mode))+
  facet_wrap(~mode)

###############
###############




neg_unique = 
  data_negative %>% 
  group_by(formula, region, site, horizon) %>% 
  mutate(n = n())

neg_unique %>% 
  left_join(meta_combined %>% dplyr::select(formula, HC, OC)) %>% 
  filter(! horizon %in% "B") %>% 
  gg_vankrev(aes(x = OC, y = HC, color = as.character(n)))+
  facet_wrap(~region + n)


neg_unique %>% 
  left_join(meta_combined %>% dplyr::select(formula, HC, OC)) %>% 
  filter(! horizon %in% "B") %>% 
  filter(n == 1) %>% 
  gg_vankrev(aes(x = OC, y = HC, color = transect_location))+
  facet_wrap(~region + transect_location)


pos_unique = 
  data_positive %>% 
  group_by(formula, region, site, horizon) %>% 
  mutate(n = n())

pos_unique %>% 
  left_join(meta_combined %>% dplyr::select(formula, HC, OC)) %>% 
  filter(! horizon %in% "B") %>% 
  gg_vankrev(aes(x = OC, y = HC, color = as.character(n)))+
  facet_wrap(~region + n)


pos_unique %>% 
  left_join(meta_combined %>% dplyr::select(formula, HC, OC)) %>% 
  filter(! horizon %in% "B") %>% 
  filter(n == 1) %>% 
  gg_vankrev(aes(x = OC, y = HC, color = transect_location))+
  facet_wrap(~region + transect_location)


###############
###############

## RELATIVE ABUNDANCE
## 
## 

#  relabund_cores_neg = compute_relabund_cores(data_longform = data_negative_longform, meta_combined)
#  relabund_cores_pos = compute_relabund_cores(data_longform = data_positive_longform, meta_combined)

#plot_relabund(relabund_cores_neg, TREATMENTS = quos(sample_label, region, transect_location, horizon))
#plot_relabund(relabund_cores_pos, TREATMENTS = quos(sample_label, region, transect_location, horizon))

# plot_relabund(icr_relabundance %>% left_join(sample_key), TREATMENTS = quos(sample_name, region, transect_location, horizon, mode))

relabund_trt = 
  icr_relabundance %>% 
  left_join(sample_key) %>% 
  group_by(sample_name, region, transect_location, horizon, mode, Class) %>% 
  dplyr::summarize(rel_abund = round(mean(relabund),2),
                   se  = round((sd(relabund/sqrt(n()))),2),
                   relative_abundance = paste(rel_abund, "\u00b1",se)) %>% 
  ungroup() %>% 
  mutate(Class = factor(Class, levels = c("aliphatic", "unsaturated/lignin", "aromatic", "condensed aromatic")))


icr_relabundance %>% 
  left_join(sample_key) %>% 
  ggplot(aes(x = sample_name, y = relabund, fill = Class))+
  geom_bar(stat = "identity")+
  facet_wrap(~mode + region + transect_location + horizon, scales = "free_x")





}