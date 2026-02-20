

# -------------------------------------------------------------------------

## spectra
process_nmr_spectra = function(nmr_spectra, sample_key){
  
  nmr_spectra_processed = 
    nmr_spectra %>% 
    filter(ppm >= 0 & ppm <= 10) %>% 
    mutate(sampleID = str_replace(sampleID, "spectra", "COMP")) %>% 
    left_join(sample_key, by = c("sampleID" = "sample_name")) %>% 
    nmr_assign_bins(binset = bins_Lynch2019)
  
}
plot_nmr_spectra_all = function(nmr_spectra_processed){
  
  plot_spectra_indiv = function(dat){

    dat %>% 
      filter(intensity < 1) %>%
      nmr_plot_spectra(binset = bins_Lynch2019,
                       mapping = aes(group = sampleID,
                                     x = ppm,
                                     y = intensity),
                       label_position = 6.5,
                       stagger = 1.2)+
      xlim(10, 0.1)+
      facet_wrap(~sample_type)+
      #  ylim(0, 1)+
      NULL    
    
  }
  
  msm_u_o = 
    nmr_spectra_processed %>% 
    filter(site == "MSM" & transect_location == "upland" & horizon == "O") %>% 
    plot_spectra_indiv(.)
  
  msm_u_b = 
    nmr_spectra_processed %>% 
    filter(site == "MSM" & transect_location == "upland" & horizon == "B") %>% 
    plot_spectra_indiv(.)
  
  msm_t_o = 
    nmr_spectra_processed %>% 
    filter(site == "MSM" & transect_location == "transition" & horizon == "O") %>% 
    plot_spectra_indiv(.)
  
  msm_t_b = 
    nmr_spectra_processed %>% 
    filter(site == "MSM" & transect_location == "transition" & horizon == "B") %>% 
    plot_spectra_indiv(.)
  
  msm_w = 
    nmr_spectra_processed %>% 
    filter(site == "MSM" & transect_location == "wetland") %>% 
    plot_spectra_indiv(.)
  
  owc_u_a = 
    nmr_spectra_processed %>% 
    filter(site == "OWC" & transect_location == "upland" & horizon == "A") %>% 
    plot_spectra_indiv(.)
  
  owc_u_b = 
    nmr_spectra_processed %>% 
    filter(site == "OWC" & transect_location == "upland" & horizon == "B") %>% 
    plot_spectra_indiv(.)
  
  owc_t_a = 
    nmr_spectra_processed %>% 
    filter(site == "OWC" & transect_location == "transition" & horizon == "A") %>% 
    plot_spectra_indiv(.)
  
  owc_t_b = 
    nmr_spectra_processed %>% 
    filter(site == "OWC" & transect_location == "transition" & horizon == "B") %>% 
    plot_spectra_indiv(.)
  
  owc_w = 
    nmr_spectra_processed %>% 
    filter(site == "OWC" & transect_location == "wetland") %>% 
    plot_spectra_indiv(.)
    
  
  
  gg_msm = cowplot::plot_grid(msm_u_o, msm_t_o, msm_w, msm_u_b, msm_t_b)
  gg_owc = cowplot::plot_grid(owc_u_a, owc_t_a, owc_w, owc_u_b, owc_t_b)
  
  list(gg_msm = gg_msm,
       gg_owc = gg_owc)
  
}
plot_nmr_spectra_subset = function(nmr_spectra_processed){
  
  subset = 
    nmr_spectra_processed %>% 
    mutate(sample_number = parse_number(sampleID)) %>% 
    group_by(site, horizon, transect_location) %>% 
    filter(sample_number == min(sample_number))
  
  
  gg_msm_subset = 
    subset %>%
    filter(site == "MSM") %>% 
    filter(intensity < 1) %>%
    nmr_plot_spectra(binset = bins_Lynch2019,
                     mapping = aes(group = sample_type,
                                   color = sample_type,
                                   x = ppm,
                                   y = intensity),
                     label_position = 3.5,
                     stagger = 0.5)+
    xlim(10, 0.1)+
    theme(legend.position = c(0.2, 0.8))+
    NULL    
  
  gg_owc_subset = 
    subset %>%
    filter(site == "OWC") %>% 
    filter(intensity < 1) %>%
    nmr_plot_spectra(binset = bins_Lynch2019,
                     mapping = aes(group = sample_type,
                                   color = sample_type,
                                   x = ppm,
                                   y = intensity),
                     label_position = 3.5,
                     stagger = 0.5)+
    xlim(10, 0.1)+
    theme(legend.position = c(0.2, 0.8))+
    NULL    

  cowplot::plot_grid(gg_msm_subset, gg_owc_subset)  
  
}    


## peaks
process_nmr_peaks = function(nmr_peaks){
  
  nmr_peaks_processed = 
    nmr_peaks %>% 
    filter(ppm >= 0 & ppm <= 10) %>% 
    mutate(sampleID = str_replace(sampleID, "peaks", "COMP")) %>% 
    filter(!`Impurity/Compound` %in% "Solvent") %>% 
    filter(`Impurity/Compound` %in% "Compound") %>% 
    nmr_assign_bins(binset = bins_Lynch2019)
  
}

## relabund
plot_nmr_relabund = function(nmr_relabundance, sample_key){
  
  nmr_relabundance %>% 
    left_join(sample_key, by = c("sampleID" = "sample_name")) %>% 
    ggplot(aes(x = sampleID, y = relabund, fill = group))+
    geom_bar(stat = "identity")+
    facet_wrap(~site + transect_location + horizon, scale = "free_x")
  
}

## stats
compute_nmr_permanova = function(nmr_relabundance, sample_key){
  
  relabund_wide = pivot_wider(nmr_relabundance, names_from = "group", values_from = "relabund") %>% 
    left_join(sample_key, by = c("sampleID" = "sample_name")) %>% 
    filter(horizon != "B")
  
  
  permanova = 
    adonis2(relabund_wide %>% dplyr::select(where(is.numeric))  ~ (site + transect_location + horizon)^2,
           data = relabund_wide)
  
  broom::tidy(permanova$aov.tab)
  
}
compute_nmr_pca = function(nmr_relabundance, sample_key){
  
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
    pivot_wider(nmr_relabundance, names_from = "group", values_from = "relabund") %>% 
    left_join(sample_key, by = c("sampleID" = "sample_name")) %>% 
    filter(horizon != "B")
  
  pca = fit_pca_function(relabund_wide)
  
  ## plot PCA
  ggbiplot(pca$pca_int, obs.scale = 1, var.scale = 1,
           groups = as.character(pca$grp$transect_location), 
           ellipse = FALSE, circle = FALSE, 
           var.axes = TRUE, alpha = 0) +
    geom_convexhull(aes(group = groups,
                        fill = groups), 
                    alpha = 0.2)+
    geom_point(size = 5,stroke = 1, #alpha = 0.5,
               aes(shape = as.character(pca$grp$site),
                   fill = groups))+
    labs(shape = "Site",
         title = "Surface horizons only")+
    scale_shape_manual(values = c(21, 23))+
    scale_fill_manual(values = pnw_palette("Bay", 3))+
    guides(fill = guide_legend(override.aes = list(shape = 21)))
    

  }