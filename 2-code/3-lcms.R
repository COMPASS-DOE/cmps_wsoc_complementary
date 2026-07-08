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

old = function(){
  
  
  combined_trt %>%
    dplyr::select(formula, HC, OC, contains("Class"), mode) %>% 
    distinct() %>% 
    ggplot(aes(x = OC, y = HC, color = mode)) +
    geom_point() +
    theme_bw()
  
  
  
  # Visualization ################################################################
  mode <- "POS"
  ORBI_LONG_MEAN = "1-data/lcms/LCMS_POS_long.csv"
  #ORBI_LONG_MEAN = paste0("1-data/lcms/JK/", mode, "/Kaizad_LCMS_", mode, "_long.csv")
  ft_long <- read.csv(ORBI_LONG_MEAN, row.names = 1)
  #plot_dir <- paste0(wd, "/", mode, "/plots")
  
  # Prelim plot-------------------------------------------------------------------
  plot <- ggplot(ft_long %>%
                   arrange(sample_name) %>%
                   mutate(Class_chem = factor(Class_chem, levels = c("lipid", "amino sugar", "protein", "carbohydrate", "unsatHC", "lignin", "tannin", "condensed HC", "other"))),
                 aes(x = OC, y = HC, fill = Class_chem, size = intensity), 
                 # color = "black", alpha = 0.2
  ) +
    geom_point(shape = 21, show.legend = c(alpha = FALSE)) +
    facet_wrap(~sample_name) +
    theme_bw() +
    scale_fill_brewer(palette = "RdYlBu") +
    labs(size = "Peak\nintensity", fill = "")
  
  plot
  
  #ggsave("chem_class_plot3.png", path = plot_dir)
  #-------------------------------------------------------------------------------
  
  # Total count-------------------------------------------------------------------
  ft_tot_count <- ft_long %>%
    group_by(sample_name, orbitrap_id, sample_type, site) %>%
    summarise(total_count = n()) %>%
    ungroup() %>%
    mutate(short_id = str_extract_all(orbitrap_id, "C\\d+")) %>%
    mutate(short_id = str_extract_all(short_id, "\\d+")) %>%
    mutate(short_id = ifelse(short_id %in% c(1:9), gsub("C", "C0", orbitrap_id), orbitrap_id)) %>%
    rename(plot_id = short_id)
  
  
  tot_count <- ggplot(ft_tot_count, aes(x = sample_type, y = total_count, fill = sample_type)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(x = "", y = "# of features") +
    facet_wrap(~site, scales = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
  
  tot_count
  
  #ggsave("tot_count.png", path = plot_dir)
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  
  # Rel. Chem.--------------------------------------------------------------------
  # relative abundance
  ra_id <- ft_long %>%
    group_by(sample_name, orbitrap_id, sample_type, site, Class_chem) %>%
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
  
  
  ra_id_plot <- ggplot(ra_id, aes(x = plot_id, y = rel_abund, fill = Class_chem)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~site, scales = "free_x", ncol = 6) +
    labs(x = "",
         y = "Relative abundance (%)")+
    scale_fill_brewer(palette = "RdYlBu") +
    theme_bw() +
    theme(legend.title = element_blank()) +
    theme(axis.text.x = element_text(angle = 90)) 
  
  ra_id_plot
  
  #ggsave("ra_chem.png", path = plot_dir)
  #-------------------------------------------------------------------------------
  
  # Rel. El.----------------------------------------------------------------------
  long_for <- ft_long %>%
    left_join(orbi_all %>% select(formula) %>% rownames_to_column(var = "feature")) %>%
    mutate(mol_cl = str_extract_all(formula, "[A-Z]")) %>%
    mutate(mol_cl = sapply(mol_cl, paste0, collapse = ""))
  
  ra_mol <- long_for %>%
    group_by(sample_name, orbitrap_id, sample_type, site, mol_cl) %>%
    dplyr::summarize(inten = sum(intensity, na.rm = TRUE), .groups = "drop") %>%
    group_by(sample_name, orbitrap_id, sample_type, site) %>%
    dplyr::mutate(total = sum(inten)) %>%
    mutate(rel_abund = (inten/total)*100) %>%
    ungroup() %>%
    mutate(short_id = str_extract_all(orbitrap_id, "C\\d+")) %>%
    mutate(short_id = str_extract_all(short_id, "\\d+")) %>%
    mutate(short_id = ifelse(short_id %in% c(1:9), gsub("C", "C0", orbitrap_id), orbitrap_id)) %>%
    rename(plot_id = short_id)
  
  ra_mol_plot <- ggplot(ra_mol, aes(x = plot_id, y = rel_abund, fill = mol_cl)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~site, scales = "free") +
    labs(x = "",
         y = "Relative abundance (%)")+
    scale_fill_manual(name = "Molecular class", values = c("CHO" = "#E16A86", "CHON" = "#909800", "CHONS" = "#00AD9A", "CHOS" = "#9183E6")) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    theme(axis.text.x = element_text(angle = 90)) 
  
  ra_mol_plot
  
  ggsave("ra_mol.png", path = plot_dir)
  #-------------------------------------------------------------------------------
  
  library(patchwork)
  tot_count / TIC / ra_id_plot / ra_mol_plot + plot_layout(heights = unit(c(1, 1, 3, 3), c('cm')))
  
  
  # Molecular properties ---------------------------------------------------------
  ft_long2 <- ft_long %>%
    mutate(short_id = str_extract_all(orbitrap_id, "C\\d+")) %>%
    mutate(short_id = str_extract_all(short_id, "\\d+")) %>%
    mutate(short_id = ifelse(short_id %in% c(1:9), gsub("C", "C0", orbitrap_id), orbitrap_id)) %>%
    rename(plot_id = short_id) %>%
    mutate(mz = str_extract(feature, "[0-9.]+$")) %>%
    mutate(mz = as.numeric(mz)) %>%
    filter(NOSC > -2 & NOSC < 2 & AImod > -0.75 & AImod < 1 & DBE > -2.5 & DBE < 20)
  
  long_dat <- pivot_longer(ft_long2, cols = c(mz, NOSC, AImod, DBE), names_to = "variable", values_to = "val") %>%
    mutate(variable = factor(variable, levels = c("mz", "NOSC", "AImod", "DBE")))
  
  
  ggplot(long_dat,
         aes(x = plot_id, y = val, fill = sample_type)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.2, fill = "lightgrey", outlier.shape = NA) +
    facet_grid(variable ~ site, scales = "free", switch = "y") +
    theme_bw() +
    labs(x = "", y = "") +
    theme(strip.placement = "outside",
          panel.spacing = unit(0, "lines"),
          axis.text.x = element_text(angle = 90))
  #-------------------------------------------------------------------------------
  
  
  
  # PCoA -------------------------------------------------------------------------
  #md_mean <- read.csv(paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_md_mean.csv"), row.names = 1)
  md_mean = read_csv("1-data/lcms/LCMS_POS_wide.csv")
  #norm_tic <- read.csv(paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_mean_norm_tic.csv"), row.names = 1)
  
  # Load libraries for PCoA
  pacman::p_load("BiocManager","ComplexHeatmap","dendextend","NbClust","cowplot",
                 "gtools", # mixedsort
                 "tidyverse",
                 "janitor", # row_to_names function
                 "factoextra", # PCA
                 "vegan") # PCoA
  
  md_sub <- md_mean
  
  # use imp0 to filter out all zero features
  
  #norm_mean2 <- norm_tic[which(rownames(norm_tic) %in% rownames(md_sub)), , drop = F]
  
  filter_all_zero_features <- function(data){
    # Remove columns where all values are NA
    data <- data[, colSums(data) > 0]
    return(data)
  }
  
  norm_mean3 <- filter_all_zero_features(md_sub)
  
  cleaned_data <- norm_mean3
  md <- md_sub %>% column_to_rownames("...1")
  
  # Verifying file consistency----------------------------------------------------
  md <- md[mixedorder(rownames(md)), ,drop=F] # ordering the md by row names
  #cleaned_data <- cleaned_data[mixedorder(rownames(cleaned_data)),, drop=F] # ordering the md by row names
  cleaned_data = md
  
  identical(rownames(cleaned_data), rownames(md))
  
  # which file names in the metadata are not in the feature table?
  setdiff(rownames(cleaned_data), rownames(md))
  
  table(rownames(md) %in% rownames(cleaned_data))
  
  ################################################################################
  
  # PCoA
  distm <- vegdist(cleaned_data, method = "bray", na.rm = FALSE) #compute distance
  
  PcoA <- cmdscale(distm, k = 2, eig = T, add = T)
  
  PcoA_points <- as.data.frame(PcoA$points) #getting the PCOs into dataframe
  variance <- round(PcoA$eig*100/sum(PcoA$eig),1) # getting the variance explained by each PCo
  names(PcoA_points)[1:2] <- paste0('PCoA', seq(1,2)) 
  
  colnames(md)
  interested_attribute_pcoa = 'sample_type'
  
  key = read.csv("1-data/lcms/orbitrap_sample_key.csv")
  
  PcoA_points = PcoA_points %>% rownames_to_column("orbitrap_id") %>% left_join(key)
  md = md %>% rownames_to_column("orbitrap_id") %>% left_join(key)
  
  
  ggplot(PcoA_points, 
         aes(x = PCoA1, 
             y = PCoA2,
             colour = as.factor(md[,interested_attribute_pcoa]), 
             label = PcoA_points$plot_id)) +
    geom_point(size=2.5, alpha = 0.5) +
    #geom_text_repel() +
    scale_colour_manual(values = c('orange','darkgreen','red','blue','black','grey','purple', 'skyblue3','magenta','green')) +
    xlab(paste('PCoA1',variance[1],'%', sep = ' ')) + 
    ylab(paste('PCoA2',variance[2],'%', sep = ' ')) + 
    labs(color = '') +
    coord_fixed() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16, face= 'bold'),
          plot.title = element_text(size = 18, face= 'bold',hjust=0.5),
          legend.title = element_text(size = 18, face= 'bold'),
          legend.text = element_text(size = 16),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
  
  #ggsave("PCoA.png", path = plot_dir)
  # ------------------------------------------------------------------------------
  
  
  
}
