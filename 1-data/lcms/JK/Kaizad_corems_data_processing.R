##### Clean up the environment and memory ######################################
rm(list = ls()) # removes all objects in the environment
gc() # Runs garbage collection to free memory
################################################################################


##### Load libraries ###########################################################
library(pacman) # install, load packages

pacman::p_load("tidyverse", "data.table", "gtools", "stringr",
               "readxl", "writexl",
               "KODAMA") # normalization
################################################################################


##### Setting up WD ############################################################
getwd()
wd <- "/Users/kimj704/Github/cmps_wsoc_complementary/1-data/lcms/JK"
setwd(wd)
################################################################################


##### Input & ouput files ######################################################
prj <- "Kaizad_LCMS"
mode <- "POS"

# INPUT FILES - metadata and HPC-processed data
META = paste0(wd, "/", mode, "/input_data/", prj, "_", mode, "_master_meta.csv")
WIDEALL <- paste0(wd, "/", mode, "/input_data/", prj, "_", mode ,"_wideall.csv")

## OUTPUT FILES
SEL_ORBI_DATA = paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_inc_formula.csv") # For formula
FT_MD_MERGED <- paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_ft_md_merged.csv")

BLK_REM <- paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_blk_rem.csv")
IMP <- paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_imp.csv")
IMP0 <- paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_imp0.csv")

ORBI_FORMULAMETA_W_FEATURE = paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_formulameta_w_features.csv") # formula meta (cals) + features
HCOC = paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_HCOC.csv")

ORBI_LONG_MEAN = paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_long_mean.csv")

META_SAMPLE <- paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_meta_sample_only.csv")
################################################################################



##### Read in data #############################################################
# Read in processed_files csv file to build metadata
md <- read.csv(paste0(wd, META), row.names = 1) %>%
  rename(sample_name = "short_id") %>%
  left_join(read_xlsx(paste0(wd, "/sample_key.xlsx"))) %>%
  mutate(sample_type = ifelse(is.na(sample_type), "blank", sample_type))

# Read in the wide form data (HPC-processed)
wideall <- read.csv(WIDEALL)
################################################################################


##### Preprocessing ############################################################
# 1) Clean up the data
# 2) Blank correction (remove all features appeared in blanks)
# 3) Replicate filter (keep only the features appeared in 2 out of 3 technical replicates)

# 1) Clean up the adata---------------------------------------------------------
# def. function
InsideLevels <- function(metatable){
  LEVELS <- c() #creating empty vector to store information
  typ <-c()
  COUNT <- c()
  for(i in 1:ncol(metatable)){ # for each metadata column
    temp <- as.data.frame(table(metatable[,i])) #table function gives the category in each column and the count of each category
    x <- temp$Var1 #getting the name of each category in every column
    if(is.double(metatable[,i])==T){ # for numeric columns in metadata table, round the category values
      x=round(as.double(x),2)
    } 
    LEVELS <- rbind(LEVELS,toString(x)) # adding all the category values in a row
    COUNT <- rbind(COUNT,toString(temp$Freq)) # getting the frequency of each level in every column
    typ <- rbind(typ,class(metatable[,i])) # getting the class of each column
  }
  out <- data.frame(INDEX=c(1:ncol(metatable)), #creating an output dataframe with 1st column as INDEX
                    ATTRIBUTES=colnames(metatable), #2nd column ATTRIBUTES will be the column name of metadata table
                    LEVELS, #3rd column LEVELS will give the different categories in each ATTRIBUTE
                    COUNT, #4th column COUNT will give the number of files present with each category
                    'ATTRIBUTE_CLASS'=typ, row.names=NULL) #Final column indicating the Class or datatype of each ATTTRIBUTE
  return(out)
}

summ <- InsideLevels(md[, 2:ncol(md)]) #excluding 1st filename
summ


# Replace the row names (features) with "XID_timeblock_mz" - 'X' is needed for row names starting with numbers
ft <- wideall %>%
  mutate(id = row_number()) %>%
  mutate(mz = rowMeans(select(., starts_with("mz_")), na.rm = TRUE),
         polarity = case_when(iontype == "protonated" ~ "POS", iontype == "de-protonated" ~ "NEG")) %>%
  select(-contains("cal_mz_"), -contains("peakarea"), contains("peakarea")) %>%
  rename_with(~ gsub("peakarea_", "", .x), starts_with("peakarea_")) 

rownames(ft) <- paste(paste0("X", ft$`id`), ft$`polarity`, ft$`block`, round(ft$`mz`,digits = 3), sep = '_')

# Save for later formula processing
write.csv(ft[,c(3:12)], SEL_ORBI_DATA, row.names = rownames(ft))


# Select relevant columns
new_ft <- ft %>%
  rownames_to_column(var = "feature") %>%
  select(feature, contains("raw")) %>%
  select(!contains("mz")) %>% # raw data columns only
  column_to_rownames(var = "feature")

new_md <- md %>%
  column_to_rownames(var = "unique.raw_all.datafile.")

# Verifying file consistency----
new_ft <- new_ft[,mixedorder(colnames(new_ft)), drop=F] #ordering the ft by its column names
new_md <- new_md[mixedorder(rownames(new_md)),, drop=F] #ordering the md by row names

# how many files in the metadata are also present in the feature table
table(rownames(new_md) %in% colnames(new_ft))

identical(rownames(new_md), colnames(new_ft))

# which file names in the metadata are not in the feature table?
setdiff(rownames(new_md),colnames(new_ft))
#-------------------------------

## Data clean-up
# Transpose
ft_t <- as.data.frame(t(new_ft)) %>%
  mutate_all(as.numeric)

rownames(ft_t) <- rownames(new_md)

# Check with metadata
identical(rownames(new_md),rownames(ft_t)) #should return TRUE now

# Merge ft with md
ft_md <- merge(new_md, ft_t, by= 0, all.x=TRUE)

# Save
write.csv(ft_md, FT_MD_MERGED)



# 2) Blank correction and replicate filters ------------------------------------
# (1) Separate sample and all blank files
md_Blank <- new_md %>%
  filter(grepl("blank|procblk", sample_type))

md_Blank1 <- md_Blank[1,]

md_Blank2 <- md_Blank[2,]

md_Samples <- new_md %>%
  filter(!grepl("Blank|blank|waterblk|procblk|QC", sample_type)) %>%
  filter(!is.na(sample_type))

Blank1 <- ft_t[which(rownames(ft_t) %in% rownames(md_Blank1)), , drop = F]
Blank2 <- ft_t[which(rownames(ft_t) %in% rownames(md_Blank2)), , drop = F]
Samples <- ft_t[which(rownames(ft_t) %in% (rownames(md_Samples))), , drop=F]

# (2) Blank correction
# first, set up a blank_mask function - this will remove any features show up in blanks.
blank_mask <- function(sample_mat, blank_mat) {
  if(!all(colnames(sample_mat) == colnames(blank_mat))) {
    stop("Column names must match and be in the same order")
  }
  
  cols_to_na <- colSums(!is.na(blank_mat)) > 0
  
  sample_mat[, cols_to_na] <- NA
  
  return(sample_mat)
}

# apply the function for blank1 and then blank2 (one blank is filter blank and the other is solution blank)
Samples_blk_cor1 <- blank_mask(Samples, Blank1)
Samples_blk_cor2 <- blank_mask(Samples_blk_cor1, Blank2)

# filter if a column is all NA across all samples
# Function
filter_all_na_features <- function(data){
  # Remove columns where all values are NA
  data <- data[, colSums(!is.na(data)) > 0]
  return(data)
}

Samples_blkcor <- filter_all_na_features(Samples_blk_cor2)


# (3) correct technical replicates (remove any features appeared in < 1/3 of technical replicates)
filter_rep <- function(data, mode){
  
  data2 <- as.data.table(data)
  original_rownames <- rownames(data)
  
  sample_ids <- regmatches(original_rownames, regexpr("C[0-9]{1,2}(_[^_]+){3}", original_rownames))
  
  unique_samples <- unique(sample_ids)
  
  filtered_data <- copy(data2)
  
  pb <- txtProgressBar(min = 0, max = length(unique_samples), style = 3) #Progress bar
  
  for (i in seq_along(unique_samples)) {
    sample <- unique_samples[i]
    
    print(paste0("processing ", sample, "..."))
    
    replicate_rows <- which(grepl(sample, original_rownames))
    
    replicate_data <- filtered_data[replicate_rows]
    valid_features <- colSums(!is.na(replicate_data)) >= 2
    invalid_cols <- names(valid_features)[!valid_features]
    
    if (length(invalid_cols) > 0) {
      filtered_data[replicate_rows, (invalid_cols) := lapply(.SD, function(x) NA), .SDcols = invalid_cols]
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Restore row names
  #filtered_data[, feature := NULL]
  class(filtered_data) <- class(as.data.frame(filtered_data))
  rownames(filtered_data) <- original_rownames
  
  return(filtered_data)
}


filtered_Samples <- filter_rep(Samples_blkcor, mode)

filtered_Samples2 <- filter_all_na_features(filtered_Samples)

filtered_Blank1 <- Blank1[, colnames(filtered_Samples2), drop = FALSE]
filtered_Blank2 <- Blank2[, colnames(filtered_Samples2), drop = FALSE]


# clean up technical replicates ------------------------------------------------
# get mean values from three technical replicates

filtered_Samples2_mean <- filtered_Samples2 %>%
  rownames_to_column(var = "file") %>%
  mutate(rep = regmatches(file, regexpr("C[0-9]{1,2}(_[^_]+){3}", file))) %>%
  column_to_rownames(var = "file") %>%
  group_by(rep) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "rep")

filtered_Samples2_mean2 <- filtered_Samples2_mean %>%
  mutate_all(~ifelse(is.nan(.), NA, .))

# (3-2) Second replicate filter-------------------------------------------------
# Five samples were collected from each site/horizon (analytical reps)
# Only keep the features that are present in all five samples

filter_analytical_rep <- function(data, mode){
  
  data2 <- as.data.table(data)
  original_rownames <- rownames(data)
  
  sample_ids <- sub(".*C\\d{1,2}_", "", original_rownames)
  
  unique_samples <- unique(sample_ids)
  
  filtered_data <- copy(data2)
  
  pb <- txtProgressBar(min = 0, max = length(unique_samples), style = 3) #Progress bar
  
  for (i in seq_along(unique_samples)) {
    sample <- unique_samples[i]
    
    print(paste0("processing ", sample, "..."))
    
    replicate_rows <- which(grepl(sample, original_rownames))
    
    replicate_data <- filtered_data[replicate_rows]
    req <- nrow(replicate_data)              # the number of replicates for each sample
    valid_features <- colSums(!is.na(replicate_data)) == req
    invalid_cols <- names(valid_features)[!valid_features]
    
    if (length(invalid_cols) > 0) {
      filtered_data[replicate_rows, (invalid_cols) := lapply(.SD, function(x) NA), .SDcols = invalid_cols]
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Restore row names
  #filtered_data[, feature := NULL]
  class(filtered_data) <- class(as.data.frame(filtered_data))
  rownames(filtered_data) <- original_rownames
  
  return(filtered_data)
}

filtered_Samples2_analytical <- filter_analytical_rep(filtered_Samples2_mean2, mode)
filtered_Samples3 <- filter_all_na_features(filtered_Samples2_analytical)

filtered_Blank1_2 <- filtered_Blank1[, colnames(filtered_Samples3), drop = FALSE]
filtered_Blank2_2 <- filtered_Blank2[, colnames(filtered_Samples3), drop = FALSE]


# Save the final data frame to blk_rem
blk_rem <- filtered_Samples3 

# Save
write.csv(blk_rem, BLK_REM, row.names = TRUE)
#blk_rem <- read.csv(BLK_REM, row.names = 1)


# Imputation (taking care of NAs)===============================================
# create random values between 0 and the minimum value in our blank-removed table and randomly replace all NA/0s with these random values

#creating bins from -1 to 10^10 using sequence function seq()
bins <- c(-1,0,(1 * 10^(seq(0,10,1)))) 

#cut function cuts the give table into its appropriate bins
scores_gapfilled <- cut(as.matrix(blk_rem),bins, labels = c('0','1',paste("1E",1:10,sep="")))

#transform function convert the tables into a column format: easy for visualization
FreqTable <- transform(table(scores_gapfilled)) #contains 2 columns: "scores_x1", "Freq"
FreqTable$Log_Freq <- log(FreqTable$Freq+1) #Log scaling the frequency values
colnames(FreqTable)[1] <- 'Range_Bins' #changing the 1st colname to 'Range Bins'

## GGPLOT2
ggplot(FreqTable, aes(x=Range_Bins, y=Log_Freq)) + 
  geom_bar(stat = "identity", position = "dodge", width=0.3) + 
  ggtitle(label = "Frequency plot - Gap Filled") +
  xlab("Range") + 
  ylab("(Log)Frequency") + 
  theme(plot.title = element_text(hjust = 0.5))

Cutoff_LOD <- min(blk_rem[blk_rem > 0], na.rm = TRUE) # lowest value except NAs
print(paste0("The limit of detection (LOD) is: ", Cutoff_LOD))

# Replacing NAs with random values
set.seed(141222) # by setting a seed, we generate the same set of random number all the time

imp <- blk_rem

for (i in 1:ncol(imp)) {
  imp[,i] <- ifelse(is.na(imp[,i]),
                    runif(nrow(imp), min = 0, max = Cutoff_LOD),
                    imp[,i])
}


# Save
write.csv(imp, IMP, row.names = TRUE)



# Replacing NAs with ZERO (0)
imp0 <- blk_rem

#data.frame(ifelse(is.na(blk_rem), 0, blk_rem))

for (i in 1:ncol(imp0)) {
  imp0[,i] <- ifelse(is.na(imp0[,i]),
                     0,
                     imp0[,i])
}


# Check if there's any zeros or NAs
sum(imp0 == 0)
sum(is.na(imp0))

# Save
write.csv(imp0, IMP0, row.names = TRUE)




# metadata - remove tech reps --------------------------------------------------
md_mean <- md_Samples %>%
  rownames_to_column(var = "sample_id") %>%
  mutate(sample_id = str_extract(sample_id, "C[0-9]{1,2}(_[^_]+){3}")) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  column_to_rownames(var = "sample_id")

write.csv(md_mean, paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_md_mean.csv"))

md_mean <- read.csv(paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_md_mean.csv"), row.names = 1)


# Formula ======================================================================
imp0_t <- as.data.frame(t(imp0)) 
ft_for <- read.csv(SEL_ORBI_DATA, row.names = 1)

# Merge the imputed data with formula
orbi_meta <- merge(imp0_t, ft_for, by.x = "row.names", by.y = "row.names", all.x = TRUE) %>%
  column_to_rownames(var="Row.names") %>%
  mutate(P = NA) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  select(-matches("MSM|OWC"), matches("MSM|OWC"))

orbi_meta


#================ FUNCTIONS adapted from Kaizad's codes ========================
compute_indices = function(dat){
  dat %>% 
    dplyr::mutate(AImod = round((1+C-(0.5*O)-S-(0.5*(N+P+H)))/(C-(0.5*O)-S-N-P),4),
                  NOSC =  round(4-(((4*C)+H-(3*N)-(2*O)-(2*S))/C),4),
                  DBE_AI = 1+C-O-S-0.5*(N+P+H),
                  DBEx =  1 + ((2*C-H + N + P))/2,
                  DBE_C = round(DBE_AI/C,4)) %>% 
    dplyr::select("C", "H", "N", "O", "S", "P", "AImod", "NOSC", "HC", "OC", "DBE_AI","DBE", "DBE_C")
}

assign_class_seidel = function(meta_clean){
  meta_clean %>%
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
                                      HC < 2.0 & HC >= 1.5 & N > 0 ~ "aliphatic+N"),
           Class_chem = case_when(OC > 0 & OC <= 0.3 & HC >= 1.5 & HC <= 2.5 ~ "lipid",
                                  OC >= 0 & OC <= 0.125 & HC >= 0.8 & HC < 1.5 ~ "unsatHC",
                                  OC > 0.3 & OC <= 0.55 & HC >= 1.5 & HC <= 2.3 ~ "protein",
                                  OC > 0.55 & OC <= 0.7 & HC >= 1.5 & HC <= 2.2 ~ "amino sugar",
                                  OC > 0.7 & OC <= 1.5 & HC >= 1.5 & HC <= 2.5 ~ "carbohydrate",
                                  OC > 0.125 & OC <= 0.65 & HC >= 0.8 & HC < 1.5 ~ "lignin",
                                  OC > 0.65 & OC <= 1.1 & HC >= 0.8 & HC < 1.5 ~ "tannin",
                                  OC >= 0 & OC <= 0.95 & HC >= 0.2 & HC < 0.8 ~ "condensed HC"),
           Class_chem = replace_na(Class_chem, "other"),) %>% 
    dplyr::select("C", "H", "N", "O", "S", "P", "AImod", "NOSC", "HC", "OC", "DBE_AI", "DBE", "DBE_C", Class, Class_detailed, Class_chem)
}
#===============================================================================

orbi_indices <- compute_indices(orbi_meta)

orbi_calculated <- assign_class_seidel(orbi_indices)

# Finally, clean out all NA, NaN, or Inf features
cleaned_orbi_meta <- orbi_calculated %>%
  filter(if_all(where(is.numeric), ~ is.finite(.) & !is.na(.)))

cleaned_orbi_meta


# Formula meta with features
orbi_all <- merge(orbi_meta %>% select(-DBE),
                  cleaned_orbi_meta %>%
                    select(-HC, -OC, -C, -H, -O, -N, -S, -P), all = FALSE, by = 0, all.y = TRUE) %>%
  column_to_rownames(var = "Row.names") %>%
  select(-matches("MSM|OWC"), matches("MSM|OWC"), -contains(".x")) # move sample columns to the end

write.csv(orbi_all, ORBI_FORMULAMETA_W_FEATURE, row.names = TRUE)


hcoc <- orbi_all[,c(1:18)]

write.csv(hcoc, HCOC, row.names = TRUE)




# Wide to long -----------------------------------------------------------------
for_to_merge <- hcoc

ft_long <- orbi_all %>%
  select(matches("MSM|OWC")) %>%
  rownames_to_column(., "feature") %>%
  reshape2::melt(id = "feature", value.name = "inten", variable.name = "sample_id") %>%
  filter(inten > 0) %>%
  inner_join(md_mean %>% rownames_to_column(var = "sample_id"), by = "sample_id") %>%
  merge(., for_to_merge %>% rownames_to_column(var = "feature"), by = "feature") %>%
  rename_with(~ sub(".y$", "", .x)) %>%
  rename(intensity = "inten") %>%
  mutate(sample_id = gsub("\\bC([12])\\b(?!\\d)", "0\\1", sample_id, perl=TRUE)) %>%
  dplyr::select(feature, formula, intensity, sample_name, sample_id, sample_type, site, transect_location, horizon,
         HC, OC, "C", "H", "N", "O", "S", "P", "AImod", "NOSC", "DBE_AI", "DBE", "DBE_C", Class, Class_detailed, Class_chem) %>%
  rename(orbitrap_id = "sample_id")

# Save
write.csv(ft_long, ORBI_LONG_MEAN)


# Look at TIC
rowSums(imp0)



# normalization ================================================================
orbi_all

dat_to_norm <- orbi_all %>%
  select(matches("MSM|OWC")) %>%
  t() %>%
  as.data.frame()

norm_tic <- normalization(dat_to_norm, 
                          method = "sum")$newXtrain

write.csv(norm_tic, paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_mean_norm_tic.csv"))
# ==============================================================================




# Visualization ################################################################
mode <- "POS"
ORBI_LONG_MEAN = paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_long_mean.csv")
ft_long <- read.csv(ORBI_LONG_MEAN, row.names = 1)
plot_dir <- paste0(wd, "/", mode, "/plots")

# Prelim plot-------------------------------------------------------------------
plot <- ggplot(ft_long %>%
                 arrange(sample_name) %>%
                 mutate(Class_chem = factor(Class_chem, levels = c("lipid", "amino sugar", "protein", "carbohydrate", "unsatHC", "lignin", "tannin", "condensed HC", "other"))),
               aes(x = OC, y = HC, fill = Class_chem, size = intensity), color = "black", alpha = 0.2) +
  geom_point(shape = 21, show.legend = c(alpha = FALSE)) +
  facet_wrap(~sample_name) +
  theme_bw() +
  scale_fill_brewer(palette = "RdYlBu") +
  labs(size = "Peak\nintensity", fill = "")

plot

ggsave("chem_class_plot3.png", path = plot_dir)
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

ggsave("tot_count.png", path = plot_dir)
#-------------------------------------------------------------------------------

# TIC---------------------------------------------------------------------------
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

TIC <- ggplot(ra_id, aes(x = sample_type, y = total, fill = sample_type)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(x = "", y = "TIC") +
  facet_wrap(~site, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))  

TIC

ggsave("TIC.png", path = plot_dir)
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

ggsave("ra_chem.png", path = plot_dir)
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
md_mean <- read.csv(paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_md_mean.csv"), row.names = 1)
norm_tic <- read.csv(paste0(wd, "/", mode, "/processed_data/", prj, "_", mode, "_mean_norm_tic.csv"), row.names = 1)

# Load libraries for PCoA
pacman::p_load("BiocManager","ComplexHeatmap","dendextend","NbClust","cowplot",
               "gtools", # mixedsort
               "tidyverse",
               "janitor", # row_to_names function
               "factoextra", # PCA
               "vegan") # PCoA

md_sub <- md_mean

# use imp0 to filter out all zero features

norm_mean2 <- norm_tic[which(rownames(norm_tic) %in% rownames(md_sub)), , drop = F]

filter_all_zero_features <- function(data){
  # Remove columns where all values are NA
  data <- data[, colSums(data) > 0]
  return(data)
}

norm_mean3 <- filter_all_zero_features(norm_mean2)

cleaned_data <- norm_mean3
md <- md_sub

# Verifying file consistency----------------------------------------------------
md <- md[mixedorder(rownames(md)), ,drop=F] # ordering the md by row names
cleaned_data <- cleaned_data[mixedorder(rownames(cleaned_data)),, drop=F] # ordering the md by row names

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

ggsave("PCoA.png", path = plot_dir)
# ------------------------------------------------------------------------------
