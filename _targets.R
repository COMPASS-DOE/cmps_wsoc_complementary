# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tidyverse", "googlesheets4", "vegan") # Packages that your targets need for their tasks.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source("2-code/0-packages.R")
tar_source("2-code/1a-fticr_functions.R")
#tar_source("2-code/1b-fticr_data.R")
tar_source("2-code/2-nmr.R")
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_target(sample_key, googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1SaA9qF4lHyOYwSkBhQZPDSV3TN8RZaCUaLHx-QlQGlk/")),
  
  # NMR -- SPE
  ## spectra
  tar_target(nmr_spe_spectra, nmr_import_spectra("1-data/nmr/processed_spe/csv_spectra", method = "mnova")),
  tar_target(nmr_spe_spectra_processed, process_nmr_spectra(nmr_spe_spectra, sample_key)),
  tar_target(gg_nmr_spe_spectra_all, plot_nmr_spectra_all(nmr_spe_spectra_processed)),
  tar_target(gg_nmr_spe_spectra_subset, plot_nmr_spectra_subset(nmr_spe_spectra_processed)),
  
  ## peaks
  tar_target(nmr_spe_peaks, nmr_import_peaks("1-data/nmr/processed_spe/csv_peaks", method = "single column")),
  tar_target(nmr_spe_peaks_processed, process_nmr_peaks(nmr_spe_peaks)),
  tar_target(nmr_spe_relabundance, nmr_relabund(nmr_spe_peaks_processed, method = "peaks")),
  tar_target(gg_nmr_spe_relabund, plot_nmr_relabund(nmr_spe_relabundance, sample_key)),
  
  ## stats
  tar_target(nmr_spe_permanova, compute_nmr_permanova(nmr_spe_relabundance, sample_key)),
  tar_target(nmr_spe_pca, compute_nmr_pca(nmr_spe_relabundance, sample_key)),
  
  
  # NMR -- WATER
  ## spectra
  tar_target(nmr_water_spectra, nmr_import_spectra("1-data/nmr/processed_water/csv_spectra", method = "mnova")),
  tar_target(nmr_water_spectra_processed, process_nmr_spectra(nmr_water_spectra, sample_key)),
  tar_target(gg_nmr_water_spectra_all, plot_nmr_spectra_all(nmr_water_spectra_processed)),
  tar_target(gg_nmr_water_spectra_subset, plot_nmr_spectra_subset(nmr_water_spectra_processed)),
  
  ## peaks
  tar_target(nmr_water_peaks, nmr_import_peaks("1-data/nmr/processed_water/csv_peaks", method = "single column")),
  tar_target(nmr_water_peaks_processed, process_nmr_peaks(nmr_water_peaks)),
  tar_target(nmr_water_relabundance, nmr_relabund(nmr_water_peaks_processed, method = "peaks")),
  tar_target(gg_nmr_water_relabund, plot_nmr_relabund(nmr_water_relabundance, sample_key)),
  
  ## stats
  tar_target(nmr_water_permanova, compute_nmr_permanova(nmr_water_relabundance, sample_key)),
  tar_target(nmr_water_pca, compute_nmr_pca(nmr_water_relabundance, sample_key)),
  
  # FTICR
  ## processing
  tar_target(icr_report_positive, read.csv("1-data/fticr/DOC 367-420 Positive Report.csv")),
  tar_target(icr_report_negative, read.csv("1-data/fticr/DOC 367-420 Negative Report.csv")),
  tar_target(icr_metadata, icr_process_metadata(icr_report_negative, icr_report_positive)),
  tar_target(icr_data_trt, icr_process_data_trt(icr_report_negative, icr_report_positive, sample_key)),
  tar_target(icr_data_cores, icr_process_data_cores(icr_report_negative, icr_report_positive, sample_key)),
  
  tar_target(icr_relabundance, compute_icr_relabund(icr_data_cores, icr_metadata)$icr_relabundance),
  tar_target(icr_relabundance_wide, compute_icr_relabund(icr_data_cores, icr_metadata)$icr_relabundance_wide)
  
#  tar_render(report, path = "3-reports/report.Rmd")

)


