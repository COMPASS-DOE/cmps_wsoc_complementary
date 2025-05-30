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
tar_source("2-code/2-nmr.R")
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_target(sample_key, googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1SaA9qF4lHyOYwSkBhQZPDSV3TN8RZaCUaLHx-QlQGlk/")),
  
  # NMR
  ## spectra
  tar_target(nmr_spectra, nmr_import_spectra("1-data/nmr/processed/csv_spectra", method = "mnova")),
  tar_target(nmr_spectra_processed, process_nmr_spectra(nmr_spectra, sample_key)),
  tar_target(gg_nmr_spectra_all, plot_nmr_spectra_all(nmr_spectra_processed)),
  tar_target(gg_nmr_spectra_subset, plot_nmr_spectra_subset(nmr_spectra_processed)),
  
  ## peaks
  tar_target(nmr_peaks, nmr_import_peaks("1-data/nmr/processed/csv_peaks", method = "single column")),
  tar_target(nmr_peaks_processed, process_nmr_peaks(nmr_peaks)),
  tar_target(nmr_relabundance, nmr_relabund(nmr_peaks_processed, method = "peaks")),
  tar_target(gg_nmr_relabund, plot_nmr_relabund(nmr_relabundance, sample_key)),
  
  ## stats
  tar_target(nmr_permanova, compute_nmr_permanova(nmr_relabundance, sample_key)),
  tar_target(nmr_pca, compute_nmr_pca(nmr_relabundance, sample_key))

)


