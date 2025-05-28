library(tidyverse)
library(nmrrr)
theme_set(theme_bw())

nmr_spectra = nmr_import_spectra("1-data/nmr/processed/csv_spectra", method = "mnova")
nmr_peaks = nmr_import_peaks("1-data/nmr/processed/csv_peaks", method = "single column")

nmr_plot_spectra(nmr_spectra, binset = bins_Lynch2019,
                 mapping = aes(group = sampleID,
                               x = ppm,
                               y = intensity))



nmr_spectra %>% 
  filter(intensity < 1) %>% 
  filter(sampleID %in% c("01", "02", "03", "04", "05")) %>% 
  nmr_plot_spectra(binset = bins_Lynch2019,
                   mapping = aes(group = sampleID,
                                 x = ppm,
                                 y = intensity),
                   label_position = 6.5,
                   stagger = 1.2)+
  xlim(10, 0.1)+
#  ylim(0, 1)+
  NULL


peaks_bins = nmr_assign_bins(nmr_peaks, binset = bins_Lynch2019)
relabund = nmr_relabund(dat = peaks_bins, method = "peaks")


relabund %>% 
  ggplot(aes(x = sampleID, y = relabund, fill = group))+
  geom_bar(stat = "identity")
