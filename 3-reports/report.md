
## FTICR-MS

### VAN KREVELEN

Unique peak assignments for + vs. - mode
![](report_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

    ## [1] "unique + vs - mode"

![](report_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

All identified peaks across the different samples

![](report_files/figure-gfm/icr_vk-1.png)<!-- -->

### PERMANOVA

|                   |  Df |  SumOfSqs |        R2 |          F | Pr(\>F) |
|:------------------|----:|----------:|----------:|-----------:|--------:|
| mode              |   1 | 2.6604864 | 0.8794898 | 2259.36464 |   0.001 |
| region            |   1 | 0.0896189 | 0.0296258 |   76.10708 |   0.001 |
| transect_location |   2 | 0.0739490 | 0.0244457 |   31.39987 |   0.001 |
| horizon           |   2 | 0.0298575 | 0.0098701 |   12.67794 |   0.001 |
| Residual          |  89 | 0.1048008 | 0.0346445 |         NA |      NA |
| Total             |  95 | 3.0250340 | 1.0000000 |         NA |      NA |

PERMANOVA across all FTICR data showed that mode (positive vs. negative)
accounted for 88% of overall variability/separation.

So we split the data by mode for all subsequent analyses.

|                   |  Df |  SumOfSqs |        R2 |         F | Pr(\>F) |
|:------------------|----:|----------:|----------:|----------:|--------:|
| region            |   1 | 0.0628111 | 0.3144840 | 144.51811 |   0.001 |
| transect_location |   2 | 0.0501636 | 0.2511601 |  57.70911 |   0.001 |
| horizon           |   2 | 0.0133861 | 0.0670219 |  15.39964 |   0.001 |
| Residual          |  42 | 0.0182542 | 0.0913957 |        NA |      NA |
| Total             |  47 | 0.1997275 | 1.0000000 |        NA |      NA |

For negative mode, region accounted for 31 % and transect 25 % of total
variability.

|                   |  Df |  SumOfSqs |        R2 |        F | Pr(\>F) |
|:------------------|----:|----------:|----------:|---------:|--------:|
| region            |   1 | 0.0298027 | 0.1808196 | 20.70815 |   0.001 |
| transect_location |   2 | 0.0381819 | 0.2316581 | 13.26518 |   0.001 |
| horizon           |   2 | 0.0196958 | 0.1194985 |  6.84271 |   0.007 |
| Residual          |  42 | 0.0604455 | 0.3667361 |       NA |      NA |
| Total             |  47 | 0.1648201 | 1.0000000 |       NA |      NA |

For positive mode, region accounted for 18 % and transect 23 % of total
variability.

### PCA

![](report_files/figure-gfm/icr_pca_biplot-1.png)<!-- -->![](report_files/figure-gfm/icr_pca_biplot-2.png)<!-- -->

### NMDS

clustering and NMDS for NEGATIVE MODE

![](report_files/figure-gfm/icr_NMDS_plot-1.png)<!-- -->![](report_files/figure-gfm/icr_NMDS_plot-2.png)<!-- -->

### RELATIVE ABUNDANCE

![](report_files/figure-gfm/icr_relabund-1.png)<!-- -->

------------------------------------------------------------------------

## NMR

We did NMR on (a) water extracts and (b) SPExtracts

### Example NMR spectra

One example spectrum per transect location

    ## [1] "WATER EXTRACTS"

![](report_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

    ## [1] "SPE EXTRACTS"

![](report_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

### PERMANOVA

|                   |  Df |  SumOfSqs |        R2 |          F | Pr(\>F) |
|:------------------|----:|----------:|----------:|-----------:|--------:|
| extraction        |   1 | 1.9789114 | 0.6932902 | 172.826234 |   0.001 |
| site              |   0 | 0.0000000 | 0.0000000 |       -Inf |      NA |
| transect_location |   2 | 0.0686294 | 0.0240436 |   2.996841 |   0.032 |
| horizon           |   0 | 0.0000000 | 0.0000000 |       -Inf |      NA |
| Residual          |  45 | 0.5152633 | 0.1805169 |         NA |      NA |
| Total             |  49 | 2.8543766 | 1.0000000 |         NA |      NA |

When analyzing WATER and SPE extracts together, the extraction type
contributed to 69% of total variability.

    ## [1] "WATER EXTRACTS"

| term              |  df |  SumOfSqs |        R2 | statistic | p.value |
|:------------------|----:|----------:|----------:|----------:|--------:|
| site              |   0 | 0.0000000 | 0.0000000 |       Inf |      NA |
| transect_location |   2 | 0.0702478 | 0.1092975 |  2.307003 |   0.086 |
| horizon           |   0 | 0.0000000 | 0.0000000 |      -Inf |      NA |
| Residual          |  17 | 0.2588235 | 0.4026995 |        NA |      NA |
| Total             |  20 | 0.6427211 | 1.0000000 |        NA |      NA |

    ## [1] "SPE EXTRACTS"

| term              |  df |  SumOfSqs |        R2 | statistic | p.value |
|:------------------|----:|----------:|----------:|----------:|--------:|
| site              |   0 | 0.0000000 | 0.0000000 |      -Inf |      NA |
| transect_location |   2 | 0.0801772 | 0.4140556 |  11.14333 |   0.001 |
| horizon           |   0 | 0.0000000 | 0.0000000 |      -Inf |      NA |
| Residual          |  25 | 0.0899386 | 0.4644656 |        NA |      NA |
| Total             |  28 | 0.1936388 | 1.0000000 |        NA |      NA |

For WATER extracts, transect location accounted for only 10% of total
variability.

For SPE extracts, transect location accounted for 41% of total
variability.

### PCA

Plotting WATER and SPE extracts together

![](report_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Separation was primarily along SPE vs. WATER extraction types.

    ## [1] "WATER EXTRACTS"

![](report_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

    ## [1] "SPE EXTRACTS"

![](report_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

### relative abundance

    ## [1] "WATER EXTRACTS"

![](report_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

    ## [1] "SPE EXTRACTS"

![](report_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

- water extracts more variable. more aliphatic
- SPE extracts less variable. more unsaturated

------------------------------------------------------------------------

## Session Info

<details>
<summary>
Session Info
</summary>

Date run: 2026-03-23

    ## R version 4.5.0 (2025-04-11)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS 26.3.2
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggh4x_0.3.1         factoextra_1.0.7    cluster_2.1.8.1    
    ##  [4] whistledown_0.1.0   ggConvexHull_0.1.0  PNWColors_0.1.0    
    ##  [7] ggbiplot_0.55       googlesheets4_1.1.1 vegan_2.7-1        
    ## [10] permute_0.9-7       nmrrr_1.0.0         lubridate_1.9.4    
    ## [13] forcats_1.0.0       stringr_1.5.1       dplyr_1.2.0        
    ## [16] purrr_1.0.4         readr_2.1.5         tidyr_1.3.2        
    ## [19] tibble_3.3.1        ggplot2_4.0.2       tidyverse_2.0.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gld_2.6.7          readxl_1.4.5       rlang_1.1.7        magrittr_2.0.3    
    ##  [5] snakecase_0.11.1   e1071_1.7-16       compiler_4.5.0     mgcv_1.9-1        
    ##  [9] callr_3.7.6        vctrs_0.7.1        pkgconfig_2.0.3    fastmap_1.2.0     
    ## [13] backports_1.5.0    labeling_0.4.3     rmarkdown_2.29     tzdb_0.5.0        
    ## [17] haven_2.5.4        ps_1.9.1           xfun_0.53          broom_1.0.12      
    ## [21] parallel_4.5.0     prettyunits_1.2.0  DescTools_0.99.60  R6_2.6.1          
    ## [25] stringi_1.8.7      RColorBrewer_1.1-3 car_3.1-3          boot_1.3-31       
    ## [29] cellranger_1.1.0   Rcpp_1.1.1         knitr_1.50         Matrix_1.7-3      
    ## [33] splines_4.5.0      igraph_2.1.4       timechange_0.3.0   tidyselect_1.2.1  
    ## [37] rstudioapi_0.17.1  abind_1.4-8        yaml_2.3.10        targets_1.11.3    
    ## [41] codetools_0.2-20   processx_3.8.6     lattice_0.22-6     plyr_1.8.9        
    ## [45] withr_3.0.2        S7_0.2.0           evaluate_1.0.3     proxy_0.4-27      
    ## [49] pillar_1.10.2      ggpubr_0.6.0       carData_3.0-5      generics_0.1.3    
    ## [53] hms_1.1.3          scales_1.4.0       rootSolve_1.8.2.4  base64url_1.4     
    ## [57] class_7.3-23       glue_1.8.0         janitor_2.2.1      lmom_3.2          
    ## [61] tools_4.5.0        data.table_1.17.0  ggsignif_0.6.4     Exact_3.3         
    ## [65] fs_1.6.6           mvtnorm_1.3-3      cowplot_1.1.3      grid_4.5.0        
    ## [69] nlme_3.1-168       googledrive_2.1.1  Formula_1.2-5      cli_3.6.5         
    ## [73] expm_1.0-0         gargle_1.5.2       gtable_0.3.6       rstatix_0.7.2     
    ## [77] digest_0.6.37      ggrepel_0.9.6      farver_2.1.2       htmltools_0.5.8.1 
    ## [81] lifecycle_1.0.5    httr_1.4.7         secretbase_1.0.5   MASS_7.3-65

</details>
