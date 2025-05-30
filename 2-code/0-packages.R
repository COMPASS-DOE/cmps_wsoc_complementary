library(tidyverse)
library(nmrrr)
library(vegan)
library(googlesheets4)
library(ggbiplot)
library(PNWColors)
library(ggConvexHull)
#devtools::install_github("cmartin/ggConvexHull")

theme_kp <- function() {  # this for all the elements common across plots
  theme_bw() %+replace%
    theme(legend.position = "top",
          legend.key=element_blank(),
          legend.title = element_text(size = 14, hjust = 0),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1.5, 'lines'),
          legend.background = element_rect(colour = NA),
          panel.border = element_rect(color="black",linewidth=1, fill = NA),
          
          plot.title = element_text(hjust = 0, size = 14),
          axis.text = element_text(size = 14, color = "black"),
          axis.title = element_text(size = 16, face = "bold", color = "black"),
          
          # formatting for facets
          panel.background = element_blank(),
          strip.background = element_rect(colour="white", fill="white"), #facet formatting
          panel.spacing.x = unit(1, "lines"), #facet spacing for x axis
          panel.spacing.y = unit(1.5, "lines"), #facet spacing for x axis
          strip.text.x = element_text(size=14, face="bold"), #facet labels
          strip.text.y = element_text(size=14, face="bold", angle = 270) #facet labels
    )
}
theme_set(theme_kp())