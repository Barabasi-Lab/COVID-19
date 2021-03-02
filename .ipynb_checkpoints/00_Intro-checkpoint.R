#######################################
### Project: Network Medicine Framework for Identifying Drug Repurposing Opportunities for COVID-19.
### Description: Descriptive analysis
### Find Disease Module, Significance and Network Visualization
### Author: Deisy Gysi
### date: March 1st 2021
#######################################

rm(list = ls())
library(tidyr)
require(magrittr)
require(data.table)
require(igraph)
require(dplyr)
require(NetSci)

PPI = fread("DatasetS2.csv") # Read PPI
COVID = fread("DatasetS1.csv") # Read COVID Genes
PPI %<>% filter(!is.na(proteinA_entrezid) & !is.na(proteinB_entrezid))
gPPI = graph_from_data_frame(PPI, directed = F) %>% 
  simplify() # Create the PPI network and remove self-loops


Target_COVID = COVID %>% 
  pull(EntrezID) %>% 
  as.character() 

gCOVID = gPPI %>% 
  induced.subgraph(., vids = Target_COVID ) # Select the subgraph from COVID

# Plot it
V(gCOVID)$label = NA
gCOVID %<>% simplify()
V(gCOVID)$size = degree(gCOVID)
V(gCOVID)$size = degree(gCOVID) %>% CoDiNA::normalize()
V(gCOVID)$size  = (V(gCOVID)$size + 0.1)*5
gCOVID %<>% delete.vertices(., degree(.) == 0)
coord = layout_with_fr(gCOVID)
V(gCOVID)$label = ifelse(V(gCOVID)$size > 2.5, V(gCOVID)$name, NA)
E(gCOVID)$curved =0.1
V(gCOVID)$color = V(gCOVID)$frame.color = "red"
V(gCOVID)$label.color = "darkred"
V(gCOVID)$label.cex = V(gCOVID)$size/5
E(gCOVID)$color = 'salmon1'
plot(gCOVID)

# Calculate the LCC significance, degree preserving
LCC = LCC_Significance(N = 1000, 
                 Targets = Target_COVID, 
                 G = gPPI)

# Plot the distribution
Histogram_LCC(LCC)



