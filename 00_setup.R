library(MSstatsPTM)
library(data.table)
library(gridExtra)
library(ggrepel)
library(MSstatsWeightedSummary) # GitHub package currently 
library(MSstatsTMT)
library(ggplot2)
library(igraph)
library(parallel)
library(splitstackshape)
library(MSstatsConvert)
library(SimulateTMT)
library(pbapply)

install.packages(c("data.table", "gridExtra", "ggrepel", "ggplot2",
                   "igraph", "parallel", "splitstackshape", "pbapply",
                   "BiocManager", "devtools", "tidyverse"))
BiocManager::install(c("MSstatsPTM", "MSstatsTMT", "MSstatsConvert"))
devtools::install_github("mstaniak/SimulateTMT")
devtools::install_github("Vitek-Lab/MSstatsWeightedSummary")