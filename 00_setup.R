# Install required packages
install.packages(c("data.table", "gridExtra", "ggrepel", "ggplot2",
                   "igraph", "parallel", "splitstackshape", "pbapply",
                   "BiocManager", "devtools", "tidyverse"))
BiocManager::install(c("MSstatsPTM", "MSstatsTMT", "MSstatsConvert"))
devtools::install_github("mstaniak/SimulateTMT")
devtools::install_github("Vitek-Lab/MSstatsWeightedSummary")