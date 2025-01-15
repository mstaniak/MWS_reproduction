# Code required to reproduce [main] Figure 3 (comparison of weighted and unique-based summaries for the PTM case study)
# Libraries ----
library(MSstatsPTM)
library(data.table)
library(gridExtra)
library(MSstatsWeightedSummary) # GitHub package currently 
library(MSstatsTMT)
library(ggplot2)
# Input data ----
load(file = "input_data/ptm/pST_MSstatsPTM_Format.rda")#shigella_peptide_pST
load(file = "input_data/ptm/global_MSstatsPTM_Format.rda")#shigella_global
load(file = "input_data/ptm/pY_MSstatsPTM_Format.rda")#shigella_peptide_pY
# Function to assign multi-site peptides to all matching modification sites
processSites = function(combined) {
  combined = as.data.table(combined)
  combined = unique(combined[, . (ProteinName, PeptideSequence, Charge, PSM, 
                                  Run, Channel, BioReplicate, Condition,
                                  Mixture, TechRepMixture,
                                  Intensity)])
  pp = unique(combined[, .(ProteinName, PeptideSequence, Run)])
  mods = as.data.table(stringr::str_extract_all(pp$ProteinName, "(_S[0-9]+)|(_T[0-9]+)|(_Y[0-9]+)", simplify = TRUE))
  
  pp = cbind(pp, mods)
  pp = melt(pp, measure.vars = c("V1", "V2"))
  pp = pp[value != ""]
  pp[order(ProteinName, Run)]
  pp[, Protein := stringr::str_replace_all(ProteinName, "(_S[0-9]+)|(_T[0-9]+)|(_Y[0-9]+)", "")]
  pp[, Protein := stringr::str_replace_all(Protein, "\\_", "")]
  pp[, ProteinNameNew := paste(Protein, value, sep = "_")]
  
  combined[, ProteinName := NULL]
  pp_merge = unique(pp[, .(ProteinName = ProteinNameNew, PeptideSequence)])
  combined = merge(combined, pp_merge, by = c("PeptideSequence"), allow.cartesian = TRUE)
  # assuming no normalization is required
  combined[, log2IntensityNormalized := log(Intensity, 2)]
  combined
}
# Pre-processing and feature / site counts ----
shigella_global = as.data.table(shigella_global)
shigella_global = shigella_global[, .(Intensity = max(Intensity)),
                                  by = .(ProteinName, PeptideSequence, 
                                         Charge, PSM, Mixture, TechRepMixture, Run, Channel, Condition, BioReplicate)]

shigella_peptide_pY = as.data.table(shigella_peptide_pY)
shigella_peptide_pST = as.data.table(shigella_peptide_pST)
shigella_peptide_pY_weight_input  = processSites(shigella_peptide_pY)
shigella_peptide_pST_weight_input = processSites(shigella_peptide_pST)

uniqueN(shigella_peptide_pY$ProteinName)
uniqueN(shigella_peptide_pY_weight_input$ProteinName)

uniqueN(shigella_peptide_pST$ProteinName)
uniqueN(shigella_peptide_pST_weight_input$ProteinName)

uniqueN(shigella_peptide_pST$PSM)
uniqueN(shigella_peptide_pST_weight_input$PSM)

graph_final = createPeptideProteinGraph(shigella_peptide_pST_weight_input)
shigella_peptide_pST_weight_input = addClusterMembership(shigella_peptide_pST_weight_input, graph_final)

as.data.table(shigella_peptide_pY)[, .(num_features_per_site = uniqueN(PSM)), by = "ProteinName"][, .(mean = mean(num_features_per_site))]
as.data.table(shigella_peptide_pY_weight_input)[, .(num_features_per_site = uniqueN(PSM)), by = "ProteinName"][, .(mean = mean(num_features_per_site))]

as.data.table(shigella_peptide_pST)[, .(num_features_per_site = uniqueN(PSM)), by = "ProteinName"][, .(mean = mean(num_features_per_site))]
as.data.table(shigella_peptide_pST_weight_input)[, .(num_features_per_site = uniqueN(PSM)), by = "ProteinName"][, .(mean = mean(num_features_per_site))]

shigella_peptide_pST_weight_input[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
shigella_peptide_pST_weight_input[, NumProteins := uniqueN(ProteinName), by = "Cluster"]

shigella_peptide_pST_weight_input[NumProteins > 1, uniqueN(Cluster)]
as.data.table(shigella_peptide_pST_weight_input)[!(IsUnique), .(num_features_per_site = uniqueN(PSM)), by = "ProteinName"][, .(mean = mean(num_features_per_site))]

unique(shigella_peptide_pST_weight_input[, .(Cluster, NumProteins)])[, .(mean = mean(NumProteins))]

as.data.table(shigella_global)[, .(num_features_per_site = uniqueN(PSM)), by = "ProteinName"][, .(mean = mean(num_features_per_site))]

NumFeatures= c(as.data.table(shigella_peptide_pST)[, .(uniqueN(PSM))],
               as.data.table(shigella_peptide_pST_weight_input)[, .(uniqueN(PSM))])

xtable::xtable(data.table(
  Processing = c("original", "proposed"),
  NumSites = c(uniqueN(shigella_peptide_pST$ProteinName), uniqueN(shigella_peptide_pST_weight_input$ProteinName)),
  NumFeaturesPerSite = c(as.data.table(shigella_peptide_pST)[, .(num_features_per_site = uniqueN(PSM)), by = "ProteinName"][, .(mean = mean(num_features_per_site))],
                         as.data.table(shigella_peptide_pST_weight_input)[, .(num_features_per_site = uniqueN(PSM)), by = "ProteinName"][, .(mean = mean(num_features_per_site))])
))

uniqueN(shigella_peptide_pST$ProteinName)
uniqueN(shigella_peptide_pST_weight_input$ProteinName)

as.data.table(shigella_peptide_pST)[, .(num_features_per_site = uniqueN(PSM)), by = "ProteinName"][, .(mean = mean(num_features_per_site))]
as.data.table(shigella_peptide_pST_weight_input)[, .(num_features_per_site = uniqueN(PSM)), by = "ProteinName"][, .(mean = mean(num_features_per_site))]

input_with_cls = copy(shigella_peptide_pST_weight_input)
input_with_cls[, NumSites := uniqueN(ProteinName), by = "Cluster"]
input_with_cls[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]

uniqueN(shigella_peptide_pST$PSM)
uniqueN(shigella_peptide_pST_weight_input$PSM)

orig_pp = unique(shigella_peptide_pST[, .(PSM, ProteinName)])
procd_pp= unique(shigella_peptide_pST_weight_input[, .(PSM, ProteinName)])

comp_pp = merge(orig_pp, procd_pp, by = c("PSM"), all.x = T, all.y = T, allow.cartesian = T)
comp_pp[, NumSitesPerPSM := uniqueN(ProteinName.y), by ="PSM"]

shigella_peptide_pST = as.data.table(shigella_peptide_pST)
# 1.98 / 1.67 = 1.185 (18% increase)

pY_input <- list(PTM = data.frame(shigella_peptide_pY), 
                 PROTEIN = data.frame(shigella_global))
pST_input <- list(PTM = data.frame(shigella_peptide_pST), 
                  PROTEIN = data.frame(shigella_global))

## Data summarization ---- (original processing / code from the MSstatsPTM paper)
# py_summarization <- dataSummarizationPTM_TMT(pY_input)
# pST_summarization <- dataSummarizationPTM_TMT(pST_input)
# saveRDS(py_summarization, "processed_data/ptm/pY_summ.RDS")
# saveRDS(pST_summarization, "processed_data/ptm/pST_summ.RDS")

py_summarization = readRDS("processed_data/ptm/pY_summ.RDS")
pST_summarization = readRDS("processed_data/ptm/pST_summ.RDS")

py_summarized_data <- py_summarization$PTM$ProteinLevelData
py_feature_data <- py_summarization$PTM$FeatureLevelData

pst_summarized_data <- pST_summarization$PTM$ProteinLevelData
pst_feature_data <- pST_summarization$PTM$FeatureLevelData

global_summarized_data <- pST_summarization$PROTEIN$ProteinLevelData
global_feature_data <- pST_summarization$PROTEIN$FeatureLevelData

## check the overlapped sites between pST and pY
shared_sites <- intersect(unique(py_summarized_data$Protein), 
                          unique(pst_summarized_data$Protein))
shared_ST <- shared_sites[grepl("_S", shared_sites)|grepl("_T", shared_sites)]
shared_Y <- shared_sites[grepl("_Y", shared_sites)]

## remove the overlapped Y sites from data_ST
data_ST <- as.data.table(pst_summarized_data)[!(Protein %in% shared_Y)]
## remove the overlapped ST sites from data_Y
data_Y <- as.data.table(py_summarized_data)[!(Protein %in% shared_ST)]

## combine pST and pY
pSTY_sum_msstats <- rbind(data_ST, data_Y)

## remove the overlapped Y sites from data_ST
data_ST_feature <- as.data.table(pst_feature_data)[!(ProteinName %in% shared_Y)]
## remove the overlapped ST sites from data_Y
data_Y_feature <- as.data.table(py_feature_data)[!(ProteinName %in% shared_ST)]

## combine pST and pY
pSTY_features <- rbind(data_ST_feature, data_Y_feature)


# global_sum_msstats_norm <- global_summarized_data %>% filter(Condition != "Norm")
# pSTY_sum_msstats_norm <- pSTY_sum_msstats %>% filter(Condition != "Norm")
# 
# pyst_sum_input <- list(PTM = list(ProteinLevelData = pSTY_sum_msstats_norm,
#                                   FeatureLevelData = pSTY_features),
#                        PROTEIN = list(ProteinLevelData = global_sum_msstats_norm,
#                                       FeatureLevelData = global_feature_data))

pyst_sum_input = readRDS("processed_data/ptm/gc_input.RDS")
ptm_quants = as.data.table(pyst_sum_input$PTM$ProteinLevelData)

# comparison <- read.delim("comparison_matrix.txt")
# rownames(comparison) <- comparison$X
# comparison <- comparison %>% dplyr::select(-X)
# comparison <- as.matrix(comparison) # this step is necessary for contrast comparison
# pSTY_model <- groupComparisonPTM(pyst_sum_input, data.type = "TMT", 
#                                  contrast.matrix = comparison)
# 
# ptm_gc = as.data.table(pSTY_model$ADJUSTED.Model)
# saveRDS(ptm_gc, "results/ptm/ptm_gc.RDS")
ptm_gc = readRDS("results/ptm/ptm_gc.RDS")

# Comparison of the summaries ----
input_with_cls[, ProteinName := stringr::str_replace_all(ProteinName, "__", "_")]
cl = input_with_cls[grepl("9Q6J5", ProteinName)][grepl("S236", ProteinName) | grepl("S240", ProteinName)]
cl[, Run := as.character(Run)]
cl[, ProteinName := as.character(ProteinName)]

summ_sh = MSstatsWeightedSummary::getWeightedProteinSummary(input_with_cls[Cluster == 5276], "Huber", norm_parameter = 1e-3, max_iter = 100)
featureWeights(summ_sh)

concat_summ = ptm_quants[grepl("9Q6J5", ProteinName)][grepl("S236", ProteinName) & grepl("S240", ProteinName)]
concat_summ_to_plot = rbind(cbind(concat_summ[, !(colnames(concat_summ) %in% c("Protein", "ProteinName")), with = F], Protein = "BD1L1_MOUSE|E9Q6J5_S240"),
                            cbind(concat_summ[, !(colnames(concat_summ) %in% c("Protein", "ProteinName")), with = F], Protein = "BD1L1_MOUSE|E9Q6J5_S236"))
concat_summ_to_plot[, Protein := stringr::str_replace_all(Protein, "_", "")]

comp = rbind(cbind(proteinData(summ_sh), Method = "proposed"),
             cbind(concat_summ_to_plot, Method = "concatenation"), use.names = T, fill = T)

concat_summ2_plot = ptm_quants[grepl("9Q6J5", ProteinName)][grepl("S236", ProteinName) | grepl("S240", ProteinName)]
concat_summ2_plot[, ProteinName := NULL]
concat_summ2_plot[, Method := "concatenation"]

comp2 = rbind(cbind(proteinData(summ_sh), Method = "proposed"),
              concat_summ2_plot, use.names = T, fill = T)
comp2[, ProteinLabel := stringr::str_replace_all(Protein, "BD1L1_?MOUSE\\|", "")]
comp2[, Site := stringr::str_replace_all(ProteinLabel, "E9Q6J5_", "")]
comp2[, Site := stringr::str_replace_all(Site, "_", ", ")]
comp2[, Site := paste("Site: ", Site)]
comp2[, RunLabel := paste("Mixture: ", ifelse(Run == "1_1", 1, 2))]

feat_data = copy(cl)
feat_data[, Protein := ProteinName]
feat_data[, ProteinLabel := stringr::str_replace_all(Protein, "BD1L1_?MOUSE\\|", "")]
feat_data[, ProteinLabel2 := ifelse(IsUnique, ProteinLabel, "E9Q6J5_S236_S240")]
feat_data[, Site2 := stringr::str_replace_all(ProteinLabel2, "E9Q6J5_", "")]
feat_data[, Site2 := stringr::str_replace_all(Site2, "_", ", ")]
feat_data[, Site2 := paste("Site: ", Site2)]
feat_data[, RunLabel := paste("Mixture: ", ifelse(Run == "1_1", 1, 2))]

feat_shared = feat_data[!(IsUnique)]
feat_shared[, Site2 := NULL]
feat_shared[, Site := NULL]

feat_shared2 = rbind(cbind(feat_shared, Site = "Site: S236"),
                     cbind(feat_shared, Site = "Site: S236, S240"),
                     cbind(feat_shared, Site = "Site: S240"))
feat_data2 = rbind(feat_data, feat_shared2, fill = T)

feat_data[, feature := ifelse(IsUnique, "unique", "shared")]
feat_data[, feature := factor(feature, levels = c("unique", "shared"))]

feature_data = feat_data2
color_choice = c("full" = "darkblue",
                 "unique" = "lightblue",
                 "concatenation" = "purple",
                 "proposed" = "red",
                 "PeCorA" = "yellow",
                 "Quantifere" = "orange",
                 "VIQoR" = "green")
full_names = c("full" = "gold standard",
               "all" = "all features",
               "shared" = "proposed",
               "unique" = "selected unique",
               "PeCorA" = "PeCorA",
               "Quantifere" = "Quantifere",
               "VIQoR" = "VIQoR")
summaries_df = comp2
colors = color_choice[unique(c(summaries_df$Method))]

chs = unique(summaries_df[, .(Run, Channel, Condition, BioReplicate)])[order(Run, BioReplicate), Channel]
chs = unique(chs)
summaries_df$ProteinName = summaries_df$Protein

summaries_df[, Channel := factor(Channel, levels = chs, ordered = T)]

feature_data[, feature := ifelse(IsUnique, "unique", "shared")]
feature_data[, feature := factor(feature, levels = c("unique", "shared"), ordered = TRUE)]
feature_data[, Channel := factor(Channel, levels = chs, ordered = T)]
feature_data[, Site := stringr::str_replace_all(Site, "  ", " ")]
comp2[, Site := stringr::str_replace_all(Site, "  ", " ")]
meths = unique(summaries_df$Method)
names(meths) = meths

feature_data[, Site := ifelse(is.na(Site), Site2, Site)]
comp2[, Site := ifelse(is.na(Site), Site2, Site)]

feature_data[, Site := stringr::str_replace(Site, "Site", "site")]
comp2[, Site := stringr::str_replace(Site, "Site", "site")]
feature_data[, Site := stringr::str_replace(Site, "  ", " ")]
comp2[, Site := stringr::str_replace(Site, "  ", " ")]

ggplot(comp2,
       aes(x = reorder(Channel, BioReplicate), y = Abundance, group = Method, color = Method)) +
  geom_line(aes(x = reorder(Channel, BioReplicate), y = log2IntensityNormalized, group = PSM, linetype = feature), 
            data = feature_data, inherit.aes = FALSE) +
  geom_point(size = 1.1) +
  geom_line(size = 1.1) +
  facet_grid(RunLabel ~ ifelse(grepl(",", Site), "site: S236_S240", Site)) +
  theme_bw() +
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 18), 
        legend.position = "bottom",
        legend.box = "vertical") +
  scale_color_manual(name = "protein-level summary based on", values = sort(colors), labels = (full_names[meths[names(sort(colors))]])) +
  scale_linetype_discrete(name = "quantitative profile of a peptide across channels") +
  ylab("log-intensity") +
  xlab("channel")
ggsave("./plots/pdf/multisite_profiles.pdf", device = "pdf", width = 10, height = 5, units = "in", scale = 1, dpi = 300)
ggsave("./plots/png/multisite_profiles.png", device = "png", width = 10, height = 5, units = "in", scale = 1, dpi = 300)
