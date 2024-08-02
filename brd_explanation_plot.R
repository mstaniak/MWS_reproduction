library(data.table)
library(MSstatsWeightedSummary)
library(ggplot2)
library(lme4)
library(gridExtra)
library(igraph)
library(parallel)
library(splitstackshape)
library(MSstatsConvert)

g00_orig = fread("./input_data/protein_degrader/G0011_OBJ0037605_msstats_input_with_protein_map.txt")
# Remove decoys
g00_orig = g00_orig[!grepl("##", ProteinName, fixed = T)]

# Extract all matching proteins for each peptide
orig_long_by_prot = cSplit(g00_orig, splitCols = "Protein", sep = ";", drop = F, direction = "long")
orig_long_by_prot = orig_long_by_prot[!grepl("##", Protein, fixed = TRUE)]
orig_long_by_prot[, ProteinName := NULL]
setnames(orig_long_by_prot, "Protein", "ProteinName")
orig_long_by_prot[, NumProteinsPerPeptide := uniqueN(ProteinName), by = "PSM"]
orig_long_by_prot[NumProteinsPerPeptide > 1]
orig_long_by_prot[grepl("BRD", ProteinName)]
peptide_protein_graph = createPeptideProteinGraph(orig_long_by_prot, "ProteinName", "PeptideSequenceNew")
quant_with_cls = addClusterMembership(orig_long_by_prot, peptide_protein_graph)
quant_with_cls[, NumProteins := uniqueN(ProteinName), by = "Cluster"]
quant_with_cls[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
quant_with_cls[, AllUnique := all(IsUnique), by = "Cluster"]

quant_with_cls[, log2Intensity := logintensity]
quant_with_cls[, log2IntensityNormalized := logintensity]

quant_no_single_shared = processIsoforms(quant_with_cls, remove_single_shared = TRUE, merge_identical = FALSE, remove_subsets = FALSE)
quant_merged_identical = processIsoforms(quant_no_single_shared, remove_single_shared = F, merge_identical = T, remove_subsets = FALSE)
quant_no_subsets = processIsoforms(quant_merged_identical, remove_single_shared = F, merge_identical = F, remove_subsets = T)

quant_with_cls[, Intensity := 2 ^ logintensity]
quant_unique_procd = MSstatsConvert::MSstatsPreprocess(quant_with_cls[, .(ProteinName, PeptideSequence, PrecursorCharge = Charge, PSM, Run, Channel, Intensity)], 
                                                       annotation = unique(quant_with_cls[, .(Run, Channel, BioReplicate, Condition, Mixture, TechRepMixture)]),
                                                       feature_columns = c("PeptideSequence", "PrecursorCharge"), 
                                                       remove_shared_peptides = TRUE, remove_single_feature_proteins = F, 
                                                       list(remove_features_with_few_measurements = FALSE,
                                                            summarize_multiple_psms = max))

quant_no_subsets[, Intensity := 2^log2Intensity]
quant_all_procd = MSstatsConvert::MSstatsPreprocess(quant_no_subsets[, .(ProteinName, PeptideSequence, PrecursorCharge = Charge, PSM, Run, Channel, Intensity)], 
                                                    annotation = unique(quant_no_subsets[, .(Run, Channel, BioReplicate, Condition, Mixture, TechRepMixture)]),
                                                    feature_columns = c("PeptideSequence", "PrecursorCharge"), 
                                                    remove_shared_peptides = FALSE, remove_single_feature_proteins = F, 
                                                    list(remove_features_with_few_measurements = FALSE,
                                                         summarize_multiple_psms = max))

quant_all_procd = as.data.table(quant_all_procd)
final_graph = createPeptideProteinGraph(quant_all_procd)
quant_all_procd = addClusterMembership(quant_all_procd, final_graph)

quant_all_procd[, NumProteins := uniqueN(ProteinName), by = "Cluster"]
quant_all_procd[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]

cm_g00 = fread("./input_data/protein_degrader/OBJ037605_G0011_comparison_matrix.txt")
group_labels = cm_g00$V1
cm_g00[, V1 := NULL]
cm_g00 = as.matrix(cm_g00)
row.names(cm_g00) = group_labels

quant_all_procd[grepl("BRD4", ProteinName), unique(Cluster)]

brd_cluster = quant_all_procd[Cluster == 529]
brd_cluster[, log2IntensityNormalized := log(Intensity, 2)]
brd_cluster[, ProteinName := stringr::str_extract(ProteinName, "BRD[2-4]")]
setnames(brd_cluster, "PrecursorCharge", "Charge")

plot_input_brd = copy(brd_cluster)
ch_order = c("126C", "127N", "127C", "128N", "128C",
             "129N", "129C", "130N", "130C", "131N")
plot_input_brd[, Channel := factor(Channel, levels = ch_order, ordered = TRUE)]


ggplot(plot_input_brd,
       aes(x = Channel, y = logintensity,
           group = PSM, color = IsUnique)) +
  geom_line() +
  facet_grid(~ProteinName) +
  theme_bw() +
  theme(legend.position = "bottom")

example_input_1 = plot_input_brd[ProteinName != "BRD3"][!(IsUnique) | PSM %in% c("K.AQAEHAEK.E_3",
                                                                                 "K.LNLPDYYK.I_2",
                                                                                 "R.KLQDVFEFR.Y_3",
                                                                                 "R.LAELQEQLR.A_2")]

full_unique_input = copy(brd_cluster)
summary_input_1 = copy(example_input_1)
summary_input_1[, Channel := as.factor(as.character(Channel))]

full_summary = getWeightedProteinSummary(full_unique_input,
                                         norm = "Huber",
                                         norm_parameter = 1e-3)
summary_1_unique = getWeightedProteinSummary(summary_input_1[(IsUnique)],
                                             norm = "Huber",
                                             norm_parameter = 1e-3)
summary_1_shared = getWeightedProteinSummary(summary_input_1,
                                             norm = "Huber",
                                             norm_parameter = 0.001)
summary_1_shared = getWeightedProteinSummary(full_unique_input,
                                             norm = "Huber",
                                             norm_parameter = 0.001,
                                             tolerance = 0.01, max_iter = 30)

weights_1 = merge(unique(summary_input_1[, .(ProteinName, PSM, IsUnique)]),
                  featureWeights(summary_1_shared),
                  by = c("ProteinName", "PSM", "IsUnique"),
                  all.x = TRUE)
weights_1[, Weight := ifelse(is.na(Weight), 0, Weight)]
weights_1[, ProteinName := stringr::str_extract(ProteinName, "BRD[2-4]")]
weights_1[!(IsUnique)][order(PSM)]

summaries_df = proteinData(summary_1_unique)
summaries_df$Channel = factor(summaries_df$Channel, levels = ch_order, 
                              ordered = TRUE)
summaries_df[, Abundance := Abundance - mean(Abundance), by = "Protein"]
summaries_df[, ProteinName := Protein]
feature_data = copy(summary_input_1[!(IsUnique)])
feature_data[, log2IntensityNormalized := log2IntensityNormalized - mean(log2IntensityNormalized),
             by = c("ProteinName", "PSM")]
feature_data$Channel = factor(feature_data$Channel, levels = ch_order, 
                              ordered = TRUE)

plot = ggplot(summaries_df, aes(x = Channel, y = Abundance, 
                                group = paste(Method, ProteinName), color = Method))
plot = plot + geom_line(aes(x = Channel, y = log2IntensityNormalized, 
                            group = PSM), data = feature_data, 
                        inherit.aes = FALSE, color = "grey", alpha = 0.5, 
                        size = 0.8)

plot = plot + geom_point(size = 1.2) + 
  geom_line(size = 1.2) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270), 
        legend.position = "bottom")
plot



feature_data2 = copy(summary_input_1)
feature_data2[, log2IntensityNormalized := log2IntensityNormalized - mean(log2IntensityNormalized),
              by = c("ProteinName", "PSM")]
feature_data2$Channel = factor(feature_data2$Channel, levels = ch_order, 
                               ordered = TRUE)
feature_data2[, ProteinName := ifelse(IsUnique, ProteinName, "shared")]
feature_data2[, feature := ifelse(IsUnique, "unique", "shared")]
feature_data2[, type := factor(ifelse(IsUnique, "unique", "shared"),
                               levels = c("unique", "shared")[c(2, 1)], ordered = TRUE)]
feature_data2[, FeatureID := as.numeric(as.factor(PSM))]
unique(feature_data2[, .(PSM, FeatureID)])
plot3 = ggplot(summaries_df[!grepl("DMSO", Condition)], aes(x = Channel, y = Abundance, 
                                                            group = paste(ProteinName), color = ProteinName))
plot3 = plot3 +
  geom_line(aes(x = Channel, y = log2IntensityNormalized, 
                group = PSM, color = ProteinName,
                linetype = type), data = feature_data2[!grepl("DMSO", Condition)], 
            inherit.aes = FALSE, alpha = 0.5, 
            size = 1.1) +
  geom_text(aes(x = Channel, y = log2IntensityNormalized, label = FeatureID),
            data = unique(feature_data2[Channel == "131N", .(PSM, log2IntensityNormalized, FeatureID, Channel)]),
            inherit.aes = FALSE)

plot3 = plot3 +
  # geom_point(size = 1.2) + 
  geom_line(size = 1.3) + 
  # facet_grid(~ProteinName) +
  scale_color_discrete(name = "feature membership") +
  scale_linetype_discrete(name = "feature") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270), 
        legend.position = "bottom",
        legend.box = "vertical")
plot3

featureWeights(summary_1_shared)[PSM %in% c("K.RQLSLDINKLPGEK.L_3", "R.LMFSNCYK.Y_2")]
featureWeights(summary_1_shared)[PSM %in% c("K.RQLSLDINKLPGEK.L_3", "R.LMFSNCYK.Y_2")]

feat_annot = unique(feature_data2[Channel == "131N", .(PSM, log2IntensityNormalized, FeatureID, Channel)][PSM %in% c("K.RQLSLDINKLPGEK.L_3", "R.LMFSNCYK.Y_2")])
feat_annot$Label2 = c("Peptide 1", "Peptide 2")

feat_annot2 = merge(feat_annot, featureWeights(summary_1_shared)[PSM %in% c("K.RQLSLDINKLPGEK.L_3", "R.LMFSNCYK.Y_2")],
                    by = c("PSM"))[ProteinName == "BRD4"]

ggplot(summaries_df[!grepl("DMSO", Condition)], 
       aes(x = Channel, y = Abundance, 
           group = paste(ProteinName), color = ProteinName)) +
  geom_line(size = 1.5) + # alpha = 0.8
  geom_line(aes(x = Channel, y = log2IntensityNormalized,
                group = PSM),
            data = feature_data2[!grepl("DMSO", Condition)][feature != "unique"][PSM %in% c("K.RQLSLDINKLPGEK.L_3", "R.LMFSNCYK.Y_2")],
            inherit.aes = FALSE, color = "#2297E6",
            size = 2, linetype = 2) +
  geom_text(aes(x = Channel, y = ifelse(Label2 == "Peptide 2", log2IntensityNormalized + 0.1, log2IntensityNormalized), label = Label2),
            data = feat_annot, size = 8,
            inherit.aes = FALSE) +
  geom_text(aes(x = Channel, y = ifelse(Label2 == "Peptide 2", log2IntensityNormalized + 0.01, log2IntensityNormalized - 0.09), label = paste("BRD4 Weight:", round(Weight, 2))),
            data = feat_annot2, size = 6,
            inherit.aes = FALSE) +
  scale_color_manual(name = "protein", values = c("#E69F00", "#56B4E9")) +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.key.spacing = unit(10, "pt"),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        legend.position = "bottom",
        legend.box = "vertical")

