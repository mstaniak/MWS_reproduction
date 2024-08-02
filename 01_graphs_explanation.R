# Libraries ----
library(data.table)
library(MSstatsWeightedSummary)
library(ggplot2)
library(lme4)
library(gridExtra)
library(igraph)
library(parallel)
library(splitstackshape)
library(MSstatsConvert)

# Data import and basic cleaning ----
g00_orig = fread("./input_data/protein_degrader/G0011_OBJ0037605_msstats_input_with_protein_map.txt")
# Remove decoys
g00_orig = g00_orig[!grepl("##", ProteinName, fixed = T)]

# Extract all matching proteins for each peptide
orig_long_by_prot = cSplit(g00_orig, splitCols = "Protein", sep = ";", drop = F, direction = "long")
orig_long_by_prot = orig_long_by_prot[!grepl("##", Protein, fixed = TRUE)]
orig_long_by_prot[, ProteinName := NULL]
setnames(orig_long_by_prot, "Protein", "ProteinName")
orig_long_by_prot[, NumProteinsPerPeptide := uniqueN(ProteinName), by = "PSM"]
peptide_protein_graph = createPeptideProteinGraph(orig_long_by_prot, "ProteinName", "PeptideSequenceNew")
quant_with_cls = addClusterMembership(orig_long_by_prot, peptide_protein_graph)
quant_with_cls[, NumProteins := uniqueN(ProteinName), by = "Cluster"]
quant_with_cls[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
quant_with_cls[, AllUnique := all(IsUnique), by = "Cluster"]

quant_with_cls[, log2Intensity := logintensity]
quant_with_cls[, log2IntensityNormalized := logintensity]

# Separated for easier comparison of filtering steps
# quant_no_single_shared = processIsoforms(quant_with_cls, remove_single_shared = TRUE, merge_identical = FALSE, remove_subsets = FALSE)
quant_merged_identical = processIsoforms(quant_with_cls, remove_single_shared = F, merge_identical = T, remove_subsets = FALSE)
# quant_no_subsets = processIsoforms(quant_merged_identical, remove_single_shared = F, merge_identical = F, remove_subsets = T)

quant_merged_identical[, Intensity := 2^log2Intensity]
quant_all_procd = MSstatsConvert::MSstatsPreprocess(quant_merged_identical[, .(ProteinName, PeptideSequence, PrecursorCharge = Charge, PSM, Run, Channel, Intensity)], 
                                                    annotation = unique(quant_merged_identical[, .(Run, Channel, BioReplicate, Condition, Mixture, TechRepMixture)]),
                                                    feature_columns = c("PeptideSequence", "PrecursorCharge"), 
                                                    remove_shared_peptides = FALSE, remove_single_feature_proteins = F, 
                                                    list(remove_features_with_few_measurements = FALSE,
                                                         summarize_multiple_psms = max))

quant_all_procd = as.data.table(quant_all_procd)
final_graph = createPeptideProteinGraph(quant_all_procd)
quant_all_procd = addClusterMembership(quant_all_procd, final_graph)

quant_all_procd[, NumProteins := uniqueN(ProteinName), by = "Cluster"]
quant_all_procd[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]

quant_all_procd[grepl("BRDT", ProteinName)]

brd_cluster = quant_all_procd[Cluster == 605]
brd_cluster[, log2IntensityNormalized := log(Intensity, 2)]
brd_cluster[, ProteinName := stringr::str_extract(ProteinName, "BRD[2-4T]")]
setnames(brd_cluster, "PrecursorCharge", "Charge")

plot_input_brd = copy(brd_cluster)
ch_order = c("126C", "127N", "127C", "128N", "128C",
             "129N", "129C", "130N", "130C", "131N")
plot_input_brd[, Channel := factor(Channel, levels = ch_order, ordered = TRUE)]

colors = palette.colors(8, "R4")[-1]
white = "#ffffff"
all_colors = c(colors, white)

pp_g_full = MSstatsWeightedSummary::createPeptideProteinGraph(unique(brd_cluster[, .(PeptideSequence, ProteinName)]))
ppdt_full = as.data.table(unique(brd_cluster[, .(PeptideSequence, ProteinName, IsUnique)]))
prots = unique(ppdt_full$ProteinName)
shared_colors = unique(ppdt_full[!(IsUnique), .(PeptideSequence)])[, .(PeptideSequence, Id = 1:.N)]
ppdt_full = merge(ppdt_full, shared_colors, by = "PeptideSequence", all.x = T)
ppdt_full[, ColorGroup := ifelse(IsUnique, 7, Id)]
colors_full = ppdt_full$ColorGroup
names(colors_full) = ppdt_full$PeptideSequence
# colors_full = colors_full + 1
colors_full = c(colors_full, rep(max(colors_full) + 1, 4))
names(colors_full)[(length(colors_full) - 3):length(colors_full)] = prots

par(mar=c(0,0,0,0)+.01)
plot(pp_g_full,
     vertex.color = all_colors[colors_full[names(V(pp_g_full))]],
     vertex.label = ifelse(grepl("BRD", names(V(pp_g_full))), names(V(pp_g_full)), ""),
     margin = -0.05, vertex.label.cex = 2)


pp_g_unique = MSstatsWeightedSummary::createPeptideProteinGraph(unique(brd_cluster[(IsUnique), .(PeptideSequence, ProteinName)]))
ppdt_unique = as.data.table(unique(brd_cluster[, .(PeptideSequence, ProteinName, IsUnique)]))
prots = unique(ppdt_unique$ProteinName)
shared_colors = unique(ppdt_unique[!(IsUnique), .(PeptideSequence)])[, .(PeptideSequence, Id = 1:.N)]
ppdt_unique = merge(ppdt_unique, shared_colors, by = "PeptideSequence", all.x = T)
ppdt_unique[, ColorGroup := ifelse(IsUnique, 7, Id)]
colors_unique = ppdt_unique$ColorGroup
names(colors_unique) = ppdt_unique$PeptideSequence
# colors_unique = colors_unique + 1
colors_unique = c(colors_unique, rep(max(colors_unique) + 1, 3))
names(colors_unique)[(length(colors_unique) - 2):length(colors_unique)] = prots

par(mar=c(0,0,0,0)+.01)
plot(pp_g_unique,
     vertex.color = all_colors[colors_unique[names(V(pp_g_unique))]],
     vertex.label = ifelse(grepl("BRD", names(V(pp_g_unique))), names(V(pp_g_unique)), ""),
     margin = -0.05, vertex.label.cex = 2)



color_vals = all_colors[colors_full[names(V(pp_g_full))]]
psms_list = names(V(pp_g_full))
names(color_vals) = psms_list

colors_full = colors_full[!grepl("BRD", colors_full)]
plot_input_2_3[, NameColor := ifelse(IsUnique, "unique", PeptideSequence)]

color_check = unique(unique(ppdt_full[, .(PeptideSequence, IsUnique, ColorGroup)])[, .(NameColor = ifelse(IsUnique, "unique", PeptideSequence),
                                                                                       ColorGroup)])
plot_input_2_3[, NameColor := factor(NameColor, levels = color_check$NameColor)]
plot_input_2_3[, Group := ifelse(grepl("DMSO", Condition), "Control", "Treatment")]
plot_input_2_3[, Time := stringr::str_extract(Condition, "_[0-9]+min")]
plot_input_2_3[, Time := stringr::str_replace(Time, "_", "")]
plot_input_2_3[, Time := stringr::str_replace(Time, "min", "")]
plot_input_2_3[, Time := stringr::str_replace(Time, "min", "")]

ggplot(plot_input_2_3[!(IsUnique)][Group == "Treatment"],
       aes(x = reorder(Time, as.numeric(Time)), y = logintensity, group = PSM, color = NameColor)) +
  geom_line(data = plot_input_2_3[(IsUnique)][Group == "Treatment"], alpha = 0.5, size = 1) +  
  geom_line(size = 1.5) +
  # geom_vline(xintercept = 5.5, color = "red", linetype = 2) +
  scale_color_manual(name = "PSM",
                     values = all_colors[color_check$ColorGroup]) +
  facet_grid(~ProteinName) +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 22),
        # axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 24), 
        legend.title = element_text(size = 20),
        legend.position = "none",
        legend.box = "vertical") +
  xlab("time") +
  ylab("log-intensity")

ggplot(plot_input_2_3[!(IsUnique)],
       aes(x = reorder(Time, as.numeric(Time)), y = logintensity, group = PSM, color = NameColor)) +
  geom_line(data = plot_input_2_3[(IsUnique)], alpha = 0.5, size = 1) +  
  geom_line(size = 1.5) +
  # geom_vline(xintercept = 5.5, color = "red", linetype = 2) +
  scale_color_manual(name = "PSM",
                     values = all_colors[color_check$ColorGroup]) +
  facet_grid(Group~ProteinName) +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 22),
        # axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 24), 
        legend.title = element_text(size = 20),
        legend.position = "none",
        legend.box = "vertical") +
  xlab("time") +
  ylab("log-intensity")







full_unique_input = copy(brd_cluster)
summary_input_1 = plot_input_brd[ProteinName != "BRD3"][!(IsUnique) | PSM %in% c("K.AQAEHAEK.E_3",
                                                                                 "K.LNLPDYYK.I_2",
                                                                                 "R.KLQDVFEFR.Y_3",
                                                                                 "R.LAELQEQLR.A_2")]
summary_input_1[, Channel := as.factor(as.character(Channel))]

full_summary = getWeightedProteinSummary(full_unique_input,
                                         norm = "Huber",
                                         norm_parameter = 1e-3)
summary_1_unique = getWeightedProteinSummary(summary_input_1[(IsUnique)],
                                             norm = "Huber",
                                             norm_parameter = 1e-3)
summary_1_shared = getWeightedProteinSummary(summary_input_1,
                                             norm = "Huber",
                                             norm_parameter = 1e-3)

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

