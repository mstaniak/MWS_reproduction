# Code required to reproduce Figure 2 (profile plot and peptide-protein graph for the protein degrader case study)
# Libraries ----
library(data.table)
library(MSstatsWeightedSummary)
library(ggplot2)
library(gridExtra)
library(igraph)
library(parallel)
library(splitstackshape)
library(MSstatsConvert)

# Data import and basic cleaning ----
g00_orig = fread("./input_data/protein_degrader/G0011_OBJ0037605_msstats_input_with_protein_map.txt")
## Remove decoys
g00_orig = g00_orig[!grepl("##", ProteinName, fixed = T)]
## Extract all matching proteins for each peptide
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

# quant_no_single_shared = processIsoforms(quant_with_cls, remove_single_shared = TRUE, merge_identical = FALSE, remove_subsets = FALSE)
quant_merged_identical = processIsoforms(quant_with_cls, remove_single_shared = F, merge_identical = T, remove_subsets = FALSE)
# quant_no_subsets = processIsoforms(quant_merged_identical, remove_single_shared = F, merge_identical = F, remove_subsets = T)
## commented: optional processing steps

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

# Peptide-protein graph plots for a single cluster ----
colors = palette.colors(8, "R4")[-1]
white = "#ffffff"
all_colors = c(colors, white)

pp_g_full = MSstatsWeightedSummary::createPeptideProteinGraph(unique(brd_cluster[, .(PeptideSequence, ProteinName)]))
# plot(pp_g_full)
brd_cluster[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
brd_cluster[!(IsUnique), unique(PSM)]
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

png("./plots/png/brd_graph_shared.png")
par(mar=c(0,0,0,0)+.01)
plot(pp_g_full,
     vertex.color = all_colors[colors_full[names(V(pp_g_full))]],
     vertex.label = ifelse(grepl("BRD", names(V(pp_g_full))), names(V(pp_g_full)), ""),
     margin = -0.05, vertex.label.cex = 2)
dev.off()

pp_g_unique = MSstatsWeightedSummary::createPeptideProteinGraph(unique(brd_cluster[(IsUnique), .(PeptideSequence, ProteinName)]))
ppdt_unique = as.data.table(unique(brd_cluster[(IsUnique), .(PeptideSequence, ProteinName, IsUnique)]))
prots = unique(ppdt_unique$ProteinName)
shared_colors = unique(ppdt_unique[!(IsUnique), .(PeptideSequence)])[, .(PeptideSequence, Id = 1:.N)]
ppdt_unique = merge(ppdt_unique, shared_colors, by = "PeptideSequence", all.x = T)
ppdt_unique[, ColorGroup := ifelse(IsUnique, 7, Id)]
colors_unique = ppdt_unique$ColorGroup
names(colors_unique) = ppdt_unique$PeptideSequence
# colors_unique = colors_unique + 1
colors_unique = c(colors_unique, rep(max(colors_unique) + 1, 3))
names(colors_unique)[(length(colors_unique) - 2):length(colors_unique)] = prots

png("./plots/png/brd_graph_noshared.png")
par(mar=c(0,0,0,0)+.01)
plot(pp_g_unique,
     vertex.color = all_colors[colors_unique[names(V(pp_g_unique))]],
     vertex.label = ifelse(grepl("BRD", names(V(pp_g_unique))), names(V(pp_g_unique)), ""),
     margin = -0.05, vertex.label.cex = 2)
dev.off()

# Profile plot for a single cluster ----
color_vals = all_colors[colors_full[names(V(pp_g_full))]]
psms_list = names(V(pp_g_full))
names(color_vals) = psms_list

colors_full = colors_full[!grepl("BRD", colors_full)]
brd_cluster[, NameColor := ifelse(IsUnique, "unique", PeptideSequence)]

color_check = unique(unique(ppdt_full[, .(PeptideSequence, IsUnique, ColorGroup)])[, .(NameColor = ifelse(IsUnique, "unique", PeptideSequence),
                                                                                       ColorGroup)])
brd_cluster[, NameColor := factor(NameColor, levels = color_check$NameColor)]
brd_cluster[, Group := ifelse(grepl("DMSO", Condition), "Control", "Treatment")]
brd_cluster[, Time := stringr::str_extract(Condition, "_[0-9]+min")]
brd_cluster[, Time := stringr::str_replace(Time, "_", "")]
brd_cluster[, Time := stringr::str_replace(Time, "min", "")]
brd_cluster[, Time := stringr::str_replace(Time, "min", "")]
brd_cluster[, Channel := factor(Channel, levels = ch_order, ordered = T)]

ggplot(brd_cluster[!(IsUnique)][ProteinName != "BRDT"],
       aes(x = reorder(Time, as.numeric(Time)), y = log2IntensityNormalized, group = PSM, color = NameColor)) +
  geom_line(data = brd_cluster[(IsUnique)][ProteinName != "BRDT"], size = 1) +  
  geom_line(size = 1.5) +
  # geom_vline(xintercept = 5.5, color = "red", linetype = 2) +
  scale_color_manual(name = "PSM",
                     values = all_colors[color_check$ColorGroup]) +
  facet_grid(Group~ProteinName) +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        # axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        legend.position = "none",
        legend.box = "vertical") +
  xlab("time") +
  ylab("log-intensity")
ggsave(filename = "plots/pdf/brd_shared_profiles_both_groups.pdf", width = 10, height = 5, scale = 1, dpi = 300)
ggsave(filename = "plots/png/brd_shared_profiles_both_groups.png", width = 10, height = 5, scale = 1, dpi = 300)
