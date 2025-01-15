# Code required to reproduce Figures:
# [main article] 3 and 4 (MSE of log-fold change estimation for the protein degrader study), 
# 7 (log-fold change estimation with outliers based on the protein degrader study)
# [supplement] 3 (MSE of log-fold change estimation for the protein degrader study by number of unique peptides),
# 4 (estimated log-fold change for the protein degrader study by number of unique peptides),
# 5 (example profile plot of the protein degrader study with outliers and estimated log-fold changes),

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
# Input data ----
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

full_summary = getWeightedProteinSummary(brd_cluster[(IsUnique)], "Huber", 1e-3)
gc_full = MSstatsTMT::groupComparisonTMT(makeMSstatsTMTInput(full_summary), cm_g00,
                                         use_log_file = FALSE)
full_summ_copy = copy(full_summary)
gc_full_triple = as.data.table(gc_full$ComparisonResult)
gc_full_triple[, Protein := stringr::str_extract(Protein, "BRD[2-4]")]
gc_full_triple = gc_full_triple[, .(Protein, Label, log2FC, SE)]

gc_gold_standard = copy(gc_full_triple)
gc_gold_standard[, TimeLabel := apply(stringr::str_extract_all(Label, "[0-9]+min", simplify = TRUE), 1, paste, sep = "-", collapse = "-")]

gs = gc_gold_standard[, .(Protein, TimeLabel, GoldStandard = log2FC)]
gs[, Time := stringr::str_extract(TimeLabel, "[0-9]+min")]
gs[, Time := stringr::str_replace(Time, "min", "")]
gs[, Time := as.numeric(Time)]

annot = unique(brd_cluster[, .(Channel, Condition, BioReplicate)])
annot[, Time := stringr::str_extract(Condition, "[0-9]+min")]
annot[, Time := stringr::str_replace(Time, "min", "")]
annot[, Time := as.numeric(Time)]
annot[, Group := ifelse(grepl("DMSO", Condition), "control", "G0011")]

# saveRDS(brd_cluster, "processed_data/protein_degrader/brd_cluster.RDS")
# saveRDS(gs, "processed_data/protein_degrader/gold_standard_tbl.RDS")
# saveRDS(annot, "processed_data/protein_degrader/brd_annot_tbl.RDS")

num_rep = 100

# sim_100 = mclapply(
#   c(1:5, 10),
#   function(num_unique) {
#     lapply(seq_len(num_rep), function(i) {
#       print(i)
#       tryCatch({
#         set.seed(i)
#         unique_peptides = unlist(lapply(unique(brd_cluster$ProteinName),
#                                         function(x) {
#                                           sample(brd_cluster[(IsUnique)][ProteinName == x, unique(PSM)], num_unique)
#                                         }), F, F)
#         summary_input_1_triple = brd_cluster[!(IsUnique) | PSM %in% unique_peptides]
#         summary_input_1_triple[, Channel := as.factor(as.character(Channel))]
# 
#         summary_1_triple_unique = getWeightedProteinSummary(summary_input_1_triple[(IsUnique)],
#                                                             norm = "Huber",
#                                                             norm_parameter = 1e-3,
#                                                             # weights_penalty = T, weights_penalty_param = 1e-5,
#                                                             initial_summary = "unique",
#                                                             tolerance = 1e-2, max_iter = 100)
# 
#         summary_1_all_triple = lapply(split(summary_input_1_triple, summary_input_1_triple[, ProteinName]),
#                                       function(x) {
#                                         x$IsUnique = TRUE
#                                         output = getWeightedProteinSummary(x,
#                                                                            norm = "Huber",
#                                                                            norm_parameter = 1e-3,
#                                                                            # weights_penalty = T, weights_penalty_param = 1e-5,
#                                                                            initial_summary = "unique",
#                                                                            tolerance = 1e-2, max_iter = 100)
#                                         output
#                                       })
# 
#         summary_1_triple_shared = getWeightedProteinSummary(summary_input_1_triple,
#                                                             norm = "Huber",
#                                                             norm_parameter = 1e-3,
#                                                             weights_penalty = T, weights_penalty_param = 1e-3,
#                                                             initial_summary = "unique",
#                                                             tolerance = 1e-2, max_iter = 100)
# 
#         list(summary_1_triple_unique,
#              summary_1_all_triple,
#              summary_1_triple_shared)
#       }, error = function(e) NULL)
#     })
#   },
#   mc.cores = 6)
# saveRDS(sim_100, "results/protein_degrader/sim_100_few.RDS")
sim_100 = readRDS("results/protein_degrader/sim_100_few.RDS")

psms_ints = brd_cluster[(IsUnique), .(ProteinName, PSM, Channel, log2IntensityNormalized)]
psms_ints_l = split(psms_ints, psms_ints$ProteinName)
psms_ints_l = lapply(psms_ints_l, function(x) {
  m = dcast(unique(x[, .(PSM, Channel, log2IntensityNormalized)]),
            PSM ~ Channel, value.var = "log2IntensityNormalized")
  res = as.matrix(m[, -1])
  row.names(res) = m$PSM
  res
})
psms_ints[, Label := ifelse(Channel == "126C", PSM, NA_character_)]

cors_1 = sort(apply(cor(t(psms_ints_l[[1]])), 1, mean, drop = T))
cors_2 = sort(apply(cor(t(psms_ints_l[[2]])), 1, mean, drop = T))
cors_3 = sort(apply(cor(t(psms_ints_l[[3]])), 1, mean, drop = T))

cors_dt = data.table(
  ProteinName = c(rep("BRD2", length(cors_1)),
                  rep("BRD3", length(cors_2)),
                  rep("BRD4", length(cors_3))),
  PSM = c(names(cors_1), names(cors_2), names(cors_3)),
  Corr = unname(c(cors_1, cors_2, cors_3))
)
cors_dt[, Rank := order(Corr), by = "ProteinName"]
cors_dt[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
cors_dt[(IsUnique)][Rank < 4]

# Example selection of features with outliers
ok_peptides = unlist(lapply(unique(cors_dt$ProteinName),
                            function(x) {
                              sample(cors_dt[(IsUnique) & Rank >= 4 & ProteinName == x, unique(PSM)], 2 - 1)
                            }), F, F)

noise_peptides = unlist(lapply(unique(cors_dt$ProteinName),
                               function(x) {
                                 sample(cors_dt[(IsUnique) & Rank < 4 & ProteinName == x, unique(PSM)], 1)
                               }), F, F)
unique_peptides = c(ok_peptides, noise_peptides)

summary_input_1_triple = brd_cluster[!(IsUnique) | PSM %in% unique_peptides]
summary_input_1_triple[, Channel := as.factor(as.character(Channel))]

summary_input_1_triple_2 = copy(summary_input_1_triple)
summary_input_1_triple_2[, Feature := factor(ifelse(IsUnique, "unique", "shared"), 
                                             levels = c("unique", "shared"),
                                             ordered = T)]
summary_input_1_triple_2[, Pattern := factor(ifelse(PSM %in% noise_peptides, "noisy", "consistent"),
                                             levels = c("noisy", "consistent"),
                                             ordered = T)]
ggplot(summary_input_1_triple_2,
       aes(x = Channel, y = log2IntensityNormalized, group = PSM,
           linetype = Feature, color = Pattern)) +
  geom_line(size = 1.2) +
  facet_grid(~ProteinName) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual(name = "pattern", values = c("#E69F00", "#56B4E9")) +
  scale_linetype_discrete(name = "peptide") +
  xlab("channel") +
  ylab("log-intensity") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 18), 
        legend.box = "vertical")
ggsave("plots/pdf/brd_noise_example.pdf", width = 10, height = 5, scale = 1, dpi = 300)
ggsave("plots/png/brd_noise_example.png", width = 10, height = 5, scale = 1, dpi = 300)


# sim_100_noise = mclapply(
#   c(1:5, 10),
#   function(num_unique) {
#     lapply(seq_len(num_rep), function(i) {
#       tryCatch({
#         set.seed(i)
#         ok_peptides = unlist(lapply(unique(cors_dt$ProteinName),
#                                     function(x) {
#                                       sample(cors_dt[(IsUnique) & Rank >= 4 & ProteinName == x, unique(PSM)], num_unique - 1)
#                                     }), F, F)
#         set.seed(i)
#         noise_peptides = unlist(lapply(unique(cors_dt$ProteinName),
#                                        function(x) {
#                                          sample(cors_dt[(IsUnique) & Rank < 4 & ProteinName == x, unique(PSM)], 1)
#                                        }), F, F)
#         unique_peptides = c(ok_peptides, noise_peptides)
#         example_input_1_triple = brd_cluster[!(IsUnique) | PSM %in% unique_peptides]
# 
#         summary_input_1_triple = copy(example_input_1_triple)
#         summary_input_1_triple[, Channel := as.factor(as.character(Channel))]
# 
#         summary_1_triple_unique = getWeightedProteinSummary(summary_input_1_triple[(IsUnique)],
#                                                             norm = "Huber",
#                                                             norm_parameter = 1e-3,
#                                                             initial_summary = "unique",
#                                                             tolerance = 1e-2, max_iter = 100)
# 
#         summary_1_all_triple = lapply(split(summary_input_1_triple, summary_input_1_triple[, ProteinName]),
#                                       function(x) {
#                                         x$IsUnique = TRUE
#                                         output = getWeightedProteinSummary(x,
#                                                                            norm = "Huber",
#                                                                            norm_parameter = 1e-3,
#                                                                            initial_summary = "unique",
#                                                                            tolerance = 1e-2, max_iter = 100)
#                                         output
#                                       })
# 
#         summary_1_triple_shared = getWeightedProteinSummary(summary_input_1_triple,
#                                                             norm = "Huber",
#                                                             norm_parameter = 1e-3,
#                                                             weights_penalty = TRUE, weights_penalty_param = 1e-4,
#                                                             initial_summary = "unique",
#                                                             tolerance = 1e-2, max_iter = 100)
# 
#         list(summary_1_triple_unique,
#              summary_1_all_triple,
#              summary_1_triple_shared)
#       }, error = function(e) NULL)
#     })
#   },
#   mc.cores = 6)
# saveRDS(sim_100_noise, "results/protein_degrader/sim_100_noise.RDS")
sim_100_noise = readRDS("results/protein_degrader/sim_100_noise.RDS")

gc_res_raw = lapply(sim_100, function(y) {
  lapply(y, function(x) {
    if (!is.null(x)) {
      all_sum = lapply(x[[2]], makeMSstatsTMTInput)
      all_sum = list(ProteinLevelData = rbindlist(lapply(all_sum, function(x) x$ProteinLevelData)),
                     FeatureLevelData = rbindlist(lapply(all_sum, function(x) x$FeatureLevelData)))
      summs = list(makeMSstatsTMTInput(x[[1]]) , 
                   all_sum,
                   makeMSstatsTMTInput(x[[3]]))
      
      gcs_comp = lapply(summs, function(x) {
        out = MSstatsTMT::groupComparisonTMT(x, cm_g00, 
                                             use_log_file = FALSE)
        out = as.data.table(out$ComparisonResult)
        out
      })
      gcs_comp = rbindlist(lapply(seq_along(gcs_comp), function(i) {
        out = gcs_comp[[i]]
        out$Method = c("unique", "all", "shared")[i]
        out
      }))
      
      gcs_comp[, TimeLabel := apply(stringr::str_extract_all(Label, "[0-9]+min", 
                                                             simplify = TRUE), 
                                    1, paste, sep = "-", collapse = "-")]
      gcs_comp
    } else {
      NULL
    }
  })
})

gc_res_raw_noise = lapply(sim_100_noise, function(y) {
  lapply(y, function(x) {
    if (!is.null(x)) {
      all_sum = lapply(x[[2]], makeMSstatsTMTInput)
      all_sum = list(ProteinLevelData = rbindlist(lapply(all_sum, function(x) x$ProteinLevelData)),
                     FeatureLevelData = rbindlist(lapply(all_sum, function(x) x$FeatureLevelData)))
      summs = list(makeMSstatsTMTInput(x[[1]]) , 
                   all_sum,
                   makeMSstatsTMTInput(x[[3]]))
      
      gcs_comp = lapply(summs, function(x) {
        out = MSstatsTMT::groupComparisonTMT(x, cm_g00, 
                                             use_log_file = FALSE)
        out = as.data.table(out$ComparisonResult)
        out
      })
      gcs_comp = rbindlist(lapply(seq_along(gcs_comp), function(i) {
        out = gcs_comp[[i]]
        out$Method = c("unique", "all", "shared")[i]
        out
      }))
      
      gcs_comp[, TimeLabel := apply(stringr::str_extract_all(Label, "[0-9]+min", 
                                                             simplify = TRUE), 
                                    1, paste, sep = "-", collapse = "-")]
      gcs_comp
    } else {
      NULL
    }
  })
})

gcs_comp_dt = rbindlist(lapply(seq_along(gc_res_raw), function(i) {
  x = gc_res_raw[[i]]
  x = rbindlist(lapply(seq_along(x), function(j) {
    y = x[[j]]
    if (!is.null(y)) {
      y$iter = j
      y
    } else {
      NULL
    }
  }))
  x$num_unique = c(1:5, 10)[i]
  x
}))
gcs_comp_dt[, Time := stringr::str_extract(TimeLabel, "[0-9]+min")]
gcs_comp_dt[, Time := stringr::str_replace(Time, "min", "")]
gcs_comp_dt[, Time := as.numeric(Time)]
gcs_comp_dt = merge(gcs_comp_dt, gs, by = c("Protein", "TimeLabel", "Time"))

gcs_comp_dt_noise = rbindlist(lapply(seq_along(gc_res_raw_noise), function(i) {
  x = gc_res_raw_noise[[i]]
  x = rbindlist(lapply(seq_along(x), function(j) {
    y = x[[j]]
    if (!is.null(y)) {
      y$iter = j
      y
    } else {
      NULL
    }
  }))
  x$num_unique = c(1:5, 10)[i]
  x
}))
gcs_comp_dt_noise[, Time := stringr::str_extract(TimeLabel, "[0-9]+min")]
gcs_comp_dt_noise[, Time := stringr::str_replace(Time, "min", "")]
gcs_comp_dt_noise[, Time := as.numeric(Time)]
gcs_comp_dt_noise = merge(gcs_comp_dt_noise, gs, 
                          by = c("Protein", "TimeLabel", "Time"))


### ----
colors = c("lightblue", "red", "purple")
colors = rev(colors)

ggplot(gcs_comp_dt[, .(SE = mean((log2FC - GoldStandard)^2)),
                   by = c("Method", "num_unique", "Protein", "iter")],
       aes(x = reorder(as.character(num_unique), num_unique), y = SE, fill = Method)) + # reorder(as.character(num_unique), num_unique)
  geom_boxplot() +
  facet_grid(~Protein, scales = "free_y") + # , nrow = 6, scales = "free_y"
  scale_fill_manual(
    name = "method",
    values = colors,
    labels = c(unique = "unique sampled peptides in each protein",
               shared = "proposed",
               all = "all sampled peptides in each protein")) +
  scale_color_manual(values = colors,
                     labels = c(unique = "unique sampled peptides in each protein",
                                shared = "proposed",
                                all = "all sampled peptides in each protein")) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18), 
        legend.direction = "vertical"
        # legend.box = "vertical"
        ) +
  xlab("number of unique peptides per protein") +
  ylab("MSE of log-fold change estimation") +
  ylim(c(0, 0.075))
ggsave("./plots/pdf/sim_brd_mse_few.pdf", width = 10, height = 5, scale = 1, dpi = 300)
ggsave("./plots/png/sim_brd_mse_few.png", width = 10, height = 5, scale = 1, dpi = 300)
colors = c(colors, "darkblue")

gt_data = unique(gcs_comp_dt[num_unique == 2, .(Time, Protein, GoldStandard, Method = "gold standard\n(all unique peptides)")])

ggplot(gcs_comp_dt[num_unique == 2],
       aes(x = reorder(as.character(Time), Time),
           y = log2FC, fill = Method)) + # , fill = Method
  geom_boxplot(linewidth = 1.2) + #, outlier.shape = NA
  geom_line(aes(x = reorder(as.character(Time), Time),
                y = GoldStandard, group = "1", color = Method),
            data = gt_data,
            inherit.aes = FALSE, size = 1.5) + # , color = "darkblue"
  facet_grid( ~ Protein) +
  scale_fill_manual(name = "method",
                    values = colors,
                    labels = c(unique = "unique sampled peptides in each protein",
                               shared = "proposed",
                               all = "all sampled peptides in each protein")) +
  scale_color_manual(values = colors[4], name = "") +
  geom_point(aes(x = reorder(as.character(Time), Time),
                 y = GoldStandard),
             inherit.aes = FALSE, color = "darkblue", size = 1.5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18), 
        legend.direction = "vertical",
        # legend.box = "vertical",
        # legend.key.spacing = unit(5, "pt"),
        # legend.spacing = unit(5, "pt"),
        legend.text.position = "right") +
  xlab("time") +
  ylab("log-fold change with respect to control") +
  ylim(c(-1.3, 0.5))
ggsave("plots/pdf/sims_brd_few.pdf", device = "pdf", width = 10, height = 5, units = "in", scale = 1, dpi = 300)
ggsave("plots/png/sims_brd_few.png", device = "png", width = 10, height = 5, units = "in", scale = 1, dpi = 300)

ggplot(gcs_comp_dt[num_unique %in% 1:3],
       aes(x = reorder(as.character(Time), Time),
           y = log2FC, fill = Method)) + # , fill = Method
  geom_boxplot(linewidth = 1.2) + #, outlier.shape = NA
  geom_line(aes(x = reorder(as.character(Time), Time),
                y = GoldStandard, group = "1", color = Method),
            data = gt_data,
            inherit.aes = FALSE, size = 1.5) + # , color = "darkblue"
  facet_grid(num_unique ~ Protein) +
  scale_fill_manual(name = "method",
                    values = colors,
                    labels = c(unique = "unique sampled peptides in each protein",
                               shared = "proposed",
                               all = "all sampled peptides in each protein")) +
  scale_color_manual(values = colors[4], name = "") +
  geom_point(aes(x = reorder(as.character(Time), Time),
                 y = GoldStandard),
             inherit.aes = FALSE, color = "darkblue", size = 1.5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18), 
        legend.direction = "vertical",
        # legend.box = "vertical",
        # legend.key.spacing = unit(5, "pt"),
        # legend.spacing = unit(5, "pt"),
        legend.text.position = "right") +
  xlab("time") +
  ylab("log-fold change with respect to control") +
  ylim(c(-1.3, 0.5))
ggsave("./plots/pdf/brd_sims_by_unique_full.pdf", width = 10, height = 5, scale = 1, dpi = 300)
ggsave("./plots/png/brd_sims_by_unique_full.png", width = 10, height = 5, scale = 1, dpi = 300)

ggplot(gcs_comp_dt_noise[, .(SE = mean((log2FC - GoldStandard)^2)),
                          by = c("Method", "num_unique", "Protein", "iter")],
       aes(x = reorder(as.character(num_unique), num_unique), y = SE, fill = Method)) + # reorder(as.character(num_unique), num_unique)
  geom_boxplot() +
  facet_wrap(~Protein, scales = "free_y") + # , nrow = 6, scales = "free_y"
  scale_fill_manual(name = "method",
                    values = colors,
                    labels = c(unique = "unique",
                               shared = "proposed",
                               all = "all")) +
  scale_color_manual(values = colors,
                     labels = c(unique = "unique",
                                shared = "proposed",
                                all = "all")) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18), 
        legend.box = "vertical") +
  xlab("number of unique peptides per protein") +
  ylab("MSE of log-fold change estimation") 
ggsave("./plots/pdf/sim_brd_mse_noise.pdf", device = "pdf", width = 10, height = 5, scale = 1, dpi = 300)
ggsave("./plots/png/sim_brd_mse_noise.png", device = "png", width = 10, height = 5, scale = 1, dpi = 300)

ggplot(gcs_comp_dt_noise[num_unique == 2],
       aes(x = reorder(as.character(Time), Time),
           y = log2FC, fill = Method)) + # , fill = Method
  geom_boxplot(linewidth = 1.2) + #, outlier.shape = NA
  geom_line(aes(x = reorder(as.character(Time), Time),
                y = GoldStandard, group = "1", color = Method),
            data = gt_data,
            inherit.aes = FALSE, size = 1.5) + # , color = "darkblue"
  facet_grid( ~ Protein) +
  scale_fill_manual(name = "method",
                    values = colors,
                    labels = c(unique = "unique sampled peptides in each protein",
                               shared = "proposed",
                               all = "all sampled peptides in each protein")) +
  scale_color_manual(values = colors[4], name = "") +
  geom_point(aes(x = reorder(as.character(Time), Time),
                 y = GoldStandard),
             inherit.aes = FALSE, color = "darkblue", size = 1.5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18), 
        # legend.box = "vertical",
        # legend.key.spacing = unit(5, "pt"),
        legend.direction = "vertical",
        # legend.spacing = unit(5, "pt"),
        legend.text.position = "right") +
  xlab("time") +
  ylab("log-fold change wrt control") +
  ylim(c(-1.3, 0.5))
ggsave("./plots/pdf/sims_brd_noise.pdf", device = "pdf", width = 10, height = 5, units = "in", scale = 1, dpi = 300)
ggsave("./plots/png/sims_brd_noise.png", device = "png", width = 10, height = 5, units = "in", scale = 1, dpi = 300)

ggplot(gcs_comp_dt_noise[num_unique %in% 1:3],
       aes(x = reorder(as.character(Time), Time),
           y = log2FC, fill = Method)) + # , fill = Method
  geom_boxplot(linewidth = 1.2) + #, outlier.shape = NA
  geom_line(aes(x = reorder(as.character(Time), Time),
                y = GoldStandard, group = "1", color = Method),
            data = gt_data,
            inherit.aes = FALSE, size = 1.5) + # , color = "darkblue"
  facet_grid(num_unique ~ Protein) +
  scale_fill_manual(name = "method",
                    values = colors,
                    labels = c(unique = "unique sampled peptides in each protein",
                               shared = "proposed",
                               all = "all sampled peptides in each protein")) +
  scale_color_manual(values = colors[4], name = "") +
  geom_point(aes(x = reorder(as.character(Time), Time),
                 y = GoldStandard),
             inherit.aes = FALSE, color = "darkblue", size = 1.5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18), 
        # legend.box = "vertical",
        # legend.key.spacing = unit(5, "pt"),
        legend.direction = "vertical",
        # legend.spacing = unit(5, "pt"),
        legend.text.position = "right") +
  xlab("time") +
  ylab("log-fold change wrt control") +
  ylim(c(-1.3, 0.5))
ggsave("./plots/pdf/brd_sims_noise_by_unique_full.pdf", device = "pdf", width = 10, height = 5, units = "in", scale = 1, dpi = 300)
ggsave("./plots/png/brd_sims_noise_by_unique_full.png", device = "png", width = 10, height = 5, units = "in", scale = 1, dpi = 300)
