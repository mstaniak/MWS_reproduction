# Code required to reproduce Figures:
# [main] 10 (comparison of statistical inference results based on different summarization methods),
# [supplemental] 7 (example cluster from the thermal profiling case study)
# Libraries ----
library(data.table)
library(MSstatsWeightedSummary)
library(ggplot2)
library(ggvenn)
library(patchwork)
# Functions ----
normalizeProteins = function(summarized_data) {
  n_runs = data.table::uniqueN(summarized_data$Run, na.rm = TRUE)
  if ((n_runs > 1)) {
    group_info = unique(summarized_data$Condition)
    if (is.element("Norm", group_info)) {
      summarized_data[!is.na(Abundance), `:=`(NumRuns, data.table::uniqueN(Run, 
                                                                           na.rm = TRUE)), by = "Protein"]
      summarized_data[!is.na(Abundance), `:=`(NumRunsWithNorm, MSstatsTMT:::.countRunsWithNorm(Run, 
                                                                                               Condition)), by = "Protein"]
      summarized_data[!is.na(Abundance), `:=`(NormalizationAbundance, 
                                              MSstatsTMT:::.getNormalizationAbundance(Abundance, Condition)), 
                      by = c("Protein", "Run")]
      summarized_data[!is.na(Abundance), `:=`(MedianNormalized, MSstatsTMT:::.getRunsMedian(.SD)), 
                      by = "Protein", .SDcols = c("Run", "NormalizationAbundance")]
      summarized_data[!is.na(Abundance), `:=`(Diff, MedianNormalized - 
                                                NormalizationAbundance)]
      summarized_data[!is.na(Abundance), `:=`(NormalizedAbundance, 
                                              Abundance + Diff)]
      summarized_data[, `:=`(Abundance, ifelse(NumRuns > 1 & NumRunsWithNorm > 
                                                 1, NormalizedAbundance, Abundance))]
      summarized_data[, `:=`(Diff, NULL)]
    } else {
      NULL
    }
  }
  summarized_data[, list(Mixture, TechRepMixture, Run, Channel, Protein, 
                         Abundance, BioReplicate, Condition)]  
}
# Input data ----
onepot_int_cls_tbl = readRDS("processed_data/onepot_tpp/onepot_sub_int_cls.RDS")
tpp_int_cls_tbl = readRDS("processed_data/onepot_tpp/tpp_sub_int_cls.RDS")
# Protein cluster processing and descriptive statistics ----
onepot_int_cls_each_uni = onepot_int_cls_tbl[(HasUnique)]
onepot_int_cls_each_uni[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
onepot_int_cls_each_uni[, NumProteins := uniqueN(ProteinName), by = "Cluster"]
onepot_int_cls_each_uni[, HasUnique := any(IsUnique), by = "ProteinName"]
onepot_int_cls_each_uni = onepot_int_cls_each_uni[(HasUnique)]
onepot_int_cls_each_uni[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
onepot_int_cls_each_uni[, NumProteins := uniqueN(ProteinName), by = "Cluster"]

tpp_int_cls_each_uni = tpp_int_cls_tbl[(HasUnique)]
tpp_int_cls_each_uni[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
tpp_int_cls_each_uni[, NumProteins := uniqueN(ProteinName), by = "Cluster"]
tpp_int_cls_each_uni[, Intensity := 2^log2IntensityNormalized]
tpp_int_cls_each_uni[, HasUnique := any(IsUnique), by = "ProteinName"]
tpp_int_cls_each_uni = tpp_int_cls_each_uni[(HasUnique)]
tpp_int_cls_each_uni[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
tpp_int_cls_each_uni[, NumProteins := uniqueN(ProteinName), by = "Cluster"]

onepot_split = split(onepot_int_cls_each_uni[NumProteins > 1], onepot_int_cls_each_uni[NumProteins > 1, Cluster])
tpp_split = split(tpp_int_cls_each_uni[NumProteins > 1], tpp_int_cls_each_uni[NumProteins > 1, Cluster])

length(onepot_split)
length(tpp_split)
# Summarization ----
# onepot_shared_summaries_int = lapply(onepot_split,
#                               function(x) {
#                                 print(unique(x$Cluster))
#                                 tryCatch({
#                                   getWeightedProteinSummary(x, "Huber", 1e-6,
#                                                             max_iter = 100, tolerance = 1e-2, initial_summary = "flat")
#                                 }, error = function(e) NULL)
#                               })
# 
# onepot_unique_summaries_int = lapply(onepot_split,
#                               function(x) {
#                                 print(unique(x$Cluster))
#                                 if (nrow(x[(IsUnique)]) > 0) {
#                                   getWeightedProteinSummary(x[(IsUnique)], "Huber", 1e-6, max_iter = 100, tolerance = 1e-2)
#                                 } else {
#                                   NULL
#                                 }
#                               })
# onepot_all_summaries_int = lapply(onepot_split,
#                            function(x) {
#                              print(unique(x$Cluster))
#                              lapply(split(x, x$ProteinName), function(y) {
#                                y$IsUnique = TRUE
#                                getWeightedProteinSummary(y, "Huber", 1e-6, max_iter = 100, tolerance = 1e-2)
#                              })
#                            })

# table(sapply(onepot_shared_summaries_int, is.null))
# table(sapply(onepot_unique_summaries_int, is.null))
# table(sapply(onepot_all_summaries_int, is.null))

# tpp_shared_summaries_int = lapply(tpp_split,
#                               function(x) {
#                                 print(unique(x$Cluster))
#                                 tryCatch({
#                                   getWeightedProteinSummary(x, "Huber", 1e-6,
#                                                             max_iter = 100, tolerance = 1e-2, initial_summary = "flat")
#                                 }, error = function(e) NULL)
#                               })

# tpp_unique_summaries_int = lapply(tpp_split,
#                               function(x) {
#                                 print(unique(x$Cluster))
#                                 if (nrow(x[(IsUnique)]) > 0) {
#                                   getWeightedProteinSummary(x[(IsUnique)], "Huber", 1e-6, max_iter = 100, tolerance = 1e-2)
#                                 } else {
#                                   NULL
#                                 }
#                               })
# tpp_all_summaries_int = lapply(tpp_split,
#                            function(x) {
#                              print(unique(x$Cluster))
#                              lapply(split(x, x$ProteinName), function(y) {
#                                y$IsUnique = TRUE
#                                getWeightedProteinSummary(y, "Huber", 1e-6, max_iter = 100, tolerance = 1e-2)
#                              })
#                            })

# saveRDS(onepot_shared_summaries_int, "processed_data/onepot_tpp/onepot_sh_summs.RDS")
# saveRDS(onepot_unique_summaries_int, "processed_data/onepot_tpp/onepot_un_summs.RDS")
# saveRDS(onepot_all_summaries_int, "processed_data/onepot_tpp/onepot_al_summs.RDS")

# saveRDS(tpp_shared_summaries_int, "processed_data/onepot_tpp/tpp_sh_summs.RDS")
# saveRDS(tpp_unique_summaries_int, "processed_data/onepot_tpp/tpp_un_summs.RDS")
# saveRDS(tpp_all_summaries_int, "processed_data/onepot_tpp/tpp_al_summs.RDS")

onepot_shared_summaries_int = readRDS("processed_data/onepot_tpp/onepot_sh_summs.RDS")
onepot_unique_summaries_int = readRDS("processed_data/onepot_tpp/onepot_un_summs.RDS")
onepot_all_summaries_int = readRDS("processed_data/onepot_tpp/onepot_al_summs.RDS")

tpp_shared_summaries_int = readRDS("processed_data/onepot_tpp/tpp_sh_summs.RDS")
tpp_unique_summaries_int = readRDS("processed_data/onepot_tpp/tpp_un_summs.RDS")
tpp_all_summaries_int = readRDS("processed_data/onepot_tpp/tpp_al_summs.RDS")


onepot_protein_data_shared = rbindlist(lapply(onepot_shared_summaries_int, proteinData))
onepot_protein_data_unique = rbindlist(lapply(onepot_unique_summaries_int, proteinData))
onepot_protein_data_all = rbindlist(lapply(onepot_all_summaries_int, function(x) rbindlist(lapply(x, proteinData))))

onepot_feat_data_shared = rbindlist(lapply(onepot_shared_summaries_int, featureData))
onepot_feat_data_unique = rbindlist(lapply(onepot_unique_summaries_int, featureData))
onepot_feat_data_all = rbindlist(lapply(onepot_all_summaries_int, function(x) rbindlist(lapply(x, featureData))))

uniqueN(onepot_feat_data_shared$ProteinName)
uniqueN(onepot_feat_data_unique$ProteinName)
uniqueN(onepot_feat_data_all$ProteinName)

# Group comparison -----
cm_onepot = readRDS("input_data/onepot_tpp/contrast_matrix.RDS")
gc_sh_onepot = MSstatsTMT::groupComparisonTMT(list(ProteinLevelData = onepot_protein_data_shared,
                                                   FeatureLevelData = onepot_feat_data_shared), cm_onepot, use_log_file = FALSE)
gc_un_onepot = MSstatsTMT::groupComparisonTMT(list(ProteinLevelData = onepot_protein_data_unique,
                                                   FeatureLevelData = onepot_feat_data_unique), cm_onepot, use_log_file = FALSE)
gc_al_onepot = MSstatsTMT::groupComparisonTMT(list(ProteinLevelData = onepot_protein_data_all,
                                                   FeatureLevelData = onepot_feat_data_all), cm_onepot, use_log_file = FALSE)

gc_sh_dt_onepot = as.data.table(gc_sh_onepot$ComparisonResult)
gc_un_dt_onepot = as.data.table(gc_un_onepot$ComparisonResult)
gc_al_dt_onepot = as.data.table(gc_al_onepot$ComparisonResult)


tpp_protein_data_shared = rbindlist(lapply(tpp_shared_summaries_int, proteinData))
tpp_protein_data_unique = rbindlist(lapply(tpp_unique_summaries_int, proteinData))
tpp_protein_data_all = rbindlist(lapply(tpp_all_summaries_int, function(x) rbindlist(lapply(x, proteinData))))

tpp_protein_data_shared = normalizeProteins(tpp_protein_data_shared)
tpp_protein_data_unique = normalizeProteins(tpp_protein_data_unique)
tpp_protein_data_all = normalizeProteins(tpp_protein_data_all)

tpp_feature_data_shared = rbindlist(lapply(tpp_shared_summaries_int, featureData))
tpp_feature_data_unique = rbindlist(lapply(tpp_unique_summaries_int, featureData))
tpp_feature_data_all = rbindlist(lapply(tpp_all_summaries_int, function(x) rbindlist(lapply(x, featureData))))

cm_curve = readRDS("./input_data/onepot_tpp/Contrast_thermal.RDS")
gc_sh_tpp = MSstatsTMT::groupComparisonTMT(list(ProteinLevelData = tpp_protein_data_shared[Condition != "Norm"],
                                            FeatureLevelData = tpp_feature_data_shared), cm_curve, use_log_file = FALSE)
gc_un_tpp = MSstatsTMT::groupComparisonTMT(list(ProteinLevelData = tpp_protein_data_unique[Condition != "Norm"],
                                            FeatureLevelData = tpp_feature_data_unique), cm_curve, use_log_file = FALSE)
gc_al_tpp = MSstatsTMT::groupComparisonTMT(list(ProteinLevelData = tpp_protein_data_all[Condition != "Norm"],
                                            FeatureLevelData = tpp_feature_data_all), cm_curve, use_log_file = FALSE)

gc_sh_dt_tpp = as.data.table(gc_sh_tpp$ComparisonResult)
gc_un_dt_tpp = as.data.table(gc_un_tpp$ComparisonResult)
gc_al_dt_tpp = as.data.table(gc_al_tpp$ComparisonResult)

interactors = readRDS("input_data/onepot_tpp/kinases_and_direct_interactors.RDS")
interactors_fixed = stringr::str_replace(interactors, "\\-1", "")

is_interact_gc_sh_onepot = rep(FALSE, nrow(gc_sh_dt_onepot))
is_interact_gc_un_onepot = rep(FALSE, nrow(gc_un_dt_onepot))
is_interact_gc_al_onepot = rep(FALSE, nrow(gc_al_dt_onepot))

is_interact_gc_sh_tpp = rep(FALSE, nrow(gc_sh_dt_tpp))
is_interact_gc_un_tpp = rep(FALSE, nrow(gc_un_dt_tpp))
is_interact_gc_al_tpp = rep(FALSE, nrow(gc_al_dt_tpp))

for (interactor in interactors_fixed) {
  is_interact_gc_sh_onepot = is_interact_gc_sh_onepot | stringr::str_detect(gc_sh_dt_onepot$Protein, interactor)
  is_interact_gc_un_onepot = is_interact_gc_un_onepot | stringr::str_detect(gc_un_dt_onepot$Protein, interactor)
  is_interact_gc_al_onepot = is_interact_gc_al_onepot | stringr::str_detect(gc_al_dt_onepot$Protein, interactor)
  
  is_interact_gc_sh_tpp = is_interact_gc_sh_tpp | stringr::str_detect(gc_sh_dt_tpp$Protein, interactor)
  is_interact_gc_un_tpp = is_interact_gc_un_tpp | stringr::str_detect(gc_un_dt_tpp$Protein, interactor)
  is_interact_gc_al_tpp = is_interact_gc_al_tpp | stringr::str_detect(gc_al_dt_tpp$Protein, interactor)
}

get_sig_prots = function(dt, lgl) dt[lgl][pvalue < 0.05, as.character(unique(Protein))]

sig_prots_sh_onepot = get_sig_prots(gc_sh_dt_onepot, is_interact_gc_sh_onepot)
sig_prots_un_onepot = get_sig_prots(gc_un_dt_onepot, is_interact_gc_un_onepot)
sig_prots_al_onepot = get_sig_prots(gc_al_dt_onepot, is_interact_gc_al_onepot)

sig_prots_sh_tpp = get_sig_prots(gc_sh_dt_tpp, is_interact_gc_sh_tpp)
sig_prots_un_tpp = get_sig_prots(gc_un_dt_tpp, is_interact_gc_un_tpp)
sig_prots_al_tpp = get_sig_prots(gc_al_dt_tpp, is_interact_gc_al_tpp)

fix_prots = function(x) unlist(stringr::str_split(x, ";"), F, F)

sig_prots_sh_onepot = fix_prots(sig_prots_sh_onepot)
sig_prots_un_onepot = fix_prots(sig_prots_un_onepot)
sig_prots_al_onepot = fix_prots(sig_prots_al_onepot)

sig_prots_sh_tpp = fix_prots(sig_prots_sh_tpp)
sig_prots_un_tpp = fix_prots(sig_prots_un_tpp)
sig_prots_al_tpp = fix_prots(sig_prots_al_tpp)

grouped_with_interact_sh_tpp = setdiff(unlist(stringr::str_split(gc_sh_dt_tpp[is_interact_gc_sh_tpp & grepl(";", Protein), unique(Protein)], ";"), F, F), interactors_fixed)
grouped_with_interact_un_tpp = setdiff(unlist(stringr::str_split(gc_un_dt_tpp[is_interact_gc_un_tpp & grepl(";", Protein), unique(Protein)], ";"), F, F), interactors_fixed)
grouped_with_interact_al_tpp = setdiff(unlist(stringr::str_split(gc_al_dt_tpp[is_interact_gc_al_tpp & grepl(";", Protein), unique(Protein)], ";"), F, F), interactors_fixed)

grouped_with_interact_sh_onepot = setdiff(unlist(stringr::str_split(gc_sh_dt_onepot[is_interact_gc_sh_onepot & grepl(";", Protein), unique(Protein)], ";"), F, F), interactors_fixed)
grouped_with_interact_un_onepot = setdiff(unlist(stringr::str_split(gc_un_dt_onepot[is_interact_gc_un_onepot & grepl(";", Protein), unique(Protein)], ";"), F, F), interactors_fixed)
grouped_with_interact_al_onepot = setdiff(unlist(stringr::str_split(gc_al_dt_onepot[is_interact_gc_al_onepot & grepl(";", Protein), unique(Protein)], ";"), F, F), interactors_fixed)

num_size = 9
lab_size = 9
tit_size = 20

# Graphical comparison of identified diferentially abundant proteins ----
ggvenn::ggvenn(list(TPP = setdiff(sig_prots_sh_tpp, grouped_with_interact_sh_tpp),
                    onePot = setdiff(sig_prots_sh_onepot, grouped_with_interact_sh_onepot),
                    known = interactors_fixed),
               show_percentage = FALSE,
               fill_color = c("red", "red", "white"),
               text_size = num_size,
               set_name_size = lab_size) +
  # ggtitle("proposed") +
  coord_cartesian(clip="off") +
  theme(plot.title = element_text(size = tit_size, hjust = 0.5)) +
  ggvenn::ggvenn(list(TPP = setdiff(sig_prots_un_tpp, grouped_with_interact_un_tpp),
                      onePot = setdiff(sig_prots_un_onepot, grouped_with_interact_un_onepot),
                      known = interactors_fixed),
                 show_percentage = FALSE,
                 fill_color = c("lightblue", "lightblue", "white"),
                 text_size = num_size,
                 set_name_size = lab_size) +
  # ggtitle("unique-only") +
  coord_cartesian(clip="off") +
  theme(plot.title = element_text(size = tit_size, hjust = 0.5)) +
  ggvenn::ggvenn(list(TPP = setdiff(sig_prots_al_tpp, grouped_with_interact_al_tpp),
                      onePot = setdiff(sig_prots_al_onepot, grouped_with_interact_al_onepot),
                      known = interactors_fixed),
                 show_percentage = FALSE,
                 fill_color = c("purple", "purple", "white"),
                 text_size = num_size,
                 set_name_size = lab_size)  +
  # ggtitle("all-peptides") +
  coord_cartesian(clip="off") +
  theme(plot.title = element_text(size = tit_size, hjust = 0.5),
        plot.margin = margin(0, 0, 0, 0))
ggsave("plots/pdf/comp3_all_rev.pdf", device = "pdf", scale = 1, width = 10, height = 5, dpi = 300)
ggsave("plots/png/comp3_all_rev.png", device = "png", scale = 1, width = 10, height = 5, dpi = 300)

# Descriptive statistics ----
tpp_int_cls_each_uni[NumProteins > 1, uniqueN(Cluster)]
tpp_int_cls_each_uni[NumProteins > 1, uniqueN(ProteinName)]

onepot_int_cls_each_uni[NumProteins > 1, uniqueN(Cluster)]
onepot_int_cls_each_uni[NumProteins > 1, uniqueN(ProteinName)]

feat_plot_1_cl = tpp_int_cls_each_uni[Cluster == 5698]
feat_plot_1_cl[, Group := ifelse(grepl("treated", Condition), "treatment", "control")]
tpp_plot_comp_1_cl = rbind(
  cbind(tpp_protein_data_shared[Protein %in% unique(feat_plot_1_cl$ProteinName)], Method = "shared"),
  cbind(tpp_protein_data_unique[Protein %in% unique(feat_plot_1_cl$ProteinName)], Method = "unique"),
  cbind(tpp_protein_data_all[Protein %in% unique(feat_plot_1_cl$ProteinName)], Method = "all")
)
tpp_plot_comp_1_cl[, Group := ifelse(grepl("treated", Condition), "treatment", "control")]
setnames(tpp_plot_comp_1_cl, "Protein", "ProteinName")
curve_annot = readRDS("input_data/onepot_tpp/Annotation_file_thermal.RDS")
feat_plot_1_cl = merge(feat_plot_1_cl, unique(as.data.table(curve_annot)[, .(Channel, temperature)]),
                       by = "Channel")
tpp_plot_comp_1_cl = merge(tpp_plot_comp_1_cl, unique(as.data.table(curve_annot)[, .(Channel, temperature)]),
                       by = "Channel")

# Visualization of summaries and profiles for a single cluster ----
ggplot(feat_plot_1_cl[Condition != "Norm"], 
       aes(x = temperature, y = log2IntensityNormalized,
           group = paste(PSM, Run), 
           color = Group,
           linetype = factor(ifelse(IsUnique, "unique", "shared"), levels = c("unique", "shared"), ordered = T))) +
  # geom_line() +
  geom_line(aes(x = temperature, y = Abundance, group = paste(Group, Run), color = Group),
            data = tpp_plot_comp_1_cl[Condition != "Norm"], size = 1.2, inherit.aes = FALSE) +
  scale_linetype_discrete(name = "peptide") +
  scale_color_manual(name = "group", values = c("#E69F00", "#56B4E9")) +
  facet_grid(ifelse(Method == "shared", "proposed", ifelse(Method == "unique", "unique", "all"))~ProteinName) +
  ylab("log-intensity") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16), 
        legend.box = "vertical")
ggsave("plots/pdf/p16591cl.pdf", device = "pdf", scale = 1, width = 10, height = 5, unit = "in", dpi = 300)
ggsave("plots/png/p16591cl.png", device = "png", scale = 1, width = 10, height = 5, unit = "in", dpi = 300)

ggplot(feat_plot_1_cl[Condition != "Norm"], 
       aes(x = temperature, y = log2IntensityNormalized,
           group = paste(PSM, Run), 
           color = Group,
           linetype = factor(ifelse(IsUnique, "unique", "shared"), levels = c("unique", "shared"), ordered = T))) +
  geom_line() +
  # geom_line(aes(x = temperature, y = Abundance, group = paste(Group, Run), color = Group),
  #           data = tpp_plot_comp_1_cl[Condition != "Norm"], size = 1.2, inherit.aes = FALSE) +
  scale_linetype_discrete(name = "peptide") +
  scale_color_manual(name = "group", values = c("#E69F00", "#56B4E9")) +
  facet_grid(~ProteinName) +
  ylab("log-intensity") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16), 
        legend.box = "vertical")
ggsave("plots/pdf/p16591cl_feature.pdf", device = "pdf", scale = 1, width = 10, height = 5, unit = "in", dpi = 300)
ggsave("plots/png/p16591cl_feature.png", device = "png", scale = 1, width = 10, height = 5, unit = "in", dpi = 300)