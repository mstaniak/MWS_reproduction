# Code required to reproduce Figures:
# [supplement] 8 (comparison of MSE of log-fold change estimation based on different starting points),
# 9 (comparison of MSE of log-fold change estimation based on different starting points in the presence of subset proteins),
# 10 (estimated log-fold changes with a subset protein based on the protein degrader study),
# 11 (influence of the loss function selection on log-fold change estimation in the presence of outliers based on the protein degrader study)
# Libraries ----
library(data.table)
library(MSstatsWeightedSummary)
library(ggplot2)
library(gridExtra)
library(parallel)
library(MSstatsConvert)
# Input data ----
brd_cluster = readRDS("processed_data/protein_degrader/brd_cluster.RDS")
gs = readRDS("processed_data/protein_degrader/gold_standard_tbl.RDS")
annot = readRDS("processed_data/protein_degrader/brd_annot_tbl.RDS")

cm_g00 = fread("./input_data/protein_degrader/OBJ037605_G0011_comparison_matrix.txt")
group_labels = cm_g00$V1
cm_g00[, V1 := NULL]
cm_g00 = as.matrix(cm_g00)
row.names(cm_g00) = group_labels

# Starting point influence with subset proteins: simulation study ----
num_rep = 100
# sim_100_starts_sub = lapply(
#   1:3 , 
#   function(num_unique) {
#     mclapply(seq_len(num_rep), function(i) {
#       print(i)
#       lapply(c("BRD2", "BRD3", "BRD4"), function(subset_protein) { #
#         tryCatch({
#           set.seed(i)
#           unique_peptides = unlist(lapply(setdiff(unique(brd_cluster$ProteinName), subset_protein),
#                                           function(x) {
#                                             sample(brd_cluster[(IsUnique)][ProteinName == x, unique(PSM)], num_unique)
#                                           }), F, F)
#           summary_input_1_triple = brd_cluster[!(IsUnique) | PSM %in% unique_peptides]
#           summary_input_1_triple[, Channel := as.factor(as.character(Channel))]
#           
#           summary_1_triple_unique = getWeightedProteinSummary(summary_input_1_triple[(IsUnique)],
#                                                               norm = "Huber",
#                                                               norm_parameter = 1e-3,
#                                                               initial_summary = "unique",
#                                                               tolerance = 1e-2, max_iter = 100)
#           
#           summary_1_all_triple = lapply(split(summary_input_1_triple, summary_input_1_triple[, ProteinName]),
#                                         function(x) {
#                                           x$IsUnique = TRUE
#                                           output = getWeightedProteinSummary(x,
#                                                                              norm = "Huber",
#                                                                              norm_parameter = 1e-3,
#                                                                              initial_summary = "unique",
#                                                                              tolerance = 1e-2, max_iter = 100)
#                                           output
#                                         })
#           
#           summary_1_triple_shared = getWeightedProteinSummary(summary_input_1_triple,
#                                                               norm = "Huber",
#                                                               norm_parameter = 1e-3,
#                                                               initial_summary = "flat",
#                                                               tolerance = 1e-2, max_iter = 100)
#           summary_1_triple_shared_mixed = getWeightedProteinSummary(summary_input_1_triple,
#                                                                     norm = "Huber",
#                                                                     norm_parameter = 1e-3,
#                                                                     initial_summary = "flat shared",
#                                                                     tolerance = 1e-2, max_iter = 100)
#           
#           list(unique = summary_1_triple_unique,
#                all = summary_1_all_triple,
#                proposed = summary_1_triple_shared,
#                proposed_mixed = summary_1_triple_shared_mixed,
#                iter =  i,
#                num_unique = num_unique,
#                subset_protein = subset_protein)
#         }, error = function(e) NULL)
#       })
#     }, mc.cores = 5)
#   })
# saveRDS(sim_100_starts_sub, "./results/protein_degrader/sim_100_starts_sub.RDS")
sim_100_starts_sub = readRDS("./results/protein_degrader/sim_100_starts_sub.RDS")

gc_res_raw_starts_sub = lapply(sim_100_starts_sub, function(y) {
  lapply(y, function(x) {
    lapply(x, function(list_of_summs) {
      if (!is.null(list_of_summs)) {
        print(x$iter)
        all_sum = lapply(list_of_summs[[2]], makeMSstatsTMTInput)
        all_sum = list(ProteinLevelData = rbindlist(lapply(all_sum, function(z) z$ProteinLevelData)),
                       FeatureLevelData = rbindlist(lapply(all_sum, function(z) z$FeatureLevelData)))
        summs = list(makeMSstatsTMTInput(list_of_summs[[1]]) , 
                     all_sum,
                     makeMSstatsTMTInput(list_of_summs[[3]]),
                     makeMSstatsTMTInput(list_of_summs[[4]]))
        
        gcs_comp = lapply(summs, function(z) {
          out = MSstatsTMT::groupComparisonTMT(z, cm_g00, 
                                               use_log_file = FALSE)
          out = as.data.table(out$ComparisonResult)
          out
        })
        gcs_comp = rbindlist(lapply(seq_along(gcs_comp), function(i) {
          out = gcs_comp[[i]]
          out$Method = c("unique", "all", "shared", "mixed")[i]
          out
        }))
        
        gcs_comp[, TimeLabel := apply(stringr::str_extract_all(Label, "[0-9]+min", 
                                                               simplify = TRUE), 
                                      1, paste, sep = "-", collapse = "-")]
        gcs_comp$num_unique = list_of_summs$num_unique
        gcs_comp$subset_protein = list_of_summs$subset_protein
        gcs_comp$iter = list_of_summs$iter
        gcs_comp
      } else {
        NULL
      }
    })
  })
})

gcs_comp_dt_starts_sub = rbindlist(unlist(unlist(gc_res_raw_starts_sub, F, F), F, F))
gcs_comp_dt_starts_sub[, Time := stringr::str_extract(TimeLabel, "[0-9]+min")]
gcs_comp_dt_starts_sub[, Time := stringr::str_replace(Time, "min", "")]
gcs_comp_dt_starts_sub[, Time := as.numeric(Time)]
gcs_comp_dt_starts_sub = merge(gcs_comp_dt_starts_sub, gs, by = c("Protein", "TimeLabel", "Time"))

# Starting points with unique peptides for every protein: simulation study ----
# sim_100_starts = lapply(
#   1:3, 
#   function(num_unique) {
#     mclapply(seq_len(num_rep), function(i) {
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
#         summary_1_triple_shared_flat = getWeightedProteinSummary(summary_input_1_triple,
#                                                                  norm = "Huber",
#                                                                  norm_parameter = 1e-3,
#                                                                  initial_summary = "flat",
#                                                                  tolerance = 1e-2, max_iter = 100)
#         summary_1_triple_shared_mixed = getWeightedProteinSummary(summary_input_1_triple,
#                                                                   norm = "Huber",
#                                                                   norm_parameter = 1e-3,
#                                                                   initial_summary = "unique",
#                                                                   tolerance = 1e-2, max_iter = 100)
#         
#         list(summary_1_triple_shared_flat,
#              summary_1_triple_shared_mixed,
#              iter =  i,
#              num_unique = num_unique)
#       }, error = function(e) NULL)
#     }, mc.cores = 6)
#   })
# saveRDS(sim_100_starts, "./results/protein_degrader/sim_100_starts.RDS")
sim_100_starts = readRDS("./results/protein_degrader/sim_100_starts.RDS")

gc_res_raw_starts = lapply(sim_100_starts, function(y) {
  lapply(y, function(list_of_summs) {
    if (!is.null(list_of_summs)) {
      summs = lapply(list_of_summs[1:2], makeMSstatsTMTInput)
      
      gcs_comp = lapply(summs, function(summ) {
        out = MSstatsTMT::groupComparisonTMT(summ, cm_g00, 
                                             use_log_file = FALSE)
        out = as.data.table(out$ComparisonResult)
        out
      })
      gcs_comp = rbindlist(lapply(seq_along(gcs_comp), function(i) {
        out = gcs_comp[[i]]
        out$Start = c("flat", "unique")[i]
        out
      }))
      gcs_comp[, iter := list_of_summs$iter]
      gcs_comp[, num_unique := list_of_summs$num_unique]
      gcs_comp[, subset_protein := list_of_summs$subset_protein]
      gcs_comp[, TimeLabel := apply(stringr::str_extract_all(Label, "[0-9]+min", 
                                                             simplify = TRUE), 
                                    1, paste, sep = "-", collapse = "-")]
      gcs_comp
    } else {
      NULL
    }
  })
})

gcs_comp_dt_starts = rbindlist(unlist(gc_res_raw_starts, F, F))
gcs_comp_dt_starts[, Time := stringr::str_extract(TimeLabel, "[0-9]+min")]
gcs_comp_dt_starts[, Time := stringr::str_replace(Time, "min", "")]
gcs_comp_dt_starts[, Time := as.numeric(Time)]
gcs_comp_dt_starts = merge(gcs_comp_dt_starts, gs, by = c("Protein", "TimeLabel", "Time"))

# Visualization of the results ----
ggplot(gcs_comp_dt_starts_sub[, .(MSE = mean((log2FC - GoldStandard)^2),
                              Var = var(log2FC)),
                          by = c("Method", "iter", "TimeLabel", "Time", "num_unique")][Method %in% c("shared", "mixed")],
       aes(x = reorder(as.character(num_unique), num_unique), y = MSE, 
           fill = ifelse(Method == "shared", "flat", "mixed"))) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(~Time) + #paste("Time =", reorder(as.character(Time), Time))
  scale_fill_manual(name = "starting point", values = c("#E69F00", "#56B4E9")) +
  scale_color_manual(name = "starting point", values = c("#E69F00", "#56B4E9")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylim(c(0, 0.2)) +
  xlab("number of unique peptides") +
  ylab("MSE of log-fold change estimation") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18), 
        legend.box = "vertical",
        legend.key.spacing = unit(5, "pt"),
        legend.spacing = unit(5, "pt"),
        legend.text.position = "right")
ggsave("./plots/pdf/start_points_brd.pdf", device = "pdf", width = 10, height = 5, units = "in", scale = 1, dpi = 300)
ggsave("./plots/png/start_points_brd.png", device = "png", width = 10, height = 5, units = "in", scale = 1, dpi = 300)

ggplot(gcs_comp_dt_starts[, .(mse = mean((log2FC - GoldStandard)^2)),
                          by = c("iter", "Start", "num_unique")],
       aes(x = as.character(num_unique), fill = Start,
           y = mse)) +
  geom_boxplot() +
  # scale_fill_manual(values = colors,
  #                   labels = c(unique = "selected unique",
  #                              shared = "proposed",
  #                              all = "all features")) +
  theme_bw() +
  scale_fill_manual(name = "starting point", values = c("#E69F00", "#56B4E9")) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16), 
        legend.box = "vertical") +
  ylab("MSE of log-fold change estimation") +
  xlab("number of unique peptides")
ggsave("./plots/pdf/starting_points.pdf", device = "pdf", width = 10, height =  5, scale = 1, units = "in", dpi = 300)
ggsave("./plots/png/starting_points.png", device = "png", width = 10, height =  5, scale = 1, units = "in", dpi = 300)

ggplot(gcs_comp_dt_starts,
       aes(x = reorder(as.character(Time), Time),
           y = log2FC, fill = Start)) + # , fill = Method
  geom_boxplot(linewidth = 1.2) + #, outlier.shape = NA
  geom_line(aes(x = reorder(as.character(Time), Time),
                y = GoldStandard, group = "1"),
            inherit.aes = FALSE, color = "darkblue", size = 1.5) +
  facet_grid(num_unique ~ Protein) +
  scale_fill_manual(name = "starting point", values = c("#E69F00", "#56B4E9")) +
  geom_point(aes(x = reorder(as.character(Time), Time),
                 y = GoldStandard),
             inherit.aes = FALSE, color = "darkblue", size = 1.5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16), 
        legend.box = "vertical",
        legend.key.spacing = unit(5, "pt"),
        legend.spacing = unit(5, "pt"),
        legend.text.position = "right") +
  xlab("time") +
  ylab("log-fold change") +
  ylim(c(-1.5, 0.5))
ggsave("./plots/pdf/starting_points_estimated.pdf", device = "pdf", width = 10, height =  5, scale = 1, units = "in", dpi = 300)
ggsave("./plots/png/starting_points_estimated.png", device = "png", width = 10, height =  5, scale = 1, units = "in", dpi = 300)

# Quantication of proteins identified only by shared peptides: result when BRD4 was the subset protein ----
colors = c("lightblue", "red", "purple")
colors = rev(colors)
ggplot(gcs_comp_dt_starts_sub[subset_protein == "BRD4"][Method != "mixed"][num_unique == 2],
       aes(x = reorder(as.character(Time), Time),
           y = log2FC, color = Method)) + # , fill = Method
  geom_boxplot(linewidth = 1.2) + #, outlier.shape = NA
  geom_line(aes(x = reorder(as.character(Time), Time),
                y = GoldStandard, group = "1"),
            inherit.aes = FALSE, color = "darkblue", size = 1.5) +
  facet_grid( ~ Protein) + # num_unique
  # scale_fill_manual(values = colors,
  #                   labels = c(unique = "unique",
  #                              shared = "proposed",
  #                              all = "all")) +
  scale_color_manual(values = colors,
                     labels = c(unique = "unique-only",
                                shared = "proposed",
                                all = "all-peptides")) +
  geom_point(aes(x = reorder(as.character(Time), Time),
                 y = GoldStandard),
             inherit.aes = FALSE, color = "darkblue", size = 1.2) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16), 
        legend.box = "vertical") +
  xlab("time") +
  ylab("log-fold change") +
  ylim(c(-1.5, 0.5))
ggsave("./plots/pdf/brd_subset_log2fcs_2uni.pdf", width = 10, height = 5, units = "in", scale = 1, dpi = 300)
ggsave("./plots/png/brd_subset_log2fcs_2uni.png", width = 10, height = 5, units = "in", scale = 1, dpi = 300)

# Influence of the loss function: simulation study ----
num_rep = 100

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

# sim_100_noise_loss_h_l2 = mclapply(
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
#         summary_1_triple_l2 = getWeightedProteinSummary(summary_input_1_triple,
#                                                             norm = "p_norm",
#                                                             norm_parameter = 2,
#                                                             initial_summary = "unique",
#                                                             tolerance = 1e-2, max_iter = 100)
# 
#         summary_1_triple_huber = getWeightedProteinSummary(summary_input_1_triple,
#                                                             norm = "Huber",
#                                                             norm_parameter = 1e-3,
#                                                             initial_summary = "unique",
#                                                             tolerance = 1e-2, max_iter = 100)
# 
#         list(summary_1_triple_l2,
#              summary_1_triple_huber)
#       }, error = function(e) NULL)
#     })
#   },
#   mc.cores = 5)
# saveRDS(sim_100_noise_loss_h_l2, "results/protein_degrader/sim_100_noise_loss_h_l2.RDS")
sim_100_noise_loss_h_l2 = readRDS("results/protein_degrader/sim_100_noise_loss_h_l2.RDS")

gc_res_raw_loss_h_l2 = lapply(sim_100_noise_loss_h_l2, function(y) {
  lapply(y, function(x) {
    if (!is.null(x)) {
      summs = lapply(x, makeMSstatsTMTInput)
      
      gcs_comp = lapply(summs, function(x) {
        out = MSstatsTMT::groupComparisonTMT(x, cm_g00, 
                                             use_log_file = FALSE)
        out = as.data.table(out$ComparisonResult)
        out
      })
      gcs_comp = rbindlist(lapply(seq_along(gcs_comp), function(i) {
        out = gcs_comp[[i]]
        out$Method = c("L2", "Huber")[i]
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

gcs_comp_dt_noise_loss_l2_h = rbindlist(lapply(seq_along(gc_res_raw_loss_h_l2), function(i) {
  x = gc_res_raw_loss_h_l2[[i]]
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
gcs_comp_dt_noise_loss_l2_h[, Time := stringr::str_extract(TimeLabel, "[0-9]+min")]
gcs_comp_dt_noise_loss_l2_h[, Time := stringr::str_replace(Time, "min", "")]
gcs_comp_dt_noise_loss_l2_h[, Time := as.numeric(Time)]
gcs_comp_dt_noise_loss_l2_h = merge(gcs_comp_dt_noise_loss_l2_h, gs, 
                                    by = c("Protein", "TimeLabel", "Time"))

### Visualization of loss function-related results ----
colors = c("#E69F00", "#56B4E9")
# colors = rev(colors)

ggplot(gcs_comp_dt_noise_loss_l2_h[num_unique %in% c(1, 2, 3, 5), .(MSE = mean((log2FC - GoldStandard)^2, na.rm = T)),
                                   by = c("Method", "num_unique", "iter")],
       aes(x = reorder(as.character(num_unique), num_unique), y = MSE, fill = Method)) + # reorder(as.character(num_unique), num_unique)
  geom_boxplot() +
  # facet_grid(, scales = "free_y") + # , nrow = 6, scales = "free_y"
  scale_fill_manual(
    name = "method",
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
        strip.text = element_text(size = 16), 
        legend.box = "vertical") +
  xlab("number of unique peptides per protein") +
  ylab("MSE of log-fold change estimation")
ggsave("./plots/pdf/sim_brd_loss.pdf", device = "pdf", width = 10, height = 5, scale = 1, dpi = 300)
ggsave("./plots/png/sim_brd_loss.png", device = "png", width = 10, height = 5, scale = 1, dpi = 300)
