# Code required to reproduces Figures:
# [main] 6, 8, 9 [log-fold change estimation and statistical inference evaluation on simulated data]
# [supplemental] 1 [example profile plots of simulated data]
# Libraries ----
library(data.table)
library(SimulateTMT)
library(MSstatsTMT)
library(MSstatsWeightedSummary)
library(ggplot2)
library(parallel)
# Useful functions ----
get_log2fc_groups = function(gc_base, log2_fcs, cm, num_proteins) {
  labels = unique(gc_base$ComparisonResult$Label)
  log2_groups = rbindlist(lapply(seq_along(log2_fcs), function(i) {
    data.table(Protein = paste("Prot", stringr::str_pad(i, pad = "0", width = 4, side = "left"), sep = "_"),
               Label = labels,
               log2Group = sapply(seq_len(nrow(cm)), function(j) {
                 log2_fcs[[i]] %*% cm[j, ]
               }))
  }))
  log2_groups = rbind(log2_groups,
                      rbindlist(lapply((length(log2_fcs) + 1):num_proteins, function(i) {
                        data.table(Protein = paste("Prot", stringr::str_pad(i, pad = "0", width = 4, side = "left"), sep = "_"),
                                   Label = labels,
                                   log2Group = 0)
                      })))  
  log2_groups
}
# Data simulation ----
combs = expand.grid(c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1))
combs = cbind(combs, sum = apply(combs, 1, sum))
combs = combs[combs$sum > 1, , drop = FALSE]
combs = combs[combs$sum == 2, ]
combs = combs[, -6, drop = FALSE]
combs = apply(combs, 2, as.logical, simplify = FALSE)
combs = as.data.frame(combs)

weights = list(c(0.7, 0.3), c(0.8, 0.2),
               c(0.5, 0.5), c(0.6, 0.4),
               c(0.9, 0.1), c(0.2, 0.8),
               c(0.3, 0.7), c(0.5, 0.5),
               c(0.4, 0.6), c(0.1, 0.9))

protein_groups = lapply(1, function(x) c(x, x + 1, x + 2, x + 3, x + 4))
protein_groups = lapply(protein_groups, function(x) stringr::str_pad(x, pad = "0", width = 4, side = "left"))
protein_groups = lapply(protein_groups, function(x) paste("Prot", x, sep = "_"))


log2_fcs = list(c(1, 0.9, 0.8, 0.5, 0.25),
                c(1, 0.5, 0.25, -0.25, -0.5),
                c(1, 0.8, 0.7, 0.3, 0)#,
                # c(1, 0.75, 0.7, 0.3, 0.25),
                # c(1, 0.7, 0.6, 0.5, 0.3)
)
tmt_design = create_TMT_design(5, length(log2_fcs), 1, 5, 2, "groupComparison")
set.seed(100)
tmt_simulated = simulate_log_abundances(tmt_design, 15, 
                                        log2_fcs, 
                                        0, list("0.8" = 0.01, "0.6" = 0.01,
                                                "0.75" = 0.01,
                                                "0.55" = 0.03,
                                                "0.5" = 0.03,
                                                "0" = 0.01,
                                                "0.45" = 0.01,
                                                "0.35" = 0.01,
                                                "0.25" = 0.05,
                                                "0.3" = 0.05,
                                                "-0.25" = 0.0,
                                                "-0.75" = 0.01,
                                                "-0.5" = 0.01,
                                                "0.9" = 0.01,
                                                "0.7" = 0.01,
                                                "1" = 0.01,
                                                "-1.25" = 0.03,
                                                "-1.5" = 0.03,
                                                "-1.75" = 0.03,
                                                "-2" = 0.03), 0.01, "constant", 0.01)
# Example profile plots
ggplot(tmt_simulated, aes(x = reorder(as.character(Channel), as.numeric(as.character(Channel))), y = Abundance, 
                          group = Protein, color = stringr::str_replace(Protein, "Prot_000", " "))) +
  geom_line(size = 1) +
  geom_point(aes(shape = Condition), size = 2) +
  theme_bw() +
  xlab("channel") +
  ylab("abundance") +
  scale_color_discrete(name = "protein") +
  scale_shape_discrete(name = "condition") +
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18), 
        legend.position = "bottom",
        legend.box = "vertical")
ggsave("plots/pdf/plot_sim_example_prots.pdf", device = "pdf", width = 10, height = 5, units = "in", scale = 1, dpi = 300)
ggsave("plots/png/plot_sim_example_prots.png", device = "png", width = 10, height = 5, units = "in", scale = 1, dpi = 300)

conditions = levels(tmt_simulated$Condition)
cm = cbind(matrix(-1, nrow = 4, ncol = 1), diag(1, 4))
colnames(cm) = conditions
row.names(cm) = c("1 vs 2", "1 vs 3", "1 vs 4", "1 vs 5")

gc_base = groupComparisonTMT(list(FeatureLevelData = data.frame(),
                                  ProteinLevelData = tmt_simulated),
                             cm, use_log_file = FALSE)
log2fc_groups = get_log2fc_groups(gc_base, log2_fcs, cm, num_proteins = 5)
# Simulation study ----
# intermediate results available on request due to data size
# for (num_biorep in 1:3) {
#   tmt_design = create_TMT_design(5, length(log2_fcs), 1, 5, num_biorep, "groupComparison")
#   tmt_design[, Channel := factor(as.character(Channel),
#                                  levels = as.character(sort(as.numeric(unique(as.character(Channel))))),
#                                  ordered = T)]
#   tmt_design[, Channel := as.character(Channel)]
#   
#   for (num_unique in c(1:3, 5)) {
#     for (num_shared in c(3, 5, 10)) {
#       for (error_sd in c(0.1, 0.2)) {
#         single_case = mclapply(1:50, function(i) {
#           tmt_simulated = simulate_log_abundances(tmt_design, 15,
#                                                   log2_fcs,
#                                                   0, list("0.8" = 0.01, "0.6" = 0.01,
#                                                           "0.75" = 0.01,
#                                                           "0.55" = 0.03,
#                                                           "0.5" = 0.03,
#                                                           "0" = 0.01,
#                                                           "0.45" = 0.01,
#                                                           "0.35" = 0.01,
#                                                           "0.25" = 0.05,
#                                                           "0.3" = 0.05,
#                                                           "-0.25" = 0.0,
#                                                           "-0.75" = 0.01,
#                                                           "-0.5" = 0.01,
#                                                           "0.9" = 0.01,
#                                                           "0.7" = 0.01,
#                                                           "1" = 0.01,
#                                                           "-1.25" = 0.03,
#                                                           "-1.5" = 0.03,
#                                                           "-1.75" = 0.03,
#                                                           "-2" = 0.03), 0.01, "constant", 0.01)
#           tmt_simulated[, AbundanceStandard := Abundance - mean(Abundance),
#                         by = "Protein"]
#           dataset = SimulateTMT::simulate_dataset(tmt_simulated, protein_groups, combs,
#                                                   num_unique, rep(num_shared, 10), weights, 0.02, error_sd, c(-1, 1), 1)
#           dataset[, log2IntensityNormalized := log(Intensity, 2)]
#           
#           summ_shared = getWeightedProteinSummary(dataset, "Huber", 1e-3, max_iter = 50)
#           summ_unique = getWeightedProteinSummary(dataset[(IsUnique)], "Huber", 1e-3, max_iter = 50)
#           summaries_all = lapply(split(dataset, dataset$ProteinName),
#                                  function(x) {
#                                    getWeightedProteinSummary(x, norm = "Huber",
#                                                              norm_parameter = 1e-3, tolerance = 1e-3, max_iter = 10, weights_mode = "contributions",
#                                                              save_weights_history = F, save_convergence_history = F)@ProteinLevelData
#                                  })
#           summaries_all = rbindlist(summaries_all)
#           
#           gc_shared = MSstatsTMT::groupComparisonTMT(makeMSstatsTMTInput(summ_shared), cm, use_log_file = F)
#           gc_unique = MSstatsTMT::groupComparisonTMT(makeMSstatsTMTInput(summ_unique), cm, use_log_file = F)
#           gc_all = MSstatsTMT::groupComparisonTMT(list(ProteinLevelData = summaries_all,
#                                                        FeatureLevelData = data.table()),
#                                                   cm, use_log_file = F)
#           
#           comp = merge(as.data.table(rbind(
#             cbind(gc_shared$ComparisonResult, Method = "shared"),
#             cbind(gc_unique$ComparisonResult, Method = "unique"),#
#             cbind(gc_all$ComparisonResult, Method = "all")
#           )), log2fc_groups, by = c("Protein", "Label"))
#           
#           list(summaries = list(summ_shared, summ_unique, summaries_all),
#                gcs = list(gc_shared, gc_unique, gc_all),
#                comparison = comp,
#                iter = i,
#                num_unique = num_unique,
#                num_biorep = num_biorep,
#                num_shared = num_shared,
#                error_sd = error_sd,
#                ground_truth = tmt_simulated,
#                feature_level = dataset
#           )
#         }, mc.cores = 5)
#         saveRDS(single_case, paste("one_case", num_shared, num_unique, num_biorep, error_sd, ".RDS", sep = "_"))
#         rm(single_case)
#         gc()
#         # single_case
#       }
#     }
#   }
# }
# raw_res_files = list.files()
# raw_res = lapply(raw_res_files, readRDS)
# res_meta_dt = rbindlist(lapply(raw_res, function(x) {
#   rbindlist(lapply(x, function(x_iter) {
#     if (!is.null(x_iter) & !inherits(x_iter, "try-error")) {
#       gc_comp = x_iter[['comparison']]
#       gc_comp[["iter"]] = x_iter[["iter"]]
#       gc_comp[["num_unique"]] = x_iter[["num_unique"]]
#       gc_comp[["num_biorep"]] = x_iter[["num_biorep"]]
#       gc_comp[["num_shared"]] = x_iter[["num_shared"]]
#       gc_comp[["error_sd"]] = x_iter[["error_sd"]]
#       gc_comp[["prop_unique_cl"]] = uniqueN(x_iter[["feature_level"]][(IsUnique), PSM]) / uniqueN(x_iter[["feature_level"]][, PSM])
#       gc_comp[["prop_unique"]] = x_iter[["num_unique"]] / (x_iter[["num_shared"]] + x_iter[["num_unique"]])
#       gc_comp
#     }
#   }))
# }))
# res_meta_dt[, PropCatCl := cut(100 * prop_unique_cl, 100 * c(0.0, 0.2, 0.4, 0.5))]
# saveRDS(res_meta_dt, "./results/simulated_data/res_meta_dt.RDS")
res_meta_dt = readRDS("./results/simulated_data/res_meta_dt.RDS")

# Visualization of the results ----
iter_mse = res_meta_dt[, .(IterMSE = mean((log2FC - log2Group) ^ 2)),
                       by = c("Method", "iter", "num_unique", "num_biorep", "num_shared", "error_sd", "PropCatCl")]
res_meta_dt2 = copy(res_meta_dt)
res_meta_dt2[, Sig := ifelse(log2Group == "0", "no difference", "difference")]
res_meta_dt2[, Method := ifelse(Method == "shared", "proposed", Method)]
iter_mse2 = res_meta_dt2[, .(IterMSE = mean((log2FC - log2Group) ^ 2)),
                         by = c("Method", "iter", "num_biorep", "error_sd", "PropCatCl")]

ggplot(res_meta_dt2[error_sd == 0.2 & num_biorep == 2 & num_unique == 2],
       aes(x = reorder(as.character(abs(log2Group)), abs(log2Group)),
           y = abs(log2FC), fill = as.character(Method))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(x = reorder(as.character(abs(log2Group)), abs(log2Group)),
                 y = abs(log2Group)),
             data = unique(res_meta_dt2[, .(log2Group = abs(log2Group))]),
             color = "white", inherit.aes = F, size = 2) +
  # facet_grid(~num_biorep) +
  scale_fill_manual(name = "method", values = c("purple", "red", "lightblue"),
                    labels = c("all simulated peptides in each protein",
                               "proposed", "unique simulated peptides in each protein")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("true log-fold change") +
  ylab("estimated log-fold change") +
  # facet_grid(paste("SD:", error_sd)~paste("#biorep:", num_biorep)) + # paste("std.dev:", error_sd)
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 18), 
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.box = "vertical") +
  ylim(c(-0.1, 1.6))
ggsave(file = "plots/pdf/plot_sim_5prot_log2fc.pdf", device = "pdf", width = 10, height = 5, scale = 1, units = "in", dpi = 300)
ggsave(file = "plots/png/plot_sim_5prot_log2fc.png", device = "png", width = 10, height = 5, scale = 1, units = "in", dpi = 300)

tnr_data = res_meta_dt[, .(tnr = sum(pvalue >= 0.05) / (sum(pvalue < 0.05 & log2Group == 0) + sum(pvalue >= 0.05 & log2Group == 0)), 
                           n = .N),
                       by = c("Method", "log2Group", "num_unique", "num_shared", "num_biorep", "error_sd", "PropCatCl")]
ggplot(tnr_data[num_biorep > 1],
       aes(x = PropCatCl, y = tnr, fill = Method)) +
  geom_boxplot() +
  facet_grid(paste("SD:", error_sd)~paste("no. biorep:", num_biorep)) + # paste("std.dev:", error_sd)
  scale_fill_manual(name = "method", values = c("purple", "red", "lightblue"),
                    labels = c("all simulated peptides in each protein",
                               "proposed", "unique simulated peptides in each protein")) +
  theme_bw() +
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 18), 
        legend.position = "bottom",
        # legend.direction = "vertical",
        legend.box = "vertical") +
  xlab("proportion of unique peptides in a cluster (%)") +
  ylab("specificity")
ggsave("plots/pdf/plot_sim_spec.pdf", device = "pdf", width = 10, height = 5, units = "in", scale = 1, dpi = 300)
ggsave("plots/png/plot_sim_spec.png", device = "png", width = 10, height = 5, units = "in", scale = 1, dpi = 300)

power_data_per_iter = res_meta_dt[, .(power = sum(adj.pvalue < 0.05 & log2Group != 0) / (sum(log2Group != 0))),
                                  by = c("Method", "iter", "num_unique", "num_biorep", "num_shared", "error_sd", "PropCatCl")]
ggplot(power_data_per_iter[num_biorep > 1],
       aes(x = PropCatCl, y = power, fill = ifelse(Method == "shared", "proposed", Method))) +
  geom_boxplot() +
  facet_grid(paste("SD:", error_sd)~paste("no. biorep:", num_biorep)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("proportion of unique peptides in a cluster (%)") +
  scale_fill_manual(name = "method", values = c("purple", "red", "lightblue"),
                    labels = c("all simulated peptides in each protein",
                               "proposed", "unique simulated peptides in each protein")) +
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 18), 
        legend.position = "bottom",
        legend.box = "vertical")
ggsave("plots/pdf/plot_sim_power_iter.pdf", device = "pdf", width = 10, height = 5, units = "in", scale = 1, dpi = 300)
ggsave("plots/png/plot_sim_power_iter.png", device = "png", width = 10, height = 5, units = "in", scale = 1, dpi = 300)

ggplot(iter_mse[error_sd < 0.5],
       aes(x = PropCatCl, y = IterMSE, fill = ifelse(Method == "shared", "proposed", Method))) +
  geom_boxplot() +
  facet_grid(paste("SD:", error_sd)~paste("no. biorep:", num_biorep)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylim(c(0, 0.12)) +
  xlab("proportion of unique peptides in a cluster (%)") +
  ylab("MSE of log-fold change estimation") +
  scale_fill_manual(name = "method", values = c("purple", "red", "lightblue")) +
  theme(legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 18),
        legend.position = "bottom",
        legend.box = "vertical")
ggsave("plots/pdf/plot_sim_mse.pdf", device = "pdf", width = 10, height = 5, units = "in", scale = 1, dpi = 300)
ggsave("plots/png/plot_sim_mse.png", device = "png", width = 10, height = 5, units = "in", scale = 1, dpi = 300)
