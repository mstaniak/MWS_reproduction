library(data.table)
library(SimulateTMT)
library(MSstatsTMT)
library(MSstatsWeightedSummary)
library(ggplot2)
library(parallel)

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
                c(1, 0.9, 0.8, 0.8, 0.7),
                c(1, 0.8, 0.7, 0.3, 0),
                # c(1, 0.75, 0.7, 0.3, 0.25),
                c(1, 0.7, 0.6, 0.5, 0.3))
tmt_design = create_TMT_design(5, length(log2_fcs), 1, 5, 2, "groupComparison")
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
                                                "1" = 0.01), 0.01, "constant", 0.01)

ggplot(tmt_simulated, aes(x = as.numeric(as.character(Channel)), y = Abundance, 
                          group = Protein, color = Protein)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom")


conditions = levels(tmt_simulated$Condition)
cm = cbind(matrix(-1, nrow = 4, ncol = 1), diag(1, 4))
colnames(cm) = conditions
row.names(cm) = c("1 vs 2", "1 vs 3", "1 vs 4", "1 vs 5")

gc_base = groupComparisonTMT(list(FeatureLevelData = data.frame(),
                                  ProteinLevelData = tmt_simulated),
                             cm, use_log_file = FALSE)
log2fc_groups = get_log2fc_groups(gc_base, log2_fcs, cm, num_proteins = 5)

# intermediate results available on request
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
#                                                           "0.5" = 0.03,
#                                                           "0.55" = 0.03,
#                                                           "0" = 0.01,
#                                                           "0.45" = 0.01,
#                                                           "0.4" = 0.01,
#                                                           "0.35" = 0.01,
#                                                           "0.3" = 0.01,
#                                                           "0.25" = 0.05,
#                                                           "-0.25" = 0.0,
#                                                           "-0.75" = 0.01,
#                                                           "-0.5" = 0.01,
#                                                           "0.9" = 0.01,
#                                                           "0.7" = 0.01,
#                                                           "1" = 0.01), 0.01, "constant", 0.01)
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
#         }, mc.cores = 10)
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
#       gc_comp[["prop_unique"]] = x_iter[["num_unique"]] / x_iter[["num_shared"]]
#       gc_comp
#     }
#   }))
# }))
# length(raw_res[[1]])
# names(raw_res[[1]][[1]])
# sapply(raw_res, length)
# res_meta_dt
# saveRDS(res_meta_dt, "../20240321_res_meta_dt.RDS")
res_meta_dt = readRDS("results/simulated_data/sim_gcs.RDS")

iter_mse = res_meta_dt[, .(IterMSE = mean((log2FC - log2Group) ^ 2)),
                       by = c("Method", "iter", "num_unique", "num_biorep", "num_shared", "error_sd", "prop_unique_cl", "prop_unique")]


ggplot(iter_mse[error_sd < 0.5],
       aes(x = reorder(as.character(round(prop_unique_cl, 2)), prop_unique_cl), y = IterMSE, fill = ifelse(Method == "shared", "proposed", Method))) +
  geom_boxplot() +
  facet_grid(paste("#biorep:", num_biorep)~paste("std.dev:", error_sd)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylim(c(0, 0.12)) +
  xlab("proportion of unique peptides") +
  ylab("MSE of log-fold change estimation") +
  scale_fill_manual(name = "method", values = c("purple", "red", "lightblue")) +
  theme(legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 14), 
        legend.position = "bottom",
        legend.box = "vertical")
ggsave("plot_sim_mse.pdf", device = "pdf", width = 10, height = 5, units = "in", scale = 1)

iter_mse_prot = res_meta_dt[, .(IterMSE = mean((log2FC - log2Group) ^ 2)),
                            by = c("Protein", "Method", "iter", "num_unique", "num_biorep", "num_shared", "error_sd", "prop_unique_cl", "prop_unique")]
ggplot(iter_mse_prot[error_sd < 0.5 & num_biorep == 2],
       aes(x = reorder(as.character(round(prop_unique_cl, 2)), prop_unique_cl), y = IterMSE, fill = Method)) +
  geom_boxplot() +
  facet_grid(Protein~error_sd) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylim(c(0, 0.15)) +
  xlab("proportion of unique peptides") +
  scale_fill_manual(name = "method", values = c("purple", "red", "lightblue")) +
  theme(legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 14), 
        legend.position = "bottom",
        legend.box = "vertical")

ppv_data = res_meta_dt[num_biorep > 1, .(ppv = sum(adj.pvalue < 0.05 & log2Group != 0) / (sum(adj.pvalue < 0.05))),
                       by = c("Method", "iter", "num_unique", "num_biorep", "num_shared", "error_sd", "prop_unique_cl", "prop_unique")]
ggplot(ppv_data[],
       aes(x = reorder(as.character(num_unique), num_unique), y = ppv, fill = Method)) +
  geom_boxplot() +
  facet_grid(num_biorep~error_sd) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylim(c(0.6, 1))

power_data = res_meta_dt[, .(power = sum(adj.pvalue < 0.05 & log2Group != 0) / (sum(log2Group != 0))),
                         by = c("Method", "iter", "num_unique", "num_biorep", "num_shared", "error_sd", "prop_unique_cl", "prop_unique")]
ggplot(power_data[num_biorep > 1],
       aes(x = reorder(as.character(round(prop_unique_cl, 2)), prop_unique_cl), y = power, fill = ifelse(Method == "shared", "proposed", Method))) +
  geom_boxplot() +
  facet_grid(paste("#biorep:", num_biorep)~paste("std.dev:", error_sd)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("proportion of unique peptides") +
  scale_fill_manual(name = "method", values = c("purple", "red", "lightblue")) +
  theme(legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 14), 
        legend.position = "bottom",
        legend.box = "vertical")
# ggsave("plot_sim_power.pdf", device = "pdf", width = 10, height = 5, units = "in", scale = 1)


fdr_data = res_meta_dt[, 
                       .(fdr = sum(adj.pvalue < 0.05 & log2Group == 0) / (sum(adj.pvalue < 0.05))),
                       by = c("Method", "iter", "num_unique", "num_biorep", "num_shared", "error_sd", "prop_unique_cl", "prop_unique")]
ggplot(fdr_data[num_biorep > 1],
       aes(x = reorder(as.character(round(prop_unique_cl, 2)), prop_unique_cl), y = fdr, fill = Method)) +
  geom_boxplot() +
  facet_grid(paste("#biorep:", num_biorep)~paste("std.dev:", error_sd)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("proportion of unique peptides") +
  ylim(c(0, 0.25))

fd_data = res_meta_dt[Protein == "Prot_0005", 
                       .(FD = sum(adj.pvalue < 0.05 & log2Group == 0),
                         n = .N),
                       by = c("Method", "iter", "num_unique", "num_biorep", "num_shared", "error_sd")]
ggplot(fd_data[num_biorep > 1 & num_shared == 10],
       aes(x = reorder(as.character(num_unique), num_unique), y = FD, fill = ifelse(Method == "shared", "proposed", Method))) +
  # geom_bar(stat = "identity", position = "dodge") +
  geom_boxplot() +
  facet_grid(paste("#biorep:", num_biorep)~paste("std.dev:", error_sd)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("proportion of unique peptides") +
  # scale_y_continuous(labels = scales::percent) +
  ylab("% of FD across replicates") +
  scale_fill_manual(name = "method", values = c("purple", "red", "lightblue")) +
  theme(legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 14), 
        legend.position = "bottom",
        legend.box = "vertical")
# ggsave("plot_sim_fd_iter.pdf", device = "pdf", width = 10, height = 5, units = "in", scale = 1)


fd_data = res_meta_dt[Protein == "Prot_0005", 
                      .(FWER = mean(any(adj.pvalue < 0.05 & log2Group == 0 )),
                        n = .N),
                      by = c("Method", "num_unique", "num_biorep", "num_shared", "error_sd")]
ggplot(fd_data[num_biorep > 1 & num_shared == 10],
       aes(x = reorder(as.character(num_unique), num_unique), y = FWER, fill = ifelse(Method == "shared", "proposed", Method))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(paste("#biorep:", num_biorep)~paste("std.dev:", error_sd)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("proportion of unique peptides") +
  # scale_y_continuous(labels = scales::percent) +
  ylab("% of FD across replicates") +
  scale_fill_manual(name = "method", values = c("purple", "red", "lightblue")) +
  theme(legend.title = element_text(size = 24),
        axis.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 22), 
        legend.position = "bottom",
        legend.box = "vertical")


ggplot(fd_data[num_biorep == 2 & Method == "unique"],
       aes(x = reorder(as.character(num_shared), num_shared), y = 100 * FD / n, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(num_unique~paste("std.dev:", error_sd)) +
  theme_bw() +
  theme(legend.position = "bottom")


res_meta_dt2 = copy(res_meta_dt)
res_meta_dt2[, Sig := ifelse(log2Group == "0", "no difference", "difference")]
res_meta_dt2[, Method := ifelse(Method == "shared", "proposed", Method)]
iter_mse2 = res_meta_dt2[, .(IterMSE = mean((log2FC - log2Group) ^ 2)),
                         by = c("Method", "iter", "num_unique", "num_biorep", "num_shared", "error_sd", "prop_unique_cl", "prop_unique")]

ggplot(res_meta_dt2[error_sd < 0.5 & num_unique == 2 & num_shared == 5 & num_biorep == 2],
       aes(x = reorder(as.character(abs(log2Group)), abs(log2Group)),
           y = abs(log2FC), fill = as.character(Method))) +
  geom_boxplot() +
  geom_point(aes(x = reorder(as.character(abs(log2Group)), abs(log2Group)),
                 y = abs(log2Group)),
             data = unique(res_meta_dt2[, .(log2Group = abs(log2Group))]),
             color = "white", inherit.aes = F, size = 2) +
  # facet_grid(~num_biorep) +
  scale_fill_manual(name = "method", values = c("purple", "red", "lightblue")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("true log-fold change") +
  ylab("estimated log-fold change") +
  theme(legend.title = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 20), 
        legend.position = "bottom",
        legend.box = "vertical") +
  ylim(c(0, 1.5))
