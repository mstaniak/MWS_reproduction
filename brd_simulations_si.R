library(data.table)
library(MSstatsWeightedSummary)
library(ggplot2)
library(gridExtra)
library(parallel)
library(MSstatsConvert)

brd_cluster = readRDS("processed_data/protein_degrader/brd_cluster.RDS")
gs = readRDS("processed_data/protein_degrader/gold_standard_tbl.RDS")
annot = readRDS("processed_data/protein_degrader/brd_annot_tbl.RDS")

cm_g00 = fread("./input_data/protein_degrader/OBJ037605_G0011_comparison_matrix.txt")
group_labels = cm_g00$V1
cm_g00[, V1 := NULL]
cm_g00 = as.matrix(cm_g00)
row.names(cm_g00) = group_labels

# Starting point influence with subset proteins

num_rep = 100

sim_100_starts = lapply(
  1:3 , 
  function(num_unique) {
    mclapply(seq_len(num_rep), function(i) {
      print(i)
      lapply(c("BRD2", "BRD3", "BRD4"), function(subset_protein) { #
        tryCatch({
          set.seed(i)
          unique_peptides = unlist(lapply(setdiff(unique(brd_cluster$ProteinName), subset_protein),
                                          function(x) {
                                            sample(brd_cluster[(IsUnique)][ProteinName == x, unique(PSM)], num_unique)
                                          }), F, F)
          summary_input_1_triple = brd_cluster[!(IsUnique) | PSM %in% unique_peptides]
          summary_input_1_triple[, Channel := as.factor(as.character(Channel))]
          
          summary_1_triple_unique = getWeightedProteinSummary(summary_input_1_triple[(IsUnique)],
                                                              norm = "Huber",
                                                              norm_parameter = 1e-3,
                                                              initial_summary = "unique",
                                                              tolerance = 1e-2, max_iter = 100)
          
          summary_1_all_triple = lapply(split(summary_input_1_triple, summary_input_1_triple[, ProteinName]),
                                        function(x) {
                                          x$IsUnique = TRUE
                                          output = getWeightedProteinSummary(x,
                                                                             norm = "Huber",
                                                                             norm_parameter = 1e-3,
                                                                             initial_summary = "unique",
                                                                             tolerance = 1e-2, max_iter = 100)
                                          output
                                        })
          
          summary_1_triple_shared = getWeightedProteinSummary(summary_input_1_triple,
                                                              norm = "Huber",
                                                              norm_parameter = 1e-3,
                                                              initial_summary = "flat",
                                                              tolerance = 1e-2, max_iter = 100)
          summary_1_triple_shared_mixed = getWeightedProteinSummary(summary_input_1_triple,
                                                                    norm = "Huber",
                                                                    norm_parameter = 1e-3,
                                                                    initial_summary = "flat shared",
                                                                    tolerance = 1e-2, max_iter = 100)
          
          list(unique = summary_1_triple_unique,
               all = summary_1_all_triple,
               proposed = summary_1_triple_shared,
               proposed_mixed = summary_1_triple_shared_mixed,
               iter =  i,
               num_unique = num_unique,
               subset_protein = subset_protein)
        }, error = function(e) NULL)
      })
    }, mc.cores = 5)
  })

gc_res_raw_starts = lapply(sim_100_starts, function(y) {
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
        gcs_comp
      } else {
        NULL
      }
    })
  })
})

gcs_comp_dt_starts = rbindlist(unlist(unlist(gc_res_raw_starts, F, F), F, F))
gcs_comp_dt_starts[, Time := stringr::str_extract(TimeLabel, "[0-9]+min")]
gcs_comp_dt_starts[, Time := stringr::str_replace(Time, "min", "")]
gcs_comp_dt_starts[, Time := as.numeric(Time)]
gcs_comp_dt_starts = merge(gcs_comp_dt_starts, gs, by = c("Protein", "TimeLabel", "Time"))

ggplot(gcs_comp_dt_starts[Protein == subset_protein][Method %in% c("mixed", "shared")],
       aes(x = reorder(as.character(Time), Time),
           y = log2FC, color = Method)) + # , fill = Method
  geom_boxplot(linewidth = 1.2) + #, outlier.shape = NA
  geom_line(aes(x = reorder(as.character(Time), Time),
                y = GoldStandard, group = "1"),
            inherit.aes = FALSE, color = "darkblue", size = 1.5) +
  facet_grid(num_unique ~ Protein) +
  scale_color_manual(name = "method", values = c("#E69F00", "#56B4E9"), labels = c("mixed", "flat")) +
  geom_point(aes(x = reorder(as.character(Time), Time),
                 y = GoldStandard),
             inherit.aes = FALSE, color = "darkblue", size = 1.2) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14), 
        legend.box = "vertical") +
  xlab("time") +
  ylab("log-fold change") +
  ylim(c(-1.5, 0.5))
# ggsave("brd_sim_subset.pdf", device = "pdf", width = 10, height = 5, scale = 1, units = "in")

# Starting points with unique peptides for every protein

sim_100_starts = lapply(
  1:3, 
  function(num_unique) {
    mclapply(seq_len(num_rep), function(i) {
      print(i)
      tryCatch({
        set.seed(i)
        unique_peptides = unlist(lapply(unique(brd_cluster$ProteinName),
                                        function(x) {
                                          sample(brd_cluster[(IsUnique)][ProteinName == x, unique(PSM)], num_unique)
                                        }), F, F)
        summary_input_1_triple = brd_cluster[!(IsUnique) | PSM %in% unique_peptides]
        summary_input_1_triple[, Channel := as.factor(as.character(Channel))]
        
        summary_1_triple_shared_flat = getWeightedProteinSummary(summary_input_1_triple,
                                                                 norm = "Huber",
                                                                 norm_parameter = 1e-3,
                                                                 initial_summary = "flat",
                                                                 tolerance = 1e-2, max_iter = 100)
        summary_1_triple_shared_mixed = getWeightedProteinSummary(summary_input_1_triple,
                                                                  norm = "Huber",
                                                                  norm_parameter = 1e-3,
                                                                  initial_summary = "unique",
                                                                  tolerance = 1e-2, max_iter = 100)
        
        list(summary_1_triple_shared_flat,
             summary_1_triple_shared_mixed,
             iter =  i,
             num_unique = num_unique)
      }, error = function(e) NULL)
    }, mc.cores = 6)
  })


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

ggplot(gcs_comp_dt_starts[, .(MSE = mean((log2FC - GoldStandard)^2),
                              Var = var(log2FC)),
                          by = c("Start", "iter", "TimeLabel", "Time", "num_unique")],
       aes(x = reorder(as.character(num_unique), num_unique), y = MSE, 
           fill = Start)) +
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
        legend.title = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 18), 
        legend.box = "vertical",
        legend.key.spacing = unit(5, "pt"),
        legend.spacing = unit(5, "pt"),
        legend.text.position = "right")
# ggsave("start_points_brd.pdf", device = "pdf", width = 10, height = 5, units = "in", scale = 1)

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
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14), 
        legend.box = "vertical") +
  ylab("MSE of log-fold change estimation") +
  xlab("number of unique peptides")
# ggsave("starting_points.pdf", device = "pdf", width = 10, height =  5, scale = 1, units = "in")

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
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14), 
        legend.box = "vertical",
        legend.key.spacing = unit(5, "pt"),
        legend.spacing = unit(5, "pt"),
        legend.text.position = "right") +
  xlab("time") +
  ylab("log-fold change") +
  ylim(c(-1.5, 0.5))
# ggsave("starting_points_estimated.pdf", device = "pdf", width = 10, height =  5, scale = 1, units = "in")

# Quantication of proteins identified only by shared peptides: result when BRD4 was the subset protein
colors = c("lightblue", "red", "purple")
colors = rev(colors)
ggplot(gcs_comp_dt_starts[subset_protein == "BRD4"][Method != "mixed"],
       aes(x = reorder(as.character(Time), Time),
           y = log2FC, color = Method)) + # , fill = Method
  geom_boxplot(linewidth = 1.2) + #, outlier.shape = NA
  geom_line(aes(x = reorder(as.character(Time), Time),
                y = GoldStandard, group = "1"),
            inherit.aes = FALSE, color = "darkblue", size = 1.5) +
  facet_grid(num_unique ~ Protein) +
  # scale_fill_manual(values = colors,
  #                   labels = c(unique = "unique",
  #                              shared = "proposed",
  #                              all = "all")) +
  scale_color_manual(values = colors,
                     labels = c(unique = "unique",
                                shared = "proposed",
                                all = "all")) +
  geom_point(aes(x = reorder(as.character(Time), Time),
                 y = GoldStandard),
             inherit.aes = FALSE, color = "darkblue", size = 1.2) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14), 
        legend.box = "vertical") +
  xlab("time") +
  ylab("log-fold change") +
  ylim(c(-1.5, 0.5))
# ggsave("brd_subset_log2fcs_2uni.pdf", width = 10, height = 5, units = "in", scale = 1)