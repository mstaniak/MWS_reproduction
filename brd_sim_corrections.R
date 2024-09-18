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

num_rep = 30

sim_100_starts = lapply(
  c(2, 5, 10), 
  function(num_unique) {
    lapply(seq_len(num_rep), function(i) {
      print(i)
      set.seed(i)
      unique_peptides = unlist(lapply(unique(brd_cluster$ProteinName),
                                      function(x) {
                                        sample(brd_cluster[(IsUnique)][ProteinName == x, unique(PSM)], num_unique)
                                      }), F, F)
      summary_input_1_triple = brd_cluster[!(IsUnique) | PSM %in% unique_peptides]
      summary_input_1_triple[, Channel := as.factor(as.character(Channel))]
      
      summary_l1 = tryCatch(getWeightedProteinSummary(summary_input_1_triple,
                                             norm = "p_norm",
                                             norm_parameter = 1,
                                             initial_summary = "unique",
                                             tolerance = 1e-2, max_iter = 100),
                            error = function(e) NULL)
      summary_l2 = tryCatch(getWeightedProteinSummary(summary_input_1_triple,
                                             norm = "p_norm",
                                             norm_parameter = 2,
                                             initial_summary = "unique",
                                             tolerance = 1e-2, max_iter = 100),
                            error = function(e) NULL) 
      summary_huber = tryCatch(getWeightedProteinSummary(summary_input_1_triple,
                                                norm = "Huber",
                                                norm_parameter = 1e-3,
                                                initial_summary = "unique",
                                                tolerance = 1e-2, max_iter = 100),
                               error = function(e) NULL)
      
          list(l1 = summary_l1,
               l2 = summary_l2,
               huber = summary_huber,
               iter =  i,
               num_unique = num_unique)
    })
  })
sim_100_starts[[1]][[2]]

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