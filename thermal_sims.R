library(data.table)
library(MSstatsWeightedSummary)
library(ggplot2)
library(splitstackshape)
library(parallel)
library(pbapply)

annot = readRDS("input_data/onepot_tpp/annotation_onePot.RDS")
psms = fread("input_data/onepot_tpp/OnePotTPPStauro_Fx_PSMs.txt")

pp_orig = unique(psms[, .(PeptideSequence = `Annotated Sequence`, Protein = `Protein Accessions`)])
pp_orig_long = cSplit(pp_orig, sep = ";", direction = "long", drop = FALSE, splitCols = "Protein")
pp_orig_long = merge(pp_orig_long, pp_orig[, .(PeptideSequence, ProteinOrig = Protein)], by = "PeptideSequence")
pp_orig_long

pp_orig[PeptideSequence == "[-R].mASGRPEELWEAVVGAAER.[F]"]
pp_orig_long[PeptideSequence == "[-R].mASGRPEELWEAVVGAAER.[F]"]
pp_orig_long[PeptideSequence == "[-R].mASGRPEELWEAVVGAAER.[F]", Protein][1]

pp_orig_long2 = copy(pp_orig_long)
pp_orig_long2[, Protein := stringr::str_replace_all(Protein, "\\-1", "")]
pp_orig_long2[, ProteinOrig := stringr::str_replace_all(ProteinOrig, "\\-1", "")]

psms2 = copy(psms)
colnames(psms2)

setnames(pp_orig_long2, "Protein", "ProteinName")
psms2 = merge(psms2, pp_orig_long2, by.x = "Annotated Sequence", by.y = "PeptideSequence", allow.cartesian = T, all.x = T, all.y = T)
psms2 = as.data.table(psms2)
psms2[, IsUnique := uniqueN(ProteinName) == 1, by = c("Annotated Sequence", "Charge")]

processed.features_sh <- MSstatsTMT::PDtoMSstatsTMTFormat(psms2,
                                                          annot,
                                                          which.proteinid = "ProteinName",
                                                          useNumProteinsColumn = FALSE,
                                                          useUniquePeptide = FALSE,
                                                          rmPSM_withfewMea_withinRun = TRUE,
                                                          rmProtein_with1Feature = TRUE,
                                                          summaryforMultipleRows = max,
                                                          use_log_file = FALSE,
                                                          append=FALSE,
                                                          verbose=TRUE)
processed.features_sh = as.data.table(processed.features_sh)
processed.features_sh[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
processed.features_sh[!(IsUnique)]
processed.features_sh[, log2Intensity := log(Intensity, 2)]

final_graph = createPeptideProteinGraph(processed.features_sh)
processed.features_sh = addClusterMembership(processed.features_sh, final_graph)
iso_withsub = processIsoforms(processed.features_sh, remove_subsets = FALSE)
iso_nosub = processIsoforms(processed.features_sh, remove_subsets = TRUE)
iso_nosub[, Intensity := 2 ^ log2Intensity]
iso_nosub = normalizeSharedPeptides(iso_nosub)
final_graph = createPeptideProteinGraph(iso_nosub)
iso_nosub = addClusterMembership(iso_nosub, final_graph)
iso_nosub[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
iso_nosub[, UniqueOnly := all(IsUnique), by = "ProteinName"]
iso_nosub[(UniqueOnly), uniqueN(ProteinName)]
iso_nosub[(UniqueOnly), uniqueN(Cluster)]

iso_nosub[, uniqueN(ProteinName)]
iso_nosub[, uniqueN(Cluster)]
# saveRDS(iso_nosub, "processed_data/onepot_tpp/iso_nosub.RDS")

ggplot(iso_nosub, aes(x = Channel, y = log2IntensityNormalized)) +
  geom_boxplot()

# shared_summaries = lapply(split(iso_nosub, iso_nosub$Cluster),
#                                function(x) {
#                                  tryCatch({
#                                    getWeightedProteinSummary(x, "Huber", 1e-3, max_iter = 100, tolerance = 1e-3)  
#                                  }, error = function(e) NULL)
#                                })
# saveRDS(shared_summaries, "processed_data/onepot_ttp/shared_summaries.RDS")
shared_summaries = readRDS("processed_data/onepot_tpp/shared_summaries.RDS")

# unique_based_nontrivial_clusters = mclapply(split(iso_nosub[!(UniqueOnly)], iso_nosub[!(UniqueOnly), Cluster]),
#                                             function(x) {
#                                               getWeightedProteinSummary(x[(IsUnique)], "Huber", 1e-3, max_iter = 100, tolerance = 1e-3)
#                                             }, mc.cores = 2)
# saveRDS(unique_based_nontrivial_clusters, "processed_data/onepot_tpp/unique_nontrivial_cls.RDS")
unique_based_nontrivial_clusters = readRDS("processed_data/onepot_tpp/unique_nontrivial_cls.RDS")

# penalty_nontrivial_clusters = mclapply(split(iso_nosub[!(UniqueOnly)], iso_nosub[!(UniqueOnly), Cluster]),
#                                             function(x) {
#                                               getWeightedProteinSummary(x, "Huber", 1e-3, max_iter = 100, tolerance = 1e-3,
#                                                                         weights_penalty = TRUE, weights_penalty_param = 1e-4)
#                                             }, mc.cores = 2)
# saveRDS(penalty_nontrivial_clusters, "processed_data/onepot_tpp/penalized_nontrivial_cls.RDS")
penalty_nontrivial_clusters = readRDS("processed_data/onepot_tpp/penalized_nontrivial_cls.RDS")
failed_ids = which(sapply(shared_summaries, is.null))

failed_list = split(iso_nosub, iso_nosub$Cluster)[failed_ids]
failed_list = lapply(failed_list, function(x) {
  x[, NumObs := sum(!is.na(log2IntensityNormalized)), by = c("PSM", "ProteinName", "Run")]
  x
})
failed_ones = lapply(failed_list, function(x) {
  tryCatch(getWeightedProteinSummary(x[NumObs > 13], "Huber", 1e-2, max_iter = 100, tolerance = 1e-3),
           error = function(e) NULL)
})

shared_all_OK = c(shared_summaries, failed_ones)
shared_all_OK = shared_all_OK[!sapply(shared_all_OK, is.null)]

table(sapply(shared_summaries, is.null))
table(sapply(penalty_nontrivial_clusters, is.null))
table(sapply(unique_based_nontrivial_clusters, is.null))
table(sapply(shared_all_OK, is.null))

protein_data_shared = rbindlist(lapply(shared_all_OK, proteinData))
protein_data_shared = merge(protein_data_shared, unique(iso_nosub[, .(ProteinName, UniqueOnly)]), by.x = "Protein", by.y = "ProteinName")

protein_unique = rbindlist(lapply(unique_based_nontrivial_clusters, proteinData))
protein_unique[, UniqueOnly := FALSE]
protein_unique = rbind(protein_unique, protein_data_shared[(UniqueOnly)])

feat_data_shared = rbindlist(lapply(shared_all_OK, featureData))
feat_data_shared = merge(feat_data_shared, unique(iso_nosub[, .(ProteinName, UniqueOnly)]), by = "ProteinName")
feat_data_shared

feat_data_unique = rbindlist(lapply(unique_based_nontrivial_clusters, featureData))
feat_data_unique[, UniqueOnly := FALSE]
feat_data_unique = rbind(feat_data_unique, feat_data_shared[(UniqueOnly)])

uniqueN(feat_data_shared$ProteinName)
uniqueN(feat_data_unique$ProteinName)

protein_data_penalty = rbindlist(lapply(penalty_nontrivial_clusters, proteinData))
protein_data_penalty[, UniqueOnly := FALSE]
protein_data_penalty = rbind(protein_data_penalty, protein_data_shared[(UniqueOnly)])

feat_data_penalty = rbindlist(lapply(penalty_nontrivial_clusters, featureData))
feat_data_penalty[, UniqueOnly := FALSE]
feat_data_penalty = rbind(feat_data_penalty, feat_data_shared[(UniqueOnly)])

uniqueN(iso_nosub$ProteinName)
uniqueN(protein_data_shared$Protein)
uniqueN(protein_unique$Protein)
uniqueN(protein_data_penalty$Protein)

cm = readRDS("input_data/onepot_tpp/contrast_matrix.RDS")

gc_sh = MSstatsTMT::groupComparisonTMT(list(ProteinLevelData = protein_data_shared,
                                            FeatureLevelData = feat_data_shared), cm, use_log_file = FALSE)
gc_un = MSstatsTMT::groupComparisonTMT(list(ProteinLevelData = protein_unique,
                                            FeatureLevelData = feat_data_unique), cm, use_log_file = FALSE)
gc_penalty = MSstatsTMT::groupComparisonTMT(list(ProteinLevelData = protein_data_penalty,
                                                 FeatureLevelData = feat_data_penalty), cm, use_log_file = FALSE)

interactors = readRDS("input_data/onepot_tpp/kinases_and_direct_interactors.RDS")
interactors_fixed = stringr::str_replace(interactors, "\\-1", "")

sh_full = as.data.table(gc_sh$ComparisonResult)[adj.pvalue < 0.05, unique(Protein)]
sh_full = stringr::str_split(sh_full, ";", simplify = T)
sh_full = as.vector(sh_full)
sh_full = sh_full[sh_full != ""]

uni_penal = as.data.table(gc_un$ComparisonResult)[adj.pvalue < 0.05, unique(Protein)]
uni_penal = stringr::str_split(uni_penal, ";", simplify = T)
uni_penal = as.vector(uni_penal)
uni_penal = uni_penal[uni_penal != ""]

sh_smallpenalty_3 = as.data.table(gc_penalty$ComparisonResult)[adj.pvalue < 0.05, unique(Protein)]
sh_smallpenalty_3 = stringr::str_split(sh_smallpenalty_3, ";", simplify = T)
sh_smallpenalty_3 = as.vector(sh_smallpenalty_3)
sh_smallpenalty_3 = sh_smallpenalty_3[sh_smallpenalty_3 != ""]

xtable::xtable(data.table("Method" = c("proposed, penalized", "proposed, no penalty", "unique"),
                          "No. discovered interactors" = c(length(intersect(sh_smallpenalty_3, interactors_fixed)),
                                                           length(intersect(sh_full, interactors_fixed)),
                                                           length(intersect(uni_penal, interactors_fixed))),
                          "No. other discoveries" = c(length(setdiff(sh_smallpenalty_3, interactors_fixed)),
                                                      length(setdiff(sh_full, interactors_fixed)),
                                                      length(setdiff(uni_penal, interactors_fixed)))))

uniqueN(gc_uni_penal$ComparisonResult$Protein)
uniqueN(gc_sh_full$ComparisonResult$Protein)
uniqueN(gc_sh_smallpenalty_3$ComparisonResult$Protein)

data.table("Method" = c("proposed, penalized", "proposed, no penalty", "unique"),
           "No. discovered interactors" = c(length(intersect(sh_smallpenalty_3, interactors_fixed)),
                                            length(intersect(sh_full, interactors_fixed)),
                                            length(intersect(uni_penal, interactors_fixed))),
           "No. other discoveries" = c(length(setdiff(sh_smallpenalty_3, interactors_fixed)),
                                       length(setdiff(sh_full, interactors_fixed)),
                                       length(setdiff(uni_penal, interactors_fixed)))) 

# iso_nosub[grepl("Q9H1B7", ProteinName), unique(Cluster)]
# unique(iso_nosub[Cluster == 1788, .(PSM, ProteinName, IsUnique)])[(IsUnique)]
# dcast(unique(iso_nosub[Cluster == 1788, .(PSM, ProteinName, IsUnique)]), ProteinName ~ ProteinName, value.var = "PSM", fun.aggregate = uniqueN)

sig_prots = intersect(union(sh_smallpenalty_3, uni_penal), interactors_fixed)
is_sig = rep(FALSE, nrow(iso_nosub))
for (prot in sig_prots) {
  is_sig = is_sig | grepl(prot, iso_nosub$ProteinName)
}
cls_sim = iso_nosub[is_sig, unique(Cluster)]
uniqueN(cls_sim)


input_sim = iso_nosub[Cluster %in% cls_sim]
input_sim[(IsUnique), .(NumUnique = uniqueN(PSM)), by = "ProteinName"]
input_sim[!(IsUnique), .(NumUnique = uniqueN(PSM)), by = "ProteinName"]
input_sim[!(IsUnique), .(NumUnique = uniqueN(PSM)), by = "Cluster"]

length(intersect(union(sh_smallpenalty_3, uni_penal), interactors_fixed))
top_3 = input_sim[(IsUnique), .(MeanIntensity = mean(log2IntensityNormalized), na.rm = TRUE),
                  by = c("ProteinName", "PSM")][order(ProteinName, -MeanIntensity)][, .(PSM, MeanIntensity, Rank = rank(-MeanIntensity)),
                                                                                    by = "ProteinName"][Rank <= 3]
uniqueN(input_sim$ProteinName)
uniqueN(top_3$ProteinName)

unique(input_sim$ProteinName)
unique(top_3$ProteinName)

uniqueN(input_sim$Cluster)
length(unique(cls_sim))

input_sim2 = input_sim[!(IsUnique) | PSM %in% unique(top_3$PSM)]

input_sim2[, uniqueN(Cluster)]
input_sim2[(IsUnique), uniqueN(Cluster)]

input_sim2[, uniqueN(ProteinName)]
input_sim2[(IsUnique), uniqueN(ProteinName)]

input_top_3_list = split(input_sim2, input_sim2$ProteinName)
input_top_3_list = input_top_3_list[sapply(input_top_3_list, function(x) nrow(x[(IsUnique)])) > 0]
reference_summaries2 = lapply(input_top_3_list, function(x) {
  getWeightedProteinSummary(x[(IsUnique)], "Huber", 1e-3, tolerance = 0.01, max_iter = 100,
                            weights_penalty = FALSE, weights_penalty_param = NULL)
})

get_gc = function(summaries) {
  mst_inputs = lapply(summaries, makeMSstatsTMTInput)
  protein_data = rbindlist(lapply(mst_inputs, function(x) x[[2]]))
  feature_data = rbindlist(lapply(mst_inputs, function(x) x[[1]]))
  
  gc_sh = MSstatsTMT::groupComparisonTMT(list(ProteinLevelData = protein_data,
                                              FeatureLevelData = feature_data),
                                         cm, use_log_file = F)
  gc_sh 
}


# ref_gc = get_gc(reference_summaries)
ref_gc2 = get_gc(reference_summaries2)
# just_shared_gc = get_gc(shared_only_summaries)

table(unique(input_sim[, .(Cluster, UniqueOnly)])[, UniqueOnly])

input_sim = input_sim[!(UniqueOnly)]
input_sim[, uniqueN(Cluster)]

cls_final_sim = unique(input_sim$Cluster)
# saveRDS(input_sim, "processed_data/onepot_tpp/input_sim.RDS")
input_sim = readRDS("processed_data/onepot_tpp/input_sim.RDS")
num_iter = 50

# sim_summs = lapply(1:2, function(num_unique) {
#   mclapply(seq_len(num_iter), function(iter) {
#     lapply(cls_final_sim, function(cl) {
#       tryCatch({
#         single_cl_data = input_sim[Cluster == cl]
#         set.seed(2 * iter)
#         random_unique_per_protein = lapply(unique(single_cl_data$ProteinName), function(prot) {
#           unique_pepts = single_cl_data[ProteinName == prot & (IsUnique), unique(PSM)]
#           sample(unique_pepts, min(num_unique, uniqueN(unique_pepts)))
#         })
#         random_unique_per_protein = unlist(random_unique_per_protein, F, F)
#         
#         summ_sh = getWeightedProteinSummary(single_cl_data[PSM %in% random_unique_per_protein | !(IsUnique)],
#                                             "Huber", 1e-3, tolerance = 0.01, max_iter = 100,
#                                             weights_penalty = TRUE, weights_penalty_param = 1e-3)
#         summ_uni = getWeightedProteinSummary(single_cl_data[PSM %in% random_unique_per_protein],
#                                              "Huber", 1e-3, tolerance = 0.01, max_iter = 100,
#                                              weights_penalty = FALSE, weights_penalty_param = 1e-3)
#         summ_sh_nopen = getWeightedProteinSummary(single_cl_data[PSM %in% random_unique_per_protein | !(IsUnique)],
#                                                   "Huber", 1e-3, tolerance = 0.01, max_iter = 100,
#                                                   weights_penalty = FALSE, weights_penalty_param = NULL)
#         summ_sh_1iter = getWeightedProteinSummary(single_cl_data[PSM %in% random_unique_per_protein | !(IsUnique)],
#                                                   "Huber", 1e-3, tolerance = 0.01, max_iter = 1,
#                                                   weights_penalty = FALSE, weights_penalty_param = NULL)
#         list(summ_sh,
#              summ_uni,
#              summ_sh_nopen,
#              summ_sh_1iter,
#              iter = iter,
#              num_unique = num_unique,
#              cl = cl)
#         
#       }, error = function(e) NULL)
#     })
#   }, mc.cores = 7)
# })
# saveRDS(sim_summs, "results/onepot_tpp/sim_res.RDS")
sim_summs = readRDS("results/onepot_tpp/sim_res.RDS")

cm = readRDS("contrast_matrix.RDS")
sim_gcs = lapply(sim_summs, function(by_num_unique) {
  lapply(by_num_unique, function(by_iter) {
    lapply(by_iter, function(by_cluster) {
      tryCatch({
        iter = by_cluster$iter
        num_unique = by_cluster$num_unique
        cl = by_cluster$cl
        
        print(paste(num_unique, iter, cl))
        
        gcs = lapply(by_cluster[1:4], function(summ) {
          MSstatsTMT::groupComparisonTMT(makeMSstatsTMTInput(summ), cm, use_log_file = FALSE)
        })
        gcs = lapply(1:4, function(i) {
          gc_dt = as.data.table(gcs[[i]]$ComparisonResult)
          gc_dt[, Method := c("penalized", "unique", "shared", "shared, 1 iter")[i]]
          gc_dt
        })
        gcs = rbindlist(gcs)
        gcs[, iter := iter]
        gcs[, num_unique := num_unique]
        gcs[, cluster := cl]
        gcs
      }, error = function(e) NULL)
    })
  })
})

sim_gcs_dt = rbindlist(unlist(unlist(sim_gcs, F, F), F, F))
problems = unique(sim_gcs_dt[, .(cluster, num_unique, iter)])[, .(num_cls = uniqueN(cluster)), by = c("num_unique", "iter")][num_cls != 14]

fixed_iters = lapply(1:2, function(issue_num_unique) {
  mclapply(problems[num_unique == issue_num_unique, unique(iter)], function(issue_iter) {
    present_clusters = sim_gcs_dt[num_unique == issue_num_unique & iter == issue_iter, unique(cluster)]
    problem_clusters = setdiff(cls_final_sim, present_clusters)
    lapply(problem_clusters, function(cl) {
      tryCatch({
        single_cl_data = input_sim[Cluster == cl]
        set.seed(2 * issue_iter)
        random_unique_per_protein = lapply(unique(single_cl_data$ProteinName), function(prot) {
          unique_pepts = single_cl_data[ProteinName == prot & (IsUnique), unique(PSM)]
          sample(unique_pepts, min(issue_num_unique, uniqueN(unique_pepts)))
        })
        random_unique_per_protein = unlist(random_unique_per_protein, F, F)
        
        summ_sh = getWeightedProteinSummary(single_cl_data[PSM %in% random_unique_per_protein | !(IsUnique)],
                                            "Huber", 1e-2, tolerance = 0.01, max_iter = 100,
                                            weights_penalty = TRUE, weights_penalty_param = 1e-3)
        summ_uni = getWeightedProteinSummary(single_cl_data[PSM %in% random_unique_per_protein],
                                             "Huber", 1e-2, tolerance = 0.01, max_iter = 100,
                                             weights_penalty = FALSE, weights_penalty_param = 1e-3)
        summ_sh_nopen = getWeightedProteinSummary(single_cl_data[PSM %in% random_unique_per_protein | !(IsUnique)],
                                                  "Huber", 1e-2, tolerance = 0.01, max_iter = 100,
                                                  weights_penalty = FALSE, weights_penalty_param = NULL)
        summ_sh_1iter = getWeightedProteinSummary(single_cl_data[PSM %in% random_unique_per_protein | !(IsUnique)],
                                                  "Huber", 1e-2, tolerance = 0.01, max_iter = 1,
                                                  weights_penalty = FALSE, weights_penalty_param = NULL)
        list(summ_sh,
             summ_uni,
             summ_sh_nopen,
             summ_sh_1iter,
             iter = issue_iter,
             num_unique = issue_num_unique,
             cl = cl)
      }, error = function(e) NULL)
    })
  }, mc.cores = 4)
})
# lapply(fixed_iters, function(x) lapply(x, function(y) sapply(y, is.null)))
sim_gcs_fixed = lapply(fixed_iters, function(by_num_unique) {
  lapply(by_num_unique, function(by_iter) {
    lapply(by_iter, function(by_cluster) {
      tryCatch({
        iter = by_cluster$iter
        num_unique = by_cluster$num_unique
        cl = by_cluster$cl
        
        print(paste(num_unique, iter, cl))
        
        gcs = lapply(by_cluster[1:4], function(summ) {
          MSstatsTMT::groupComparisonTMT(makeMSstatsTMTInput(summ), cm, use_log_file = FALSE)
        })
        gcs = lapply(1:4, function(i) {
          gc_dt = as.data.table(gcs[[i]]$ComparisonResult)
          gc_dt[, Method := c("penalized", "unique", "shared", "shared, 1 iter")[i]]
          gc_dt
        })
        gcs = rbindlist(gcs)
        gcs[, iter := iter]
        gcs[, num_unique := num_unique]
        gcs[, cluster := cl]
        gcs
      }, error = function(e) NULL)
    })
  })
})
sim_gcs_fixed_dt = rbindlist(unlist(unlist(sim_gcs_fixed, F, F), F, F))

sim_gcs_dt = rbind(sim_gcs_dt, sim_gcs_fixed_dt)

top_3 = input_sim[(IsUnique), .(MeanIntensity = mean(log2IntensityNormalized), na.rm = TRUE),
                  by = c("ProteinName", "PSM")][order(ProteinName, -MeanIntensity)][, .(PSM, MeanIntensity, Rank = rank(-MeanIntensity)),
                                                                                    by = "ProteinName"][Rank <= 3]
input_sim2 = input_sim[!(IsUnique) | PSM %in% unique(top_3$PSM)]
input_top_3_list = split(input_sim2, input_sim2$ProteinName)
input_top_3_list = input_top_3_list[sapply(input_top_3_list, function(x) nrow(x[(IsUnique)])) > 0]
reference_summaries2 = lapply(input_top_3_list, function(x) {
  getWeightedProteinSummary(x[(IsUnique)], "Huber", 1e-3, tolerance = 0.01, max_iter = 100,
                            weights_penalty = FALSE, weights_penalty_param = NULL)
})

get_gc = function(summaries) {
  mst_inputs = lapply(summaries, makeMSstatsTMTInput)
  protein_data = rbindlist(lapply(mst_inputs, function(x) x[[2]]))
  feature_data = rbindlist(lapply(mst_inputs, function(x) x[[1]]))
  
  gc_sh = MSstatsTMT::groupComparisonTMT(list(ProteinLevelData = protein_data,
                                              FeatureLevelData = feature_data),
                                         cm, use_log_file = F)
  gc_sh 
}

ref_gc2 = get_gc(reference_summaries2)

sim_gcs_dt = merge(sim_gcs_dt, as.data.table(ref_gc2$ComparisonResult)[, .(Protein, Label, adj.pvalue, log2FC)],
                   by = c("Protein", "Label"), all.x = T, suffixes = c("", "_ref"))

sim_gcs_dt_itermse = sim_gcs_dt[, .(IterMSE = mean((log2FC - log2FC_ref)^2, na.rm = T)),
                                by = c("iter", "Method", "num_unique")]
sim_gcs_dt_itermse[, NumUnique2 := ifelse(num_unique == 1, "#unique = 1", "#unique = 2")]
sim_gcs_dt_itermse[, Method2 := ifelse(Method == "penalized", "shared, penalized", Method)]
sim_gcs_dt_itermse[, Method3 := stringr::str_replace_all(Method2, "shared", "proposed")]
sim_gcs_dt_itermse[, Method3 := stringr::str_replace_all(Method3, "iter", "iteration")]
unique(sim_gcs_dt_itermse$Method2)
unique(sim_gcs_dt_itermse$Method3)

ggplot(sim_gcs_dt_itermse[!is.na(Method)],
       aes(x = Method3, y = IterMSE, fill = as.character(num_unique))) +
  geom_boxplot() +
  # facet_grid(~NumUnique2) +
  theme_bw() +
  # theme(legend.position = "bottom",
  #       axis.text = element_text(size = 22),
  #       axis.) +
  scale_fill_manual(name = "number of unique peptides", values = c("#E69F00", "#56B4E9")) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14), 
        legend.position = "bottom",
        legend.box = "vertical") +
  coord_flip() +
  xlab("method") +
  ylab("MSE of log-fold change estimation")
ggsave("resampling_mse.pdf", device = "pdf", width = 10, height = 5, units = "in")
