# Libraries ----
library(data.table)
library(MSstatsWeightedSummary)
library(ggplot2)
library(splitstackshape)
library(parallel)
library(ggVennDiagram)
library(patchwork)
# Raw data ----
# ## onePot
# psms_onepot = fread("./input_data/onepot_tpp/OnePotTPPStauro_Fx_PSMs.txt")
# annot_onepot = readRDS("./input_data/onepot_tpp/annotation_onePot.RDS")
# ## TPP
# curve_input = fread("input_data/onepot_tpp/XuY_ACSChemBio_2021_FullCurve_PSMs (2).txt")
# curve_annot = readRDS("input_data/onepot_tpp/Annotation_file_thermal.RDS")
# curve_annot$Channel = stringr::str_replace_all(curve_annot$Channel, "Abundance", "Abundance:")
# colnames(curve_input)[grepl("Abundance", colnames(curve_input))]
# colnames(curve_input) = stringr::str_replace_all(colnames(curve_input), "Abundance", "Abundance:")
# colnames(curve_input)[grepl("Abundance", colnames(curve_input))]
# colnames(curve_input) = stringr::str_replace(colnames(curve_input), "131N", "131")
# curve_input[, `Abundance: 131C` := NULL]
# Original pre-processing ----
# ## onePot
# orig_onepot <- MSstatsTMT::PDtoMSstatsTMTFormat(psms_onepot,
#                                                 annot_onepot,
#                                                 which.proteinid = "Master.Protein.Accessions",
#                                                 useNumProteinsColumn = FALSE,
#                                                 useUniquePeptide = TRUE,
#                                                 rmPSM_withfewMea_withinRun = TRUE,
#                                                 rmProtein_with1Feature = TRUE,
#                                                 summaryforMultipleRows = max,
#                                                 use_log_file = FALSE,
#                                                 append=FALSE,
#                                                 verbose=TRUE)
# 
# ## TPP
# orig_tpp <- MSstatsTMT::PDtoMSstatsTMTFormat(curve_input,
#                                              curve_annot,
#                                              which.proteinid = "Master.Protein.Accessions",
#                                              useNumProteinsColumn = FALSE,
#                                              useUniquePeptide = TRUE,
#                                              rmPSM_withfewMea_withinRun = TRUE,
#                                              rmProtein_with1Feature = TRUE,
#                                              summaryforMultipleRows = max,
#                                              use_log_file = FALSE,
#                                              append=FALSE,
#                                              verbose=TRUE)
# saveRDS(orig_onepot, "processed_data/onepot_tpp/orig_onepot.RDS")
# saveRDS(orig_tpp, "processed_data/onepot_tpp/orig_tpp.RDS")
orig_onepot = readRDS("processed_data/onepot_tpp/orig_onepot.RDS")
orig_tpp = readRDS("processed_data/onepot_tpp/orig_tpp.RDS")
# Full peptide-protein graphs ----
# ## onePot
# onepot_pp_orig = unique(psms_onepot[, .(PeptideSequence = `Annotated Sequence`, ProteinName = `Protein Accessions`)])
# onepot_pp = cSplit(onepot_pp_orig, sep = ";", direction = "long", drop = FALSE, splitCols = "ProteinName")
# onepot_pp[, ProteinName := stringr::str_replace_all(ProteinName, "\\-1", "")]
# onepot_pp_graph = createPeptideProteinGraph(onepot_pp)
# ## TPP
# curve_pp_raw = unique(curve_input[, .(PeptideSequence = `Annotated Sequence`, ProteinName = `Protein Accessions`)])
# curve_pp = cSplit(curve_pp_raw, splitCols = "ProteinName", sep = "; ", direction = "long", drop = FALSE)
# curve_pp[, ProteinName := stringr::str_replace_all(ProteinName, "\\-1", "")]
# curve_pp_graph = createPeptideProteinGraph(curve_pp)
# # Cluster identification ----
# ## onePot
# onepot_pp = addClusterMembership(onepot_pp, onepot_pp_graph)
# onepot_input_all_prots = merge(psms_onepot, onepot_pp, by.x = "Annotated Sequence", by.y = "PeptideSequence", allow.cartesian = T, all.x = T, all.y = T)
# onepot_input_all_prots[, `Master Protein Accessions` := NULL]
# onepot_input_all_prots[, `Protein Accessions` := NULL]
# ## TPP
# curve_data_pp = addClusterMembership(curve_pp, curve_pp_graph)
# curve_input_all_prots = merge(curve_input, curve_data_pp,
#                               by.x = "Annotated Sequence", by.y = "PeptideSequence",
#                               all.x = T, all.y = T, allow.cartesian = TRUE)
# curve_input_all_prots[, `Master Protein Accessions` := NULL]
# curve_input_all_prots[, `Protein Accessions` := NULL]
# # MSstatsTMT pre-processing ----
# rm(onepot_pp_graph, curve_pp_graph, onepot_pp_orig, curve_pp_raw)
# gc()
## onePot
# onepot_procd = MSstatsTMT::PDtoMSstatsTMTFormat(onepot_input_all_prots, annot_onepot,
#                                                 which.proteinid = "ProteinName",
#                                                 useNumProteinsColumn = FALSE,
#                                                 useUniquePeptide = FALSE,
#                                                 rmPSM_withfewMea_withinRun = TRUE,
#                                                 rmProtein_with1Feature = TRUE,
#                                                 summaryforMultipleRows = max,
#                                                 use_log_file = FALSE,
#                                                 append=FALSE,
#                                                 verbose=TRUE)
# saveRDS(onepot_procd, "processed_data/onepot_tpp/onepot_procd.RDS")
# onepot_procd = readRDS("processed_data/onepot_tpp/onepot_procd.RDS")
# onepot_procd = as.data.table(onepot_procd)
# rm(psms_onepot, onepot_pp, onepot_input_all_prots, curve_data_pp, curve_pp, curve_input)
# gc()
## TPP
# curve_procd = MSstatsTMT::PDtoMSstatsTMTFormat(curve_input_all_prots, curve_annot,
#                                                which.proteinid = "ProteinName",
#                                                useNumProteinsColumn = FALSE,
#                                                useUniquePeptide = FALSE,
#                                                rmPSM_withfewMea_withinRun = TRUE,
#                                                rmProtein_with1Feature = TRUE,
#                                                summaryforMultipleRows = max,
#                                                use_log_file = FALSE,
#                                                append=FALSE,
#                                                verbose=TRUE)
# saveRDS(curve_procd, "processed_data/onepot_tpp/curve_procd.RDS")
# curve_procd = readRDS("processed_data/onepot_tpp/curve_procd.RDS")
# curve_procd = as.data.table(curve_procd)
# Process isoforms ---
# ## onePot
# onepot_procd[, log2Intensity := log(Intensity, 2)]
# onepot_graph_procd = createPeptideProteinGraph(onepot_procd)
# onepot_procd = addClusterMembership(onepot_procd, onepot_graph_procd)
# onepot_procd_iso = processIsoforms(onepot_procd, T, T, F)
# saveRDS(onepot_procd_iso, "processed_data/onepot_tpp/onepot_procd_iso.RDS")
onepot_procd_iso = readRDS("processed_data/onepot_tpp/onepot_procd_iso.RDS")
# ## TPP
# curve_graph_procd = createPeptideProteinGraph(curve_procd)
# curve_procd = addClusterMembership(curve_procd, curve_graph_procd)
# curve_procd[, log2Intensity := log(Intensity, 2)]
# curve_procd_iso = processIsoforms(curve_procd, T, T, F)
# saveRDS(curve_procd_iso, "processed_data/onepot_tpp/curve_procd_iso.RDS")
curve_procd_iso = readRDS("processed_data/onepot_tpp/curve_procd_iso.RDS")
# Normalization ----
onepot_procd_iso[, Intensity := 2 ^ log2Intensity]
onepot_procd_iso = normalizeSharedPeptides(onepot_procd_iso)
onepot_procd_iso_graph = createPeptideProteinGraph(onepot_procd_iso)
onepot_procd_iso = addClusterMembership(onepot_procd_iso, onepot_procd_iso_graph)
rm(onepot_procd_iso_graph)
curve_procd_iso[, log2IntensityNormalized := log2Intensity]
# Cluster statistics ----
onepot_procd_iso = getClusterStatistics(onepot_procd_iso, TRUE)
curve_procd_iso = getClusterStatistics(curve_procd_iso, TRUE)
# Known interactors ----
interactors = readRDS("input_data/onepot_tpp/kinases_and_direct_interactors.RDS")
interactors_fixed = stringr::str_replace(interactors, "\\-1", "")

is_interact_onepot = rep(FALSE, nrow(onepot_procd_iso))
is_interact_tpp = rep(FALSE, nrow(curve_procd_iso))

for (interactor in interactors_fixed) {
  is_interact_onepot = is_interact_onepot | grepl(interactor, onepot_procd_iso$ProteinName)
  is_interact_tpp = is_interact_tpp | grepl(interactor, curve_procd_iso$ProteinName)
}

onepot_cls_int = onepot_procd_iso[is_interact_onepot][NumProteins > 1, unique(Cluster)]
tpp_cls_int = curve_procd_iso[is_interact_tpp][NumProteins > 1, unique(Cluster)]

uniqueN(onepot_cls_int)
uniqueN(tpp_cls_int)

# saveRDS(onepot_procd_iso[Cluster %in% onepot_cls_int], "processed_data/onepot_tpp/onepot_sub_int_cls.RDS")
# saveRDS(curve_procd_iso[Cluster %in% tpp_cls_int], "processed_data/onepot_tpp/tpp_sub_int_cls.RDS")

# Statistics for these datasets ---- 
## Cluster of interests - uniqueness
unique(onepot_procd_iso[Cluster %in% onepot_cls_int][, .(Cluster, EachHasUnique, NumProteins)])
onepot_procd_iso[Cluster %in% onepot_cls_int][(HasUnique), .(NumProteins = uniqueN(ProteinName)), by = "Cluster"][NumProteins > 1]
curve_procd_iso[Cluster %in% tpp_cls_int][(HasUnique), .(NumProteins = uniqueN(ProteinName)), by = "Cluster"][NumProteins > 1]
# Counts of features and proteins
uniqueN(orig_onepot$ProteinName)
uniqueN(onepot_procd_iso$ProteinName)

uniqueN(orig_tpp$ProteinName)
uniqueN(curve_procd_iso$ProteinName)

uniqueN(orig_onepot$PSM)
uniqueN(onepot_procd_iso$PSM)

uniqueN(orig_tpp$PSM)
uniqueN(curve_procd_iso$PSM)

curve_procd_iso[(HasUnique), uniqueN(ProteinName)]
curve_procd_iso[(HasUnique), uniqueN(PSM)]

uniqueN(onepot_procd_iso$Cluster)
uniqueN(curve_procd_iso$Cluster)

unique(onepot_procd_iso[, .(Cluster, NumProteins)])[, mean(NumProteins)]
unique(curve_procd_iso[, .(Cluster, NumProteins)])[, mean(NumProteins)]
unique(curve_procd_iso[(HasUnique), .(Cluster, NumProteins)])[, mean(NumProteins)]

unique(onepot_procd_iso[][!(IsUnique), .(NumPSM = uniqueN(PSM)), by = "Cluster"])[, mean(NumPSM)]
unique(curve_procd_iso[][!(IsUnique), .(NumPSM = uniqueN(PSM)), by = "Cluster"])[, mean(NumPSM)]

unique(onepot_procd_iso[(HasUnique)][!(IsUnique), .(NumPSM = uniqueN(PSM)), by = "Cluster"])[, mean(NumPSM)]
unique(curve_procd_iso[(HasUnique)][!(IsUnique), .(NumPSM = uniqueN(PSM)), by = "Cluster"])[, mean(NumPSM)]

# Example cluster (Table 2)
onepot_procd_iso[ProteinName == "Q7Z5L9", uniqueN(PSM)]
onepot_procd_iso[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
unique(onepot_procd_iso[ProteinName == "Q7Z5L9", .(PSM, IsUnique)])

onepot_procd_iso[ProteinName == "Q7Z5L9-2", uniqueN(PSM)]
unique(onepot_procd_iso[ProteinName == "Q7Z5L9-2", .(PSM, IsUnique)])

onepot_procd_iso[ProteinName == "Q9H1B7", uniqueN(PSM)]
unique(onepot_procd_iso[ProteinName == "Q9H1B7", .(PSM, IsUnique)])
