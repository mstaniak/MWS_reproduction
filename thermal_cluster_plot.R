library(data.table)
library(MSstatsWeightedSummary)
library(ggplot2)

iso_nosub = readRDS("processed_data/onepot_tpp/iso_nosub.RDS")
iso_nosub[Cluster == 1788]

cl_plot_dt = iso_nosub[Cluster == 1788]
cl_plot_dt[, XLabel := stringr::str_replace(BioReplicate, "_rep", "_")]
cl_plot_dt[, XLabel := stringr::str_replace(XLabel, "stauro", "")]
cl_plot_dt[, XLabel := stringr::str_replace(XLabel, "_", ", ")]
xlab_order = c("DMSO, 1", "DMSO, 2", "DMSO, 3",
               "1x, 1", "1x, 2", "1x, 3",
               "5x, 1", "5x, 2", "5x, 3",
               "10x, 1", "10x, 2", "10x, 3",
               "25x, 1", "25x, 2", "25x, 3")
xlab_order2 = c("0x, 1", "0x, 2", "0x, 3",
                "1x, 1", "1x, 2", "1x, 3",
                "5x, 1", "5x, 2", "5x, 3",
                "10x, 1", "10x, 2", "10x, 3",
                "25x, 1", "25x, 2", "25x, 3")
cl_plot_dt[, XLabel2 := stringr::str_replace(XLabel, "DMSO", "0x")]
cl_plot_dt[, XLabel := factor(XLabel, labels = xlab_order, ordered = TRUE)]
cl_plot_dt[, XLabel2 := factor(XLabel2, labels = xlab_order2, ordered = TRUE)]
cl_plot_dt[, feature := ifelse(IsUnique, "unique", "shared")]
cl_plot_dt[, feature := factor(feature, levels = c("unique", "shared"), ordered = TRUE)]
cl_plot_dt[, AbundanceStandard := log2IntensityNormalized - mean(log2IntensityNormalized), by = c("ProteinName", "PSM")]
cl_plot_dt[, AbundanceStandard2 := log2IntensityNormalized / sd(log2IntensityNormalized), by = c("ProteinName", "PSM")]
unique(cl_plot_dt$ProteinName)
# color = ifelse(PSM %in% most_outl, "outlying", "consistent")

summ_sh = getWeightedProteinSummary(cl_plot_dt, norm = "Huber", norm_parameter = 1e-3, tolerance = 1e-2, max_iter = 100, 
                                    weights_penalty = T, weights_penalty_param = 1e-3)
summ_un = getWeightedProteinSummary(cl_plot_dt[(IsUnique)], norm = "Huber", norm_parameter = 1e-3, tolerance = 1e-2, max_iter = 100)

cm_onepot = readRDS("input_data/onepot_tpp/contrast_matrix.RDS")

gc_sh = MSstatsTMT::groupComparisonTMT(makeMSstatsTMTInput(summ_sh), cm_onepot)
gc_un = MSstatsTMT::groupComparisonTMT(makeMSstatsTMTInput(summ_un), cm_onepot)

featureWeights(summ_sh)[, uniqueN(PSM)]
featureWeights(summ_sh)[ProteinName == "Q7Z5L9-2", uniqueN(PSM)]

mean_cors_3 = apply(cor(as.matrix(dcast(split(cl_plot_dt, cl_plot_dt$ProteinName)[[3]], Run + Channel ~ PSM, value.var = "log2IntensityNormalized")[, -(1:2), with = F]),
                        method = "spearman"), 1, function(x) mean(x[x != 1]))
names(sort(mean_cors_3)[1:2])

mean_cors_2 = apply(cor(as.matrix(dcast(split(cl_plot_dt, cl_plot_dt$ProteinName)[[2]], Run + Channel ~ PSM, value.var = "log2IntensityNormalized")[, -(1:2), with = F]),
                        method = "spearman"), 1, function(x) mean(x[x != 1]))
names(sort(mean_cors_2)[1:2])

mean_cors_1 = apply(cor(as.matrix(dcast(split(cl_plot_dt, cl_plot_dt$ProteinName)[[1]], Run + Channel ~ PSM, value.var = "log2IntensityNormalized")[, -(1:2), with = F]),
                        method = "spearman"), 1, function(x) mean(x[x != 1]))
names(sort(mean_cors_1)[1:2])

most_outl = c(names(sort(mean_cors_1)[1]), names(sort(mean_cors_2)[1]), names(sort(mean_cors_3)[1]))
most_outl = c(names(sort(mean_cors_1)[1:2]), names(sort(mean_cors_2)[1:2]), names(sort(mean_cors_3)[1:2]))

summ_comp = rbind(
  cbind(proteinData(summ_sh), Method = "proposed"),
  cbind(proteinData(summ_un), Method = "unique"),
  fill = T, use.names = T
)
summ_comp[, ProteinName := Protein]
summ_comp[, XLabel := stringr::str_replace(BioReplicate, "_rep", "_")]
summ_comp[, XLabel := stringr::str_replace(XLabel, "stauro", "")]
summ_comp[, XLabel := stringr::str_replace(XLabel, "_", ", ")]
summ_comp[, XLabel2 := stringr::str_replace(XLabel, "DMSO", "0x")]
summ_comp[, XLabel := factor(XLabel, labels = xlab_order, ordered = TRUE)]
summ_comp[, XLabel2 := factor(XLabel2, labels = xlab_order2, ordered = TRUE)]

# Profile plot for all proteins in the example cluster
ggplot(cl_plot_dt,#[PSM %in% names(sort(mean_cors_1)[1:2])],#[ProteinName == "Q9H1B7"][PSM == "[K].lGEEQQR.[Q]_2"],#[grepl("EEWAS", PSM)],
       aes(x = XLabel2, y = log2IntensityNormalized, 
           group = PSM, linetype = feature)) +
  geom_line(size = 0.9, color = "grey") + # 
  geom_line(aes(x = XLabel2, y = log2IntensityNormalized, group = PSM, linetype = feature),
            color = "purple", data = cl_plot_dt[PSM %in% most_outl], size = 0.9) +
  geom_line(aes(x = XLabel2, y = Abundance, group = Method, color = Method),
            data = summ_comp, size = 1.2, inherit.aes = FALSE) +
  scale_color_manual(name = "method", values = c("red", "lightblue")) +
  facet_grid( ~ ProteinName) +
  # scale_color_manual(name = "profile", values = c("grey", "red")) +
  # guides(color = "none") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 270)) +
  ylab("log-intensity") +
  xlab("biological replicate") +  
  theme(legend.position = "bottom",
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 18),
        legend.box = "vertical")

# Features with the small correlation to other features
apply(cor(as.matrix(dcast(split(cl_plot_dt, cl_plot_dt$ProteinName)[[3]], Run + Channel ~ PSM, value.var = "log2IntensityNormalized")[, -(1:2), with = F])), 1, mean)
apply(cor(as.matrix(dcast(split(cl_plot_dt, cl_plot_dt$ProteinName)[[3]], Run + Channel ~ PSM, value.var = "log2IntensityNormalized")[, -(1:2), with = F])), 1, function(x) max(x[x != 1]))
apply(cor(as.matrix(dcast(split(cl_plot_dt, cl_plot_dt$ProteinName)[[3]], Run + Channel ~ PSM, value.var = "log2IntensityNormalized")[, -(1:2), with = F])), 1, function(x) min(x[x != 1]))

ggplot(cl_plot_dt,#[PSM %in% names(sort(mean_cors_1)[1:2])],#[ProteinName == "Q9H1B7"][PSM == "[K].lGEEQQR.[Q]_2"],#[grepl("EEWAS", PSM)],
       aes(x = XLabel2, y = AbundanceStandard, 
           group = PSM, linetype = feature)) +
  geom_line(size = 0.9, color = "grey") + # 
  geom_line(aes(x = XLabel2, y = AbundanceStandard, group = PSM, linetype = feature),
            color = "purple", data = cl_plot_dt[PSM %in% most_outl], size = 0.9) +
  # geom_line(aes(x = XLabel2, y = Abundance, group = Method, color = Method),
  #           data = summ_comp, size = 1.5, inherit.aes = FALSE) +
  scale_color_manual(name = "method", values = c("red", "lightblue")) +
  facet_grid( ~ ProteinName) +
  # scale_color_manual(name = "profile", values = c("grey", "red")) +
  # guides(color = "none") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 270)) +
  ylab("log-intensity") +
  xlab("biological replicate") +  
  theme(legend.position = "bottom",
        legend.title = element_text(size = 22),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 22),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 20), 
        legend.box = "vertical")
