library(data.table)
library(ggplot2)

get_optim_problem = function(feature_data, channel_order) {
  feature_obs = dcast(unique(feature_data[, .(PSM, Channel, log2IntensityNormalized)]),
                      PSM ~ Channel, value.var = "log2IntensityNormalized")
  channel_names = colnames(feature_obs)[-1]
  feature_obs = as.matrix(feature_obs[, -1, with = FALSE])
  feature_obs = feature_obs[, channel_order]
  
  delta_matrix = dcast(unique(feature_data[, .(PSM, ProteinName, Present = 1)]),
                       PSM ~ ProteinName, value.var = "Present", fill = 0)
  features = delta_matrix[["PSM"]]
  delta_matrix = as.matrix(delta_matrix[, -1, with = FALSE])
  delta_sums = rowSums(delta_matrix)
  delta_sums_shared = delta_sums > 1
  delta_num_params = delta_sums - 1
  names(delta_num_params) = paste0("w", seq_along(delta_num_params))
  weights_counts = delta_num_params[delta_num_params > 0]
  
  prot_matrix = matrix(1, nrow = 1, ncol = uniqueN(feature_data$Channel))
  
  num_prot = uniqueN(feature_data$ProteinName)
  num_prot_params = num_prot - 1
  num_channels = uniqueN(feature_data$Channel)
  num_channel_params = num_channels - 1
  num_features = uniqueN(feature_data$PSM)
  num_weights = sum(weights_counts)
  
  function(parameters) {
    mu = parameters[1]
    prot_channels_params = parameters[2:(1 + num_channel_params * num_prot)]
    
    prot_channels = matrix(prot_channels_params,
                           nrow = num_prot, ncol = num_channel_params, byrow = TRUE)
    prot_channels = cbind(prot_channels, -rowSums(prot_channels))
    weights_params = parameters[(2 + num_channel_params * num_prot):(1 + num_channel_params * num_prot + num_weights)]
    weights_params_list = lapply(seq_along(weights_counts),
                                 function(i) {
                                   num_params = weights_counts[i]
                                   if (i == 1) {
                                     weights = weights_params[1:num_params]
                                   } else {
                                     weights = weights_params[(sum(weights_counts[1:(i - 1)]) + 1):((sum(weights_counts[1:(i - 1)]) + weights_counts[i]))]
                                   }
                                   weights = c(weights, 1 - sum(weights))
                                   weights
                                 })
    names(weights_params_list) = names(weights_counts)
    weight_matrix = delta_matrix
    
    for (i in seq_len(nrow(weight_matrix))) {
      if (delta_num_params[i] > 0) {
        weight_matrix[i, weight_matrix[i, ] != 0] = weights_params_list[[names(delta_num_params)[i]]]
      }
    }
    
    feature_params = parameters[(2 + num_channel_params * num_prot + num_weights):(length(parameters))]
    feature_params = c(feature_params, -sum(feature_params))
    
    feature_matrix = matrix(feature_params, ncol = 1) %*% matrix(1, ncol = num_channels)
    
    predicted = mu + (delta_matrix * weight_matrix) %*% prot_channels + feature_matrix
    diff = feature_obs - predicted
    
    sum(diag(t(diff) %*% diff))
  }
}

simulateFeatureModelNoiseless = function(num_channels,
                                         weights_matrix,
                                         intercept,
                                         protein_values,
                                         channel_values, 
                                         feature_values) {
  num_features = nrow(weights_matrix)
  num_proteins = ncol(weights_matrix)
  
  data_matrix = matrix(intercept, nrow = num_features,
                       ncol = num_channels)
  simulated_data = data_matrix + weights_matrix %*% (c(protein_values, -sum(protein_values)) %*% matrix(1, ncol = num_channels) + cbind(channel_values, -rowSums(channel_values))) + c(feature_values, -sum(feature_values)) %*% matrix(1, ncol = num_channels)
  row.names(simulated_data) = paste("Feature", 1:num_features, sep = "_")
  colnames(simulated_data) = paste("Channel", 1:num_channels, sep = "_")
  simulated_data
}

addNoise = function(noiseless_data, sd) {
  noisy_data = noiseless_data + matrix(rnorm(nrow(noiseless_data) * ncol(noiseless_data), sd = sd),
                                       nrow = nrow(noiseless_data))
  noisy_data
}

makeLongDataFormat = function(noisy_data, weights_matrix) {
  feature_ids = row.names(noisy_data)
  noisy_data = as.data.table(noisy_data)
  noisy_data[, PSM := feature_ids]
  noisy_data = melt(noisy_data, variable.name = "Channel", 
                    value.name = "log2IntensityNormalized",
                    variable.factor = FALSE)
  
  weights_matrix = as.data.table(weights_matrix)
  weights_matrix[, PSM := feature_ids]
  pp = melt(weights_matrix, variable.name = "ProteinName",
            value.name = "Weight", variable.factor = FALSE,
            id.vars = "PSM")
  pp = pp[Weight > 0]
  pp = unique(pp[, .(PSM, ProteinName, Present = ceiling(Weight))])
  pp[, Present := NULL]
  noisy_data = merge(noisy_data, pp, by = "PSM", allow.cartesian = TRUE)
  noisy_data
}

getWeightsMatrix = function(feature_data) {
  delta_matrix = dcast(unique(feature_data[, .(PSM, ProteinName, Present = 1)]),
                       PSM ~ ProteinName, value.var = "Present", fill = 0)
  features = delta_matrix[["PSM"]]
  delta_matrix = as.matrix(delta_matrix[, -1, with = FALSE])
  
  delta_sums = rowSums(delta_matrix)
  delta_sums_shared = delta_sums > 1
  delta_num_params = delta_sums - 1
  names(delta_num_params) = paste0("w", seq_along(delta_num_params))
  weights_counts = delta_num_params[delta_num_params > 0]
  
  num_prot = uniqueN(feature_data$ProteinName)
  num_prot_params = num_prot - 1
  num_channels = uniqueN(feature_data$Channel)
  num_channel_params = num_channels - 1
  num_features = uniqueN(feature_data$PSM)
  num_weights = sum(weights_counts)
  
  analytical_gradient = function(parameters) {
    mu = parameters[1]
    prot_channels_params = parameters[2:(1 + num_channel_params * num_prot)]
    
    prot_channels = matrix(prot_channels_params,
                           nrow = num_prot, ncol = num_channel_params, byrow = TRUE)
    prot_channels = cbind(prot_channels, -rowSums(prot_channels))
    
    feature_params = parameters[(2 + num_channel_params * num_prot + num_weights):(length(parameters))]
    feature_params = c(feature_params, -sum(feature_params))
    
    feature_matrix = matrix(feature_params, ncol = 1) %*% matrix(1, ncol = num_channels)
    
    weights_params = parameters[(2 + num_channel_params * num_prot):(1 + num_channel_params * num_prot + num_weights)]
    weights_params_list = lapply(seq_along(weights_counts),
                                 function(i) {
                                   num_params = weights_counts[i]
                                   if (i == 1) {
                                     weights = weights_params[1:num_params]
                                   } else {
                                     weights = weights_params[(sum(weights_counts[1:(i - 1)]) + 1):((sum(weights_counts[1:(i - 1)]) + weights_counts[i]))]
                                   }
                                   weights = c(weights, 1 - sum(weights))
                                   weights
                                 })
    names(weights_params_list) = names(weights_counts)
    weight_matrix = delta_matrix
    
    for (i in seq_len(nrow(weight_matrix))) {
      if (delta_num_params[i] > 0) {
        weight_matrix[i, weight_matrix[i, ] != 0] = weights_params_list[[names(delta_num_params)[i]]]
      }
    }
    
    row.names(weight_matrix) = features
    weight_matrix
  }
}

weights_matrix = matrix(c(1, 0,
                          1, 0,
                          0, 1,
                          0, 1,
                          0.7, 0.3,
                          0.3, 0.7,
                          0.5,0.5,
                          0.6, 0.4,
                          0.4, 0.6),
                        ncol = 2, byrow = T)
colnames(weights_matrix) = paste("Protein", 1:ncol(weights_matrix), 
                                 sep = "_")

protein_values = 0.5
intercept = 10
channel_values = matrix(c(seq(0.5, 0.01, length.out = 9),
                          seq(0.5, 0.1, length.out = 9)),
                        nrow = 2, byrow = T)
feature_values = runif(nrow(weights_matrix) - 1, -1, 1)
num_channels = 10

simulated_small = simulateFeatureModelNoiseless(num_channels,
                                                weights_matrix,
                                                intercept,
                                                protein_values,
                                                channel_values,
                                                feature_values)

simulated_small_noise = addNoise(simulated_small, 0.05)

feature_data = makeLongDataFormat(simulated_small_noise, weights_matrix)

order_ch = unique(feature_data$Channel)
feature_data_plot = copy(feature_data)
feature_data_plot[, Channel := factor(Channel, levels = order_ch, ordered = T)]
feature_data_plot[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]

ggplot(feature_data_plot, 
       aes(x = Channel, y = log2IntensityNormalized, group = PSM, color = IsUnique)) +
  geom_line() +
  facet_wrap(~ProteinName) +
  theme_bw() +
  theme(legend.position = "bottom")

channel_values_sets = lapply(1:5, function(i) {
  matrix(c(seq(0.5, 0.01, length.out = 9),
           seq(0.5, c(0.3, 0.2, 0.1, 0.05, 0.01)[i], length.out = 9)),
         nrow = 2, byrow = T)
  
})

corrs = sapply(channel_values_sets, function(x) cor(x[1, ], x[2, ]))
mses = sapply(channel_values_sets, function(x) sum(((x[1, ] - mean(x[1, ])) - (x[2, ] - mean(x[2, ])))^2))

simulated_nonoise = lapply(channel_values_sets, function(set) {
  simulated_small = simulateFeatureModelNoiseless(num_channels,
                                                  weights_matrix,
                                                  intercept,
                                                  protein_values,
                                                  set,
                                                  feature_values)
})

simulated_noise = lapply(simulated_nonoise, function(x) {
  lapply(1:50, function(i) {
    simulated_small_noise = addNoise(x, 0.05)
    feature_data = makeLongDataFormat(simulated_small_noise, weights_matrix)
    feature_data
  })
})

simulated_1example = lapply(simulated_noise, function(x) x[[1]])
simulated_for_plot = lapply(simulated_1example, function(x) {
  feature_data_plot = copy(x)
  feature_data_plot[, Channel := factor(Channel, levels = order_ch, ordered = T)]
  feature_data_plot[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
  feature_data_plot  
})
simulated_for_plot = lapply(seq_along(simulated_for_plot), function(i) {
  x = simulated_for_plot[[i]]
  x$set = i
  x
})
simulated_for_plot = rbindlist(simulated_for_plot)
ggplot(simulated_for_plot[set == 1], 
       aes(x = Channel, y = log2IntensityNormalized, group = PSM, color = IsUnique)) +
  geom_line() +
  facet_grid(set~ProteinName) +
  theme_bw() +
  theme(legend.position = "bottom")

# sims_30_by_diff = lapply(simulated_noise, function(x) {
#   lapply(x, function(y) {
#     simulated_dataset = copy(y)
#     
#     feature_obs = dcast(unique(simulated_dataset[, .(PSM, Channel, log2IntensityNormalized)]),
#                         PSM ~ Channel, value.var = "log2IntensityNormalized")
#     channel_names = colnames(feature_obs)[-1]
#     feature_obs = as.matrix(feature_obs[, -1, with = FALSE])
#     feature_obs = feature_obs[, order_ch]
#     
#     delta_matrix = dcast(unique(simulated_dataset[, .(PSM, ProteinName, Present = 1)]),
#                          PSM ~ ProteinName, value.var = "Present", fill = 0)
#     features = delta_matrix[["PSM"]]
#     delta_matrix = as.matrix(delta_matrix[, -1, with = FALSE])
#     delta_sums = rowSums(delta_matrix)
#     delta_sums_shared = delta_sums > 1
#     delta_num_params = delta_sums - 1
#     names(delta_num_params) = paste0("w", seq_along(delta_num_params))
#     weights_counts = delta_num_params[delta_num_params > 0]
#     
#     prot_matrix = matrix(1, nrow = 1, ncol = uniqueN(simulated_dataset$Channel))
#     
#     num_prot = uniqueN(simulated_dataset$ProteinName)
#     num_prot_params = num_prot - 1
#     num_channels = uniqueN(simulated_dataset$Channel)
#     num_channel_params = num_channels - 1
#     num_features = uniqueN(simulated_dataset$PSM)
#     num_weights = sum(weights_counts)
#     
#     parameters = c(1, rep(0.5, num_prot * num_channel_params), 
#                    unlist(sapply(weights_counts, function(x) rep(1 / (x + 1), times = x), simplify = TRUE)),
#                    rep(1, num_features - 1))
#     
#     to_optimize = get_optim_problem(simulated_dataset, order_ch)
#     
#     lower_bounds = rep(-Inf, length(parameters))
#     upper_bounds = rep(Inf, length(parameters))
#     
#     lower_bounds[(2 + num_channel_params * num_prot):(1 + num_channel_params * num_prot + num_weights)] = 0
#     upper_bounds[(2 + num_channel_params * num_prot):(1 + num_channel_params * num_prot + num_weights)] = 1
#     
#     test_optim_sim = optim(parameters, to_optimize, control = list(maxit = 100000),
#                            method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds)
#     test_optim_sim$par
#     
#     weights_matrix_full = weights_matrix
#     
#     list(solution = test_optim_sim,
#          data = y,
#          optimized = to_optimize,
#          hessian = numDeriv::hessian(to_optimize, test_optim_sim$par))
#   })
# })

# saveRDS(sims_30_by_diff, "20230104_sims30bydiff.RDS")

sims_30_by_diff = readRDS("results/simulated_data/covariance_sims.RDS")

get_weights = getWeightsMatrix(feature_data)

weights = rbindlist(lapply(1:5, function(i) {
  # min_effect = c(0.3, 0.2, 0.1, 0.05, 0.01)[i]
  min_effect = mses[i]
  data = rbindlist(lapply(sims_30_by_diff[[i]], function(x) {
    data.table(WeightID = paste("Weight", 1:5, sep = "_"),
               Weight = x$solution$par[paste0("w", 5:9)])
  })) 
  data$min_effect = min_effect
  data
}))

weights

ggplot(weights, aes(x = reorder(as.character(min_effect), min_effect), y = Weight)) +
  geom_boxplot() +
  facet_wrap(~WeightID, ncol = 5) +
  theme_bw() +
  xlab("Minimum difference within protein profile") +
  ylab("Estimated weight")


ggplot(weights[, .(weight_var = var(Weight)),
               by = c("min_effect")], aes(x = reorder(as.character(min_effect), min_effect), y = weight_var, group = 1)) +
  geom_point() +
  geom_line() +
  # facet_wrap(~WeightID, ncol = 5) +
  theme_bw() +
  xlab("Minimum difference in protein profile") +
  ylab("Variance of weight parameter") +
  theme(legend.title = element_text(size = 24),
        axis.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        axis.text = element_text(size = 22),
        # axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 24), 
        legend.position = "bottom",
        legend.box = "vertical") 
ggplot(weights[, .(weight_var = var(Weight)),
               by = c("min_effect")], aes(x = min_effect, y = weight_var, group = 1)) +
  geom_point() +
  geom_line() +
  # facet_wrap(~WeightID, ncol = 5) +
  theme_bw() +
  xlab("Mean squared distance between protein profiles") +
  ylab("Variance of weight parameter") +
  theme(legend.title = element_text(size = 24),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 24),
        axis.text = element_text(size = 22),
        # axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 24), 
        legend.position = "bottom",
        legend.box = "vertical") 

channels = rbindlist(lapply(1:5, function(i) {
  # min_effect = c(0.3, 0.2, 0.1, 0.05, 0.01)[i]
  min_effect = mses[i]
  data = rbindlist(lapply(sims_30_by_diff[[i]], function(x) {
    data.table(ProteinName = rep(c("Protein_1", "Protein_2"), each = 9),
               ChannelID = rep(paste("Channel", 1:9, sep = "_"), times = 2),
               Channel = x$solution$par[2:19])
  })) 
  data$min_effect = min_effect
  data
}))

ggplot(channels, aes(y = Channel, x = as.character(min_effect), fill = ProteinName)) +
  geom_boxplot() +
  theme_bw()
ggplot(channels[, .(Channel = sd(Channel)), 
                by = c("ProteinName", "min_effect")],
       aes(y = Channel, x = min_effect, color = ProteinName, group = ProteinName)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom")


ggplot(channels[, .(sd = sd(Channel, na.rm = T)),
                by = c("ProteinName", "min_effect", "ChannelID")], 
       aes(x = reorder(as.character(min_effect), min_effect), y = sd, color = ProteinName)) +
  geom_point() +
  facet_grid(~ChannelID) +
  theme_bw() +
  xlab("Minimum difference within protein profile") +
  theme(legend.position = "bottom")
ggplot(channels[, .(Channel = Channel - mean(Channel)), 
                by = c("min_effect", "ChannelID", "ProteinName")],
       aes(x = reorder(as.character(min_effect), min_effect), y = Channel, fill = ProteinName)) +
  geom_boxplot() +
  facet_grid(~ChannelID) +
  theme_bw() +
  xlab("Minimum difference within protein profile") +
  theme(legend.position = "bottom")
ggplot(channels[, .(Channel = Channel - mean(Channel)), 
                by = c("min_effect", "ChannelID")],
       aes(x = reorder(as.character(min_effect), min_effect), y = Channel)) +
  geom_boxplot() +
  facet_grid(~ChannelID) +
  theme_bw() +
  xlab("Minimum difference within protein profile") +
  theme(legend.position = "bottom")

ggplot(channels[, .(Channel = Channel - mean(Channel)), 
                by = c("min_effect", "ChannelID", "ProteinName")], aes(x = reorder(as.character(min_effect), min_effect), y = Channel)) +
  geom_boxplot() +
  facet_grid(ProteinName~ChannelID) +
  theme_bw() +
  xlab("Minimum difference within protein profile")
ggplot(channels, aes(x = reorder(as.character(min_effect), min_effect), y = Channel)) +
  geom_boxplot() +
  facet_grid(~ChannelID) +
  theme_bw() +
  xlab("Minimum difference within protein profile")


channels[, .(mean = mean(Channel),
             var = var(Channel)),
         by = c("ProteinName", "ChannelID", "min_effect")]

ggplot(channels[, .(mean = mean(Channel),
                    var = var(Channel)),
                by = c("ProteinName", "ChannelID", "min_effect")][, .(var = mean(var)), by = "min_effect"],
       aes(x = min_effect, y = var)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  ylab("parameter variance") +
  xlab("minimum difference in protein profile") +
  theme(legend.title = element_text(size = 24),
        axis.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        axis.text = element_text(size = 22),
        # axis.text.x = element_text(angle = 270),
        strip.text = element_text(size = 24), 
        legend.position = "bottom",
        legend.box = "vertical") 
