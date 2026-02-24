library(ggplot2)
library(ggpubr)
library(randomcoloR)
# load dataframe using load_output_xlsx_as_dataframe.R

# get max and mean for vector measurements - capillary
for (i in 1:96) {
  df$avg_tortuosity[i] = mean(df$tortuosity[[i]])
  df$max_tortuosity[i] = max(df$tortuosity[[i]])
  df$avg_wtort[i] = mean(df$weighted_tortuosity[[i]])
  df$max_wtort[i] = max(df$weighted_tortuosity[[i]])
  df$avg_seg_leng[i] = mean(df$segment_length[[i]])
  df$max_seg_leng[i] = max(df$segment_length[[i]])
}
temp = data.frame(df$avg_tortuosity,
                  df$max_tortuosity,
                  df$avg_wtort, 
                  df$max_wtort,
                  df$avg_seg_leng,
                  df$max_seg_leng)

# get max and mean for vector measurements - neuron
for (i in 1:24) {
  df_neuron$avg_nvdist[i] = mean(df_neuron$NV_distance[[i]])
  df_neuron$max_nvdist[i] = max(df_neuron$NV_distance[[i]])
  df_neuron$avg_nndist[i] = mean(df_neuron$NN_distance[[i]])
  df_neuron$max_nndist[i] = max(df_neuron$NN_distance[[i]])
}
temp = data.frame(df_neuron$avg_nvdist,
                  df_neuron$max_nvdist,
                  df_neuron$avg_nndist,
                  df_neuron$max_nndist)

# get subject-specific mean path scores
df$mean_path_score = NA
for(i in 0:23) {
  score = sum(df[((4*i)+1):((4*i)+4),6:8])
  df$mean_path_score[((4*i)+1):((4*i)+4)] = score # H&E, tau, abeta for all 4 CA subregions
}
df$mean_path_score_CA = NA
df$mean_path_score_CA = (df$H_E + df$tau + df$abeta) / 3


