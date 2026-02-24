library(ggplot2)
library(ggpubr)
library(randomcoloR)
# load dataframe using load_output_xlsx_as_dataframe.R

# quality assessment: % neurons with NVdist>40um
df_neuron$withNVdist <- NA
for(i in 1:24){
  df_neuron$withNVdist[i] = length(df_neuron$NV_distance[[i]]) / length(df_neuron$NN_distance[[i]])
}
cor.test(df_neuron$withNVdist, df_neuron$mean_path_score_CA, method="kendall", nan.action=na.exclude)
cor.test(df_neuron$withNVdist, df_neuron$age, method="pearson", nan.action=na.exclude)

# relationship between path scores and NV features
# stats: pearson correlation test, two sided p value, exclude NA values
neuron_vessel_correlation <- function(df_curr){
  print(cor.test(df_curr$neuron_density, df_curr$normalized_vess_vol, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$avg_nvdist, df_curr$normalized_vess_vol, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$avg_nndist, df_curr$normalized_vess_vol, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$max_nndist, df_curr$normalized_vess_vol, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$withNVdist, df_curr$normalized_vess_vol, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$neuron_density, df_curr$normalized_vess_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$avg_nvdist, df_curr$normalized_vess_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$avg_nndist, df_curr$normalized_vess_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$max_nndist, df_curr$normalized_vess_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$withNVdist, df_curr$normalized_vess_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$neuron_density, df_curr$node_per_vol, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$avg_nvdist, df_curr$node_per_vol, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$avg_nndist, df_curr$node_per_vol, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$max_nndist, df_curr$node_per_vol, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$withNVdist, df_curr$node_per_vol, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$neuron_density, df_curr$node_per_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$avg_nvdist, df_curr$node_per_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$avg_nndist, df_curr$node_per_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$max_nndist, df_curr$node_per_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$withNVdist, df_curr$node_per_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$neuron_density, df_curr$avg_tortuosity, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$avg_nvdist, df_curr$avg_tortuosity, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$avg_nndist, df_curr$avg_tortuosity, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$max_nndist, df_curr$avg_tortuosity, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$withNVdist, df_curr$avg_tortuosity, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$neuron_density, df_curr$avg_seg_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$avg_nvdist, df_curr$avg_seg_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$avg_nndist, df_curr$avg_seg_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$max_nndist, df_curr$avg_seg_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$withNVdist, df_curr$avg_seg_leng, method="pearson", na.action=na.exclude))
}

neuron_vessel_correlation <- function(df_curr){
  print(cor.test(df_curr$withNVdist, df_curr$normalized_vess_vol, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$withNVdist, df_curr$normalized_vess_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$withNVdist, df_curr$node_per_vol, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$withNVdist, df_curr$node_per_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$withNVdist, df_curr$avg_tortuosity, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$withNVdist, df_curr$avg_seg_leng, method="pearson", na.action=na.exclude))
  
}
# pooling from all subjects
neuron_vessel_correlation(df) # pooling from all CA samples
neuron_vessel_correlation(df[df$region=="CA1",]) # pooling from CA1 only
neuron_vessel_correlation(df[df$region=="CA2",]) # pooling from CA2 only
neuron_vessel_correlation(df[df$region=="CA3",]) # pooling from CA3 only
neuron_vessel_correlation(df[df$region=="CA4",]) # pooling from CA4 only
