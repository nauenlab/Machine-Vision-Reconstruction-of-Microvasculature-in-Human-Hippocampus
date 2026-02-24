library(ggplot2)
library(ggpubr)
library(randomcoloR)
# load dataframe using load_output_xlsx_as_dataframe.R
df_low = df[df$lowAD_selection=="y",]
df_high = df[df$lowAD_selection=="n",]

# relationship between path scores and NV features
# stats: Kendall rank correlation test, two sided p value, exclude NA values
path_nv_analysis <- function(df_curr){
  print(cor.test(df_curr$H_E, df_curr$normalized_vess_vol, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$H_E, df_curr$normalized_vess_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$H_E, df_curr$avg_vess_diameter, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$H_E, df_curr$node_per_vol, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$H_E, df_curr$node_per_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$H_E, df_curr$avg_wtort, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$H_E, df_curr$max_wtort, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$H_E, df_curr$avg_seg_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$H_E, df_curr$max_seg_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$normalized_vess_vol, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$normalized_vess_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$avg_vess_diameter, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$node_per_vol, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$node_per_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$avg_wtort, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$max_wtort, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$avg_seg_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$max_seg_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$normalized_vess_vol, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$normalized_vess_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$avg_vess_diameter, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$node_per_vol, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$node_per_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$avg_wtort, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$max_wtort, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$avg_seg_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$max_seg_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$normalized_vess_vol, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$normalized_vess_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$avg_vess_diameter, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$node_per_vol, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$node_per_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$avg_wtort, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$max_wtort, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$avg_seg_leng, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$max_seg_leng, method="kendall", na.action=na.exclude))
}
# pooling from all subjects
path_nv_analysis(df) # pooling from all CA samples
path_nv_analysis(df[df$region=="CA1",]) # pooling from CA1 only
path_nv_analysis(df[df$region=="CA2",]) # pooling from CA2 only
path_nv_analysis(df[df$region=="CA3",]) # pooling from CA3 only
path_nv_analysis(df[df$region=="CA4",]) # pooling from CA4 only
# pooling from only low AD subjects
path_nv_analysis(df_low) # pooling from all CA samples
path_nv_analysis(df_low[df_low$region=="CA1",]) # pooling from CA1 only
path_nv_analysis(df_low[df_low$region=="CA2",]) # pooling from CA2 only
path_nv_analysis(df_low[df_low$region=="CA3",]) # pooling from CA3 only
path_nv_analysis(df_low[df_low$region=="CA4",]) # pooling from CA4 only
# pooling from only high AD subjects
path_nv_analysis(df_high) # pooling from all CA samples
path_nv_analysis(df_high[df_high$region=="CA1",]) # pooling from CA1 only
path_nv_analysis(df_high[df_high$region=="CA2",]) # pooling from CA2 only
path_nv_analysis(df_high[df_high$region=="CA3",]) # pooling from CA3 only
path_nv_analysis(df_high[df_high$region=="CA4",]) # pooling from CA4 only


# relationship between path scores and NV features (neuron analysis only)
# stats: Kendall rank correlation test, two sided p value, exclude NA values
path_neuron_analysis <- function(df_curr){
  print(cor.test(df_curr$H_E, df_curr$neuron_density, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$H_E, df_curr$avg_nvdist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$H_E, df_curr$max_nvdist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$H_E, df_curr$avg_nndist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$H_E, df_curr$max_nndist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$H_E, df_curr$withNVdist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$neuron_density, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$avg_nvdist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$max_nvdist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$avg_nndist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$max_nndist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$withNVdist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$neuron_density, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$avg_nvdist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$max_nvdist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$avg_nndist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$max_nndist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$withNVdist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$neuron_density, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$avg_nvdist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$max_nvdist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$avg_nndist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$max_nndist, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$mean_path_score_CA, df_curr$withNVdist, method="kendall", na.action=na.exclude))
}
path_neuron_analysis(df_neuron) # pooling from all CA samples
path_neuron_analysis(df_neuron[df_neuron$region=="CA1",]) # pooling from CA1 only
path_neuron_analysis(df_neuron[df_neuron$region=="CA2",]) # pooling from CA2 only
path_neuron_analysis(df_neuron[df_neuron$region=="CA3",]) # pooling from CA3 only
path_neuron_analysis(df_neuron[df_neuron$region=="CA4",]) # pooling from CA4 only
