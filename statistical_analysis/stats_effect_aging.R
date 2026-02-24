library(ggplot2)
library(ggpubr)
library(randomcoloR)
# load dataframe using load_output_xlsx_as_dataframe.R
df_low = df[df$lowAD_selection=="y",]
df_high = df[df$lowAD_selection=="n",]

# relationship between path scores and age
# stats: Kendall rank correlation test, two sided p value, exclude NA values
path_age_analysis <- function(df_curr){
  print(cor.test(df_curr$H_E, df_curr$age, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$tau, df_curr$age, method="kendall", na.action=na.exclude))
  print(cor.test(df_curr$abeta, df_curr$age, method="kendall", na.action=na.exclude))
}
path_age_analysis(df_low) # pooling from all CA samples
path_age_analysis(df_low[df_low$region=="CA1",]) # pooling from CA1 only
path_age_analysis(df_low[df_low$region=="CA2",]) # pooling from CA2 only
path_age_analysis(df_low[df_low$region=="CA3",]) # pooling from CA3 only
path_age_analysis(df_low[df_low$region=="CA4",]) # pooling from CA4 only


# relationship between age and NV features
# stats: Pearson correlation test, two sided p value, exclude NA values
age_nv_analysis <- function(df_curr){
  print(cor.test(df_curr$age, df_curr$normalized_vess_vol, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$age, df_curr$normalized_vess_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$age, df_curr$avg_vess_diameter, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$age, df_curr$node_per_vol, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$age, df_curr$node_per_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$age, df_curr$avg_wtort, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$age, df_curr$max_wtort, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$age, df_curr$avg_seg_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$age, df_curr$max_seg_leng, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$age, df_curr$neuron_density, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$age, df_curr$avg_nvdist, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$age, df_curr$max_nvdist, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$age, df_curr$avg_nndist, method="pearson", na.action=na.exclude))
  print(cor.test(df_curr$age, df_curr$max_nndist, method="pearson", na.action=na.exclude))
}
# pool from low AD subjects
age_nv_analysis(df_low) # pooling from all CA samples
age_nv_analysis(df_low[df_low$region=="CA1",]) # pooling from CA1 only
age_nv_analysis(df_low[df_low$region=="CA2",]) # pooling from CA2 only
age_nv_analysis(df_low[df_low$region=="CA3",]) # pooling from CA3 only
age_nv_analysis(df_low[df_low$region=="CA4",]) # pooling from CA4 only
# neurons
age_nv_analysis(df_neuron)
age_nv_analysis(df_neuron[df_neuron$region=="CA1",])
age_nv_analysis(df_neuron[df_neuron$region=="CA2",])
age_nv_analysis(df_neuron[df_neuron$region=="CA3",])
age_nv_analysis(df_neuron[df_neuron$region=="CA4",])
