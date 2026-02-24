library(ggplot2)
library(ggpubr)
library(randomcoloR)
# load dataframe using load_output_xlsx_as_dataframe.R
df_low = df[df$lowAD_selection=="y",]
df_high = df[df$lowAD_selection=="n",]

# comparing path scores and NV params in low VS high AD subjects
low_vs_high_path_analysis <- function(col_name, low, high) {
  low1 = low[low$region=="CA1",]
  low2 = low[low$region=="CA2",]
  low3 = low[low$region=="CA3",]
  low4 = low[low$region=="CA4",]
  high1 = high[high$region=="CA1",]
  high2 = high[high$region=="CA2",]
  high3 = high[high$region=="CA3",]
  high4 = high[high$region=="CA4",]
  print(wilcox.test(as.numeric(unlist(low[,col_name])), as.numeric(unlist(high[,col_name])), na.action=na.exclude))
  print(wilcox.test(as.numeric(unlist(low1[,col_name])), as.numeric(unlist(high1[,col_name])), na.action=na.exclude))
  print(wilcox.test(as.numeric(unlist(low2[,col_name])), as.numeric(unlist(high2[,col_name])), na.action=na.exclude))
  print(wilcox.test(as.numeric(unlist(low3[,col_name])), as.numeric(unlist(high3[,col_name])), na.action=na.exclude))
  print(wilcox.test(as.numeric(unlist(low4[,col_name])), as.numeric(unlist(high4[,col_name])), na.action=na.exclude))
}
low_vs_high_path_analysis("H_E", df_low, df_high)
low_vs_high_path_analysis("tau", df_low, df_high)
low_vs_high_path_analysis("abeta", df_low, df_high)

low_vs_high_nv_analysis <- function(col_name, low, high) {
  low1 = low[low$region=="CA1",]
  low2 = low[low$region=="CA2",]
  low3 = low[low$region=="CA3",]
  low4 = low[low$region=="CA4",]
  high1 = high[high$region=="CA1",]
  high2 = high[high$region=="CA2",]
  high3 = high[high$region=="CA3",]
  high4 = high[high$region=="CA4",]
  print(t.test(low[,col_name], high[,col_name], na.action=na.exclude))
  print(t.test(low1[,col_name], high1[,col_name], na.action=na.exclude))
  print(t.test(low2[,col_name], high2[,col_name], na.action=na.exclude))
  print(t.test(low3[,col_name], high3[,col_name], na.action=na.exclude))
  print(t.test(low4[,col_name], high4[,col_name], na.action=na.exclude))
}
low_vs_high_nv_analysis("normalized_vess_vol", df_low, df_high)
low_vs_high_nv_analysis("normalized_vess_leng", df_low, df_high)
low_vs_high_nv_analysis("avg_vess_diameter", df_low, df_high)
low_vs_high_nv_analysis("node_per_vol", df_low, df_high)
low_vs_high_nv_analysis("node_per_leng", df_low, df_high)
low_vs_high_nv_analysis("avg_wtort", df_low, df_high)
low_vs_high_nv_analysis("max_wtort", df_low, df_high)
low_vs_high_nv_analysis("avg_seg_leng", df_low, df_high)
low_vs_high_nv_analysis("max_seg_leng", df_low, df_high)
