# laod libary
library(readxl)
folder_name = "/Users/yunonnebai/NauenLab/2023_data_analysis/dataset/"
file_name = "outputs_numerical_NEWNAME.xlsx"
file_name2 = "XtraData_total_vessel_voxel.xlsx"
dataset_file_path = paste(folder_name, file_name, sep="")
dataset_file_path2 = paste(folder_name, file_name2, sep="")

# load dataset numerical outputs
df = read_xlsx(dataset_file_path)
df2 = read_xlsx(dataset_file_path2)
df_neuron = df[df$NeuN_selection=='y',]
df_neuron = df_neuron[!is.na(df_neuron$research_num),]

# load dataset vector outputs
df['tortuosity'] <- NA
df_tort = read_xlsx(dataset_file_path, sheet=2)
for(i in 1:96) {
  curr = as.vector(t(df_tort[i,-(1:4)]))
  curr = list(curr[!is.na(curr)])
  df$tortuosity[i] = curr
}
df['weighted_tortuosity'] <- NA
df_wtort = read_xlsx(dataset_file_path, sheet=3)
for(i in 1:96) {
  curr = as.vector(t(df_wtort[i,-(1:4)]))
  curr = list(curr[!is.na(curr)])
  df$weighted_tortuosity[i] = curr
}
df['segment_length'] <- NA
df_sl = read_xlsx(dataset_file_path, sheet=4)
for(i in 1:96) {
  curr = as.vector(t(df_sl[i,-(1:4)]))
  curr = list(curr[!is.na(curr)])
  df$segment_length[i] = curr
}
df_neuron['NV_distance'] <- NA
df_nvdist = read_xlsx(dataset_file_path, sheet=5)
for(i in 1:24) {
  curr = as.vector(t(df_nvdist[i,-(1:4)]))
  curr = list(curr[!is.na(curr)])
  df_neuron$NV_distance[i] = curr
}
df_neuron['NN_distance'] <- NA
df_nndist = read_xlsx(dataset_file_path, sheet=6)
for(i in 1:24) {
  curr = as.vector(t(df_nndist[i,-(1:4)]))
  curr = list(curr[!is.na(curr)])
  df_neuron$NN_distance[i] = curr
}

