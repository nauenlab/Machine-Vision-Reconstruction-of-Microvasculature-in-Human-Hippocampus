library(ggplot2)
library(ggpubr)
library(randomcoloR)
# load dataframe using load_output_xlsx_as_dataframe.R
df_low = df[df$lowAD_selection=="y",]
df_high = df[df$lowAD_selection=="n",]

# prep: calculate RELATIVE values (capillary measurements)
# rVALUE = CA_value / subj avg across all 4 CA
get_subj_relative <- function(df){
  df$rH_E <- NA
  df$rTau <- NA
  df$rAbeta <- NA
  df$rMean_path <- NA
  df$rNVV <- NA
  df$rNVL <- NA
  df$rAVD <- NA
  df$rNPV <- NA
  df$rNPL <- NA
  df$rSPL <- NA
  df$rTavg <- NA
  df$rTmax <- NA
  df$rWTavg <- NA
  df$rWTmax <- NA
  df$rSLavg <- NA
  df$rSLmax <- NA
  for (i in 1:(nrow(df)/4)){
    for (j in 6:21){
      curr_avg = mean(as.numeric(unlist(df[(((i-1)*4)+1):(i*4),j])), na.rm=TRUE)
      df[(((i-1)*4)+1),j+27] = df[(((i-1)*4)+1),j] / curr_avg
      df[(((i-1)*4)+2),j+27] = df[(((i-1)*4)+2),j] / curr_avg
      df[(((i-1)*4)+3),j+27] = df[(((i-1)*4)+3),j] / curr_avg
      df[(i*4),j+27] = df[(i*4),j] / curr_avg
    }
  }
  return(df)
}
df = get_subj_relative(df)
df_low = get_subj_relative(df_low)
df_high = get_subj_relative(df_high)

# prep: calculate RELATIVE values (neuron measurements)
# rVALUE = CA_value / subj avg across all 4 CA
df_neuron$rNeuronD <- NA
df_neuron$rNVdistavg <- NA
df_neuron$rNVdistmax <- NA
df_neuron$rNNdistavg <- NA
df_neuron$rNNdistmax <- NA
df_neuron$rwithNVdist <- NA
for (i in 1:(nrow(df_neuron)/4)){
  for (j in 22:27){
    curr_avg = mean(as.numeric(unlist(df_neuron[(((i-1)*4)+1):(i*4),j])), na.rm=TRUE)
    df_neuron[(((i-1)*4)+1),j+11] = df_neuron[(((i-1)*4)+1),j] / curr_avg
    df_neuron[(((i-1)*4)+2),j+11] = df_neuron[(((i-1)*4)+2),j] / curr_avg
    df_neuron[(((i-1)*4)+3),j+11] = df_neuron[(((i-1)*4)+3),j] / curr_avg
    df_neuron[(i*4),j+11] = df_neuron[(i*4),j] / curr_avg
  }
}

# Variation in different CA regions
# stat: T test, two sided p value, exclude NA values
loc_analysis <- function(col_name, df_curr) {
  df1 = df_curr[df_curr$region=="CA1",]
  df2 = df_curr[df_curr$region=="CA2",]
  df3 = df_curr[df_curr$region=="CA3",]
  df4 = df_curr[df_curr$region=="CA4",]
  print(t.test(df1[,col_name], df2[,col_name], na.action=na.exclude))
  print(t.test(df1[,col_name], df3[,col_name], na.action=na.exclude))
  print(t.test(df1[,col_name], df4[,col_name], na.action=na.exclude))
  print(t.test(df2[,col_name], df3[,col_name], na.action=na.exclude))
  print(t.test(df2[,col_name], df4[,col_name], na.action=na.exclude))
  print(t.test(df3[,col_name], df4[,col_name], na.action=na.exclude))
}
# pool from low AD subjects (relative values)
loc_analysis("rH_E", df_low)
loc_analysis("rTau", df_low)
loc_analysis("rAbeta", df_low)
loc_analysis("rNVV", df_low)
loc_analysis("rNVL", df_low)
loc_analysis("rAVD", df_low)
loc_analysis("rNPV", df_low)
loc_analysis("rNPL", df_low)
loc_analysis("rWTavg", df_low)
loc_analysis("rWTmax", df_low)
loc_analysis("rSLavg", df_low)
loc_analysis("rSLmax", df_low)
# pool from high AD subjects (relative values)
loc_analysis("rH_E", df_high)
loc_analysis("rTau", df_high)
loc_analysis("rAbeta", df_high)
loc_analysis("rNVV", df_high)
loc_analysis("rNVL", df_high)
loc_analysis("rAVD", df_high)
loc_analysis("rNPV", df_high)
loc_analysis("rNPL", df_high)
loc_analysis("rWTavg", df_high)
loc_analysis("rWTmax", df_high)
loc_analysis("rSLavg", df_high)
loc_analysis("rSLmax", df_high)
# pool from all subjects (relative values)
loc_analysis("rTau", df)
loc_analysis("rH_E", df)
loc_analysis("rAbeta", df)
# pool from subjects with neuron analysis (relative values)
loc_analysis("rNeuronD", df_neuron)
loc_analysis("rNVdistavg", df_neuron)
loc_analysis("rNVdistmax", df_neuron)
loc_analysis("rNNdistavg", df_neuron)
loc_analysis("rNNdistmax", df_neuron)
loc_analysis("rwithNVdist", df_neuron)

