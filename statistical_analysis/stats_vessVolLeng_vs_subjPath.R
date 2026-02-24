library(ggplot2)
library(ggpubr)
library(randomcoloR)
# load df: outputs_numerical_NEWNAME
# load df2: XtraData_total_vessel_voxel

#####
# fisher transformation and plotting
fisher_transformation <- function(r, n) { # 95% confidence interval
  z = log((1+r)/(1-r)) / 2
  l = z - (1.96 / sqrt(n-3))
  u = z + (1.96 / sqrt(n-3))
  CI_l = (exp(2*l)-1) / (exp(2*l)+1)
  CI_u = (exp(2*u)-1) / (exp(2*u)+1)
  out <- data.frame(z, CI_l, CI_u)
  return(out)
}
# separate AD from non-AD
df2_low = df2[df2$lowAD_selection=="y",]
df2_high = df2[df2$lowAD_selection=="n",]
# get CA1
df_CA1 = df2[df2$CA=="CA1",]
df_CA1_low = df2_low[df2_low$CA=="CA1",]
df_CA1_high = df2_high[df2_high$CA=="CA1",]


##### 
# ALL SUBJECTS ANALYZED TOGETHER
#####
# subject NVV and mean path score
stats_pear <- cor.test(df2$mean_path_score, df2$sum_vess_vol/df2$sum_image_vol, na.action="na.omit")
stats_fish <- fisher_transformation(stats_pear$estimate, 24)
# CA1 NVV and mean path score
stats_pear <- cor.test(df_CA1$mean_path_score, df_CA1$normalized_vess_vol, na.action="na.omit")
stats_fish <- fisher_transformation(stats_pear$estimate, 24)
# CA1 NVV and CA1 mean path score
stats_pear <- cor.test(df_CA1$CA1_mean_path, df_CA1$normalized_vess_vol, na.action="na.omit")
stats_fish <- fisher_transformation(stats_pear$estimate, 24)

# subject average vess length and mean path score
stats_pear <- cor.test(df2$mean_path_score, df2$sum_vess_leng, na.action="na.omit")
stats_fish <- fisher_transformation(stats_pear$estimate, 24)
# CA1 vess length and mean path score
stats_pear <- cor.test(df_CA1$mean_path_score, df_CA1$normalized_vess_leng, na.action="na.omit")
stats_fish <- fisher_transformation(stats_pear$estimate, 24)


######
# NON AD
######
# subj vol & sub mean path score
stats_pear <- cor.test(df2_low$mean_path_score, df2_low$sum_vess_vol/df2_low$sum_image_vol, na.action="na.omit")
stats_fish <- fisher_transformation(stats_pear$estimate, 11)
# CA1 vol & sub mean path score
stats_pear <- cor.test(df_CA1_low$mean_path_score, df_CA1_low$normalized_vess_vol, na.action="na.omit")
stats_fish <- fisher_transformation(stats_pear$estimate, 11)
# subj length & subj mean path score
stats_pear <- cor.test(df2_low$mean_path_score, df2_low$sum_vess_leng, na.action="na.omit")
stats_fish <- fisher_transformation(stats_pear$estimate, 11)
# CA1 vess length and mean path score
stats_pear <- cor.test(df_CA1_low$mean_path_score, df_CA1_low$normalized_vess_leng, na.action="na.omit")
stats_fish <- fisher_transformation(stats_pear$estimate, 11)


######
# AD
######
# subj vol & sub mean path score
stats_pear <- cor.test(df2_high$mean_path_score, df2_high$sum_vess_vol/df2_high$sum_image_vol, na.action="na.omit")
stats_fish <- fisher_transformation(stats_pear$estimate, 13)
# CA1 vol & sub mean path score
stats_pear <- cor.test(df_CA1_high$mean_path_score, df_CA1_high$normalized_vess_vol, na.action="na.omit")
stats_fish <- fisher_transformation(stats_pear$estimate, 13)
# subj length & subj mean path score
stats_pear <- cor.test(df2_high$mean_path_score, df2_high$sum_vess_leng, na.action="na.omit")
stats_fish <- fisher_transformation(stats_pear$estimate, 13)
# CA1 vess length and mean path score
stats_pear <- cor.test(df_CA1_high$mean_path_score, df_CA1_high$normalized_vess_leng, na.action="na.omit")
stats_fish <- fisher_transformation(stats_pear$estimate, 13)



