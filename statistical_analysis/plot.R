library(ggplot2)
library(ggpubr)
library(randomcoloR)
# load dataframe using load_dataset_as_dataframe.R
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
      df[(((i-1)*4)+1),j+28] = df[(((i-1)*4)+1),j] / curr_avg
      df[(((i-1)*4)+2),j+28] = df[(((i-1)*4)+2),j] / curr_avg
      df[(((i-1)*4)+3),j+28] = df[(((i-1)*4)+3),j] / curr_avg
      df[(i*4),j+28] = df[(i*4),j] / curr_avg
    }
  }
  return(df)
}
df = get_subj_relative(df)
df_low = get_subj_relative(df_low)
df_high = get_subj_relative(df_high)



#############
# plot 1: vessel branching
#############
#fig1a-d
df_sl = rbind(data.frame(sl=c(unlist(df_low[df_low$region=="CA1",]$segment_length)),
                         group="Non-AD", region="CA1"),
              data.frame(sl=c(unlist(df_low[df_low$region=="CA2",]$segment_length)),
                         group="Non-AD", region="CA2"),
              data.frame(sl=c(unlist(df_low[df_low$region=="CA3",]$segment_length)),
                         group="Non-AD", region="CA3"),
              data.frame(sl=c(unlist(df_low[df_low$region=="CA4",]$segment_length)),
                         group="Non-AD", region="CA4"),
              data.frame(sl=c(unlist(df_high[df_high$region=="CA1",]$segment_length)),
                         group="AD", region="CA1"),
              data.frame(sl=c(unlist(df_high[df_high$region=="CA2",]$segment_length)),
                         group="AD", region="CA2"),
              data.frame(sl=c(unlist(df_high[df_high$region=="CA3",]$segment_length)),
                         group="AD", region="CA3"),
              data.frame(sl=c(unlist(df_high[df_high$region=="CA4",]$segment_length)),
                         group="AD", region="CA4"))
ggplot(df_sl[df_sl$region=="CA1",], aes(sl, color=group)) + 
  stat_ecdf() + theme_bw() + theme(aspect.ratio=1) +
  xlab("Segment length (um)") + ylab("eCDF") + ggtitle("CA1") +
  labs(color="Group") + scale_color_manual(values=c("forestgreen", "magenta")) +
  xlim(0,500)
ggplot(df_sl[df_sl$region=="CA2",], aes(sl, color=group)) + 
  stat_ecdf() + theme_bw() + theme(aspect.ratio=1) +
  xlab("Segment length (um)") + ylab("eCDF") + ggtitle("CA2") +
  labs(color="Group") + scale_color_manual(values=c("forestgreen", "magenta")) +
  xlim(0,500)
ggplot(df_sl[df_sl$region=="CA3",], aes(sl, color=group)) + 
  stat_ecdf() + theme_bw() + theme(aspect.ratio=1) +
  xlab("Segment length (um)") + ylab("eCDF") + ggtitle("CA3") +
  labs(color="Group") + scale_color_manual(values=c("forestgreen", "magenta")) +
  xlim(0,500)
ggplot(df_sl[df_sl$region=="CA4",], aes(sl, color=group)) + 
  stat_ecdf() + theme_bw() + theme(aspect.ratio=1) +
  xlab("Segment length (um)") + ylab("eCDF") + ggtitle("CA4") +
  labs(color="Group") + scale_color_manual(values=c("forestgreen", "magenta")) +
  xlim(0,500)
#1e
ggplot(df, aes(x=mean_path_score_CA, y=node_per_vol)) +
  geom_jitter(width=0.03) +
  theme_bw() + ggtitle("(p>0.1)") + theme(aspect.ratio=1) +
  xlab("Mean pathology score") + ylab("Branch point density (1/um^3)") +
  scale_y_continuous(trans='log10')
#1f: No capillary atrophy/change in vessel diameter
ggplot(df, aes(x=mean_path_score_CA, y=avg_vess_diameter)) +
  geom_jitter(width=0.03) +
  theme_bw() + ggtitle("(p>0.1)") + theme(aspect.ratio=1) +
  xlab("Mean pathology score") + ylab("Average vessel diameter (um)")



###############
# plot 2: vessel volume
###############
#2a
ggplot(data.frame(name=c("Non-AD", "AD"),
                             mean=c(0.01662066, 0.0185058)*100,
                             std=c(0.004080402, 0.004236102)*100),
                  aes(x=reorder(name,+mean), y=mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=0.2, position=position_dodge(.9)) +
  geom_signif(stat="identity", 
              data=data.frame(x=1,xend=2, y=2.4, yend=2.4, annotation="*"),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  theme_bw() + ylab("Vessel volume [%]") + xlab("") + ylim(0.0, 2.5) + theme(aspect.ratio = 3)
#2b
ggplot(df, aes(y=normalized_vess_vol*100, x=tau)) +
  geom_point() +
  geom_smooth(method="lm", formula=y~x) +
  scale_color_manual(values=colors) +
  theme_bw() + xlab("Tau") + xlim(0,3) + ylab("Vessel volume [%]") +
  theme(aspect.ratio = 1) +
  ggtitle("T=0.176 (p=0.028)")
#2c
ggplot(df, aes(y=normalized_vess_vol*100, x=mean_path_score_CA)) +
  geom_point() +
  geom_smooth(method="lm", formula=y~x) +
  scale_color_manual(values=colors) +
  theme_bw() + xlab("Mean pathology score") + xlim(0,3) + ylab("Vessel volume [%]") +
  theme(aspect.ratio = 1) +
  ggtitle("T=0.174 (p=0.025)")
#2d
ggplot(df2, aes(x=mean_path_score, y=sum_vess_vol/sum_image_vol*100)) + 
  geom_point() +
  geom_smooth(method="lm", formula=y~x) +
  theme_bw() + xlab("Subject mean path score") + ylab("Subject vessel volume [%]") +
  ylim(1.0, 2.9) + theme(aspect.ratio = 1) +
  ggtitle("r=0.442 (p=0.031)")
#2e
df_CA1 = df2[df2$CA=="CA1",]
ggplot(df_CA1, aes(x=mean_path_score, y=normalized_vess_vol*100)) +
  geom_point() +
  geom_smooth(method="lm", formula=y~x) +
  theme_bw() + xlab("Subject mean path score") + ylab("CA1 vessel volume [%]") +
  ylim(1.0, 2.9) + theme(aspect.ratio = 1) +
  ggtitle("r=0.529 (p=0.0078)")



###############
# plot 3: tau
###############
#3a: regional variation
df_plot = data.frame(matrix(ncol=3, nrow=8))
colnames(df_plot) <- c("region", "mean", "std")
df_plot$region = rep(c("CA1", "CA2", "CA3", "CA4"), 2)
make_mean_std_df <- function(df_plot, df_data, row_num, col_num) {
  df_plot$mean[row_num] = mean(as.numeric(unlist(df_data[df_data$region=="CA1",col_num])), na.rm=TRUE)
  df_plot$mean[row_num+1] = mean(as.numeric(unlist(df_data[df_data$region=="CA2",col_num])), na.rm=TRUE)
  df_plot$mean[row_num+2] = mean(as.numeric(unlist(df_data[df_data$region=="CA3",col_num])), na.rm=TRUE)
  df_plot$mean[row_num+3] = mean(as.numeric(unlist(df_data[df_data$region=="CA4",col_num])), na.rm=TRUE)
  df_plot$std[row_num] = sd(as.numeric(unlist(df_data[df_data$region=="CA1",col_num])), na.rm=TRUE)
  df_plot$std[row_num+1] = sd(as.numeric(unlist(df_data[df_data$region=="CA2",col_num])), na.rm=TRUE)
  df_plot$std[row_num+2] = sd(as.numeric(unlist(df_data[df_data$region=="CA3",col_num])), na.rm=TRUE)
  df_plot$std[row_num+3] = sd(as.numeric(unlist(df_data[df_data$region=="CA4",col_num])), na.rm=TRUE)
  return(df_plot)
}
df_plot = make_mean_std_df(df_plot, df, row_num=1, col_num=35)  #relative tau
colors <- c("CA1"="red", "CA2"="gold", "CA3"="cyan", "CA4"="blue")
ggplot(df_plot, aes(x=factor(region), y=mean)) +
  geom_bar(aes(fill=region), stat="identity", position=position_dodge()) +
  geom_errorbar(aes(fill=region, ymin=mean-std, ymax=mean+std), width=0.2, position=position_dodge(.9)) +
  scale_fill_manual(values=colors) +
  theme_bw() + xlab("Subregion") + ylab("Relative tau") + ylim(0,1.8) + 
  theme(aspect.ratio=2) + theme(legend.position="none") +
  geom_signif(stat="identity",
              data=data.frame(x=c(1, 1, 1), 
                              xend=c(2, 3, 4),
                              y=c(1.47, 1.61, 1.76), 
                              annotation=c("*", "***", " *** ")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))
#3B: tau doesn't change with age
df_low$subj_avg_vess_vol = NA
df_low$tau_subject = NA
for (i in 0:10){
  df_low$subj_avg_vess_vol[((4*i)+1):((4*i)+4)] = mean(df_low$normalized_vess_vol[((4*i)+1):((4*i)+4)])
  df_low$tau_subject[((4*i)+1):((4*i)+4)] = mean(df_low$tau[((4*i)+1):((4*i)+4)])
}
coeff = 135
ggplot(df_low, aes(x=age)) + 
  geom_point(aes(y=tau_subject, color="Tau"), size=2) + 
  geom_point(data=df_low, aes(x=age, y=subj_avg_vess_vol * coeff, color="Vessel volume"), size=2) +
  scale_y_continuous(name = "Average tau pathology",
                     sec.axis = sec_axis(~./coeff, name="Average vessel volume fraction"),
                     limits = c(0,3)) +
  geom_smooth(aes(x=age, y=tau_subject), method="lm", formula=y~x) +
  theme_bw() + theme(aspect.ratio=0.5) + 
  ggtitle("Tau: T=0.368 (p=0.0037) \n Vessel volume: (p>0.05)") + xlab("Age (years)") +
  scale_color_manual(values = c("Tau"="blue", "Vessel volume"="red")) +
  theme(legend.position="bottom",
        legend.text = element_text(size=11)) + 
  guides(color = guide_legend(title="")) 


#############
# supp1: plot_subject_path_age
# some data manipulations first
df$sum_path_score = NA
df$sum_path_score_plotting = NA
for(i in 0:23) {
  score = sum(df[((4*i)+1):((4*i)+4),6:8])
  df$sum_path_score[((4*i)+1):((4*i)+4)] = score # H&E, tau, abeta for all 4 CA subregions
  df$sum_path_score_plotting[(4*i)+1] = score
}
df$age_for_plotting = NA
for(i in 0:23) {
  df$age_for_plotting[(4*i)+1] = df$age[(4*i)+1]
}
# plotting
shapes <- c("H&E"=15, "Tau"=17, "Beta-amyloid"=16, "Scaled NIA-AA"=8)
colors <- c("CA1"="red", "CA2"="gold", "CA3"="cyan", "CA4"="blue")
fills <- c("Mean pathology \nscore"="grey60", "Age"="white")
ggplot(df, aes(x=reorder(research_num, sum_path_score))) + 
  geom_col(aes(y=age_for_plotting/100*3, fill="Age"), color="grey60", width=0.4) +
  geom_col(aes(y=sum_path_score_plotting/12, fill="Mean pathology \nscore"), alpha=0.8) + 
  geom_jitter(aes(y=H_E, shape="H&E", color=region), size=2.5, width=0.15, height=0.1) +
  geom_jitter(aes(y=tau, shape="Tau", color=region), size=2.5, width=0.15, height=0.1) + 
  geom_jitter(aes(y=abeta, shape="Beta-amyloid", color=region), size=2.5, width=0.15, height=0.1) + 
  geom_point(aes(y=AD_score*3, shape="Scaled NIA-AA"), size=2.5) +
  xlab("Subjects ranked by mean pathology score") +
  labs(shape="Pathology score", color="Region", fill="") +
  scale_shape_manual(breaks=c("H&E", "Beta-amyloid", "Tau", "Scaled NIA-AA"), values=shapes) +
  scale_color_manual(values=colors) + 
  scale_fill_manual(values=fills) +
  theme_bw() + 
  theme(axis.text.x=element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name="Pathology score", 
                     sec.axis = sec_axis(~.*100/3, name="Age")) +
  guides(fill=guide_legend(order=1), 
         color=guide_legend(order=0), 
         shape=guide_legend(order=2))

#############
# supp 2: age association with path scores (all subjects)
df$he_subject = NA
df$tau_subject = NA
df$abeta_subject = NA
for(i in 0:23){
  df$he_subject[((4*i)+1):((4*i)+4)] = mean(df$H_E[((4*i)+1):((4*i)+4)])
  df$tau_subject[((4*i)+1):((4*i)+4)] = mean(df$tau[((4*i)+1):((4*i)+4)])
  df$abeta_subject[((4*i)+1):((4*i)+4)] = mean(df$abeta[((4*i)+1):((4*i)+4)])
}
he <- ggplot(df) + 
  geom_point(aes(x=age, y=he_subject)) + 
  theme_bw() + theme(aspect.ratio=0.5) + ylab("Average H&E") + ylim(0,3) +
  xlab("Age (years)") + ggtitle("(p>0.1)")
tau <- ggplot(df, aes(x=age, y=tau_subject)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=0.5) + ylab("Average tau") + ylim(0,3) +
  xlab("Age (years)") + ggtitle("(p>0.1)")
abeta <- ggplot(df) + 
  geom_point(aes(x=age, y=abeta_subject)) + 
  theme_bw() + theme(aspect.ratio=0.5) + ylab("Average abeta") + ylim(0,3) +
  xlab("Age (years)") + ggtitle("(p>0.1)")
meanpath <- ggplot(df) + 
  geom_point(aes(x=age, y=mean_path_score)) + 
  theme_bw() + theme(aspect.ratio=0.5) + ylab("Mean pathology score") + ylim(0,3) +
  xlab("Age (years)") + ggtitle("(p>0.1)")
ggarrange(he, tau, abeta, meanpath, ncol=2, nrow=2,widths = c(0.5,0.5,0.5,0.5))











####################################################################################
####################################################################################
## ARCHIVE ##
####################################################################################
####################################################################################
# figure 4G: relative vessel volume in VA regions nonAD vs AD
df_plot = data.frame(matrix(ncol=4, nrow=8))
colnames(df_plot) <- c("var_name", "region", "mean", "std")
df_plot$region = rep(c("CA1", "CA2", "CA3", "CA4"), 2)
df_plot$var_name = c(rep(c("Non-AD"), 4),
                     rep(c("AD"), 4))
df_plot = make_mean_std_df(df_plot, df_low, row_num=1, col_num=38) # relative % vess volume in low
df_plot = make_mean_std_df(df_plot, df_high, row_num=5, col_num=38) # relative % vess volume in high
colors <- c("CA1"="red", "CA2"="gold", "CA3"="cyan", "CA4"="blue")
level_order <- c("Non-AD", "AD")
rPVV <- ggplot(df_plot, aes(x=factor(var_name, level=level_order), y=mean)) +
  geom_bar(aes(fill=region), stat="identity", position=position_dodge()) +
  geom_errorbar(aes(fill=region, ymin=mean-std, ymax=mean+std), width=0.2, position=position_dodge(.9)) +
  scale_fill_manual(values=colors) +
  theme_bw() + xlab("") + ylab("Relative vessel volume") + ylim(0,1.75) + 
  theme(legend.position="top", 
        legend.key.size=unit(0.2,'cm'), 
        legend.title = element_blank()) +
  geom_signif(stat="identity",
              data=data.frame(x=c(0.9, 0.65, 1.9), 
                              xend=c(1.4, 0.9, 2.4),
                              y=c(1.6, 1.4, 1.6), 
                              annotation=c("***", " * ", " *** ")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))


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

###############
# plot 5: vessel length
# figure 5A
highlow <- ggplot(data.frame(name=c("Non-AD", "AD"),
                             mean=c(0.00536357, 0.006287909),
                             std=c(0.001146992, 0.001533897)),
                  aes(x=reorder(name,+mean), y=mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=0.2, position=position_dodge(.9)) +
  geom_signif(stat="identity", 
              data=data.frame(x=1,xend=2, y=0.008, yend=0.008, annotation="**"),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  theme_bw() + ylab("Vessel length [1/um^2]") + xlab("") + ylim(0,0.0085)
# figure 5B
tau <- ggplot(df, aes(y=normalized_vess_leng, x=tau)) +
  geom_point() +
  geom_smooth(method="lm", formula=y~x) +
  scale_color_manual(values=colors) +
  theme_bw() + xlab("Tau") + xlim(0,3) + ylab("Vessel length [1/um^2]") +
  theme(aspect.ratio = 1) +
  ggtitle("T=0.207 (p=0.001)")
# figure 5C
# abeta <- ggplot(df, aes(y=normalized_vess_leng, x=abeta)) +
#   geom_point() +
#   geom_smooth(method="lm", formula=y~x) +
#   scale_color_manual(values=colors) +
#   theme(aspect.ratio = 1) +
#   theme_bw() + xlab("Abeta") + xlim(0,3) + ylab("Vessel length [1/um^2]") +
#   ggtitle("T=0.197 (p=0.019)")
# figure 5D
meanpath <- ggplot(df, aes(y=normalized_vess_leng, x=mean_path_score_CA)) +
  geom_point() +
  geom_smooth(method="lm", formula=y~x) +
  scale_color_manual(values=colors) + 
  theme_bw() + xlab("Mean pathology score") + xlim(0,3) + ylab("Vessel length [1/um^2]") +
  theme(aspect.ratio = 1) +
  ggtitle("T=0.171 (p=0.028)")
# figure 5E
SUBJleng_SUBJpath <- ggplot(df2, aes(x=mean_path_score, y=sum_vess_leng)) + 
  geom_point() +
  geom_smooth(method="lm", formula=y~x) +
  theme_bw() + xlab("Subject mean path score") + ylab("Subject vessel length [1/um^2]") +
  ylim(0.004, 0.012) + theme(aspect.ratio = 1) +
  ggtitle("r=0.674 (p=0.00031)")
# figure 5F
df_CA1 = df2[df2$CA=="CA1",]
CA1leng_SUBJpath <- ggplot(df_CA1, aes(x=mean_path_score, y=normalized_vess_leng)) +
  geom_point() +
  geom_smooth(method="lm", formula=y~x) +
  theme_bw() + xlab("Subject mean path score") + ylab("CA1 vessel length [1/um^2]") +
  ylim(0.004, 0.012) + theme(aspect.ratio = 1) +
  ggtitle("r=0.631 (p=0.00096)")
# figure 5G
df_plot = data.frame(matrix(ncol=4, nrow=8))
colnames(df_plot) <- c("var_name", "region", "mean", "std")
df_plot$region = rep(c("CA1", "CA2", "CA3", "CA4"), 2)
df_plot$var_name = c(rep(c("Non-AD"), 4),
                     rep(c("AD"), 4))
df_plot = make_mean_std_df(df_plot, df_low, row_num=1, col_num=39) # relative vess length in low
df_plot = make_mean_std_df(df_plot, df_high, row_num=5, col_num=39) # relative vess length in high
colors <- c("CA1"="red", "CA2"="gold", "CA3"="cyan", "CA4"="blue")
level_order <- c("Non-AD", "AD")
rNVL <- ggplot(df_plot, aes(x=factor(var_name, level=level_order), y=mean)) +
  geom_bar(aes(fill=region), stat="identity", position=position_dodge()) +
  geom_errorbar(aes(fill=region, ymin=mean-std, ymax=mean+std), width=0.2, position=position_dodge(.9)) +
  scale_fill_manual(values=colors) +
  theme_bw() + xlab("") + ylab("Relative vessel length") + ylim(0,1.6) + 
  theme(legend.position="top", 
        legend.key.size=unit(0.2,'cm'), 
        legend.title = element_blank()) +
  geom_signif(stat="identity",
              data=data.frame(x=c(0.65, 1.1, 2.1), 
                              xend=c(0.9, 1.3, 2.35),
                              y=c(1.4, 1.4, 1.3), 
                              annotation=c("***", " *** ", "*")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))
# put everything together
# ggarrange(highlow, 
#           ggarrange(tau, abeta, meanpath, SUBJleng_SUBJpath, CA1leng_SUBJpath, rNVL,
#                     ncol=3, nrow=2),
#           widths=c(0.4,2), ncol=2, nrow=1)
ggarrange(highlow, 
          ggarrange(tau, meanpath, SUBJleng_SUBJpath, CA1leng_SUBJpath,
                    ncol=2, nrow=2),
          rNVL,
          widths=c(0.4,1,0.6), ncol=3, nrow=1)

# neuron-vessel correlation
# figure A
df_plot = data.frame(matrix(ncol=3, nrow=8))
colnames(df_plot) <- c("region", "mean", "std")
df_plot$region = rep(c("CA1", "CA2", "CA3", "CA4"), 2)
df_plot = make_mean_std_df(df_plot, df_neuron, row_num=1, col_num=36)
colors <- c("CA1"="red", "CA2"="gold", "CA3"="cyan", "CA4"="blue")
rNNDavg <- ggplot(df_plot, aes(x=factor(region), y=mean)) +
  geom_bar(aes(fill=region), stat="identity", position=position_dodge()) +
  geom_errorbar(aes(fill=region, ymin=mean-std, ymax=mean+std), width=0.2, position=position_dodge(.9)) +
  scale_fill_manual(values=colors) +
  theme_bw() + xlab("Subregion") + ylab("Relative average neuron-neuron distance") + ylim(0,1.45) + 
  theme(aspect.ratio=2) + theme(legend.position="none") +
  geom_signif(stat="identity",
              data=data.frame(x=c(1, 2, 3, 2), 
                              xend=c(2, 3, 4, 4),
                              y=c(1.1, 1.15, 1.29, 1.4), 
                              annotation=c("**", "*", " ** ", "***")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))
# figure B
df_plot = make_mean_std_df(df_plot, df_neuron, row_num=1, col_num=34)
rNVDavg <- ggplot(df_plot, aes(x=factor(region), y=mean)) +
  geom_bar(aes(fill=region), stat="identity", position=position_dodge()) +
  geom_errorbar(aes(fill=region, ymin=mean-std, ymax=mean+std), width=0.2, position=position_dodge(.9)) +
  scale_fill_manual(values=colors) +
  theme_bw() + xlab("Subregion") + ylab("Relative average neuron-vessel distance") + ylim(0,1.27) + 
  theme(aspect.ratio=2) + theme(legend.position="none") +
  geom_signif(stat="identity",
              data=data.frame(x=c(1, 2), 
                              xend=c(2, 4),
                              y=c(1.2, 1.15), 
                              annotation=c("**", "*")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))
# figure C-E
HE <- ggplot(df_neuron, aes(x=H_E, y=withNVdist*100)) +
  geom_point() + 
  geom_smooth(method="lm", formula=y~x) +
  xlab("H&E") + ylab("Neurons with NVD<40um [%]") +
  ggtitle("T=0.524 (p=0.00125)") +
  theme_bw() + theme(aspect.ratio = 1)
abeta <- ggplot(df_neuron, aes(x=abeta, y=withNVdist*100)) +
  geom_point() + 
  geom_smooth(method="lm", formula=y~x) +
  xlab("Abeta") + ylab("Neurons with NVD<40um [%]") +
  ggtitle("T=0.378 (p=0.0245)") +
  theme_bw() + theme(aspect.ratio = 1)
meanpath <- ggplot(df_neuron, aes(x=mean_path_score_CA, y=withNVdist*100)) +
  geom_point() + 
  geom_smooth(method="lm", formula=y~x) +
  xlab("Mean pathology score") + ylab("Neurons with NVD<40um [%]") +
  ggtitle("T=0.415 (p=0.00719)") +
  theme_bw() + theme(aspect.ratio = 1)
# figure 6F-K
vessvol_nv <- ggplot(df_neuron, aes(x=normalized_vess_vol, y=avg_nvdist)) +
  geom_point(aes(color=mean_path_score_CA)) + 
  geom_smooth(method="lm", formula=y~x) +
  xlab("Vessel volume [%]") + ylab("Avg NV distance [um]") +
  ggtitle("r=-0.769 (p=1.11e-5)") + labs(color="Mean pathology score") +
  scale_color_gradient2() + theme_bw() + theme(aspect.ratio = 1)
vessleng_nv <- ggplot(df_neuron, aes(x=normalized_vess_leng, y=avg_nvdist)) +
  geom_point(aes(color=mean_path_score_CA)) +
  geom_smooth(method="lm", formula=y~x) +
  xlab("Vessel length [1/um^2]") + ylab("Avg NV distance [um]") +
  ggtitle("r=-0.780 (p=7.10e-6)") + labs(color="Mean pathology score") +
  scale_color_gradient2() + theme_bw() + theme(aspect.ratio = 1)
node_nv <- ggplot(df_neuron, aes(x=node_per_vol, y=avg_nvdist)) +
  geom_point(aes(color=mean_path_score_CA)) +
  geom_smooth(method="lm", formula=y~x) +
  xlab("Branching density [1/um^3]") + ylab("Avg NV distance [um]") +
  ggtitle("r=-0.556 (p=0.00474)") + labs(color="Mean pathology score") +
  scale_color_gradient2() + theme_bw() + theme(aspect.ratio = 1)
vessvol_nn <- ggplot(df_neuron, aes(x=normalized_vess_vol, y=avg_nndist)) +
  geom_point(aes(color=mean_path_score_CA)) +
  geom_smooth(method="lm", formula=y~x) +
  xlab("Vessel volume [%]") + ylab("Avg NN distance [um]") +
  ggtitle("r=-0.450 (p=0.0274)") + labs(color="Mean pathology score") +
  scale_color_gradient2() + theme_bw() + theme(aspect.ratio = 1)
vessleng_nn <- ggplot(df_neuron, aes(x=normalized_vess_leng, y=avg_nndist)) +
  geom_point(aes(color=mean_path_score_CA)) +
  geom_smooth(method="lm", formula=y~x) +
  xlab("Vessel length [1/um^2]") + ylab("Avg NN distance [um]") +
  ggtitle("r=-0.668 (p=0.00036)") + labs(color="Mean pathology score") +
  scale_color_gradient2() + theme_bw() + theme(aspect.ratio = 1)
vessleng_nd <- ggplot(df_neuron, aes(x=normalized_vess_leng, y=neuron_density)) +
  geom_point(aes(color=mean_path_score_CA)) +
  geom_smooth(method="lm", formula=y~x) +
  xlab("Vessel length [1/um^2]") + ylab("Neuron density [1/um^3]") +
  ggtitle("r=0.475 (p=0.0191)") + labs(color="Mean pathology score") +
  scale_color_gradient2() + theme_bw() + theme(aspect.ratio = 1)
# putting everything together
ggarrange(
  ggarrange(rNNDavg, rNVDavg, ncol=1, nrow=2),
  ggarrange(HE, abeta, meanpath,
            vessvol_nv, vessleng_nv, node_nv, 
            vessvol_nn, vessleng_nn, vessleng_nd,
            ncol=3, nrow=3, legend="bottom", common.legend=TRUE),
  ncol=2, nrow=1, widths=c(0.5,2))

