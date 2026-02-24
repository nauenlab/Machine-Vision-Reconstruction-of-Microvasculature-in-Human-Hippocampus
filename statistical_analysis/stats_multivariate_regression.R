# multivariate multiple predicting model
# inputs: age, vasculature features (normalized vessel volume and length,
#   node per volume and length, avg/max weighted tortuosity, avg/max segment length)
# prediction outputs: H&E, tau, abeta patho scores
# ========
# multivariate regression code: https://library.virginia.edu/data/articles/getting-started-with-multivariate-multiple-regression
# interpret outputs: https://towardsdatascience.com/understanding-linear-regression-output-in-r-7a9cbda948b3

# multivariate regression on all samples
multivariate_regression <- function(df) {
  m_he <- lm(H_E ~ age + normalized_vess_vol + normalized_vess_leng +
               node_per_vol + node_per_leng + avg_wtort + max_wtort +
               avg_seg_leng + max_seg_leng + avg_vess_diameter,
             data = df)
  m_tau <- lm(tau ~ age + normalized_vess_vol + normalized_vess_leng +
                node_per_vol + node_per_leng + avg_wtort + max_wtort +
                avg_seg_leng + max_seg_leng + avg_vess_diameter,
              data = df)
  m_abeta <- lm(abeta ~ age + normalized_vess_vol + normalized_vess_leng +
                  node_per_vol + node_per_leng + avg_wtort + max_wtort +
                  avg_seg_leng + max_seg_leng + avg_vess_diameter,
                data = df)
  m_meanpath <- lm(mean_path_score ~ age + normalized_vess_vol + normalized_vess_leng +
                     node_per_vol + node_per_leng + avg_wtort + max_wtort +
                     avg_seg_leng + max_seg_leng + avg_vess_diameter,
                   data = df)
  print(summary(m_he))
  print(summary(m_tau))
  print(summary(m_abeta))
  print(summary(m_meanpath))
}
df_low = df[df$lowAD_selection=='y',]
df_high = df[df$lowAD_selection=='n',]
multivariate_regression(df)
multivariate_regression(df[df$region=='CA1',])
multivariate_regression(df[df$region=='CA2',])
multivariate_regression(df[df$region=='CA3',])
multivariate_regression(df[df$region=='CA4',])
multivariate_regression(df_low)
multivariate_regression(df_high)



