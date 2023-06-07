rm(list = ls())
methods_str <- c("ordDP","ordPYP","stdPYP","lsX1","lsK")
tmp = load("results/crossval_BRCA_FB.Rdata")
Ks_perr_BRCA <- Ks_perr
W1_perr_BRCA <- W1_perr
tmp = load("results/crossval_EDAR_FB.Rdata")
Ks_perr_EDAR <- Ks_perr
W1_perr_EDAR <- W1_perr

library(dplyr)
library(tidyr)
library(ggplot2)
tmp_Ks_BRCA = pivot_longer(as.data.frame(Ks_perr_BRCA), cols = colnames(Ks_perr_BRCA), names_to = "method")
tmp_Ks_BRCA$variable = "K"
tmp_Ks_BRCA$data = "BRCA"
tmp_Ks_EDAR = pivot_longer(as.data.frame(Ks_perr_EDAR), cols = colnames(Ks_perr_EDAR), names_to = "method")
tmp_Ks_EDAR$variable = "K"
tmp_Ks_EDAR$data = "EDAR"
tmp_W1s_BRCA = pivot_longer(as.data.frame(W1_perr_BRCA), cols = colnames(W1_perr_BRCA), names_to = "method")
tmp_W1s_BRCA$variable = "W1"
tmp_W1s_BRCA$data = "BRCA"
tmp_W1s_EDAR = pivot_longer(as.data.frame(W1_perr_EDAR), cols = colnames(W1_perr_EDAR), names_to = "method")
tmp_W1s_EDAR$variable = "W1"
tmp_W1s_EDAR$data = "EDAR"

results <- rbind(tmp_Ks_BRCA,tmp_Ks_EDAR,tmp_W1s_BRCA,tmp_W1s_EDAR)
results$method[results$method == "lsX1"] <- "lsM1"
results$method <- factor(results$method, levels = c("FB","ordPYP","ordDP","stdPYP","lsM1","lsK"))

p = ggplot(results) + geom_boxplot(aes(x = method, y = value)) +
  facet_nested( ~ variable+data, scale = "free") + scale_y_log10() + ylab("Percentage error")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Predictive performance on real data")
p
ggsave(filename = "genetics_crossval.png", plot = p, device = "png", path = "figures", 
       width = 10, height = 3)

############### Table with results


library(ggh4x)
library(kableExtra)

tmp1b = results %>% group_by(data,variable,method) %>% summarize(RMSE = median(value, na.rm = T)) 
temp1b = pivot_wider(tmp1b, 
                     names_from = c("variable","data"),
                     values_from = "RMSE")
temp1b <- temp1b[,1+c(0,1,3,2,4)]

kbl(temp1b, booktabs = T, format = "latex", digits = 3) %>% kable_styling(latex_options = "striped")
