## this uses results3_tidy
rm(list = ls())
library(dplyr)
library(ggplot2)
library(ggh4x)

tmp = load("results/results_model.Rdata")
results_model <- results3_tidy
data_model <- results_all
tmp = load("results/results_DP.Rdata")
results_DP <- results3_tidy
data_DP <- results_comb
tmp = load("results/results_PYP.Rdata")
results_PYP <- results3_tidy
data_PYP <- results_comb
tmp = load("results/results_zipf.Rdata") ## note: this is infinite zipf (zeta distr), not finzipf (finite zipf)
results_zipf <- results3_tidy
data_zipf <- results_comb

results_model$order <- "model"
results_DP$true_alpha <- NA
results_zipf$true_theta <- NA
names(results_zipf)[4] <- "true_alpha"

results_zipf$order[results_zipf$order == "arrival_weight"] <- "arrival-weighted"
results_zipf$order[results_zipf$order == "astable"] <- "alpha-stable"
results_DP$order[results_DP$order == "arrival_weight"] <- "arrival-weighted"
results_DP$order[results_DP$order == "astable"] <- "alpha-stable"
results_PYP$order[results_PYP$order == "arrival_weight"] <- "arrival-weighted"
results_PYP$order[results_PYP$order == "astable"] <- "alpha-stable"

results_zipf$distr <- "zipf"
results_DP$distr <- "DP"
results_PYP$distr <- "PYP"
results_model$distr <- "model"

results_model <- results_model[,c("sim","true_theta","true_alpha","distr","order","param","method","RMSE","error")]
results_DP <- results_DP[,c("sim","true_theta","true_alpha","distr","order","param","method","RMSE","error")]
results_PYP <- results_PYP[,c("sim","true_theta","true_alpha","distr","order","param","method","RMSE","error")]
results_zipf <- results_zipf[,c("sim","true_theta","true_alpha","distr","order","param","method","RMSE","error")]

combined <- rbind(results_model,results_DP, results_PYP,results_zipf)
combined$param[combined$param == "W1_A1"] <- "W1|A1"
combined$param[combined$param == "W1_B1"] <- "W1|B1"

combined$method[combined$method == "ord"] <- "ordPYP"
combined$method[combined$method == "std"] <- "stdPYP"

combined$method <- factor(combined$method, levels = c("stdPYP","ordPYP","ordDP","lsX1","lsK"))





##########################################################################################
##########################################################################################
p1 = ggplot(combined %>% filter(error == "percentage") %>% 
              filter(param %in% c("W1","K","W1|A1","W1|B1")) %>%
              group_by(method,param,order,distr,sim) %>% 
              summarize(RMSE = median(RMSE, na.rm = T))) + 
  geom_boxplot(aes(x = method, y = RMSE), outlier.size = 0.1) +scale_y_log10() +
  facet_nested(param ~ order + distr , scale = "free") + 
  ylab("Percentage error")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Predictive performance under misspecification")
p1
ggsave(filename = "misspec_simulation_group_median.png", plot = p1, device = "png", path = "figures", 
       width = 8, height = 6)
##########################################################################################
##########################################################################################
p2 = ggplot(combined %>% filter(error == "percentage",distr  == "model") %>% 
              filter(param %in% c("W1","K","W1|A1","W1|B1")) %>%
              group_by(method,param,order,distr,sim) %>% 
              summarize(RMSE = median(RMSE, na.rm = T))) + 
  geom_boxplot(aes(x = method, y = RMSE), outlier.size = 0.1) + scale_y_log10() +
  # facet_grid( ~ param, scales = "free") + 
  facet_wrap(~param, ncol = 4) + ylab("Percentage error")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Predictive performance under correct specification")
p2
ggsave(filename = "model_simulation_group_median.png", plot = p2, device = "png", path = "figures", 
       width = 8, height = 2.5)







##########################################################################################
######################### table with numeric results #####################################
##########################################################################################
tmp1 = combined %>% filter(error == "percentage") %>% 
  filter(param %in% c("W1","K","W1|A1","W1|B1")) %>%
  group_by(method,param,order,distr,sim) %>% 
  summarize(RMSE = median(RMSE, na.rm = T)) 
tmp1b = tmp1 %>% group_by(order,distr,param,method,) %>% summarize(RMSE = median(RMSE, na.rm = T)) 


library(tidyr)
library(kableExtra)
temp1b = pivot_wider(tmp1b, 
            names_from = c("order","distr"),
            values_from = "RMSE")
temp1b2 <- rbind(temp1b[1,],
                 temp1b[1:5,])
for(i in 2:4){
  temp1b2 <- rbind(temp1b2, 
                   temp1b[(i-1)*5+1,],
                   temp1b[(i-1)*5+1:5,])
}
dim(temp1b2)
temp1b2 <- temp1b2[,-1]

kbl(temp1b2, booktabs = T, format = "latex") %>% kable_styling(latex_options = "striped")
