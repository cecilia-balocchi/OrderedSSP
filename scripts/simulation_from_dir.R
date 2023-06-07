rm(list = ls())
library(VGAM)

source("scripts/funs.R")
Rcpp::sourceCpp("scripts/funs.cpp")

Nsim <- 100
sim_x_batch <- 10
n_batches <- 10
Nsim_post <- 25

n <- 500
m <- 5000
Ns <- 1:(n/10)*10
Ms <- 1:(m/10)*10

mh_niter <- 50500
mh_burnin <- 500
mh_thin <- 50
post_subset <- seq(mh_burnin + 1, mh_niter, by = mh_thin)

set.seed(20210916)
thetas <- runif(Nsim, 0,10)
alphas <- rep(0, Nsim)
Ktot <- 10 # for finite dir (which was not included in paper)

fun_ordEPPF_PY <- function(par,ms){
  -ordEPPF_PY(ms, par[1], par[2])
}
fun_ordEPPF_DP <- function(par,ms){
  -ordEPPF_PY(ms, par, 0)
}
fun_EPPF_PY <- function(par,ms){
  -EPPF_PY(ms, par[1], par[2])
}
fun_EPPF_DP <- function(par,ms){
  -EPPF_PY(ms, par, 0)
}
fun_X1LS <- function(par, X1s, Ns){
  if(par[1] < 0) return(NA)
  if(par[2] < 0) return(NA)
  if(par[2] > 1) return(NA)
  sum((X1s - EX1(Ns, par[1], par[2]))^2)
}
fun_Kls <- function(par, Ks, ns){
  if(par[1] < 0) return(NA)
  if(par[2] < 0) return(NA)
  if(par[2] > 1) return(NA)
  sum((EK_PY(ns, par[1], par[2])-Ks)^2)
}
method_str <- c("_ord", "_std", "_ordDP", "_lsK", "_lsX1", "_ordFB")

set.seed(20210916)

# "dir" , "DP" (DP stands for Dirichlet Process, dir for the finite-support dirichlet distribution)
# "dir" was not used in the paper
distr_str <- "DP" 

for(which_order in c("arrival_weight", "astable")){
  for(batch in 1:n_batches){
    rstudioapi::jobRunScript("scripts/simulate_dir.R", name = paste0("simulate_dir_batch",batch), importEnv = TRUE)
  }
}
for(which_order in c("arrival_weight", "astable")){
  for(batch in 1:n_batches){
    rstudioapi::jobRunScript("scripts/simulate_dir_FB.R", name = paste0("simulate_dir_batch",batch), importEnv = TRUE)
  }
}

## analysis of the results

percentage_error <- TRUE
distr_str <- "DP" # "dir" , "DP"

iter <- 1
for(order_str in c("arrival_weight", "astable")){
  which_order = order_str
  results_all <- data.frame(distr = NA, 
                            order = NA,
                            batch = rep(1:n_batches, each = Nsim_post*sim_x_batch),
                            sim = rep(1:Nsim, each = Nsim_post), 
                            true_theta = NA,# true_alpha = NA, 
                            # true_par = NA,
                            theta_ord = NA,alpha_ord = NA, 
                            theta_std = NA,alpha_std = NA, 
                            theta_ordDP = NA,alpha_ordDP = NA, 
                            theta_lsK = NA,alpha_lsK = NA, 
                            theta_lsX1 = NA,alpha_lsX1 = NA, 
                            K = NA, X1 = NA, K2 = NA, A1 = NA, W1 = NA)
  
  for(batch in 1:n_batches){
    sims = (batch-1)*sim_x_batch*Nsim_post + 1:(sim_x_batch*Nsim_post)
    filestr <- paste0("results/",ifelse(distr_str %in% c("dir","DP"),"simdir/","simzipf/"))
    filestr <- paste0(filestr, "sim",distr_str)
    filestr <- paste0(filestr, "_",which_order,"order")
    filestr <- paste0(filestr, "_n500_m5000_batch",batch,".rdata")
    tmpload = load(filestr)
    results_all[sims, ] <- cbind(distr_str,which_order,results)
  }
  if(iter == 1){
    results_comb <- results_all
  } else {
    results_comb <- rbind(results_comb, results_all)
  }
  iter <- iter + 1
}


################################################################
###############  here we do not take the mean ################
################################################################
percentage_error <- TRUE
results3_tidy <- data.frame(sim = rep(1:Nsim, each = Nsim_post*6*length(method_str)),
                            sim_post = rep(rep(1:Nsim_post, each = 6), Nsim*length(method_str)),
                            method = rep(rep(substr(method_str,2,10), each = Nsim_post*6),Nsim),
                            true_theta = NA,
                            param = NA,  RMSE = NA)
for(order_str in c("arrival_weight", "astable")){
  which_order = order_str
  if(order_str == "arrival_weight"){
    results_all <- results_comb[1:2500,]
  } else {
    results_all <- results_comb[2500+1:2500,]
  }
  for(batch in 1:n_batches){
    sims = (batch-1)*sim_x_batch + 1:sim_x_batch
    filestr <- paste0("results/",ifelse(distr_str %in% c("dir","DP"),"simdir/","simzipf/"))
    filestr <- paste0(filestr, "sim",distr_str)
    filestr <- paste0(filestr, "_",which_order,"order")
    filestrA <- paste0(filestr, "_n500_m5000_batch",batch,".rdata")
    tmpload = load(filestrA)
    filestrB <- paste0(filestr, "_n500_m5000_FB_batch",batch,".rdata")
    tmpload = load(filestrB)
    
    W1_curves <- get(paste0("W1_curves_",batch))
    res_fs <- get(paste0("res_fs_",batch))
    W1curve_fs <- get(paste0("W1curve_fs_",batch))
    
    W1curve_fs_post <- get(paste0("W1curve_fs_post_",batch))
    res_fs_post <- get(paste0("res_fs_post_",batch))
    
    # let's incorporate W1curve_fs_post into W1curve_fs
    for(ind in 1:length(W1curve_fs_post)){
      W1curve_fs[[paste0(sims[ind],"_ordFB")]] <- W1curve_fs_post[[ind]]
    }
    for(ind in 1:length(res_fs_post)){
      res_fs[[paste0(sims[ind],"_ordFB")]] <- res_fs_post[[ind]]
    }
    thetas_post_list <- get(paste0("thetas_post_list_",batch))
    alphas_post_list <- get(paste0("alphas_post_list_",batch))

    for(sim in sims){
      for(i_method in 1:length(method_str)){
        res_f <- res_fs[[paste0(sim,method_str[i_method])]]
        EK_f <- res_f[1]; EW1_f <- res_f[2]
        pA1_f <- res_f[3]
        EK_nm_f <- res_f[6]; Eunseen_f <- res_f[7]
        K <- results_all[(sim-1)*Nsim_post + 1, "K"]
        EKunseen_nm_f <- Eunseen_f + K
        EW1curve <- W1curve_fs[[paste0(sim,method_str[i_method])]]
        if(i_method < 6){
          est_theta <- results_all[(sim-1)*Nsim_post + 1,paste0("theta",method_str[i_method])]
          est_alpha <- results_all[(sim-1)*Nsim_post + 1,paste0("alpha",method_str[i_method])]
          pB <- pB1(n,m,est_theta, est_alpha)
        } else {
          est_thetas <- thetas_post_list[[sim-(batch-1)*sim_x_batch]][post_subset]
          est_alphas <- alphas_post_list[[sim-(batch-1)*sim_x_batch]][post_subset]
          pB <- mean(pB1(n,m,est_thetas, est_alphas))
        }
        EW1gA1 <- res_f[4]/(1-pB)
        EW1gB1 <- res_f[5]/pB
        
        method_str_tmp = substr(method_str,2,10)[i_method]
        true_theta <- results_comb[(sim-1)*Nsim_post + 1, "true_theta"]
        
        for(sim_post in 1:Nsim_post){
          index_sim_post <- (sim-1)*Nsim_post + sim_post
          W1 <- results_all[index_sim_post, "W1"]
          A1 <- results_all[index_sim_post, "A1"]
          W1curve <- W1_curves[[paste0(sim,"_",sim_post)]]
          K2 <- results_all[index_sim_post, "K2"]
          
          if(percentage_error){
            index = which(results3_tidy[,1] == sim &
                            results3_tidy[,2] == sim_post &
                            results3_tidy[,3] == method_str_tmp)
            results3_tidy[index[1],4:6] <- c(true_theta,"W1curve",
                                             sqrt(mean(( (W1curve - EW1curve)/W1curve )^2)) )
            results3_tidy[index[2],4:6] <- c(true_theta,"W1",
                                             abs(W1-EW1_f)/W1 )
            results3_tidy[index[3],4:6] <- c(true_theta,"A1",
                                             abs(A1-pA1_f) )
            results3_tidy[index[4],4:6] <- c(true_theta,"K",
                                             abs(K2-EKunseen_nm_f)/K2 )
            
            if(A1){
              # RMSE_W1_A1_tmp <- abs(W1-res_f[4])/W1
              RMSE_W1_A1_tmp <- abs(W1-EW1gA1)/W1
              RMSE_W1_B1_tmp <- NA
            } else {
              RMSE_W1_A1_tmp <- NA
              # RMSE_W1_B1_tmp <- abs(W1-res_f[5])/W1
              RMSE_W1_B1_tmp <- abs(W1-EW1gB1)/W1
            }
            results3_tidy[index[5],4:6] <- c(true_theta,"W1_A1",
                                             RMSE_W1_A1_tmp )
            results3_tidy[index[6],4:6] <- c(true_theta,"W1_B1",
                                             RMSE_W1_B1_tmp )
          } else {
            warning("normal error not implemented!")
          }
        }
      }
    }
  }
  results3_tidy$error <- "percentage"
  results3_tidy$true_theta <- as.numeric(results3_tidy$true_theta)
  results3_tidy$RMSE <- as.numeric(results3_tidy$RMSE)
  results3_tidy$order <- order_str
  assign(paste0("results3_tidy_",order_str),results3_tidy)
}
results3_tidy <- rbind(results3_tidy_astable, results3_tidy_arrival_weight)

######################################################
######################################################

save(list = c("results3_tidy","results_comb"), file = "results/results_DP_FB.Rdata")

