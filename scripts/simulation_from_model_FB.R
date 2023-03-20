rm(list = ls())

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

mh_niter <- 1000
mh_burnin <- 150
mh_thin <- 10
post_subset <- seq(mh_burnin + 1, mh_niter, by = mh_thin)

set.seed(20210916)
thetas <- runif(Nsim, 0,10)
alphas <- rbeta(Nsim, 2,4)

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
method_str <- c("_ord", "_std", "_ordDP", "_lsK", "_lsX1")

set.seed(20210916)
for(batch in 1:n_batches){
  rstudioapi::jobRunScript("scripts/simulate_model_FB.R", name = paste0("simulate_model_FB_batch",batch), importEnv = TRUE)
}

## analysis of the results

percentage_error <- TRUE

results_all <- data.frame(batch = rep(1:n_batches, each = Nsim_post*sim_x_batch),
                          sim = rep(1:Nsim, each = Nsim_post), 
                          sim_post = rep(1:Nsim_post, times = Nsim), 
                          true_theta = NA,true_alpha = NA, 
                          theta_ord = NA,alpha_ord = NA, 
                          theta_std = NA,alpha_std = NA, 
                          theta_ordDP = NA,alpha_ordDP = NA, 
                          theta_lsK = NA,alpha_lsK = NA, 
                          theta_lsX1 = NA,alpha_lsX1 = NA, 
                          K = NA, X1 = NA, A1 = NA, W1 = NA)

for(batch in 1:n_batches){
  sims = (batch-1)*sim_x_batch*Nsim_post + 1:(sim_x_batch*Nsim_post)
  tmpload = load(paste0("results/simmodel/simmodel_n500_m5000_batch",batch,".rdata"))
  tmpload = load(paste0("results/simmodel/simmodel_FB_n500_m5000_batch",batch,".rdata"))
  results_all[sims, ] <- results
}

## NOTE: results_tidy is not used for plots!
results_tidy <- data.frame(sim = NA, true_theta = NA, true_alpha = NA, 
                            param = NA, method = NA, RMSE = NA)
for(batch in 1:n_batches){
  sims = (batch-1)*sim_x_batch + 1:sim_x_batch
  W1_curves <- get(paste0("W1_curves_",batch)) # this is from the sim_post data
  res_fs <- get(paste0("res_fs_",batch)) # these are from formulas and parameters
  W1curve_fs <- get(paste0("W1curve_fs_",batch)) # these are from formulas and parameters
  
  for(sim in sims){
    for(i_method in 1:5){
      RMSE_W1curve_tmp <- 0
      RMSE_W1_tmp <- 0
      RMSE_A1_tmp <- 0
      RMSE_K_tmp <- 0
      RMSE_W1_A1_tmp <- c()
      RMSE_W1_B1_tmp <- c()
      
      W1_curve_tmp <- rep(0, 500) # this depends on the length of Ms
      W1_tmp <- 0
      A1_tmp <- 0
      K_tmp <- 0
      W1_A1_tmp <- c()
      W1_B1_tmp <- c()
      
      res_f <- res_fs[[paste0(sim,method_str[i_method])]] 
      EK_f <- res_f[1]; EW1_f <- res_f[2]
      pA1_f <- res_f[3]
      EW1curve <- W1curve_fs[[paste0(sim,method_str[i_method])]]
      
      for(sim_post in 1:Nsim_post){
        index_sim_post <- (sim-1)*Nsim_post + sim_post
        
        W1 <- results_all[index_sim_post, "W1"]
        A1 <- results_all[index_sim_post, "A1"]
        W1curve <- W1_curves[[paste0(sim,"_",sim_post)]]
        K <- results_all[index_sim_post, "K"]
        
        if(percentage_error){
          RMSE_W1curve_tmp <- RMSE_W1curve_tmp +
            sqrt(mean(( (W1curve - EW1curve)/W1curve )^2))
          RMSE_W1_tmp <- RMSE_W1_tmp + abs(W1-EW1_f)/W1
          RMSE_A1_tmp <- RMSE_A1_tmp + abs(A1-pA1_f)
          RMSE_K_tmp <- RMSE_K_tmp + abs(K-EK_f)/K
          if(A1){
            RMSE_W1_A1_tmp <- c(RMSE_W1_A1_tmp,abs(W1-res_f[4])/W1)
            RMSE_W1_B1_tmp <- c(RMSE_W1_B1_tmp,NA)
          } else {
            RMSE_W1_A1_tmp <- c(RMSE_W1_A1_tmp,NA)
            RMSE_W1_B1_tmp <- c(RMSE_W1_B1_tmp,abs(W1-res_f[5])/W1)
          }
        } else {
          RMSE_W1curve_tmp <- RMSE_W1curve_tmp +
            sqrt(mean(( W1curve - EW1curve )^2))
          RMSE_W1_tmp <- RMSE_W1_tmp + abs(W1-EW1_f)
          RMSE_A1_tmp < + abs(A1-pA1_f)
          RMSE_K_tmp <- RMSE_K_tmp + abs(K-EK_f)
          if(A1){
            RMSE_W1_A1_tmp <- c(RMSE_W1_A1_tmp,abs(W1-res_f[4]))
            RMSE_W1_B1_tmp <- c(RMSE_W1_B1_tmp,NA)
          } else {
            RMSE_W1_A1_tmp <- c(RMSE_W1_A1_tmp,NA)
            RMSE_W1_B1_tmp <- c(RMSE_W1_B1_tmp,abs(W1-res_f[5]))
          }
        }
        
        W1_curve_tmp <- W1_curve_tmp + W1_curves[[paste0(sim,"_",sim_post)]]
        W1_tmp <- W1_tmp + W1
        A1_tmp <- A1_tmp + A1
        K_tmp <- K_tmp + results_all[index_sim_post, "K"]
        if(A1){
          W1_A1_tmp <- c(W1_A1_tmp,W1)
          W1_B1_tmp <- c(W1_B1_tmp,NA)
        } else {
          W1_A1_tmp <- c(W1_A1_tmp,NA)
          W1_B1_tmp <- c(W1_B1_tmp,W1)
        }
        
      }
      
      method_str_tmp = substr(method_str,2,10)[i_method]
      true_theta <- results_all[index_sim_post, "true_theta"]
      true_alpha <- results_all[index_sim_post, "true_alpha"]
      
      tmp = list(sim,true_theta,true_alpha,"W1curve",method_str_tmp,RMSE_W1curve_tmp/Nsim_post)
      names(tmp) <- colnames(results_tidy)
      results_tidy <- rbind(results_tidy, tmp)
      
      tmp = list(sim,true_theta,true_alpha,"W1",method_str_tmp,RMSE_W1_tmp/Nsim_post)
      names(tmp) <- colnames(results_tidy)
      results_tidy <- rbind(results_tidy, tmp)
      
      tmp = list(sim,true_theta,true_alpha,"A1",method_str_tmp,RMSE_A1_tmp/Nsim_post)
      names(tmp) <- colnames(results_tidy)
      results_tidy <- rbind(results_tidy, tmp)
      
      tmp = list(sim,true_theta,true_alpha,"K",method_str_tmp,RMSE_K_tmp/Nsim_post)
      names(tmp) <- colnames(results_tidy)
      results_tidy <- rbind(results_tidy, tmp)
      
      tmp = list(sim,true_theta,true_alpha,"W1_A1",method_str_tmp,mean(RMSE_W1_A1_tmp, na.rm = T))
      names(tmp) <- colnames(results_tidy)
      results_tidy <- rbind(results_tidy, tmp)
      
      tmp = list(sim,true_theta,true_alpha,"W1_B1",method_str_tmp,mean(RMSE_W1_B1_tmp, na.rm = T))
      names(tmp) <- colnames(results_tidy)
      results_tidy <- rbind(results_tidy, tmp)
      
      ### 
      if(percentage_error){
        tmp1 = list(sim,true_theta,true_alpha,"W1curve2",method_str_tmp,
                    sqrt(mean(( (W1_curve_tmp/Nsim_post - W1curve_fs[[paste0(sim,method_str[i_method])]])/(W1_curve_tmp/Nsim_post) )^2)) )
        tmp2 = list(sim,true_theta,true_alpha,"W12",method_str_tmp,abs(W1_tmp/Nsim_post-EW1_f)/(W1_tmp/Nsim_post))
        tmp3 = list(sim,true_theta,true_alpha,"K2",method_str_tmp,abs(K_tmp/Nsim_post-EK_f)/(K_tmp/Nsim_post))
        tmp4 = list(sim,true_theta,true_alpha,"W1_A12",method_str_tmp,abs(mean(W1_A1_tmp, na.rm = TRUE) - res_f[4])/mean(W1_A1_tmp, na.rm = TRUE))
        tmp5 = list(sim,true_theta,true_alpha,"W1_B12",method_str_tmp,abs(mean(W1_B1_tmp, na.rm = T) - res_f[5])/mean(W1_B1_tmp, na.rm = TRUE))
        tmp6 = list(sim,true_theta,true_alpha,"A12",method_str_tmp,abs(A1_tmp/Nsim_post-pA1_f)/(A1_tmp/Nsim_post))
      } else {
        tmp1 = list(sim,true_theta,true_alpha,"W1curve2",method_str_tmp,
                    sqrt(mean(( W1_curve_tmp/Nsim_post - W1curve_fs[[paste0(sim,method_str[i_method])]] )^2)))
        tmp2 = list(sim,true_theta,true_alpha,"W12",method_str_tmp,abs(W1_tmp/Nsim_post-EW1_f))
        tmp3 = list(sim,true_theta,true_alpha,"K2",method_str_tmp,abs(K_tmp/Nsim_post-EK_f))
        tmp4 = list(sim,true_theta,true_alpha,"W1_A12",method_str_tmp,abs(mean(W1_A1_tmp, na.rm = TRUE) - res_f[4]))
        tmp5 = list(sim,true_theta,true_alpha,"W1_B12",method_str_tmp,abs(mean(W1_B1_tmp, na.rm = T) - res_f[5]))
        tmp6 = list(sim,true_theta,true_alpha,"A12",method_str_tmp,abs(A1_tmp/Nsim_post-pA1_f))
      }
      
      
      tmp = tmp1; names(tmp) <- colnames(results_tidy)
      results_tidy <- rbind(results_tidy, tmp)
      tmp = tmp2; names(tmp) <- colnames(results_tidy)
      results_tidy <- rbind(results_tidy, tmp)
      tmp = tmp3; names(tmp) <- colnames(results_tidy)
      results_tidy <- rbind(results_tidy, tmp)
      tmp = tmp4; names(tmp) <- colnames(results_tidy)
      results_tidy <- rbind(results_tidy, tmp)
      tmp = tmp5; names(tmp) <- colnames(results_tidy)
      results_tidy <- rbind(results_tidy, tmp)
      tmp = tmp6; names(tmp) <- colnames(results_tidy)
      results_tidy <- rbind(results_tidy, tmp)
    }
  }
}
results_tidy <- results_tidy[-1,]
results_tidy$error <- "percentage"

## let's add the results from Full-Bayes (FB)
for(batch in 1:n_batches){
  sims = (batch-1)*sim_x_batch + 1:sim_x_batch
  W1_curves <- get(paste0("W1_curves_",batch)) # this is from the sim_post data
  res_fs <- get(paste0("res_fs_",batch)) # these are from formulas and parameters
  W1curve_fs <- get(paste0("W1curve_fs_",batch)) # these are from formulas and parameters
  
  # these are all posterior formulas (for using MH)
  # thetas_post_list <- get(paste0("thetas_post_list_",batch))
  # alphas_post_list <- get(paste0("alphas_post_list_",batch))
  EK_list <- get(paste0("EK_list_",batch))
  pA1_list <- get(paste0("pA1_list_",batch))
  EW1_list <- get(paste0("EW1_list_",batch))
  EW1_A1_list <- get(paste0("EW1_A1_list_",batch))
  EW1_B1_list <- get(paste0("EW1_B1_list_",batch))
  W1curve_f_list <- get(paste0("W1curve_f_list_",batch))
  # this contains the posterior means.
  res_fs_post <- get(paste0("res_fs_post_",batch))
  
  for(sim in sims){
    ## FullBayes
    RMSE_W1curve_tmp <- 0
    RMSE_W1_tmp <- 0
    RMSE_A1_tmp <- 0
    RMSE_K_tmp <- 0
    RMSE_W1_A1_tmp <- c()
    RMSE_W1_B1_tmp <- c()
    
    W1_curve_tmp <- rep(0, 500) # this depends on the length of Ms
    W1_tmp <- 0
    A1_tmp <- 0
    K_tmp <- 0
    W1_A1_tmp <- c()
    W1_B1_tmp <- c()
    
    res_f <- res_fs_post[[paste0(sim)]] 
    EK_f <- res_f[1]; EW1_f <- res_f[2]; pA1_f <- res_f[3]
    
    # EK_post <- EK_list[[paste0(sim)]]
    # pA1_post <- pA1_list[[paste0(sim)]]
    # EW1_post <- EW1_list[[paste0(sim)]]
    # EW1_A1_post <- EW1_A1_list[[paste0(sim)]]
    # EW1_B1_post <- EW1_B1_list[[paste0(sim)]]
    # EK_f <- mean(EK_post); EW1_f <- mean(EW1_post); pA1_f <- mean(pA1_post)
    
    EW1curve_post <- W1curve_f_list[[paste0(sim)]]
    EW1curve <- colMeans(EW1curve_post)
    
    for(sim_post in 1:Nsim_post){
      index_sim_post <- (sim-1)*Nsim_post + sim_post
      
      W1 <- results_all[index_sim_post, "W1"]
      A1 <- results_all[index_sim_post, "A1"]
      W1curve <- W1_curves[[paste0(sim,"_",sim_post)]]
      K <- results_all[index_sim_post, "K"]
      
      if(percentage_error){
        RMSE_W1curve_tmp <- RMSE_W1curve_tmp +
          sqrt(mean(( (W1curve - EW1curve)/W1curve )^2))
        # apply(EW1curve_post, MARGIN = 1, function(x) sqrt(mean(( (W1curve - x)/W1curve )^2))  )
        
        RMSE_W1_tmp <- RMSE_W1_tmp + abs(W1-EW1_f)/W1
        RMSE_A1_tmp <- RMSE_A1_tmp + abs(A1-pA1_f)
        RMSE_K_tmp <- RMSE_K_tmp + abs(K-EK_f)/K
        if(A1){
          RMSE_W1_A1_tmp <- c(RMSE_W1_A1_tmp,abs(W1-res_f[4])/W1)
          RMSE_W1_B1_tmp <- c(RMSE_W1_B1_tmp,NA)
        } else {
          RMSE_W1_A1_tmp <- c(RMSE_W1_A1_tmp,NA)
          RMSE_W1_B1_tmp <- c(RMSE_W1_B1_tmp,abs(W1-res_f[5])/W1)
        }
      } else {
        RMSE_W1curve_tmp <- RMSE_W1curve_tmp +
          sqrt(mean(( W1curve - EW1curve )^2))
        RMSE_W1_tmp <- RMSE_W1_tmp + abs(W1-EW1_f)
        RMSE_A1_tmp < + abs(A1-pA1_f)
        RMSE_K_tmp <- RMSE_K_tmp + abs(K-EK_f)
        if(A1){
          RMSE_W1_A1_tmp <- c(RMSE_W1_A1_tmp,abs(W1-res_f[4]))
          RMSE_W1_B1_tmp <- c(RMSE_W1_B1_tmp,NA)
        } else {
          RMSE_W1_A1_tmp <- c(RMSE_W1_A1_tmp,NA)
          RMSE_W1_B1_tmp <- c(RMSE_W1_B1_tmp,abs(W1-res_f[5]))
        }
      }
      
      W1_curve_tmp <- W1_curve_tmp + W1_curves[[paste0(sim,"_",sim_post)]]
      W1_tmp <- W1_tmp + W1
      A1_tmp <- A1_tmp + A1
      K_tmp <- K_tmp + results_all[index_sim_post, "K"]
      if(A1){
        W1_A1_tmp <- c(W1_A1_tmp,W1)
        W1_B1_tmp <- c(W1_B1_tmp,NA)
      } else {
        W1_A1_tmp <- c(W1_A1_tmp,NA)
        W1_B1_tmp <- c(W1_B1_tmp,W1)
      }
      
    }
    
    method_str_tmp = "ordFB"
    true_theta <- results_all[index_sim_post, "true_theta"]
    true_alpha <- results_all[index_sim_post, "true_alpha"]
    
    tmp = list(sim,true_theta,true_alpha,"W1curve",method_str_tmp,RMSE_W1curve_tmp/Nsim_post)
    names(tmp) <- colnames(results_tidy)
    results_tidy <- rbind(results_tidy, tmp)
    
    tmp = list(sim,true_theta,true_alpha,"W1",method_str_tmp,RMSE_W1_tmp/Nsim_post)
    names(tmp) <- colnames(results_tidy)
    results_tidy <- rbind(results_tidy, tmp)
    
    tmp = list(sim,true_theta,true_alpha,"A1",method_str_tmp,RMSE_A1_tmp/Nsim_post)
    names(tmp) <- colnames(results_tidy)
    results_tidy <- rbind(results_tidy, tmp)
    
    tmp = list(sim,true_theta,true_alpha,"K",method_str_tmp,RMSE_K_tmp/Nsim_post)
    names(tmp) <- colnames(results_tidy)
    results_tidy <- rbind(results_tidy, tmp)
    
    tmp = list(sim,true_theta,true_alpha,"W1_A1",method_str_tmp,mean(RMSE_W1_A1_tmp, na.rm = T))
    names(tmp) <- colnames(results_tidy)
    results_tidy <- rbind(results_tidy, tmp)
    
    tmp = list(sim,true_theta,true_alpha,"W1_B1",method_str_tmp,mean(RMSE_W1_B1_tmp, na.rm = T))
    names(tmp) <- colnames(results_tidy)
    results_tidy <- rbind(results_tidy, tmp)
    
    ### 
    if(percentage_error){
      tmp1 = list(sim,true_theta,true_alpha,"W1curve2",method_str_tmp,
                  sqrt(mean(( (W1_curve_tmp/Nsim_post - EW1curve)/(W1_curve_tmp/Nsim_post) )^2)) )
      tmp2 = list(sim,true_theta,true_alpha,"W12",method_str_tmp,abs(W1_tmp/Nsim_post-EW1_f)/(W1_tmp/Nsim_post))
      tmp3 = list(sim,true_theta,true_alpha,"K2",method_str_tmp,abs(K_tmp/Nsim_post-EK_f)/(K_tmp/Nsim_post))
      tmp4 = list(sim,true_theta,true_alpha,"W1_A12",method_str_tmp,abs(mean(W1_A1_tmp, na.rm = TRUE) - res_f[4])/mean(W1_A1_tmp, na.rm = TRUE))
      tmp5 = list(sim,true_theta,true_alpha,"W1_B12",method_str_tmp,abs(mean(W1_B1_tmp, na.rm = T) - res_f[5])/mean(W1_B1_tmp, na.rm = TRUE))
      tmp6 = list(sim,true_theta,true_alpha,"A12",method_str_tmp,abs(A1_tmp/Nsim_post-pA1_f)/(A1_tmp/Nsim_post))
    } else {
      tmp1 = list(sim,true_theta,true_alpha,"W1curve2",method_str_tmp,
                  sqrt(mean(( W1_curve_tmp/Nsim_post - EW1curve )^2)))
      tmp2 = list(sim,true_theta,true_alpha,"W12",method_str_tmp,abs(W1_tmp/Nsim_post-EW1_f))
      tmp3 = list(sim,true_theta,true_alpha,"K2",method_str_tmp,abs(K_tmp/Nsim_post-EK_f))
      tmp4 = list(sim,true_theta,true_alpha,"W1_A12",method_str_tmp,abs(mean(W1_A1_tmp, na.rm = TRUE) - res_f[4]))
      tmp5 = list(sim,true_theta,true_alpha,"W1_B12",method_str_tmp,abs(mean(W1_B1_tmp, na.rm = T) - res_f[5]))
      tmp6 = list(sim,true_theta,true_alpha,"A12",method_str_tmp,abs(A1_tmp/Nsim_post-pA1_f))
    }
    
    
    tmp = tmp1; names(tmp) <- colnames(results_tidy)
    results_tidy <- rbind(results_tidy, tmp)
    tmp = tmp2; names(tmp) <- colnames(results_tidy)
    results_tidy <- rbind(results_tidy, tmp)
    tmp = tmp3; names(tmp) <- colnames(results_tidy)
    results_tidy <- rbind(results_tidy, tmp)
    tmp = tmp4; names(tmp) <- colnames(results_tidy)
    results_tidy <- rbind(results_tidy, tmp)
    tmp = tmp5; names(tmp) <- colnames(results_tidy)
    results_tidy <- rbind(results_tidy, tmp)
    tmp = tmp6; names(tmp) <- colnames(results_tidy)
    results_tidy <- rbind(results_tidy, tmp)
  }
}

################################################################
##### let's try something different: do not take the mean ######
################################################################
method_str <- c("_ord", "_std", "_ordDP", "_lsK", "_lsX1", "_ordFB")

for(batch in 1:n_batches){
  sims = (batch-1)*sim_x_batch*Nsim_post + 1:(sim_x_batch*Nsim_post)
  tmpload = load(paste0("results/simmodel/simmodel_n500_m5000_batch",batch,".rdata"))
  tmpload = load(paste0("results/simmodel/simmodel_FB_n500_m5000_batch",batch,".rdata"))
}

results3_tidy <- data.frame(sim = rep(1:Nsim, each = Nsim_post*6*length(method_str)), 
                            sim_post = rep(rep(1:Nsim_post, each = 6), Nsim*length(method_str)), 
                            method = rep(rep(substr(method_str,2,10), each = Nsim_post*6),Nsim),
                            true_theta = NA, true_alpha = NA, 
                            param = NA,  RMSE = NA)
for(batch in 1:n_batches){
  sims = (batch-1)*sim_x_batch + 1:sim_x_batch
  W1_curves <- get(paste0("W1_curves_",batch))
  res_fs <- get(paste0("res_fs_",batch))
  W1curve_fs <- get(paste0("W1curve_fs_",batch))
  
  W1curve_f_list <- get(paste0("W1curve_f_list_",batch))
  res_fs_post <- get(paste0("res_fs_post_",batch))
  
  
  for(ind in 1:length(W1curve_f_list)){
    W1curve_fs[[paste0(sims[ind],"_ordFB")]] <- colMeans(W1curve_f_list[[ind]])
  }
  for(ind in 1:length(res_fs_post)){
    res_fs[[paste0(sims[ind],"_ordFB")]] <- res_fs_post[[ind]]
  }
  
  for(sim in sims){
    for(i_method in 1:length(method_str)){
      res_f <- res_fs[[paste0(sim,method_str[i_method])]] 
      EK_f <- res_f[1]; EW1_f <- res_f[2]
      pA1_f <- res_f[3]
      EW1curve <- W1curve_fs[[paste0(sim,method_str[i_method])]]
      
      method_str_tmp = substr(method_str,2,10)[i_method]
      true_theta <- results_all[(sim-1)*Nsim_post + 1, "true_theta"]
      true_alpha <- results_all[(sim-1)*Nsim_post + 1, "true_alpha"]
      
      for(sim_post in 1:Nsim_post){
        index_sim_post <- (sim-1)*Nsim_post + sim_post
        W1 <- results_all[index_sim_post, "W1"]
        A1 <- results_all[index_sim_post, "A1"]
        W1curve <- W1_curves[[paste0(sim,"_",sim_post)]]
        K <- results_all[index_sim_post, "K"]
        
        if(percentage_error){
          index = which(results3_tidy[,1] == sim & 
                          results3_tidy[,2] == sim_post &
                          results3_tidy[,3] == method_str_tmp)
          results3_tidy[index[1],4:7] <- c(true_theta,true_alpha,"W1curve",
                                           sqrt(mean(( (W1curve - EW1curve)/W1curve )^2)) )
          results3_tidy[index[2],4:7] <- c(true_theta,true_alpha,"W1",
                                           abs(W1-EW1_f)/W1 )
          results3_tidy[index[3],4:7] <- c(true_theta,true_alpha,"A1",
                                           abs(A1-pA1_f) )
          results3_tidy[index[4],4:7] <- c(true_theta,true_alpha,"K",
                                           abs(K-EK_f)/K )
          

          if(A1){
            RMSE_W1_A1_tmp <- abs(W1-res_f[4])/W1
            RMSE_W1_B1_tmp <- NA
          } else {
            RMSE_W1_A1_tmp <- NA
            RMSE_W1_B1_tmp <- abs(W1-res_f[5])/W1
          }
          results3_tidy[index[5],4:7] <- c(true_theta,true_alpha,"W1_A1",
                                           RMSE_W1_A1_tmp )
          results3_tidy[index[6],4:7] <- c(true_theta,true_alpha,"W1_B1",
                                           RMSE_W1_B1_tmp )
        } else {
          warning("normal error not implemented!")
        }
      }
    }
  }
}
results3_tidy$error <- "percentage"
results3_tidy$true_alpha <- as.numeric(results3_tidy$true_alpha)
results3_tidy$true_theta <- as.numeric(results3_tidy$true_theta)
results3_tidy$RMSE <- as.numeric(results3_tidy$RMSE)

######################################################
######################################################

save(list = c("results3_tidy","results_tidy","results_all"), file = "results/results_model.Rdata")
