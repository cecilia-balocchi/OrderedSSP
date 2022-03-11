rm(list = ls())

source("scripts/funs.R")
Rcpp::sourceCpp("scripts/funs.cpp")
library(VGAM)

Nsim <- 100
sim_x_batch <- 10
n_batches <- 10
Nsim_post <- 25

n <- 500
m <- 5000
Ns <- 1:(n/10)*10
Ms <- 1:(m/10)*10

set.seed(20210916)
params <- rbeta(Nsim, 4,2)
Ktot <- m*10 # for finite zipf (which was not included in paper)



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

# "zipf" , "finzipf" (zipf is actually the zeta distribution, finzipf is the finite-support zipf distribution)
# "finzipf" was not used in the paper
distr_str <- "zipf" # "zipf", "finzipf"
which_order <- "arrival_weight" # "arrival_weight", "astable"
# we ran both for "arrival_weight" and "astable"

for(batch in 1:n_batches){
  rstudioapi::jobRunScript("scripts/simulate_zipf.R", name = paste0("simulate_zipf_batch",batch), importEnv = TRUE)
}

##

percentage_error <- TRUE
distr_str <- "zipf" 


iter <- 1
for(order_str in c("arrival_weight", "astable")){
  which_order = order_str
  results_all <- data.frame(distr = NA, 
                            order = NA,
                            batch = rep(1:n_batches, each = Nsim_post*sim_x_batch),
                            sim = rep(1:Nsim, each = Nsim_post), 
                            # true_theta = NA,true_alpha = NA, 
                            true_par = NA,
                            theta_ord = NA,alpha_ord = NA, 
                            theta_std = NA,alpha_std = NA, 
                            theta_ordDP = NA,alpha_ordDP = NA, 
                            theta_lsK = NA,alpha_lsK = NA, 
                            theta_lsX1 = NA,alpha_lsX1 = NA, 
                            K = NA, X1 = NA, Kpost = NA, A1 = NA, W1 = NA)
  
  for(batch in 1:n_batches){
    sims = (batch-1)*sim_x_batch*Nsim_post + 1:(sim_x_batch*Nsim_post)
    filestr <- paste0("results/",ifelse(distr_str %in% c("zipf","finzipf"),"simzipf/","simdir/"))
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

percentage_error <- TRUE
results_tidy <- data.frame(order = NA, sim = NA, true_par = NA, 
                            param = NA, method = NA, RMSE = NA)
order_tmp <- 0
for(order_str in c("arrival_weight", "astable", "inv", "invsqrt")){
  which_order = order_str
  for(batch in 1:n_batches){
    sims = (batch-1)*sim_x_batch + 1:sim_x_batch
    filestr <- paste0("results/",ifelse(distr_str %in% c("zipf","finzipf"),"simzipf/","simdir/"))
    filestr <- paste0(filestr, "sim",distr_str)
    filestr <- paste0(filestr, "_",which_order,"order")
    filestr <- paste0(filestr, "_n500_m5000_batch",batch,".rdata")
    tmpload = load(filestr)
    
    W1_curves <- get(paste0("W1_curves_",batch))
    res_fs <- get(paste0("res_fs_",batch))
    W1curve_fs <- get(paste0("W1curve_fs_",batch))
    
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
        
        res_f <- res_fs[[paste0(sim,method_str[i_method])]] ## it's the same for all sim_post
        EK_f <- res_f[1]; EW1_f <- res_f[2]
        pA1_f <- res_f[3]
        EW1curve <- W1curve_fs[[paste0(sim,method_str[i_method])]]
        
        if(i_method == 1) results_comb$W1_pred[order_tmp*Nsim*Nsim_post+ (sim-1)*Nsim_post + 1:Nsim_post] <- EW1_f
        if(i_method == 1) results_comb$K_pred[order_tmp*Nsim*Nsim_post+ (sim-1)*Nsim_post + 1:Nsim_post] <- EK_f
        
        for(sim_post in 1:Nsim_post){
          index_sim_post <- order_tmp*Nsim*Nsim_post+ (sim-1)*Nsim_post + sim_post
          W1 <- results_comb[index_sim_post, "W1"]
          A1 <- results_comb[index_sim_post, "A1"]
          W1curve <- W1_curves[[paste0(sim,"_",sim_post)]]
          K <- results_comb[index_sim_post, "K"]
          if(percentage_error){
            RMSE_W1curve_tmp <- RMSE_W1curve_tmp +
              sqrt(mean(( (W1curve - EW1curve)/W1curve )^2))
            RMSE_W1_tmp <- RMSE_W1_tmp + abs(W1-EW1_f)/W1
            RMSE_A1_tmp <- RMSE_A1_tmp + abs(A1-pA1_f)
            RMSE_K_tmp <- RMSE_K_tmp + abs(K-EK_f)/K
            if(A1){
              RMSE_W1_A1_tmp <- c(RMSE_W1_A1_tmp,abs(W1-res_f[4]))/W1
              RMSE_W1_B1_tmp <- c(RMSE_W1_B1_tmp,NA)
            } else {
              RMSE_W1_A1_tmp <- c(RMSE_W1_A1_tmp,NA)
              RMSE_W1_B1_tmp <- c(RMSE_W1_B1_tmp,abs(W1-res_f[5]))/W1
            }
          } else {
            RMSE_W1curve_tmp <- RMSE_W1curve_tmp +
              sqrt(mean(( W1curve - EW1curve )^2))
            RMSE_W1_tmp <- RMSE_W1_tmp + abs(W1-EW1_f)
            RMSE_A1_tmp <- RMSE_A1_tmp + abs(A1-pA1_f)
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
          K_tmp <- K_tmp + results_comb[index_sim_post, "K"]
          if(A1){
            W1_A1_tmp <- c(W1_A1_tmp,W1)
            W1_B1_tmp <- c(W1_B1_tmp,NA)
          } else {
            W1_A1_tmp <- c(W1_A1_tmp,NA)
            W1_B1_tmp <- c(W1_B1_tmp,W1)
          }
          
        }
        
        tmp = list(order_str,sim,results_comb$true_par[index_sim_post],"W1curve",substr(method_str,2,10)[i_method],RMSE_W1curve_tmp/Nsim_post)
        names(tmp) <- colnames(results_tidy)
        results_tidy <- rbind(results_tidy, tmp)
        
        tmp = list(order_str,sim,results_comb$true_par[index_sim_post],"W1",substr(method_str,2,10)[i_method],RMSE_W1_tmp/Nsim_post)
        names(tmp) <- colnames(results_tidy)
        results_tidy <- rbind(results_tidy, tmp)
        
        tmp = list(order_str,sim,results_comb$true_par[index_sim_post],"A1",substr(method_str,2,10)[i_method],RMSE_A1_tmp/Nsim_post)
        names(tmp) <- colnames(results_tidy)
        results_tidy <- rbind(results_tidy, tmp)
        
        tmp = list(order_str,sim,results_comb$true_par[index_sim_post],"K",substr(method_str,2,10)[i_method],RMSE_K_tmp/Nsim_post)
        names(tmp) <- colnames(results_tidy)
        results_tidy <- rbind(results_tidy, tmp)
        
        tmp = list(order_str,sim,results_comb$true_par[index_sim_post],"W1_A1",substr(method_str,2,10)[i_method],mean(RMSE_W1_A1_tmp, na.rm = T))
        names(tmp) <- colnames(results_tidy)
        results_tidy <- rbind(results_tidy, tmp)
        
        tmp = list(order_str,sim,results_comb$true_par[index_sim_post],"W1_B1",substr(method_str,2,10)[i_method],mean(RMSE_W1_B1_tmp, na.rm = T))
        names(tmp) <- colnames(results_tidy)
        results_tidy <- rbind(results_tidy, tmp)
        
        ### 
        if(percentage_error){
          tmp1 = list(order_str,sim,results_comb$true_par[index_sim_post],"W1curve2",substr(method_str,2,10)[i_method],
                      sqrt(mean(( (W1_curve_tmp/Nsim_post - W1curve_fs[[paste0(sim,method_str[i_method])]])/(W1_curve_tmp/Nsim_post) )^2)) )
          tmp2 = list(order_str,sim,results_comb$true_par[index_sim_post],"W12",substr(method_str,2,10)[i_method],abs(W1_tmp/Nsim_post-EW1_f)/(W1_tmp/Nsim_post))
          tmp3 = list(order_str,sim,results_comb$true_par[index_sim_post],"K2",substr(method_str,2,10)[i_method],abs(K_tmp/Nsim_post-EK_f)/(K_tmp/Nsim_post))
          tmp4 = list(order_str,sim,results_comb$true_par[index_sim_post],"W1_A12",substr(method_str,2,10)[i_method],abs(mean(W1_A1_tmp, na.rm = TRUE) - res_f[4])/mean(W1_A1_tmp, na.rm = TRUE))
          tmp5 = list(order_str,sim,results_comb$true_par[index_sim_post],"W1_B12",substr(method_str,2,10)[i_method],abs(mean(W1_B1_tmp, na.rm = T) - res_f[5])/mean(W1_B1_tmp, na.rm = TRUE))
          tmp6 = list(order_str,sim,results_comb$true_par[index_sim_post],"A12",substr(method_str,2,10)[i_method],abs(A1_tmp/Nsim_post-pA1_f)/(A1_tmp/Nsim_post))
          
        } else {
          tmp1 = list(order_str,sim,results_comb$true_par[index_sim_post],"W1curve2",substr(method_str,2,10)[i_method],
                      sqrt(mean(( W1_curve_tmp/Nsim_post - W1curve_fs[[paste0(sim,method_str[i_method])]] )^2)))
          tmp2 = list(order_str,sim,results_comb$true_par[index_sim_post],"W12",substr(method_str,2,10)[i_method],abs(W1_tmp/Nsim_post-EW1_f))
          tmp3 = list(order_str,sim,results_comb$true_par[index_sim_post],"K2",substr(method_str,2,10)[i_method],abs(K_tmp/Nsim_post-EK_f))
          tmp4 = list(order_str,sim,results_comb$true_par[index_sim_post],"W1_A12",substr(method_str,2,10)[i_method],abs(mean(W1_A1_tmp, na.rm = TRUE) - res_f[4]))
          tmp5 = list(order_str,sim,results_comb$true_par[index_sim_post],"W1_B12",substr(method_str,2,10)[i_method],abs(mean(W1_B1_tmp, na.rm = T) - res_f[5]))
          tmp6 = list(order_str,sim,results_comb$true_par[index_sim_post],"A12",substr(method_str,2,10)[i_method],abs(A1_tmp/Nsim_post-pA1_f))
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
  order_tmp <- order_tmp + 1
}
results_tidy <- results_tidy[-1,]
results_tidy$error = "percentage"



################################################################
##### let's try something different: do not take the mean ######
################################################################
percentage_error <- TRUE
results3_tidy <- data.frame(sim = rep(1:Nsim, each = Nsim_post*6*length(method_str)),
                            sim_post = rep(rep(1:Nsim_post, each = 6), Nsim*length(method_str)),
                            method = rep(rep(substr(method_str,2,10), each = Nsim_post*6),Nsim),
                            true_par = NA,
                            param = NA,  RMSE = NA)
for(order_str in c("arrival_weight", "astable")){
  which_order = order_str
  for(batch in 1:n_batches){
    sims = (batch-1)*sim_x_batch + 1:sim_x_batch
    filestr <- paste0("results/",ifelse(distr_str %in% c("zipf","finzipf"),"simzipf/","simdir/"))
    filestr <- paste0(filestr, "sim",distr_str)
    filestr <- paste0(filestr, "_",which_order,"order")
    filestr <- paste0(filestr, "_n500_m5000_batch",batch,".rdata")
    tmpload = load(filestr)
    
    W1_curves <- get(paste0("W1_curves_",batch))
    res_fs <- get(paste0("res_fs_",batch))
    W1curve_fs <- get(paste0("W1curve_fs_",batch))
    
    for(sim in sims){
      for(i_method in 1:5){
        res_f <- res_fs[[paste0(sim,method_str[i_method])]]
        EK_f <- res_f[1]; EW1_f <- res_f[2]
        pA1_f <- res_f[3]
        EW1curve <- W1curve_fs[[paste0(sim,method_str[i_method])]]
        
        method_str_tmp = substr(method_str,2,10)[i_method]
        true_par <- results_comb[(sim-1)*Nsim_post + 1, "true_par"]
        
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
            results3_tidy[index[1],4:6] <- c(true_par,"W1curve",
                                             sqrt(mean(( (W1curve - EW1curve)/W1curve )^2)) )
            results3_tidy[index[2],4:6] <- c(true_par,"W1",
                                             abs(W1-EW1_f)/W1 )
            results3_tidy[index[3],4:6] <- c(true_par,"A1",
                                             abs(A1-pA1_f) )
            results3_tidy[index[4],4:6] <- c(true_par,"K",
                                             abs(K-EK_f)/K )
            
            if(A1){
              RMSE_W1_A1_tmp <- abs(W1-res_f[4])/W1
              RMSE_W1_B1_tmp <- NA
            } else {
              RMSE_W1_A1_tmp <- NA
              RMSE_W1_B1_tmp <- abs(W1-res_f[5])/W1
            }
            results3_tidy[index[5],4:6] <- c(true_par,"W1_A1",
                                             RMSE_W1_A1_tmp )
            results3_tidy[index[6],4:6] <- c(true_par,"W1_B1",
                                             RMSE_W1_B1_tmp )
          } else {
            warning("normal error not implemented!")
          }
        }
      }
    }
  }
  results3_tidy$error <- "percentage"
  results3_tidy$true_par <- as.numeric(results3_tidy$true_par)
  results3_tidy$RMSE <- as.numeric(results3_tidy$RMSE)
  # results3_tidy$order <- order_str
  assign(paste0("results3_tidy_",order_str),results3_tidy)
}

results3_tidy_astable$order <- "astable"
results3_tidy_arrival_weight$order <- "arrival_weight"
results3_tidy <- rbind(results3_tidy_astable, results3_tidy_arrival_weight)


save(list = c("results3_tidy","results_tidy","results_comb"), file = "results/results_zipf.Rdata")
