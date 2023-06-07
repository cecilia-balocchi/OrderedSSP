rm(list = ls())
library(optimization)

source("scripts/funs.R")
tmp = load("data/Gen1000_cleaned_BRCA.RData")

methods_str <- c("ordDP","ordPYP","stdPYP","lsX1","lsK","FB")
methods_str2 <- c("ordDP","ordPYP","stdPYP","lsM1","lsK","FB")

set.seed(123)
seeds <- sample(10000,size = 100)

n_obs <- sum(Gen1000_gt_keep01)
mut_index_array <- which(Gen1000_gt_keep01 == 1, arr.ind = T)

####

train_size <- round(n_obs/20) 
test_size <-  n_obs - train_size 
N = train_size
M = test_size
X1ns <- floor(seq(1, N, length = 100))
Kns <- floor(seq(1, N, length = 100))
Kns_post <- floor(seq(1, M, length = 25))
W1ns <- floor(seq(1, M, length = 25))

for(i_seed in 1:length(seeds)){
  seed <- seeds[i_seed]; set.seed(seed)
  train_index <- sample(n_obs, size = train_size, replace = F)
  test_index <- sample(setdiff(1:n_obs,train_index))
  
  mut_index_train = mut_index_array[train_index,]
  Gen1000_gt_train01 <- Gen1000_gt_keep01*0
  Gen1000_gt_train01[mut_index_train] <- 1
  
  mut_index_test = mut_index_array[test_index,]
  totcounts_train <- rowSums(Gen1000_gt_train01)
  
  
  variants_train <- mut_index_train[,1]
  rank_variants_train <- Gen1000_fix_keep$rank[variants_train]
  Kcurve2 <- sapply(Kns, function(m) length(unique(variants_train[1:m])))
  K.pre <- length(unique(variants_train))
  X1 <- sum(rank_variants_train == max(rank_variants_train))
  
  X1s <- sapply(X1ns, function(x) sum(rank_variants_train[1:x] == max(rank_variants_train[1:x])))
  
  ms <- totcounts_train[order(Gen1000_fix_keep$AGE_MODE, decreasing = FALSE)]
  ms <- ms[ms > 0]
  
  variants_test <- mut_index_test[,1]
  rank_variants_test <- Gen1000_fix_keep$rank[variants_test]
  
  Kcurve2_post <- sapply(Kns_post, function(m) length(unique(c(variants_test[1:m])))) 
  K.post <- Kcurve2_post[length(Kcurve2_post)]
  
  max_pre <- max(rank_variants_train)
  W1 <- ifelse(max(rank_variants_test) > max_pre,
               sum(rank_variants_test == max(rank_variants_test)), 
               sum(rank_variants_test == max_pre)+X1 )
  A1 = ifelse(max(rank_variants_test) > max_pre, TRUE, FALSE)
  W1curve2 <- sapply(W1ns, function(x) ifelse(max(rank_variants_test[1:x]) > max_pre,
                                              sum(rank_variants_test[1:x] == max(rank_variants_test[1:x])), 
                                              sum(rank_variants_test[1:x] == max_pre)+X1 ))
  
  assign(paste0("ms","_seed",i_seed),ms)
  assign(paste0("X1s","_seed",i_seed),X1s)
  assign(paste0("Kcurve2","_seed",i_seed),Kcurve2)
  assign(paste0("K.post","_seed",i_seed),K.post)
  assign(paste0("Kcurve2_post","_seed",i_seed),Kcurve2_post)
  assign(paste0("W1","_seed",i_seed),W1)
  assign(paste0("A1","_seed",i_seed),A1)
  assign(paste0("W1curve2","_seed",i_seed),W1curve2)
}

####

fun_ordEPPF_PY <- function(par, ms){
  -ordEPPF_PY(ms, par[1], par[2])
}
fun_ordEPPF_DP <- function(par, ms){
  -ordEPPF_PY(ms, par, 0)
}
fun_EPPF_PY <- function(par, ms){
  -EPPF_PY(ms, par[1], par[2])
}
fun_EPPF_DP <- function(par, ms){
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

thetas <- array(NA, dim = c(length(seeds),5))
alphas <- array(NA, dim = c(length(seeds),5))
thetas_list <- list()
alphas_list <- list()
times_post <- numeric(length(seeds))

mh_niter <- 100500
mh_burnin <- 500
mh_thin <- 100
post_subset <- seq(mh_burnin + 1, mh_niter, by = mh_thin)

for(i_seed in 1:length(seeds)){
  ms <- get(paste0("ms","_seed",i_seed))
  X1s <- get(paste0("X1s","_seed",i_seed))
  Kcurve2 <- get(paste0("Kcurve2","_seed",i_seed))
  
  tmp_opt <- optimize(f = fun_ordEPPF_DP, ms = ms, interval = c(0,1000))
  theta <- tmp_opt$minimum
  alpha <- 0
  i = 1; assign(paste0("theta_",methods_str[i],"_seed",i_seed),theta); assign(paste0("alpha_",methods_str[i],"_seed",i_seed),alpha)
  thetas[i_seed,i] <- theta; alphas[i_seed,i] <- alpha
  
  tmp_nm <- optim(par = c(1,0.5), fn = fun_ordEPPF_PY, ms = ms, method = "Nelder-Mead")
  theta <- tmp_nm$par[1]
  alpha <- tmp_nm$par[2]
  i = 2; assign(paste0("theta_",methods_str[i],"_seed",i_seed),theta); assign(paste0("alpha_",methods_str[i],"_seed",i_seed),alpha)
  thetas[i_seed,i] <- theta; alphas[i_seed,i] <- alpha
  
  tmp_nm <- optim(par = c(1,0.5), fn = fun_EPPF_PY, ms = ms, method = "Nelder-Mead")
  theta <- tmp_nm$par[1]
  alpha <- tmp_nm$par[2]
  i = 3; assign(paste0("theta_",methods_str[i],"_seed",i_seed),theta); assign(paste0("alpha_",methods_str[i],"_seed",i_seed),alpha)
  thetas[i_seed,i] <- theta; alphas[i_seed,i] <- alpha
  
  tmp_nm <- optim(par = c(1,0.5), fn = fun_X1LS, X1s = X1s, Ns = X1ns, method = "Nelder-Mead")
  theta <- tmp_nm$par[1]
  alpha <- tmp_nm$par[2]
  i = 4; assign(paste0("theta_",methods_str[i],"_seed",i_seed),theta); assign(paste0("alpha_",methods_str[i],"_seed",i_seed),alpha)
  thetas[i_seed,i] <- theta; alphas[i_seed,i] <- alpha
  
  tmp_nm <- optim(par = c(1,0.5), fn = fun_Kls, Ks = Kcurve2, ns = Kns, method = "Nelder-Mead")
  theta <- tmp_nm$par[1]
  alpha <- tmp_nm$par[2]
  i = 5; assign(paste0("theta_",methods_str[i],"_seed",i_seed),theta); assign(paste0("alpha_",methods_str[i],"_seed",i_seed),alpha)
  thetas[i_seed,i] <- theta; alphas[i_seed,i] <- alpha
  
  time2 <- Sys.time()
  post <- posterior_MH_hyper(ms, niter = mh_niter)
  times_post[i_seed] <- Sys.time()-time2
  
  thetas_post <- post$thetas
  alphas_post <- post$alphas
  thetas_thin <- thetas_post[post_subset]
  alphas_thin <- alphas_post[post_subset]
  
  thetas_list[[i_seed]] <- thetas_thin
  alphas_list[[i_seed]] <- alphas_thin
}

#### this computes the number of distinct in the additional sample (not conditional on the first sample)
#### we don't need to compute the CI for every seed, given that it takes so long.

Ks_perr <- array(0, dim = c(length(seeds), length(methods_str)))
Ks_curve_perr <- array(0, dim = c(length(seeds), length(methods_str)))
colnames(Ks_perr) <- methods_str
colnames(Ks_curve_perr) <- methods_str

Kcurve_true <- array(0, dim = c(length(seeds), length(Kns_post)))
Kcurve_est <- array(0, dim = c(length(seeds), length(Kns_post)))
for(i in 1:length(methods_str)){
  assign(paste0("Kcurve_",methods_str[i]),Kcurve_est)
}
# KcurveCIs_est <- array(0, dim = c(length(seeds), length(Kns_post),2))
# for(i in 1:length(methods_str)){
#   assign(paste0("KcurveCIs_",methods_str[i]),KcurveCIs_est)
# }

for(i_seed in 1:length(seeds)){
  for(i in 1:5){
    theta = get(paste0("theta_",methods_str[i],"_seed",i_seed))
    alpha = get(paste0("alpha_",methods_str[i],"_seed",i_seed))

    Kpre <- get(paste0("Kcurve2","_seed",i_seed))[length(Kns)]
    
    EKcurve <- EK_PY(Kns_post, theta, alpha)
    if(alpha < 0.01){
      EKcurve <- theta * log(Kns_post)
    }
    EK <- EKcurve[length(EKcurve)]
    
    assign(paste0("EK_",methods_str[i]),EK)
    assign(paste0("EKcurve_",methods_str[i]),EKcurve)

    Kcurve_est <- get(paste0("Kcurve_",methods_str[i]))
    Kcurve_est[i_seed, ] <- EKcurve
    assign(paste0("Kcurve_",methods_str[i]),Kcurve_est)
  }
  i = 6
    
  theta = thetas_list[[i_seed]]
  alpha = alphas_list[[i_seed]]
  Kpre <- get(paste0("Kcurve2","_seed",i_seed))[length(Kns)]
  EKcurve_post <- sapply(Kns_post, function(ns) EK_PY(ns, theta, alpha))
  # if(mean(alpha) < 0.01){
  #   EKcurve_post <- sapply(Kns_post, function(ns) theta * log(ns))
  # }
  EKcurve <- colMeans(EKcurve_post)
  EK <- EKcurve[length(EKcurve)]
  assign(paste0("EK_",methods_str[i]),EK)
  assign(paste0("EKcurve_",methods_str[i]),EKcurve)
  Kcurve_est <- get(paste0("Kcurve_",methods_str[i]))
  Kcurve_est[i_seed, ] <- EKcurve
  assign(paste0("Kcurve_",methods_str[i]),Kcurve_est)
  
  ## end of i = 6
  
  K.post <- get(paste0("K.post","_seed",i_seed))
  Ks_perr[i_seed,] <- sapply(sapply(paste0("EK_",methods_str), function(x) get(x)),
                             function(x) sqrt(mean((K.post-x)^2/K.post^2)))
  Kcurve2_post <- get(paste0("Kcurve2_post","_seed",i_seed))
  Ks_curve_perr[i_seed,] <- sapply(lapply(paste0("EKcurve_",methods_str), function(x) get(x)),
                                   function(x) sqrt(mean((Kcurve2_post-x)^2/Kcurve2_post^2)))
  Kcurve_true[i_seed,] <- Kcurve2_post
}
Kcurve_avg <- colMeans(Kcurve_true)


#####

W1_perr <- array(0, dim = c(length(seeds), length(methods_str)))
# W1A1_perr <- array(NA, dim = c(length(seeds), length(methods_str)))
# W1B1_perr <- array(NA, dim = c(length(seeds), length(methods_str)))
W1curve_perr <- array(0, dim = c(length(seeds), length(methods_str)))
colnames(W1_perr) <- methods_str
colnames(W1curve_perr) <- methods_str

W1curve_avg <- W1ns*0
W1curve_true <- array(0, dim = c(length(seeds), length(W1ns)))
W1curve_est <- array(0, dim = c(length(seeds), length(W1ns)))
for(i in 1:length(methods_str)){
  assign(paste0("W1curve_",methods_str[i]),W1curve_est)
}

for(i_seed in 1:length(seeds)){
  X1s <- get(paste0("X1s","_seed",i_seed))
  X1 <- X1s[length(X1s)]
  for(i in 1:5){
    theta = get(paste0("theta_",methods_str[i],"_seed",i_seed))
    alpha = get(paste0("alpha_",methods_str[i],"_seed",i_seed))
    # pW1_M <- sapply(W1ns, function(w) exp(lpW1_comb(w,X1,N,M, theta, alpha))) # NOT USED?
    
    EW1tmp <- EW1(X1,N,M, theta, alpha); EW1curve <- EW1(X1,N,W1ns, theta, alpha)
    
    # EW1A1 <- EW1_A1(N,M,theta, alpha); EW1B1 <- EW1_B1(X1,N,M, theta, alpha)
    assign(paste0("EW1_",methods_str[i]),EW1tmp); assign(paste0("EW1curve_",methods_str[i]),EW1curve)
    # assign(paste0("EW1A1_",methods_str[i]),EW1A1); assign(paste0("EW1B1_",methods_str[i]),EW1B1)
    
    W1curve_est <- get(paste0("W1curve_",methods_str[i]))
    W1curve_est[i_seed, ] <- EW1curve
    assign(paste0("W1curve_",methods_str[i]),W1curve_est)
  }
  
  ## i = 6
  i = 6
  theta = thetas_list[[i_seed]]
  alpha = alphas_list[[i_seed]]
  
  EW1tmp <- mean(EW1(X1,N,M, theta, alpha)); 
  EW1curve_post <- sapply(W1ns, function(ns) EW1(X1,N,ns, theta, alpha))
  EW1curve <- colMeans(EW1curve_post)
  
  assign(paste0("EW1_",methods_str[i]),EW1tmp); assign(paste0("EW1curve_",methods_str[i]),EW1curve)
  W1curve_est <- get(paste0("W1curve_",methods_str[i]))
  W1curve_est[i_seed, ] <- EW1curve
  assign(paste0("W1curve_",methods_str[i]),W1curve_est)
  
  ## end of i = 6
  
  W1 <- get(paste0("W1","_seed",i_seed))
  W1curve2 <- get(paste0("W1curve2","_seed",i_seed))
  
  W1_perr[i_seed,] <- sapply(sapply(paste0("EW1_",methods_str), function(x) get(x)), 
                             function(x) abs((W1-x)^2/W1^2))
  W1curve_perr[i_seed,] <- sapply(lapply(paste0("EW1curve_",methods_str), function(x) get(x)), 
                                  function(x) sqrt(mean((W1curve2-x)^2/W1curve2^2))) 
  W1curve_avg <- W1curve_avg + W1curve2
  W1curve_true[i_seed,] <- W1curve2
}
W1curve_avg <- W1curve_avg / length(seeds)


## compute CI for the needed seeds
## this is actually not used, we use the empirical quantiles

median_seeds_W1curve <- numeric(length(methods_str))
median_seeds_Kcurve <- numeric(length(methods_str))
for(i_method in 1:length(methods_str)){
  median_seeds_W1curve[i_method] = which.min(abs(W1curve_perr[,i_method] - median(W1curve_perr[,i_method])))
  median_seeds_Kcurve[i_method] = which.min(abs(Ks_curve_perr[,i_method] - median(Ks_curve_perr[,i_method])))
}

KcurveCIs_est <- array(0, dim = c(length(methods_str), length(Kns_post),2))
W1curveCIs_est <- array(0, dim = c(length(methods_str), length(W1ns),2))

for(i_method in 1:5){
  i_seed = median_seeds_W1curve[i_method]
  theta = get(paste0("theta_",methods_str[i_method],"_seed",i_seed))
  alpha = get(paste0("alpha_",methods_str[i_method],"_seed",i_seed))
  X1s <- get(paste0("X1s","_seed",i_seed))
  X1 <- X1s[length(X1s)]
  
  for(j in 1:length(W1ns)){
    Mi <- W1ns[j]+X1
    if(j == 1){
      lower_bound_temp = 1
      Mi_temp = Mi
    } else {
      lower_bound_temp = max(1,W1curveCIs_est[i_method,j-1,1] - 100)
      Mi_temp = W1curveCIs_est[i_method,j-1,2] + 100
    }
    pi <- exp(sapply(lower_bound_temp:Mi_temp, function(x) lpW1_comb(x,X1,N,Mi,theta,alpha)))
    while(sum(pi) < 0.95){
      lower_bound_temp = max(1,lower_bound_temp - 100)
      Mi_temp = Mi_temp + 100
      pi <- exp(sapply(lower_bound_temp:Mi_temp, function(x) lpW1_comb(x,X1,N,Mi,theta,alpha)))
    }
    pi_sorted = sort(pi, decreasing = T)
    threshold = min(pi_sorted[1:min(which(cumsum(pi_sorted) > 0.95))])
    index_tmp = which(pi>=threshold)
    index = (lower_bound_temp:Mi_temp)[index_tmp]
    W1curveCIs_est[i_method,j,] <- c(min(which(index[-1]-index[-length(index)] == 1)),max(index))
  }
}
###
i_method = 6 
i_seed = median_seeds_W1curve[i_method]
theta = thetas_list[[i_seed]]
alpha = alphas_list[[i_seed]]
X1s <- get(paste0("X1s","_seed",i_seed))
X1 <- X1s[length(X1s)]
for(j in 1:length(W1ns)){
  Mi <- W1ns[j]+X1
  if(j == 1){
    lower_bound_temp = 1
    Mi_temp = Mi
  } else {
    lower_bound_temp = max(1,W1curveCIs_est[i_method,j-1,1] - 100)
    Mi_temp = W1curveCIs_est[i_method,j-1,2] + 100
  }
  pi <- colMeans(  exp(sapply(lower_bound_temp:Mi_temp, function(x) lpW1_comb(x,X1,N,Mi,theta,alpha)))  ) ## maybe?
  while(sum(pi) < 0.95){
    lower_bound_temp = max(1,lower_bound_temp - 100)
    Mi_temp = Mi_temp + 100
    pi <- colMeans(  exp(sapply(lower_bound_temp:Mi_temp, function(x) lpW1_comb(x,X1,N,Mi,theta,alpha)))  ) ## maybe?
  }
  pi_sorted = sort(pi, decreasing = T)
  threshold = min(pi_sorted[1:min(which(cumsum(pi_sorted) > 0.95))])
  index_tmp = which(pi>=threshold)
  index = (lower_bound_temp:Mi_temp)[index_tmp]
  W1curveCIs_est[i_method,j,] <- c(min(which(index[-1]-index[-length(index)] == 1)),max(index))
}
### end i_method= 6

save(list = c("Ks_perr","W1_perr","Kns_post","W1ns","N","M","alphas","thetas",
              "W1curve_perr",paste0("W1curve_",methods_str),#paste0("CIs_",methods_str), 
              "W1curve_true", "W1curve_avg",
              "KcurveCIs_est","W1curveCIs_est",
              "median_seeds_W1curve","median_seeds_Kcurve",
              "Ks_curve_perr",paste0("Kcurve_",methods_str),#paste0("KcurveCIs_",methods_str),
              "Kcurve_true","Kcurve_avg"), 
     file = "results/crossval_BRCA_FB.Rdata")


## we don't use these for the analysis
# for(i_method in 1:length(methods_str)){
#   i_seed = median_seeds_Kcurve[i_method]
#   theta = get(paste0("theta_",methods_str[i_method],"_seed",i_seed))
#   alpha = get(paste0("alpha_",methods_str[i_method],"_seed",i_seed))
#   
#   i_m = 1
#   KcurveCIs_est[i_method,i_m,] <- get_CI_Kcurve(1)
#   if(i_method == 1){ 
#     # does not work for large n
#     # for(i_m in 2:length(Kns_post)){
#     #   lpi <- Kn_DP_distribution(Kns_post[i_m],theta)
#     #   KcurveCIs_est[i_method,i_m,] <- get_CI_Kcurve(exp(lpi))
#     # }
#   } else {
#     for(i_m in 2:length(Kns_post)){
#       lpi <- Kn_emp_distribution(Kns_post[i_m],theta,alpha,nsamples = 10000)
#       # lpi <- Kn_distribution(Kns_post[i_m],theta, alpha)
#       KcurveCIs_est[i_method,i_m,] <- get_CI_Kcurve(exp(lpi))
#     }
#   }
# }
# 
# 
# save(list = c("Ks_perr","W1_perr","Kns_post","W1ns","N","M","alphas","thetas",
#               "W1curve_perr",paste0("W1curve_",methods_str),#paste0("CIs_",methods_str), 
#               "W1curve_true", "W1curve_avg",
#               "KcurveCIs_est","W1curveCIs_est",
#               "median_seeds_W1curve","median_seeds_Kcurve",
#               "Ks_curve_perr",paste0("Kcurve_",methods_str),#paste0("KcurveCIs_",methods_str),
#               "Kcurve_true","Kcurve_avg"), 
#      file = "results/crossval_BRCA.Rdata")
# # tmp = load("results/crossval_BRCA.Rdata")


## using empirical quantiles
png("figures/median_KcurveQ_BRCA.png", width = 10, height = 2.5, units = "in", res = 300)
par(mfrow = c(1,5), oma = c(0, 0, 2, 0), mar= c(3.1,4.1,2.1,2.1))
for(i_method in c(3,2,1,4,5)){
  i_seed = which.min(abs(Ks_curve_perr[,i_method] - median(Ks_curve_perr[,i_method])))
  Kcurve_est <- get(paste0("Kcurve_",methods_str[i_method]))
  Kcurve_mean <- Kcurve_est[i_seed,]
  lCI <- apply(Kcurve_est, MARGIN = 2, quantile, probs = 0.025)
  uCI <- apply(Kcurve_est, MARGIN = 2, quantile, probs = 0.975)
  plot(Kns_post, Kcurve_true[i_seed,], type ="l", ylim = range(c(Kcurve_est,Kcurve_true[i_seed,])),
       main = methods_str2[i_method], xlab = "", ylab = "K")
  mtext("samples", side = 1, line = 2, cex = 0.8,outer = FALSE)
  polygon(c(rev(Kns_post), Kns_post),c(rev(lCI), uCI), col = rgb(1,0,0,0.3), border = NA)
  lines(Kns_post, Kcurve_mean, lwd = 2, col = "red")
  lines(Kns_post, Kcurve_true[i_seed,], lwd = 2)
}
mtext("Prediction of K as a function of the number of samples", side = 3, line = -0, cex = 1.2,outer = TRUE)
dev.off()

png("figures/median_W1curveQ_BRCA.png", width = 10, height = 2.5, units = "in", res = 300)
par(mfrow = c(1,5), oma = c(0, 0, 2, 0), mar= c(3.1,4.1,2.1,2.1))
for(i_method in c(3,2,1,4,5)){
  i_seed = which.min(abs(W1curve_perr[,i_method] - median(W1curve_perr[,i_method])))
  W1curve_est <- get(paste0("W1curve_",methods_str[i_method]))
  W1curve_mean <- W1curve_est[i_seed,]
  # sds <- apply(W1curve_est, MARGIN = 2, sd)
  lCI <- apply(W1curve_est, MARGIN = 2, quantile, probs = 0.025)
  uCI <- apply(W1curve_est, MARGIN = 2, quantile, probs = 0.975)
  plot(W1ns, W1curve_true[i_seed,], type ="l", ylim = range(c(W1curve_est,W1curve_true[i_seed,])),
       main = methods_str2[i_method], xlab = "", ylab = "W1")
  mtext("samples", side = 1, line = 2, cex = 0.8,outer = FALSE)
  polygon(c(rev(W1ns), W1ns),c(rev(lCI), uCI), col = rgb(1,0,0,0.3), border = NA)
  lines(W1ns, W1curve_mean, lwd = 2, col = "red")
  lines(W1ns, W1curve_true[i_seed,], lwd = 2)
}
mtext("Prediction of W1 as a function of the number of samples", side = 3, line = -0, cex = 1.2,outer = TRUE)
dev.off()

## FB
# tmp = load("results/crossval_BRCA_FB.Rdata")

png("figures/median_KcurveQ_BRCA_FB.png", width = 10, height = 2.5, units = "in", res = 300)
par(mfrow = c(1,6), oma = c(0, 0, 2, 0), mar= c(3.1,4.1,2.1,2.1))
for(i_method in c(6,2,1,3,4,5)){
  i_seed = which.min(abs(Ks_curve_perr[,i_method] - median(Ks_curve_perr[,i_method])))
  Kcurve_est <- get(paste0("Kcurve_",methods_str[i_method]))
  Kcurve_mean <- Kcurve_est[i_seed,]
  lCI <- apply(Kcurve_est, MARGIN = 2, quantile, probs = 0.025)
  uCI <- apply(Kcurve_est, MARGIN = 2, quantile, probs = 0.975)
  plot(Kns_post, Kcurve_true[i_seed,], type ="l", ylim = range(c(Kcurve_est,Kcurve_true[i_seed,])),
       main = methods_str2[i_method], xlab = "", ylab = "K")
  mtext("samples", side = 1, line = 2, cex = 0.8,outer = FALSE)
  polygon(c(rev(Kns_post), Kns_post),c(rev(lCI), uCI), col = rgb(1,0,0,0.3), border = NA)
  lines(Kns_post, Kcurve_mean, lwd = 2, col = "red")
  lines(Kns_post, Kcurve_true[i_seed,], lwd = 2)
}
mtext("Prediction of K as a function of the number of samples", side = 3, line = -0, cex = 1.2,outer = TRUE)
dev.off()


png("figures/median_W1curveQ_BRCA_FB.png", width = 10, height = 2.5, units = "in", res = 300)
par(mfrow = c(1,6), oma = c(0, 0, 2, 0), mar= c(3.1,4.1,2.1,2.1))
for(i_method in c(6,2,1,3,4,5)){
  i_seed = which.min(abs(W1curve_perr[,i_method] - median(W1curve_perr[,i_method])))
  W1curve_est <- get(paste0("W1curve_",methods_str[i_method]))
  W1curve_mean <- W1curve_est[i_seed,]
  # sds <- apply(W1curve_est, MARGIN = 2, sd)
  lCI <- apply(W1curve_est, MARGIN = 2, quantile, probs = 0.025)
  uCI <- apply(W1curve_est, MARGIN = 2, quantile, probs = 0.975)
  plot(W1ns, W1curve_true[i_seed,], type ="l", ylim = range(c(W1curve_est,W1curve_true[i_seed,])),
       main = methods_str2[i_method], xlab = "", ylab = "W1")
  mtext("samples", side = 1, line = 2, cex = 0.8,outer = FALSE)
  polygon(c(rev(W1ns), W1ns),c(rev(lCI), uCI), col = rgb(1,0,0,0.3), border = NA)
  lines(W1ns, W1curve_mean, lwd = 2, col = "red")
  lines(W1ns, W1curve_true[i_seed,], lwd = 2)
}
mtext("Prediction of W1 as a function of the number of samples", side = 3, line = -0, cex = 1.2,outer = TRUE)
dev.off()

