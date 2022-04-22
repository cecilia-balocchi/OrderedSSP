rm(list = ls())
library(optimization)
wdstr <- ""
source(paste0(wdstr,"scripts/funs.R"))

methods_str <- c("ordDP","ordPYP","stdPYP","lsX1","lsK")
methods_str2 <- c("ordDP","ordPYP","stdPYP","lsM1","lsK")

tmp = load(paste0(wdstr,"data/bibtex_clean.rdata"))           # bibtex_data
tmp1 = load(paste0(wdstr,"data/ordered_cit_data.rdata"))      # ordered_data
path = "data/SCC2016/"
paperList = read.table(paste0(wdstr,path, "paperList.txt"), header = T, sep = ",")
paperCitAdj = read.table(paste0(wdstr,path, "paperCitAdj.txt"), header = F, sep = " ")
journals <- unique(bibtex_data$JOURNAL)

set.seed(123)
seeds <- sample(10000,size = 100)

tmp_order <- order(bibtex_data$ranking, decreasing = T)
Adj_ord <- paperCitAdj[tmp_order,tmp_order]
List_ord <- paperList[tmp_order,]

data <- ordered_data
n_obs <- nrow(data) # the observations are the citations (== sum(Adj_ord))

####

train_size <- round(n_obs/10) 
test_size <-  n_obs - train_size 
N = train_size
M = test_size
X1ns <- floor(seq(1, N, length = 100))
Kns <- floor(seq(1, N, length = 100))
Kns_post <- floor(seq(1, M, length = 100))
Kns_post_CI <- floor(seq(1, M, length = 25))
W1ns <- floor(seq(1, M, length = 100))
W1ns_CI <- floor(seq(1, M, length = 25))

# i_seed = 1
for(i_seed in 1:length(seeds)){
  seed <- seeds[i_seed]; set.seed(seed)
  train_index <- sample(n_obs, size = train_size, replace = F)
  test_index <- sample(setdiff(1:n_obs,train_index))
  
  data_train = data[train_index,]  
  data_test = data[test_index,]    
  
  ## remember: oldest corresponds to largest rank
  
  Kcurve2 <- sapply(Kns, function(m) 
    length(unique(data_train$rank_cited[1:m])))
  K.pre <- Kcurve2[length(Kcurve2)]
  
  X1s <- sapply(X1ns, function(x) 
    sum(data_train$rank_cited[1:x] == max(data_train$rank_cited[1:x])))
  X1 <- X1s[length(X1s)]
  
  clusters <- data.frame(rank_cited = sort(unique(data_train$rank_cited), decreasing = FALSE))
  clusters$nk <- NA
  for(i in 1:K.pre){
    clusters$nk[i] <- sum(data_train$rank_cited == clusters$rank_cited[i])
  }
  ms <- clusters$nk
  if(sum(ms == 0)> 0){
    warning("some zero in ms!")
  }
  
  Kcurve2_post <- sapply(Kns_post, function(m) 
    length(unique(data_test$rank_cited[1:m]))) 
  K.post <- Kcurve2_post[length(Kcurve2_post)]
  
  max_pre <- max(data_train$rank_cited)
  W1curve2 <- sapply(W1ns, function(x) ifelse(max(data_test$rank_cited[1:x]) > max_pre,
                                              sum(data_test$rank_cited[1:x] == max(data_test$rank_cited[1:x])), 
                                              sum(data_test$rank_cited[1:x] == max_pre)+X1 ))
  W1 <- W1curve2[length(W1curve2)]
  A1 = ifelse(max(data_test$rank_cited) > max_pre, TRUE, FALSE)
  
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
}

#### this computes the number of distinct in the additional sample (not conditional on the first sample)
#### we don't need to compute the CI for every seed, given that it takes so long.

Ks_perr <- array(0, dim = c(length(seeds), 5))
Ks_curve_perr <- array(0, dim = c(length(seeds), 5))
colnames(Ks_perr) <- methods_str
colnames(Ks_curve_perr) <- methods_str

Kcurve_true <- array(0, dim = c(length(seeds), length(Kns_post)))
Kcurve_est <- array(0, dim = c(length(seeds), length(Kns_post)))
for(i in 1:5){
  assign(paste0("Kcurve_",methods_str[i]),Kcurve_est)
}
# KcurveCIs_est <- array(0, dim = c(length(seeds), length(Kns_post_CI),2))
# for(i in 1:5){
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
    
    ## do not compute the CIs
    # KcurveCIs_est <- get(paste0("KcurveCIs_",methods_str[i]))
    # i_m = 1
    # KcurveCIs_est[i_seed,i_m,] <- get_CI_Kcurve(1)
    # if(i == 1){ # need a different function for DP
    #   for(i_m in 2:length(Kns_post_CI)){
    #     lpi <- Kn_DP_distribution(Kns_post_CI[i_m],theta)
    #     KcurveCIs_est[i_seed,i_m,] <- get_CI_Kcurve(exp(lpi))
    #   }
    # } else {
    #   for(i_m in 2:length(Kns_post_CI)){
    #     lpi <- Kn_distribution(Kns_post_CI[i_m],theta, alpha)
    #     KcurveCIs_est[i_seed,i_m,] <- get_CI_Kcurve(exp(lpi))
    #   }
    # }
    # assign(paste0("KcurveCIs_",methods_str[i]),KcurveCIs_est)
    
    assign(paste0("EK_",methods_str[i]),EK)
    assign(paste0("EKcurve_",methods_str[i]),EKcurve)
    
    Kcurve_est <- get(paste0("Kcurve_",methods_str[i]))
    Kcurve_est[i_seed, ] <- EKcurve
    assign(paste0("Kcurve_",methods_str[i]),Kcurve_est)
  }
  K.post <- get(paste0("K.post","_seed",i_seed))
  Ks_perr[i_seed,] <- sapply(sapply(paste0("EK_",methods_str), function(x) get(x)), 
                             function(x) sqrt(mean((K.post-x)^2/K.post^2)))
  Kcurve2_post <- get(paste0("Kcurve2_post","_seed",i_seed))
  Ks_curve_perr[i_seed,] <- sapply(lapply(paste0("EKcurve_",methods_str), function(x) get(x)), 
                                   function(x) sqrt(mean((Kcurve2_post-x)^2/Kcurve2_post^2)))
  Kcurve_true[i_seed,] <- Kcurve2_post
}
Kcurve_avg <- colMeans(Kcurve_true)


###

W1_perr <- array(0, dim = c(length(seeds), 5))
# W1A1_perr <- array(NA, dim = c(length(seeds), 5))
# W1B1_perr <- array(NA, dim = c(length(seeds), 5))
W1curve_perr <- array(0, dim = c(length(seeds), 5))
colnames(W1_perr) <- methods_str
colnames(W1curve_perr) <- methods_str

W1curve_avg <- W1ns*0
W1curve_true <- array(0, dim = c(length(seeds), length(W1ns)))
W1curve_est <- array(0, dim = c(length(seeds), length(W1ns)))
for(i in 1:5){
  assign(paste0("W1curve_",methods_str[i]),W1curve_est)
}

# CIs_est <- array(0, dim = c(length(seeds), length(W1ns_CI),2))
# for(i in 1:5){
#   assign(paste0("CIs_",methods_str[i]),CIs_est)
# }

for(i_seed in 1:length(seeds)){
  # cat(i_seed, " ")
  X1s <- get(paste0("X1s","_seed",i_seed))
  X1 <- X1s[length(X1s)]
  for(i in 1:5){
    theta = get(paste0("theta_",methods_str[i],"_seed",i_seed))
    alpha = get(paste0("alpha_",methods_str[i],"_seed",i_seed))
    # pW1_M <- sapply(W1ns, function(w) exp(lpW1_comb(w,X1,N,M, theta, alpha))) # NOT USED?
    
    EW1tmp <- EW1(X1,N,M, theta, alpha); EW1curve <- EW1(X1,N,W1ns, theta, alpha)
    
    # CIs_est <- get(paste0("CIs_",methods_str[i]))
    # for(j in 1:length(W1ns_CI)){
    #   Mi <- W1ns_CI[j]+X1
    #   if(j == 1){
    #     lower_bound_temp = 1
    #     Mi_temp = Mi
    #   } else {
    #     lower_bound_temp = max(1,CIs_est[i_seed,j-1,1] - 100)
    #     Mi_temp = CI2s_est[i_seed,j-1] + 100
    #   }
    #   pi <- exp(sapply(lower_bound_temp:Mi_temp, function(x) lpW1_comb(x,X1,N,Mi,theta,alpha)))
    #   # if(pi[1] < pi[2]) warning(paste("increasing pi",i_seed,i,j))
    #   while(sum(pi) < 0.95){
    #     lower_bound_temp = max(1,lower_bound_temp - 100)
    #     Mi_temp = Mi_temp + 100
    #     pi <- exp(sapply(lower_bound_temp:Mi_temp, function(x) lpW1_comb(x,X1,N,Mi,theta,alpha)))
    #   }
    #   pi_sorted = sort(pi, decreasing = T)
    #   threshold = min(pi_sorted[1:min(which(cumsum(pi_sorted) > 0.95))])
    #   index_tmp = which(pi>=threshold)
    #   index = (lower_bound_temp:Mi_temp)[index_tmp]
    #   CIs_est[i_seed,j,] <- c(min(which(index[-1]-index[-length(index)] == 1)),max(index))
    # }
    # assign(paste0("CIs_",methods_str[i]),CIs_est)
    
    # EW1A1 <- EW1_A1(N,M,theta, alpha); EW1B1 <- EW1_B1(X1,N,M, theta, alpha)
    # pA1(N,M,theta, alpha)
    assign(paste0("EW1_",methods_str[i]),EW1tmp); assign(paste0("EW1curve_",methods_str[i]),EW1curve)
    # assign(paste0("EW1A1_",methods_str[i]),EW1A1); assign(paste0("EW1B1_",methods_str[i]),EW1B1)
    
    W1curve_est <- get(paste0("W1curve_",methods_str[i]))
    W1curve_est[i_seed, ] <- EW1curve
    assign(paste0("W1curve_",methods_str[i]),W1curve_est)
  }
  W1 <- get(paste0("W1","_seed",i_seed))
  W1curve2 <- get(paste0("W1curve2","_seed",i_seed))
  
  W1_perr[i_seed,] <- sapply(sapply(paste0("EW1_",methods_str), function(x) get(x)), 
                             function(x) abs((W1-x)^2/W1^2))
  # if(get(paste0("A1","_seed",i_seed))){
  #   W1A1_perr[i_seed,] <- sapply(sapply(paste0("EW1A1_",methods_str), function(x) get(x)), 
  #                                function(x) abs((W1-x)^2/W1^2))
  # } else {
  #   W1B1_perr[i_seed,] <- sapply(sapply(paste0("EW1B1_",methods_str), function(x) get(x)), 
  #                                 function(x) abs((W1-x)^2/W1^2))
  # }
  W1curve_perr[i_seed,] <- sapply(lapply(paste0("EW1curve_",methods_str), function(x) get(x)), 
                                  function(x) sqrt(mean((W1curve2-x)^2/W1curve2^2))) 
  W1curve_avg <- W1curve_avg + W1curve2
  W1curve_true[i_seed,] <- W1curve2
}
W1curve_avg <- W1curve_avg / length(seeds)

## compute CI for the needed seeds


median_seeds_W1curve <- numeric(5)
median_seeds_Kcurve <- numeric(5)
for(i_method in 1:5){
  median_seeds_W1curve[i_method] = which.min(abs(W1curve_perr[,i_method] - median(W1curve_perr[,i_method])))
  if(i_method == 4){
    median_seeds_W1curve[i_method] = sort(abs(W1curve_perr[,i_method] - median(W1curve_perr[,i_method])), decreasing = F, index.return=T)$ix[2]
  }
  median_seeds_Kcurve[i_method] = which.min(abs(Ks_curve_perr[,i_method] - median(Ks_curve_perr[,i_method])))
}

KcurveCIs_est <- array(0, dim = c(5, length(Kns_post_CI),2))
W1curveCIs_est <- array(0, dim = c(5, length(W1ns_CI),2))

for(i_method in 1:5){
  i_seed = median_seeds_W1curve[i_method]
  theta = get(paste0("theta_",methods_str[i_method],"_seed",i_seed))
  alpha = get(paste0("alpha_",methods_str[i_method],"_seed",i_seed))
  X1s <- get(paste0("X1s","_seed",i_seed))
  X1 <- X1s[length(X1s)]
  
  for(j in 1:length(W1ns_CI)){
    Mi <- W1ns_CI[j]+X1
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
    # pi_sorted = sort(pi, decreasing = T)
    # threshold = min(pi_sorted[1:min(which(cumsum(pi_sorted) > 0.95))])
    # index_tmp = which(pi>=threshold)
    # index = (lower_bound_temp:Mi_temp)[index_tmp]
    # W1curveCIs_est[i_method,j,] <- c(min(which(index[-1]-index[-length(index)] == 1)),max(index))
    W1curveCIs_est[i_method,j,] <- get_CI_general(pi, (lower_bound_temp:Mi_temp))
  }
}

save(list = c("Ks_perr","W1_perr","Kns_post","W1ns","Kns_post_CI","W1ns_CI","N","M","alphas","thetas",
              "W1curve_perr",paste0("W1curve_",methods_str),#paste0("CIs_",methods_str), 
              "W1curve_true", "W1curve_avg",
              "KcurveCIs_est","W1curveCIs_est",
              "median_seeds_W1curve","median_seeds_Kcurve",
              "Ks_curve_perr",paste0("Kcurve_",methods_str),#paste0("KcurveCIs_",methods_str),
              "Kcurve_true","Kcurve_avg"), 
     file = "results/crossval_bibtex_CI.Rdata") # or without CI?



for(i_method in 1:5){
  i_seed = median_seeds_Kcurve[i_method]
  theta = get(paste0("theta_",methods_str[i_method],"_seed",i_seed))
  alpha = get(paste0("alpha_",methods_str[i_method],"_seed",i_seed))
  
  i_m = 1
  KcurveCIs_est[i_method,i_m,] <- get_CI_Kcurve(1)
  if(i_method == 1){ # need a different function for DP
    # for(i_m in 2:length(Kns_post_CI)){
    #   lpi <- Kn_DP_distribution(Kns_post_CI[i_m],theta)
    #   KcurveCIs_est[i_method,i_m,] <- get_CI_general(exp(lpi),1:length(lpi))
    # }
  } else {
    for(i_m in 2:length(Kns_post_CI)){
      lpi <- Kn_distribution(Kns_post_CI[i_m],theta, alpha)
      KcurveCIs_est[i_method,i_m,] <- get_CI_general(exp(lpi),1:length(lpi))
    }
  }
}


save(list = c("Ks_perr","W1_perr","Kns_post","W1ns","Kns_post_CI","W1ns_CI","N","M","alphas","thetas",
              "W1curve_perr",paste0("W1curve_",methods_str),#paste0("CIs_",methods_str), 
              "W1curve_true", "W1curve_avg",
              "KcurveCIs_est","W1curveCIs_est",
              "median_seeds_W1curve","median_seeds_Kcurve",
              "Ks_curve_perr",paste0("Kcurve_",methods_str),#paste0("KcurveCIs_",methods_str),
              "Kcurve_true","Kcurve_avg"), 
     file = "results/crossval_bibtex.Rdata")

tmp = load("results/crossval_bibtex.Rdata")
# using quantiles
png("figures/median_KcurveQ_bibtex.png", width = 10, height = 2.5, units = "in", res = 300)
par(mfrow = c(1,5), oma = c(0, 0, 2, 0), mar= c(3.1,4.1,2.1,2.1))
for(i_method in c(3,2,1,4,5)){
  i_seed = which.min(abs(Ks_curve_perr[,i_method] - median(Ks_curve_perr[,i_method])))
  Kcurve_est <- get(paste0("Kcurve_",methods_str[i_method]))
  Kcurve_mean <- Kcurve_est[i_seed,]
  # sds <- apply(Kcurve_est, MARGIN = 2, sd)
  lCI <- apply(Kcurve_est, MARGIN = 2, quantile, probs = 0.025)
  uCI <- apply(Kcurve_est, MARGIN = 2, quantile, probs = 0.975)
  plot(Kns_post, Kcurve_true[i_seed,], type ="l", ylim = range(c(Kcurve_est,Kcurve_true[i_seed,])),
       main = methods_str2[i_method], xlab = "", ylab = "W1")
  mtext("samples", side = 1, line = 2, cex = 0.8,outer = FALSE)
  polygon(c(rev(Kns_post), Kns_post),c(rev(lCI), uCI), col = rgb(1,0,0,0.3), border = NA)
  lines(Kns_post, Kcurve_mean, lwd = 2, col = "red")
  lines(Kns_post, Kcurve_true[i_seed,], lwd = 2)
}
mtext("Prediction of K as a function of the number of samples", side = 3, line = -0, cex = 1.2,outer = TRUE)
dev.off()


png("figures/median_W1curveQ_bibtex.png", width = 10, height = 2.5, units = "in", res = 300)
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

####

tmp = load("results/crossval_bibtex.Rdata")
Ks_perr_bibtex <- Ks_perr
W1_perr_bibtex <- W1_perr

library(dplyr)
library(ggplot2)
tmp_Ks_bibtex = pivot_longer(as.data.frame(Ks_perr_bibtex), cols = colnames(Ks_perr_bibtex), names_to = "method")
tmp_Ks_bibtex$variable = "K"
tmp_Ks_bibtex$data = "bibtex"
tmp_W1s_bibtex = pivot_longer(as.data.frame(W1_perr_bibtex), cols = colnames(W1_perr_bibtex), names_to = "method")
tmp_W1s_bibtex$variable = "W1"
tmp_W1s_bibtex$data = "bibtex"

results <- rbind(tmp_Ks_bibtex,tmp_W1s_bibtex)
results$method[results$method == "lsX1"] <- "lsM1"
results$method <- factor(results$method, levels = c("stdPYP","ordPYP","ordDP","lsM1","lsK"))

p = ggplot(results) + geom_boxplot(aes(x = method, y = value)) +
  facet_nested( ~ variable, scale = "free") + scale_y_log10() + ylab("Percentage error")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Predictive performance on real data")
p
ggsave(filename = "bibtex_crossval.png", plot = p, device = "png", path = "figures",
       width = 5, height = 3)


