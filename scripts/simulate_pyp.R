library(VGAM)

sims <- (batch-1)*sim_x_batch + 1:sim_x_batch

results <- data.frame(batch = rep(batch, Nsim_post*sim_x_batch),
                      sim = rep(sims, each = Nsim_post), 
                      true_theta = NA, true_alpha = NA, 
                      # true_par = NA,
                      theta_ord = NA, alpha_ord = NA, 
                      theta_std = NA, alpha_std = NA, 
                      theta_ordDP = NA, alpha_ordDP = NA, 
                      theta_lsK = NA, alpha_lsK = NA, 
                      theta_lsX1 = NA, alpha_lsX1 = NA, 
                      K = NA, X1 = NA, K2 = NA, A1 = NA, W1 = NA)
W1_curves <- list()
W1_curve_avgs <- list()
res_fs <- list()
W1curve_fs <- list()


index <- 1
for(sim in sims){
  set.seed(1609+sim)
  cat("\n** ",sim, " **")
  
  theta <- thetas[sim]
  alpha <- alphas[sim]
  if(distr_str == "PYP"){
    tmp_pre <- PYP(theta, alpha, n)
    x_pre <- tmp_pre
  }
  pre_clusters <- data.frame(unique_clusters = unique(tmp_pre))
  pre_clusters$ms <- sapply(pre_clusters$unique_clusters, function(x) sum(tmp_pre == x))
  if(which_order == "arrival_weight"){
    pre_clusters$weights <- rexp(nrow(pre_clusters), rate = 1/(1:nrow(pre_clusters)))
    pre_clusters$rank <- match(pre_clusters$weights, sort(pre_clusters$weights, decreasing = T))
  }
  if(which_order == "astable"){
    pre_clusters <- get_rank_alphastable(tmp_pre, pre_clusters)
  }
  if(which_order == "inv"){
    pre_clusters <- get_rank_test1(tmp_pre, pre_clusters)
  }
  if(which_order == "invsqrt"){
    pre_clusters <- get_rank_test2(tmp_pre, pre_clusters)
  }
  
  x_pre <- pre_clusters$rank[match(x_pre, pre_clusters$unique_clusters)]
  pre_clusters$ms <- sapply(pre_clusters$unique_clusters, function(x) sum(tmp_pre == x))
  ms <- pre_clusters$ms[order(pre_clusters$rank)]
  
  X1s <- sapply(Ns, function(ni) {xx = x_pre[1:ni]; sum(xx == max(xx))})
  X1 <- X1s[length(X1s)]
  Ks <- sapply(Ns, function(ni) {length(unique(x_pre[1:ni]))})
  K <- Ks[length(Ks)]
  
  tmp_nm <- error_safe(optim(par = c(1,0.5), fn = fun_ordEPPF_PY, ms = ms, method = "Nelder-Mead"))
  parameters_ord <- c(tmp_nm$par[1],tmp_nm$par[2])
  
  tmp_nm2 <- error_safe(optim(par = c(1,0.5), fn = fun_EPPF_PY, ms = ms, method = "Nelder-Mead"))
  parameters_std <- c(tmp_nm2$par[1],tmp_nm2$par[2])
  
  tmp_nm3 <- error_safe(optimize(f = fun_ordEPPF_DP, ms = ms, interval = c(0,1000)))
  parameters_ordDP <- c(tmp_nm3$minimum,0)
  
  tmp_nm4 <- optim(par = c(1,0.5), fn = fun_Kls, Ks = Ks, ns = Ns, method = "Nelder-Mead")
  parameters_lsK <- c(tmp_nm4$par[1],tmp_nm4$par[2])
  
  tmp_nm5 <- optim(par = c(1,0.5), fn = fun_X1LS, X1s = X1s, Ns = Ns, method = "Nelder-Mead")
  parameters_lsX1 <- c(tmp_nm5$par[1],tmp_nm5$par[2])
  
  results[(index-1)*Nsim_post + 1:Nsim_post, c("true_theta","true_alpha")] <- rep(c(theta,alpha),each = Nsim_post)
  results[(index-1)*Nsim_post + 1:Nsim_post, c("theta_ord", "alpha_ord")] <- rep(parameters_ord, each = Nsim_post)
  results[(index-1)*Nsim_post + 1:Nsim_post, c("theta_std", "alpha_std")] <- rep(parameters_std, each = Nsim_post)
  results[(index-1)*Nsim_post + 1:Nsim_post, c("theta_ordDP", "alpha_ordDP")] <- rep(parameters_ordDP, each = Nsim_post)
  results[(index-1)*Nsim_post + 1:Nsim_post, c("theta_lsK", "alpha_lsK")] <- rep(parameters_lsK, each = Nsim_post)
  results[(index-1)*Nsim_post + 1:Nsim_post, c("theta_lsX1", "alpha_lsX1")] <- rep(parameters_lsX1, each = Nsim_post)
  results[(index-1)*Nsim_post + 1:Nsim_post, c("K","X1")] <- rep(c(K, X1), each = Nsim_post)
  
  
  for(i_method in 1:5){
    str <- method_str[i_method]
    parameters <- get(paste0("parameters",str))
    EK_f <- EK_PY(n, parameters[1], parameters[2])
    EK_nm_f <- EK_PY(n+m, parameters[1], parameters[2])
    Eunseen_f <- Eunseen_PY(n,m,K, parameters[1], parameters[2])
    
    pA1_f <- pA1(n,m,parameters[1], parameters[2])
    EW1_f <- EW1(X1, n, m, parameters[1], parameters[2])
    EW1_A1_f <- EW1_A1(n,m,parameters[1], parameters[2])
    EW1_B1_f <- EW1_B1(X1,n,m,parameters[1], parameters[2])
    W1curve_f <- sapply(Ms, function(mi) EW1(X1, n,mi,parameters[1], parameters[2]))
    res_f <- c(EK_f, EW1_f, pA1_f, EW1_A1_f, EW1_B1_f,EK_nm_f,Eunseen_f) # at the end so they don't mess up the rest (for now)
    res_fs[[paste0(sim,str)]] <- res_f
    W1curve_fs[[paste0(sim,str)]] <- W1curve_f
  }
  W1_curve_avg_tmp <- 0*Ms
  for(sim_post in 1:Nsim_post){
    cat(sim_post, " ")
    if(distr_str == "PYP"){
      tmp_post <- PYP_post(tmp_pre, theta, alpha, m)
    }
    tmp <- c(tmp_pre, tmp_post)
    clusters <- data.frame(unique_clusters = unique(tmp))
    clusters$ms <- sapply(clusters$unique_clusters, function(x) sum(tmp == x))
    
    if(which_order == "arrival_weight"){
      clusters$weights <- NA
      clusters$weights[match(pre_clusters$unique_clusters,clusters$unique_clusters)] <- pre_clusters$weights
      nw <- nrow(clusters)-nrow(pre_clusters)
      clusters$weights[is.na(clusters$weights)] <- rexp(nw, rate = 1/(1:nw))
      clusters$rank <- match(clusters$weights, sort(clusters$weights, decreasing = T))
    }
    if(which_order == "astable"){
      clusters <- get_rank_post_alphastable(tmp, clusters, pre_clusters)
    }
    if(which_order == "inv"){
      clusters <- get_rank_post_test1(tmp, clusters, pre_clusters)
    }
    if(which_order == "invsqrt"){
      clusters <- get_rank_post_test2(tmp, clusters, pre_clusters)
    }
    x_pre2 <- clusters$rank[match(tmp[1:n], clusters$unique_clusters)]
    x_post <- clusters$rank[match(tmp[n+(1:m)], clusters$unique_clusters)]
    
    
    
    ## find W1 (empirical and formulas)
    max_pre <- max(x_pre2)
    A1 = ifelse(max(x_post)>max_pre, TRUE, FALSE)
    W1 = ifelse(A1, sum(x_post == max(x_post)), sum(c(x_pre2,x_post) == max_pre))
    W1curve <- sapply(Ms, function(x) ifelse(max(x_post[1:x])>max_pre,
                                             sum(x_post[1:x] == max(x_post[1:x])),
                                             sum(c(x_pre2,x_post[1:x]) == max_pre)))
    
    results[(index-1)*Nsim_post + sim_post,"K2"] <- length(clusters$unique_clusters)
    results[(index-1)*Nsim_post + sim_post,"A1"] <- A1
    results[(index-1)*Nsim_post + sim_post,"W1"] <- W1
    
    W1_curves[[paste0(sim,"_",sim_post)]] <- W1curve
    W1_curve_avg_tmp <- W1_curve_avg_tmp + W1curve
  }
  W1_curve_avgs[[sim]] <- W1_curve_avg_tmp/Nsim_post
  index <- index + 1
}

assign(paste0("W1_curves_",batch),W1_curves)
assign(paste0("W1_curve_avgs_",batch),W1_curve_avgs)
assign(paste0("res_fs_",batch),res_fs)
assign(paste0("W1curve_fs_",batch),W1curve_fs)



save_list <- c("results",paste0(c("W1_curves_","W1_curve_avgs_","res_fs_","W1curve_fs_"),batch),
               "Ns","Ms","distr_str","which_order")
filestr <- paste0("../results/",ifelse(distr_str %in% c("PYP"),"simpyp/","simzipf/"))
filestr <- paste0(filestr, "sim",distr_str)
filestr <- paste0(filestr, "_",which_order,"order")

filestr <- paste0(filestr, "_n",n,"_m",m,"_batch",batch,".rdata")
save(list= save_list, 
     file = filestr)




