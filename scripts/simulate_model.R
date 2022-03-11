
sims <- (batch-1)*sim_x_batch + 1:sim_x_batch

results <- data.frame(batch = rep(batch, Nsim_post*sim_x_batch),
                      sim = rep(sims, each = Nsim_post), 
                      sim_post = rep(1:Nsim_post, times = sim_x_batch), 
                      true_theta = NA, true_alpha = NA, 
                      theta_ord = NA, alpha_ord = NA, 
                      theta_std = NA, alpha_std = NA, 
                      theta_ordDP = NA, alpha_ordDP = NA, 
                      theta_lsK = NA, alpha_lsK = NA, 
                      theta_lsX1 = NA, alpha_lsX1 = NA, 
                      K = NA, X1 = NA, A1 = NA, W1 = NA)
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
  tmp_pre = error_safe(generate_ordPY_Cpp(n, theta, alpha))
  
  ms <- tmp_pre$ms
  X1s <- sapply(Ns, function(ni) {xx = tmp_pre$x[1:ni]; sum(xx == max(xx))})
  X1 <- X1s[length(X1s)]
  Ks <- sapply(Ns, function(ni) {length(unique(tmp_pre$x[1:ni]))})
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
    if(i_method == 3){
      EK_f <- parameters[1] * log(n)
    } else {
      EK_f <- EK_PY(n, parameters[1], parameters[2])
    }

    pA1_f <- pA1(n,m,parameters[1], parameters[2])
    EW1_f <- EW1(X1, n, m, parameters[1], parameters[2])
    EW1_A1_f <- EW1_A1(n,m,parameters[1], parameters[2])
    EW1_B1_f <- EW1_B1(X1,n,m,parameters[1], parameters[2])
    W1curve_f <- sapply(Ms, function(mi) EW1(X1, n,mi,parameters[1], parameters[2]))
    res_f <- c(EK_f, EW1_f, pA1_f, EW1_A1_f, EW1_B1_f)
    res_fs[[paste0(sim,str)]] <- res_f
    W1curve_fs[[paste0(sim,str)]] <- W1curve_f
  }
  W1_curve_avg_tmp <- 0*Ms
  for(sim_post in 1:Nsim_post){
    cat(sim_post, " ")
    tmp_post = generate_ordPY_post_Cpp(tmp_pre, m, theta, alpha)
    x_pre <- tmp_pre$x
    x_pre2 <- tmp_post$x[1:n]
    x_post <- tmp_post$x[n+(1:m)]

    ## find W1 (empirical and formulas)
    max_pre <- max(x_pre2)
    A1 = ifelse(max(x_post)>max_pre, TRUE, FALSE)
    W1 = ifelse(A1, sum(x_post == max(x_post)), sum(c(x_pre2,x_post) == max_pre))
    W1curve <- sapply(Ms, function(x) ifelse(max(x_post[1:x])>max_pre,
                                          sum(x_post[1:x] == max(x_post[1:x])),
                                          sum(c(x_pre2,x_post[1:x]) == max_pre)))
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

folder_str = "simmodel/"
save_list <- c("results",paste0(c("W1_curves_","res_fs_","W1curve_fs_"),batch))
save(list= save_list, 
     file = paste0("../results/",folder_str,"simmodel_n500_m5000_batch",batch,".rdata"))




