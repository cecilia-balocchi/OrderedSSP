Rcpp::sourceCpp("../scripts/funs.cpp")

sims <- (batch-1)*sim_x_batch + 1:sim_x_batch

EK_list <- list()
pA1_list <- list()
EW1_list <- list()
EW1_A1_list <- list()
EW1_B1_list <- list()
W1curve_f_list <- list()
res_fs_post <- list()

thetas_post_list <- list()
alphas_post_list <- list()

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
  
  post <- posterior_MH_hyper(ms, niter = mh_niter)
  thetas_post <- post$thetas
  alphas_post <- post$alphas
  
  # plot(thetas_post, type = "l")
  # plot(alphas_post, type = "l")
  
  thetas_thin <- thetas_post[post_subset]
  alphas_thin <- alphas_post[post_subset]
  
  # hist(thetas_thin)
  # summary(thetas_thin)
  
  # a_theta <- 0.1; b_theta <- 0.1
  # prior_theta <- rgamma(1000, shape = a_theta, rate = b_theta)
  # hist(prior_theta)
  # mean(prior_theta)
  
  EK_post <- EK_PY(n, thetas_thin, alphas_thin)
  pA1_post <- pA1(n,m, thetas_thin, alphas_thin)
  EW1_post <- EW1(X1, n, m, thetas_thin, alphas_thin)
  EW1_A1_post <- EW1_A1(n,m,thetas_thin, alphas_thin)
  EW1_B1_post <- EW1_B1(X1,n,m,thetas_thin, alphas_thin)
  W1curve_f_post <- sapply(Ms, function(mi) EW1(X1, n,mi,thetas_thin, alphas_thin))
  # each row corresponds to a posterior sample, the columns to the different Ms
  
  EK_f = mean(EK_post)
  EW1_f = mean(EW1_post)
  pA1_f = mean(pA1_post)
  EW1_A1_f = mean(EW1_A1_post)
  EW1_B1_f = mean(EW1_B1_post)
  res_f <- c(EK_f, EW1_f, pA1_f, EW1_A1_f, EW1_B1_f)
  res_fs_post[[paste0(sim)]] <- res_f

  thetas_post_list[[paste0(sim)]] <- thetas_post
  alphas_post_list[[paste0(sim)]] <- alphas_post
  
  EK_list[[paste0(sim)]] <- EK_post
  pA1_list[[paste0(sim)]] <- pA1_post
  EW1_list[[paste0(sim)]] <- EW1_post
  EW1_A1_list[[paste0(sim)]] <- EW1_A1_post
  EW1_B1_list[[paste0(sim)]] <- EW1_B1_post
  W1curve_f_list[[paste0(sim)]] <- W1curve_f_post
}


assign(paste0("res_fs_post_",batch),res_fs_post)
assign(paste0("thetas_post_list_",batch),thetas_post_list)
assign(paste0("alphas_post_list_",batch),alphas_post_list)
assign(paste0("EK_list_",batch),EK_list)
assign(paste0("pA1_list_",batch),pA1_list)
assign(paste0("EW1_list_",batch),EW1_list)
assign(paste0("EW1_A1_list_",batch),EW1_A1_list)
assign(paste0("EW1_B1_list_",batch),EW1_B1_list)
assign(paste0("W1curve_f_list_",batch),W1curve_f_list)

folder_str = "simmodel/"
save_list <- c(paste0(c("thetas_post_list_","alphas_post_list_","res_fs_post_",
                        "EK_list_","pA1_list_","EW1_list_",
                        "EW1_A1_list_","EW1_B1_list_","W1curve_f_list_"),batch))
save(list= save_list, 
     file = paste0("../results/",folder_str,"simmodel_FB_n500_m5000_batch",batch,".rdata"))




