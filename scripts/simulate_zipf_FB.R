library(VGAM)

sims <- (batch-1)*sim_x_batch + 1:sim_x_batch

EK_post_list <- list()
pA1_post_list <- list()
EW1_post_list <- list()
EW1_A1_post_list <- list()
EW1_B1_post_list <- list()
W1curve_f_post_list <- list()

W1curve_fs_post <- list()
res_fs_post <- list()

thetas_post_list <- list()
alphas_post_list <- list()

times_opt <- c()
times_post <- c()

for(sim in sims){
  set.seed(1609+sim)
  cat("\n** ",sim, " **")
  
  par <- params[sim]
  if(distr_str == "zipf"){
    tmp_pre <- rzeta(n = n, shape = par)
    x_pre <- tmp_pre
  } 
  if(distr_str == "finzipf"){
    tmp_pre <- rzipf(n = n, N = Ktot, shape = par)
    x_pre <- tmp_pre
  }
  pre_clusters <- data.frame(unique_clusters = unique(tmp_pre))
  pre_clusters$ms <- sapply(pre_clusters$unique_clusters, function(x) sum(tmp_pre == x))
  if(which_order == "int"){
    pre_clusters$rank <- match(pre_clusters$unique_clusters, sort(pre_clusters$unique_clusters, decreasing = T))
  }
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
  
  time2 <- Sys.time()
  post <- posterior_MH_hyper(ms, niter = mh_niter)
  time_post <- Sys.time()-time2
  times_post <- c(times_post, time_post)

  thetas_post <- post$thetas
  alphas_post <- post$alphas
  thetas_thin <- thetas_post[post_subset]
  alphas_thin <- alphas_post[post_subset]

  EK_post <- sapply(1:length(post_subset), function(i) EK_PY(n, thetas_thin[i], alphas_thin[i]))
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
  EK_nm_f <- mean(sapply(1:length(post_subset), function(i) EK_PY(n+m, thetas_thin[i], alphas_thin[i])))
  Eunseen_f <- mean(sapply(1:length(post_subset), function(i) Eunseen_PY(n,m,K, thetas_thin[i], alphas_thin[i]))) 
  res_f <- c(EK_f, EW1_f, pA1_f, EW1_A1_f, EW1_B1_f,EK_nm_f,Eunseen_f)

  res_fs_post[[paste0(sim)]] <- res_f
  W1curve_fs_post[[paste0(sim)]] <- colMeans(W1curve_f_post)

  thetas_post_list[[paste0(sim)]] <- thetas_post
  alphas_post_list[[paste0(sim)]] <- alphas_post
  
  EK_post_list[[paste0(sim)]] <- EK_post
  pA1_post_list[[paste0(sim)]] <- pA1_post
  EW1_post_list[[paste0(sim)]] <- EW1_post
  EW1_A1_post_list[[paste0(sim)]] <- EW1_A1_post
  EW1_B1_post_list[[paste0(sim)]] <- EW1_B1_post
  W1curve_f_post_list[[paste0(sim)]] <- W1curve_f_post
}

assign(paste0("res_fs_post_",batch),res_fs_post)
assign(paste0("W1curve_fs_post_",batch),W1curve_fs_post)
assign(paste0("thetas_post_list_",batch),thetas_post_list)
assign(paste0("alphas_post_list_",batch),alphas_post_list)
assign(paste0("EK_post_list_",batch),EK_post_list)
assign(paste0("pA1_post_list_",batch),pA1_post_list)
assign(paste0("EW1_post_list_",batch),EW1_post_list)
assign(paste0("EW1_A1_post_list_",batch),EW1_A1_post_list)
assign(paste0("EW1_B1_post_list_",batch),EW1_B1_post_list)
assign(paste0("W1curve_f_post_list_",batch),W1curve_f_post_list)

save_list <- c("distr_str","which_order",
              paste0(c("res_fs_post_","W1curve_fs_post_",
                        "thetas_post_list_","alphas_post_list_",
                        "EK_post_list_","pA1_post_list_","EW1_post_list_",
                        "EW1_A1_post_list_","EW1_B1_post_list_","W1curve_f_post_list_"),batch),
              "times_post","times_opt")
filestr <- paste0("../results/",ifelse(distr_str %in% c("zipf","finzipf"),"simzipf/","simdir/"))
filestr <- paste0(filestr, "sim",distr_str)
filestr <- paste0(filestr, "_",which_order,"order")

filestr <- paste0(filestr, "_n",n,"_m",m,"_FB_batch",batch,".rdata")
save(list= save_list, 
     file = filestr)
