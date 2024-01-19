rm(list = ls())
library(ggplot2)
library(dplyr)

# this script computes the order distribution for a new cluster
# and creates the corresponding figure in the manuscript

############################# FUNCTIONS #############################
risfact <- function(theta, x) {
  # (theta)_{(x)} (usually x is integer)
  lgamma(theta+x)-lgamma(theta)
}
ordEPPF_PY <- function(ms, theta, alpha){
  if((alpha < 0) || (alpha >= 1)) return(NA)
  if(theta < 0) return(NA)
  K <- length(ms)
  # rs <- cumsum(ms)
  rs <- c(rev(cumsum(rev(ms))),0)
  
  if(all.equal(theta,0)==TRUE){
    res <- sum(risfact(1-alpha, ms-1)) - lfactorial(sum(ms)-1)
    res <- res + sum(log(alpha*rs[-c(1,K+1)] + theta*ms[-K])) # I removed the last term of the product
    res <- res + log(ms[K]) # this is what remains when theta is cancelled out
    res <- res -sum(log(rs[-(K+1)]))
    return(res)
  }
  
  res <- sum(risfact(1-alpha, ms-1)) - risfact(theta,sum(ms)) 
  res <- res + sum(log(alpha*rs[-1] + theta*ms))-sum(log(rs[-(K+1)]))
  if(is.nan(res)) cat(theta, " ", alpha, "\n")
  res
}

# this function creates "all" possible configurations for a varying number of elements
# divided into k clusters. Note that the way this is written will work only for k <= 10
# some of the generated configurations can be repeated (so we will consider only the unique ones later)
rec <- function(ms, k){
  # k up to 10!
  if(sum(ms)<n_max){
    # the first entry of switch is for EXPR, the (i+1)th one tells what to do if k == i.
    switch (k,
            rbind(ms, rec(ms+1,k)),
            rbind(ms, rec(ms+diag(k)[1,],k),rec(ms+diag(k)[2,],k)),
            rbind(ms, rec(ms+diag(k)[1,],k),rec(ms+diag(k)[2,],k),
                  rec(ms+diag(k)[3,],k)),
            rbind(ms, rec(ms+diag(k)[1,],k),rec(ms+diag(k)[2,],k),
                  rec(ms+diag(k)[3,],k),rec(ms+diag(k)[4,],k)),
            rbind(ms, rec(ms+diag(k)[1,],k),rec(ms+diag(k)[2,],k),
                  rec(ms+diag(k)[3,],k),rec(ms+diag(k)[4,],k),
                  rec(ms+diag(k)[5,],k)),
            rbind(ms, rec(ms+diag(k)[1,],k),rec(ms+diag(k)[2,],k),
                  rec(ms+diag(k)[3,],k),rec(ms+diag(k)[4,],k),
                  rec(ms+diag(k)[5,],k),rec(ms+diag(k)[6,],k)),
            rbind(ms, rec(ms+diag(k)[1,],k),rec(ms+diag(k)[2,],k),
                  rec(ms+diag(k)[3,],k),rec(ms+diag(k)[4,],k),
                  rec(ms+diag(k)[5,],k),rec(ms+diag(k)[6,],k),
                  rec(ms+diag(k)[7,],k)),
            rbind(ms, rec(ms+diag(k)[1,],k),rec(ms+diag(k)[2,],k),
                  rec(ms+diag(k)[3,],k),rec(ms+diag(k)[4,],k),
                  rec(ms+diag(k)[5,],k),rec(ms+diag(k)[6,],k),
                  rec(ms+diag(k)[7,],k),rec(ms+diag(k)[8,],k)),
            rbind(ms, rec(ms+diag(k)[1,],k),rec(ms+diag(k)[2,],k),
                  rec(ms+diag(k)[3,],k),rec(ms+diag(k)[4,],k),
                  rec(ms+diag(k)[5,],k),rec(ms+diag(k)[6,],k),
                  rec(ms+diag(k)[7,],k),rec(ms+diag(k)[8,],k),
                  rec(ms+diag(k)[9,],k)),
            rbind(ms, rec(ms+diag(k)[1,],k),rec(ms+diag(k)[2,],k),
                  rec(ms+diag(k)[3,],k),rec(ms+diag(k)[4,],k),
                  rec(ms+diag(k)[5,],k),rec(ms+diag(k)[6,],k),
                  rec(ms+diag(k)[7,],k),rec(ms+diag(k)[8,],k),
                  rec(ms+diag(k)[9,],k),rec(ms+diag(k)[10,],k))
    )
  }
}
# this computes the probability that a new observation goes into a new cluster with order j
qj_new <- function(j, ms, theta, alpha, rs = NULL){
  # j can range from 1 to k+1!
  k = length(ms)
  n = sum(ms)
  if(is.null(rs)){
    rs = c(rev(cumsum(rev(ms))),0)
  }
  if((all.equal(theta,0)==TRUE) && (j == k+1)){
    # logpr = log(alpha)-log(n)+log(rs[1])-log(1+rs[1]) # computed by hand
    logpr = -log(theta+n) -log(1+rs[j]) # first term put later: log(theta+alpha*rs[j])
    ind1 = 1:(j-2) # removed term j-1
    ind2 = 1+ind1
    if(j>2){
      logpr = logpr +
        sum(log(rs[ind1])-log(1+rs[ind1])) +
        sum(log(alpha * rs[ind2] + alpha + theta*ms[ind1])) -
        sum(log(alpha * rs[ind2] + theta*ms[ind1]))
    }
    # this is term j-1
    logpr = logpr +
      sum(log(rs[j-1])-log(1+rs[j-1])) +
      sum(log(alpha * rs[j] + alpha)) #-
    # sum(log(alpha * rs[j])) +
    # log(alpha*rs[j]) # this was the initial term
    # !! removed both terms as they cancel out !!
  } else {
    logpr = log(theta+alpha*rs[j]) - log(theta+n) - log(1+rs[j])
    if(j > 1){
      ind1 = 1:(j-1)
      ind2 = 1+ind1
      logpr = logpr + 
        sum(log(rs[ind1])-log(1+rs[ind1])) + 
        sum(log(alpha * rs[ind2] + alpha + theta*ms[ind1])) - 
        sum(log(alpha * rs[ind2] + theta*ms[ind1]))
    }
  }
  exp(logpr)
}
# this is needed to adjust the probabilities of the configurations so that they sum to 1
multiplier2 <- function(conf){
  n = sum(conf)
  exp(lfactorial(n) - sum(lfactorial(conf)))
}



###################### Create the configurations ######################
# here we remove the duplicates
for(k in 2:10){
  n_max = 11
  configs <- rec(rep(1,k),k)
  # get unique configs
  configs_unique <- configs[1,,drop = FALSE]
  if(nrow(configs) >=2){
    for(j in 2:nrow(configs)){
      new_bool <- TRUE
      for(l in 1:nrow(configs_unique)){
        if(all.equal(configs[j,],configs_unique[l,])==TRUE){
          new_bool <- FALSE
          break
        }
      }
      if(new_bool == TRUE){
        configs_unique <- rbind(configs_unique, configs[j,])
      }
    }
  }
  assign(paste0("config",k),configs_unique)
}
# for k = 1
k= 1
configs_unique <- array(1:10,dim = c(10,1))
assign(paste0("config",k),configs_unique)


###################### Compute the probabilities ######################
for(temp in 1:5){
  theta = switch (temp,
                  0.5,
                  0.5,
                  0.5,
                  0,
                  0.5)
  alpha = switch (temp,
                  0,
                  0.3,
                  0.5,
                  0.5,
                  0.7)
  cat(theta, " ", alpha, "\n")
  
  combined_all <- data.frame(j = NULL, 
                             pr_j = NULL, 
                             i = NULL,
                             pr_i = NULL,
                             k = NULL, 
                             n = NULL,
                             type = NULL,
                             theta = NULL,
                             alpha = NULL)
  for(k in 1:10){
    configs = get(paste0("config",k))
    probs = array(NA, dim=c(nrow(configs),k+1))
    probs_i = numeric(nrow(configs))
    for(i in 1:nrow(configs)){
      ms = configs[i,,drop=FALSE]
      qs = sapply(1:(k+1), FUN = function(j) qj_new(j, ms, theta, alpha))
      probs[i,] <- qs/sum(qs)
      probs_i[i] <- exp(ordEPPF_PY(ms,theta,alpha)) * multiplier2(as.numeric(ms))
    }
    combined <- data.frame(j = rep(1:(k+1),nrow(configs)), # j is the order
                           pr_j = NA, 
                           i = rep(1:nrow(configs), each = k+1), # configuration identifier
                           pr_i = rep(probs_i, each = k+1), # prob of configuration
                           k = rep(k,(k+1)*nrow(configs)), 
                           n = rep(rowSums(configs),each = k+1),
                           type = "individual",
                           theta = theta,
                           alpha = alpha)
    for(i in 1:nrow(probs)){
      combined$pr_j[(i-1)*(k+1) + 1:(k+1)] <- probs[i,] # the i-th row in probs corresponds to the i-th block of k+1 elements in combined
    }
    combined_all <- rbind(combined_all, combined) # combines all k
  }
  
  # this computes the average across configurations for a given k and given n
  average_n <- combined_all %>% group_by(n,k,j) %>% summarize(pr_j = sum(pr_j*pr_i)/sum(pr_i),
                                                              pr_i = sum(pr_i),
                                                              .groups="drop")
  combined_all <- rbind(combined_all,
                        data.frame(
                          j = average_n$j,
                          pr_j = average_n$pr_j,
                          i = 0,
                          pr_i = average_n$pr_i, # this is actually pr(k)
                          k = average_n$k,
                          n = average_n$n,
                          type = "average_n",
                          theta = theta,
                          alpha = alpha)
  )
  
  # this computes the average across k for a given n
  average_k <- combined_all %>% filter(type == "average_n") %>% group_by(n,j) %>% summarize(pr_j = sum(pr_j*pr_i)/sum(pr_i),
                                                                                            pr_i = sum(pr_i), .groups="drop")
  combined_all <- rbind(combined_all, 
                        data.frame(
                          j = average_k$j, 
                          pr_j = average_k$pr_j, 
                          i = 0,
                          pr_i = average_k$pr_i,
                          k = 0, 
                          n = average_k$n,
                          type = "average_k",
                          theta = theta,
                          alpha = alpha)
  )
  
  if(temp == 1){
    combined_all_n10 <- combined_all %>% filter(n == 10)
  } else {
    combined_all_n10 <- rbind(combined_all_n10,
                              combined_all %>% filter(n == 10))
  }
}


combined_all_n10_selected <- combined_all_n10 #%>% filter(theta != 0)
combined_all_n10_selected$alpha_theta = interaction(combined_all_n10_selected$alpha,combined_all_n10_selected$theta,sep = ", ")
combined_all_n10_selected$alpha_theta = factor(combined_all_n10_selected$alpha_theta, 
                                               levels = c("0, 0.5", "0.3, 0.5","0.5, 0.5","0.7, 0.5", "0.5, 0"))

p <- ggplot(subset(combined_all_n10_selected, type == "individual"), aes(x = j, y = pr_j)) + 
  geom_point(aes(col = factor(k), group = i), alpha = 0.05) + 
  # facet_grid( ~ alpha+theta, labeller = label_both)+ 
  facet_grid( ~ alpha_theta, labeller = label_both)+ 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + scale_color_brewer(palette = 9, type = "div") + ylab("Probability")+ xlab("Order of new cluster") +#xlab("Rank of new cluster") +
  labs(color = "Number of\nprevious\nclusters")+
  # ggtitle("Rank distribution for new cluster", subtitle = "After observing n = 10 observations")+
  ggtitle("Order distribution for new cluster", subtitle = "After observing n = 10 observations")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))
p <- p + geom_line(data = subset(combined_all_n10_selected, type == "average_n"), 
                   mapping = aes(x = j, y = pr_j,col = factor(k)), alpha= 1,linewidth=1)
p <- p + geom_line(data = subset(combined_all_n10_selected, type == "average_k"),
                   mapping = aes(x = j, y = pr_j), alpha= 1, col= "black", linetype = "dashed")
p

ggsave(filename = "ordering_distr_n10_order.png",plot = p, device = "png", path = "figures", width = 12, height = 5)
