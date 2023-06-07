generate_ordPY <- function(n, theta, alpha, print = FALSE){
  x <- numeric(n)
  k <- 0
  ms <- c()
  rs <- c()
  rs1 <- c()
  for(i in 1:n){
    if(i == 1){
      x[1] <- 1
      k <- 1
      ms <- 1
      rs <- 1
      rs1 <- 0
      if(print) cat(i, " ", k, " ", x, " ",ms,  "\n")
    } else {
      # log prob for clusters j = 1..k
      m <- i-1
      lp1 <- log(ms[1]-alpha) + log(ms[1]+1) - log(theta+m) - log(ms[1])  +log(rs[1]/(rs[1]+1))
      if(k>1){
        lp1 <- lp1 + sum(log(rs[2:k]/(rs[2:k]+1)) + log( alpha*(rs1[2:k]+1) + theta*ms[2:k] ) - log( alpha*rs1[2:k] + theta*ms[2:k] ))
        lps <- log(ms[-1]-alpha) + log(alpha*rs1[-1] + theta*(ms[-1]+1)) - log(theta+m) - log(alpha*rs1[-1] + theta*ms[-1]) +log(rs[-1]/(rs[-1]+1))
        lps <- lps + sapply(2:k, function(j) {
          if(j == k) return(0)
          ind <- (j+1):k
          sum(log(rs[ind]/(rs[ind]+1)) + log( alpha*(rs1[ind]+1) + theta*ms[ind] ) - log( alpha*rs1[ind] + theta*ms[ind] ))
        })
        lps <- c(lp1, lps)
      } else {
        lps <- lp1
      }
      
      # log prob for NEW clusters with RANK j = 1..k (and k+1, but we'll add that later)
      ## I separated j=1 from the rest, added the first term of the product over l
      lq1 <- log(alpha + theta*ms[1]) - log(theta+m) - log(ms[1]+1)
      ## now we add the second part of the product for l = 2..k if k>=2.
      if(k>1){
        lq1 <- lq1 + 
          sum(log(rs[2:k]/(rs[2:k]+1)) + log( alpha*(rs1[2:k]+1) + theta*ms[2:k] ) - log( alpha*rs1[2:k] + theta*ms[2:k] ))
        
        lqs <- log(theta + alpha * rs1[-1]) - log(1+rs1[-1]) - log(theta + m)
        lqs <- lqs + sapply(2:k, function(j){
          ind <- j:k
          sum(log(rs[ind]/(rs[ind]+1)) + log( alpha*(rs1[ind]+1) + theta*ms[ind] ) - log( alpha*rs1[ind] + theta*ms[ind] ))
        })
        lqs <- c(lq1, lqs)
      } else {
        lqs <- lq1
      } 
      lqs <- c(lqs, log(theta + alpha*rs[k])-log(1+rs[k])-log(theta+m) ) # add k+1
      # exp(lqs)/sum(exp(lqs))
      # sum(c(exp(lps),exp(lqs)))
      temp <- sample(2*k+1, size = 1, prob = c(exp(lps),exp(lqs)))
      if(temp <= k){
        x[i] <- temp
        ms[temp] <- ms[temp] + 1
        rs[temp:k] <- rs[temp:k]+1
        if(temp < k){
          rs1[(temp+1):k] <- rs1[(temp+1):k] + 1
        }
      } else {
        if(temp <= 2*k){
          # casino, bisogna shiftare tutto
          new_cl <- temp-k
          x[which(x >= new_cl)] <- x[which(x >= new_cl)] + 1
          x[i] <- new_cl
          if(new_cl == 1){
            newms <- c(1, ms)
            newrs <- c(1, rs+1)
            if(k == 1){
              newrs1 <- c(rs1[1:new_cl], rs1[new_cl]+1)
            } else {
              newrs1 <- c(rs1[1:new_cl], rs1[new_cl]+1, rs1[(new_cl+1):k]+1)
            }
          } else if(new_cl == k) {
            newms <- c(ms[1:(new_cl-1)], 1, ms[new_cl:k])
            newrs <- c(rs[1:(new_cl-1)], rs[new_cl-1]+1, rs[new_cl:k]+1)
            newrs1 <- c(rs1[1:new_cl],rs1[new_cl]+1)
          } else {
            newms <- c(ms[1:(new_cl-1)], 1, ms[new_cl:k])
            newrs <- c(rs[1:(new_cl-1)], rs[new_cl-1]+1, rs[new_cl:k]+1)
            newrs1 <- c(rs1[1:new_cl], rs1[new_cl]+1, rs1[(new_cl+1):k]+1)
          }
          ms <- newms
          rs <- newrs
          rs1 <- newrs1
          k <- k + 1
        } else {
          x[i] <- k+1
          ms <- c(ms, 1)
          rs <- c(rs, rs[k]+1)
          rs1 <- c(rs1, rs[k])
          k <- k + 1
        }
      }
      if(print) cat(i, " ", k, " ", x, " ",ms,  "\n")
    }
  }
  return(list(x = x, ms = ms, rs = rs, rs1 = rs1))
}



risfact <- function(theta, x) {
  lgamma(theta+x)-lgamma(theta)
}
ordEPPF_PY <- function(ms, theta, alpha){
  if((alpha < 0) || (alpha >= 1)) return(NA)
  if(theta < 0) return(NA)
  K <- length(ms)
  rs <- cumsum(ms)
  res <- sum(risfact(1-alpha, ms-1)) - risfact(theta+1,sum(ms)-1) 
  res <- res + sum(log(alpha*rs[-K] + theta*ms[-1]))-sum(log(rs[-1]))
  if(is.nan(res)) cat(theta, " ", alpha, "\n")
  res
}

# \frac{\theta}{\alpha}[\frac{(\theta+\alpha)_{(n)}}{(\theta)_{(n)}}-1]
EK_PY <- function(n, theta, alpha){
  if(alpha < 1e-10) alpha <- 1e-8
  # compute the exact expectation of the number of clusters under the PY with a sample size of n
  theta/alpha * ( exp(risfact(theta+alpha,n) - risfact(theta, n))-1 )
}
# (k+\frac{\theta}{\alpha}) \left(\frac{(\theta+n+\alpha)_{(m)}}{(\theta+n)_{(m)}}-1\right).
# k is the number of clusters in the first n. (x)_{m} = \prod_{i=0,..m-1} (x+i)

Eunseen_PY <- function(n, m, K, theta, alpha){
  if(alpha < 1e-10) alpha <- 1e-8
  (K+theta/alpha) * ( exp(lgamma(theta+n+alpha+m)-lgamma(theta+n+alpha) - lgamma(theta+n+m)+lgamma(theta+n)) -1)
}


EPPF_PY <- function(es, theta, alpha){
  if((alpha < 0) || (alpha >= 1)) return(NA)
  if(theta < 0) return(NA)
    # if(theta <= -alpha) return(NA)
  K <- length(es)
  res <- sum(risfact(1-alpha, es-1)) - risfact(theta+1,sum(es)-1)
  res <- res + ifelse(K>1, sum(log( theta + alpha *(1:(K-1)))), 0)
  res
}


error_safe <- function(expr){
  tryCatch(expr,
           error = function(e){
             message("An error occurred:\n", e)
             NA
           })
}

pX1 <- function(x, n, theta, alpha){
  # formula (22) [13-07-21]
  logp <- log((alpha*(n-x)+theta*x)/n)
  # if x = 1, we want log(rising fact) to be = 0, 
  logp <- logp + lgamma(x-alpha)-lgamma(1-alpha) - lfactorial(x)
  logp <- logp + lfactorial(n) - lgamma(theta+n) + lgamma(theta+1)
  logp <- logp + lgamma(theta+n-x)-lgamma(theta+1)- lfactorial(n-x)
  logp
}

pX1_checks <- function(x, n, theta, alpha){
  if(x == n){
    # logp <- log((alpha*(n-x)+theta*x)/n)
    # logp <- logp + sum(log(1-alpha + 1:(x-1)-1)) 
    # logp <- logp - sum(log(theta+1 + 1:(n-1)-1))
    ## the previous was missing a "-log(theta)"
    logp <- ordEPPF_PY(x, theta, alpha)
  } else {
    logp <- pX1(x,n, theta, alpha)
  }
  return(logp)
}

post_sum <- function(n,m, theta, alpha){
  tmp <- lgamma(n+m)+lgamma(1-alpha+m)+lgamma(theta+n)
  tmp <- tmp - lgamma(1+m) - lgamma(n) - lgamma(1-alpha) - lgamma(n+m+theta)
  1-exp(tmp)
}
### correct
# pX1B <- function(x, n, theta, alpha){
#   lpr <- log(alpha * (n-x) + theta*x) - log(n) 
#   lpr <- lpr + lgamma(x - alpha) - lgamma(1 - alpha) 
#   lpr <- lpr + lfactorial(n) - lfactorial(x) - lfactorial(n-x) 
#   lpr <- lpr - ( lgamma(n + theta) - lgamma(1 + theta) )
#   lpr <- lpr + ifelse(n-x>0, lgamma(n-x + theta) - lgamma(1 + theta), -log(theta)) 
#   return(lpr)
# }


EX1 <- function(n, theta, alpha){
  lE <- lgamma(theta + 1 -alpha +n)- lgamma(theta + 2 -alpha) - lgamma(theta + n) + lgamma(theta+1)
  exp(lE)
}

E2X1 <- function(n, theta, alpha){
  lE2 <- log(2*n*(1-alpha)+theta+alpha) + lgamma(theta-alpha+1+n)-lgamma(theta-alpha+3) -lgamma(theta+n)+lgamma(theta+1)
  exp(lE2)
}

VarX1 <- function(n, theta, alpha){
  E2X1(n, theta, alpha)-EX1(n, theta, alpha)^2
}

pW1 <- function(w,n,m, theta, alpha){
  # formula (25) [13-07-21]
  # equal to lpW1_A1
  logp <- log( alpha*(n+m-w) + theta*w ) - log(n+m)
  logp <- logp + lgamma(1-alpha +w-1)-lgamma(1-alpha) - lfactorial(w)
  logp <- logp + lfactorial(m) - lgamma(theta+1+n+m-1) + lgamma(theta+1)
  logp <- logp + lgamma(theta+1+n+m-w-1)-lgamma(theta+1)- lfactorial(m-w)
  logp
}
lpW1_A1 <- function(w,n,m, theta, alpha){
  # formula (25) [13-07-21] # same as pW1
  logp <- log( alpha*(n+m-w) + theta*w ) - log(n+m)
  logp <- logp + lfactorial(m) - lfactorial(w) - lfactorial(m-w)
  logp <- logp + lgamma(1-alpha +w-1)-lgamma(1-alpha) 
  logp <- logp - lgamma(theta+1+n+m-1) #+ lgamma(theta+1) # si cancellano
  logp <- logp + lgamma(theta+1+n+m-w-1) #-lgamma(theta+1)
  logp
}
lpW1_B1 <- function(w,x,n,m, theta, alpha){
  # formula (27) [13-07-21]
  logp <- log( alpha*(n+m-w-x) + theta*(w+x) ) - log(n+m) 
  logp <- logp + log(n) - log( alpha*(n-x) + theta*x )
  logp <- logp + lfactorial(m) - lfactorial(w) - lfactorial(m-w)
  logp <- logp + lgamma(1-alpha +w+x-1) #-lgamma(1-alpha) # si cancellano
  logp <- logp - lgamma(1-alpha +x-1) #+lgamma(1-alpha)
  logp <- logp + lgamma(theta+1+n-1) #- lgamma(theta+1) # si cancellano
  logp <- logp + lgamma(theta+1+n+m-w-x-1) #- lgamma(theta+1)
  logp <- logp - lgamma(theta+1+n+m-1) #+lgamma(theta+1)
  logp <- logp - lgamma(theta+1+n-x-1) #+lgamma(theta+1)
  logp
}
lpW1_B1_bis <- function(w,x,n,m, theta, alpha){
  # formula (28) [13-07-21] # simpler formula, should be equal to (27)
  logp <- log( alpha*(n+m-w-x) + theta*(w+x) ) - log( alpha*(n-x) + theta*x )
  logp <- logp + log(n) - log(n+m) 
  logp <- logp + lfactorial(m) - lfactorial(w) - lfactorial(m-w)
  logp <- logp + lgamma(theta + n - x + m-w)-lgamma(theta + n - x) 
  logp <- logp + lgamma(x-alpha+w) - lgamma(x-alpha) 
  logp <- logp - lgamma(theta+n+m)+lgamma(theta+n) 
  logp
}
lpW1_comb <- function(w,x,n,m, theta, alpha){
  # formula (29) [13-07-21]
  logp <- log( alpha*(n+m-w) + theta*w ) - log( n+m )
  logp <- logp + lgamma(1-alpha +w-1)-lgamma(1-alpha) 
  logp <- logp + lgamma(theta+1+n+m-w-1) #- lgamma(theta+1)
  logp <- logp - lgamma(theta+1+n+m-1) #+ lgamma(theta+1)
  bool1 <- (w <= m)
  bool2 <- ( (w-x >= 0)&(w-x <= m) )
  tmp1 <- tmp2 <- -Inf
  if(bool1){
    tmp1 <- lfactorial(m) - lfactorial(w) - lfactorial(m-w)
  }
  if(bool2){
    tmp2 <- lfactorial(m) - lfactorial(w-x) - lfactorial(m-w+x)
    tmp2 <- tmp2 + log(n)-log( alpha*(n-x) + theta*x )
    tmp2 <- tmp2 + lgamma(theta+1+n-1)-lgamma(theta+1+n-x-1)-lgamma(1-alpha +x-1)+lgamma(1-alpha)
  }
  logp <- logp + ifelse(bool1 & !bool2, tmp1, 0)
  logp <- logp + ifelse(!bool1 & bool2, tmp2, 0)
  if(bool1 & bool2){
    tmp1 <- logp + tmp1
    tmp2 <- logp + tmp2
    logp <- log(exp(tmp1)+exp(tmp2))
  }
  logp
}
EW1_A1 <- function(n,m, theta, alpha){
  # formula (30) [13-07-21]
  logp <- log(m/(n+m)) + log(theta+alpha*n)-log(theta+n+1-alpha)
  logp <- logp + lgamma(theta+n+1-alpha+m) - lgamma(theta+n+1-alpha)
  logp <- logp - lgamma(theta+n+m) + lgamma(theta+n)
  exp(logp)
}
EW1_B1 <- function(x,n,m, theta, alpha){
  # formula (31) [13-07-21]
  logp <- log(n/(n+m)) - log(alpha*(n-x)+theta*x) + log(Const(n,m,theta,alpha,x))
  logp <- logp + lgamma(theta+n-alpha+m) - lgamma(theta+n-alpha)
  logp <- logp - lgamma(theta+n+m) + lgamma(theta+n)
  exp(logp)
}
EW1 <- function(x,n,m, theta, alpha){
  # formula (32) [13-07-21]
  tmp1 <- log(m/(n+m)) + log(theta+n*alpha)-log(theta+n+1-alpha)
  tmp1 <- tmp1 + lgamma(theta+n+1-alpha+m) - lgamma(theta+n+1-alpha)
  tmp1 <- tmp1 - lgamma(theta+n+m) + lgamma(theta+n)
  tmp2 <- log(n/(n+m)) + lgamma(theta+n-alpha+m)-lgamma(theta+n-alpha)
  tmp2 <- tmp2 - lgamma(theta+n+m) + lgamma(theta+n)
  tmp2 <- tmp2 + log(Const(n,m,theta,alpha,x))
  tmp2 <- tmp2 - log(alpha*(n-x)+theta*x)
  exp(tmp1)+exp(tmp2)
}
Const <- function(n, m, theta,alpha, x){
  C <- (x + m*(x-alpha)/(theta+n-alpha)) * (alpha*(n+m-x) + theta*x)
  C <- C + (x-alpha)/(theta+n-alpha)*m*(theta-alpha)*(x+1)
  C <- C + (x-alpha)*(x+1-alpha)/(theta+n-alpha)/(theta+n+1-alpha)*m*(theta-alpha)*(m-1)
  C
}
pB1 <- function(n, m, theta,alpha){
  # formula (34) 
  logp <- log(n/(n+m)) + lgamma(theta+n+1-alpha+m) - lgamma(theta+n+1-alpha)
  logp <- logp - lgamma(theta+n+m) + lgamma(theta+n)
  exp(logp)
}
pA1 <- function(n, m, theta,alpha){
  1-pB1(n, m, theta,alpha)
}

generate_ordPY_post <- function(res, M, theta, alpha, print = FALSE){
  old_x <- res$x
  n <- length(old_x)
  x <- numeric(n+M)
  x[1:n] <- old_x
  k <- length(res$ms)
  ms <- res$ms
  rs <- res$rs
  rs1 <- res$rs1
  for(i in 1:M){
    # log prob for clusters j = 1..k
    n_temp <- n+i-1
    
    lp1 <- log(ms[1]-alpha) + log(ms[1]+1) - log(theta+n_temp) - log(ms[1])  +log(rs[1]/(rs[1]+1))
    if(k>1){
      lp1 <- lp1 + sum(log(rs[2:k]/(rs[2:k]+1)) + log( alpha*(rs1[2:k]+1) + theta*ms[2:k] ) - log( alpha*rs1[2:k] + theta*ms[2:k] ))
      lps <- log(ms[-1]-alpha) + log(alpha*rs1[-1] + theta*(ms[-1]+1)) - log(theta+n_temp) - log(alpha*rs1[-1] + theta*ms[-1]) +log(rs[-1]/(rs[-1]+1))
      lps <- lps + sapply(2:k, function(j) {
        if(j == k) return(0)
        ind <- (j+1):k
        sum(log(rs[ind]/(rs[ind]+1)) + log( alpha*(rs1[ind]+1) + theta*ms[ind] ) - log( alpha*rs1[ind] + theta*ms[ind] ))
      })
      lps <- c(lp1, lps)
    } else {
      lps <- lp1
    }
    # log prob for NEW clusters with RANK j = 1..k (and k+1, but we'll add that later)
    ## I separated j=1 from the rest, added the first term of the product over l
    lq1 <- log(alpha + theta*ms[1]) - log(theta+n_temp) - log(ms[1]+1)
    ## now we add the second part of the product for l = 2..k if k>=2.
    if(k>1){
      lq1 <- lq1 + 
        sum(log(rs[2:k]/(rs[2:k]+1)) + log( alpha*(rs1[2:k]+1) + theta*ms[2:k] ) - log( alpha*rs1[2:k] + theta*ms[2:k] ))
      
      lqs <- log(theta + alpha * rs1[-1]) - log(1+rs1[-1]) - log(theta + n_temp)
      lqs <- lqs + sapply(2:k, function(j){
        ind <- j:k
        sum(log(rs[ind]/(rs[ind]+1)) + log( alpha*(rs1[ind]+1) + theta*ms[ind] ) - log( alpha*rs1[ind] + theta*ms[ind] ))
      })
      lqs <- c(lq1, lqs)
    } else {
      lqs <- lq1
    } 
    lqs <- c(lqs, log(theta + alpha*rs[k])-log(1+rs[k])-log(theta+n_temp) ) # add k+1
    # sum(c(exp(lps),exp(lqs)))
    temp <- sample(2*k+1, size = 1, prob = c(exp(lps),exp(lqs)))
    if(temp <= k){
      x[n+i] <- temp
      ms[temp] <- ms[temp] + 1
      rs[temp:k] <- rs[temp:k]+1
      if(temp < k){
        rs1[(temp+1):k] <- rs1[(temp+1):k] + 1
      }
    } else {
      if(temp <= 2*k){
        # casino, bisogna shiftare tutto
        new_cl <- temp-k
        x[which(x >= new_cl)] <- x[which(x >= new_cl)] + 1
        x[n+i] <- new_cl
        if(new_cl == 1){
          newms <- c(1, ms)
          newrs <- c(1, rs+1)
          if(k == 1){
            newrs1 <- c(rs1[1:new_cl], rs1[new_cl]+1)
          } else {
            newrs1 <- c(rs1[1:new_cl], rs1[new_cl]+1, rs1[(new_cl+1):k]+1)
          }
        } else if(new_cl == k) {
          newms <- c(ms[1:(new_cl-1)], 1, ms[new_cl:k])
          newrs <- c(rs[1:(new_cl-1)], rs[new_cl-1]+1, rs[new_cl:k]+1)
          newrs1 <- c(rs1[1:new_cl],rs1[new_cl]+1)
        } else {
          newms <- c(ms[1:(new_cl-1)], 1, ms[new_cl:k])
          newrs <- c(rs[1:(new_cl-1)], rs[new_cl-1]+1, rs[new_cl:k]+1)
          newrs1 <- c(rs1[1:new_cl], rs1[new_cl]+1, rs1[(new_cl+1):k]+1)
        }
        ms <- newms
        rs <- newrs
        rs1 <- newrs1
        k <- k + 1
      } else {
        x[n+i] <- k+1
        ms <- c(ms, 1)
        rs <- c(rs, rs[k]+1)
        rs1 <- c(rs1, rs[k])
        k <- k + 1
      }
    }
    if(print) cat(i, " ", k, " ", x, " ",ms,  "\n")
  }
  return(list(x = x, ms = ms, rs = rs, rs1 = rs1))
}



PYP <- function(theta, alpha, n){
  # Sample from the distribution on partitions from 
  # the Pitman Yor process by choosing the first, 
  # then the second conditional on the first, and so on.
  ids <- numeric(n)
  cl_sizes <- rep(NA, n)
  pr <- rep(NA, n+1)
  k <- 0
  ids[1] <- 1
  cl_sizes[k+1] <- 1
  k <- k + 1
  pr[1] <- 1-alpha
  pr[k + 1] <- theta+alpha*k
  for(i in 2:n){
    cluster <- sample(x = k+1, size = 1, prob = pr[1:(k+1)])
    if(cluster <= k){
      ids[i] <- cluster
      cl_sizes[cluster] <- cl_sizes[cluster] + 1
      pr[cluster] <- pr[cluster] + 1 
    } else {
      ids[i] <- k + 1
      cl_sizes[k + 1] <- 1
      k <- k + 1
      pr[k] <- 1-alpha
      pr[k + 1] <- theta+alpha*k
    }
  }
  return(ids)
}

PYP_post <- function(res, theta, alpha, m){
  n <- length(res)
  
  ids <- numeric(n+m)
  ids[1:n] <- res
  old_cl <- as.numeric(table(res))
  k <- length(old_cl)
  
  cl_sizes <- rep(NA, n)
  cl_sizes[1:k] <- old_cl-alpha
  
  pr <- rep(NA, n+1)
  pr[1:k] <- old_cl
  pr[k + 1] <- theta+alpha*k
  
  for(i in (1:m)+n){
    cluster <- sample(x = k+1, size = 1, prob = pr[1:(k+1)])
    if(cluster <= k){
      ids[i] <- cluster
      cl_sizes[cluster] <- cl_sizes[cluster] + 1
      pr[cluster] <- pr[cluster] + 1 
    } else {
      ids[i] <- k + 1
      cl_sizes[k + 1] <- 1
      k <- k + 1
      pr[k] <- 1-alpha
      pr[k + 1] <- theta+alpha*k
    }
  }
  return(ids[n+(1:m)])
}



### rank functions


get_rank_alphastable <- function(tmp_pre1, pre_clusters){
  ## tmp_pre1 is the vector with the cluster assignments (length = n obs)
  ## pre_clusters is a data.frame with columns "unique_clusters"
  ## and as many rows as the number of clusters
  pre_clusters$rank <- NA
  pre_clusters$rank[1] <- 1
  c1 <- pre_clusters$unique_clusters[1]
  if(nrow(pre_clusters)>1){
    for(i in 2:nrow(pre_clusters)){
      K <- i-1
      when_i <- min(which(tmp_pre1 == pre_clusters$unique_clusters[i]))
      m1 <- sum(tmp_pre1[1:(when_i-1)] == c1)
      
      pr <- c(1/m1,rep(1, K))
      r = sample(K+1, size = 1, prob = pr)
      if(r == K+1){
        pre_clusters$rank[i] <- r
      } else {
        index <- pre_clusters$rank[1:K] >= r
        pre_clusters$rank[index] <- pre_clusters$rank[index] + 1
        pre_clusters$rank[i] <- r
        if(r == 1){
          c1 = pre_clusters$unique_clusters[i]
        }
      }
    }
  }
  return(pre_clusters)
}

get_rank_post_alphastable <- function(tmp, clusters, pre_clusters){
  clusters$rank <- NA
  clusters$rank[1:nrow(pre_clusters)] <- pre_clusters$rank
  
  c1 = clusters$unique_clusters[which(clusters$rank == 1)]
  if(nrow(clusters)>nrow(pre_clusters)){
    for(i in (nrow(pre_clusters)+1):nrow(clusters) ){
      K <- i-1
      when_i <- min(which(tmp == clusters$unique_clusters[i]))
      m1 <- sum(tmp[1:(when_i-1)] == c1)
      
      pr <- c(1/m1,rep(1, K))
      r = sample(K+1, size = 1, prob = pr)
      if(r == K+1){
        clusters$rank[i] <- r
      } else {
        index <- clusters$rank[1:K] >= r
        clusters$rank[index] <- clusters$rank[index] + 1
        clusters$rank[i] <- r
        if(r == 1){
          c1 = clusters$unique_clusters[i]
        }
      }
    }
  }
  
  return(clusters)
}


get_rank_test1 <- function(tmp_pre1, pre_clusters){
  ## tmp_pre1 is the vector with the cluster assignments (length = n obs)
  ## pre_clusters is a data.frame with columns "unique_clusters" and "ms"
  ## and as many rows as the number of clusters
  pre_clusters$rank <- NA
  pre_clusters$rank[1] <- 1
  # c1 <- pre_clusters$unique_clusters[1]
  for(i in 2:nrow(pre_clusters)){
    K <- i-1
    
    pr <- 1/(1:(K+1))
    pr <- pr/sum(pr)
    
    r = sample(K+1, size = 1, prob = pr)
    if(r == K+1){
      pre_clusters$rank[i] <- r
    } else {
      index <- pre_clusters$rank[1:K] >= r
      pre_clusters$rank[index] <- pre_clusters$rank[index] + 1
      pre_clusters$rank[i] <- r
      if(r == 1){
        c1 = pre_clusters$unique_clusters[i]
      }
    }
  }
  return(pre_clusters)
}
get_rank_post_test1 <- function(tmp, clusters, pre_clusters){
  clusters$rank <- NA
  clusters$rank[1:nrow(pre_clusters)] <- pre_clusters$rank
  
  if(nrow(clusters) > nrow(pre_clusters)){
    for(i in (nrow(pre_clusters)+1):nrow(clusters) ){
      K <- i-1
      
      pr <- 1/(1:(K+1))
      pr <- pr/sum(pr)
      r = sample(K+1, size = 1, prob = pr)
      if(r == K+1){
        clusters$rank[i] <- r
      } else {
        index <- clusters$rank[1:K] >= r
        clusters$rank[index] <- clusters$rank[index] + 1
        clusters$rank[i] <- r
        if(r == 1){
          c1 = clusters$unique_clusters[i]
        }
      }
    }
  }
  
  return(clusters)
}


get_rank_test2 <- function(tmp_pre1, pre_clusters){
  ## tmp_pre1 is the vector with the cluster assignments (length = n obs)
  ## pre_clusters is a data.frame with columns "unique_clusters" and "ms"
  ## and as many rows as the number of clusters
  pre_clusters$rank <- NA
  pre_clusters$rank[1] <- 1
  # c1 <- pre_clusters$unique_clusters[1]
  for(i in 2:nrow(pre_clusters)){
    K <- i-1
    
    pr <- 1/sqrt(1:(K+1))
    pr <- pr/sum(pr)
    
    r = sample(K+1, size = 1, prob = pr)
    if(r == K+1){
      pre_clusters$rank[i] <- r
    } else {
      index <- pre_clusters$rank[1:K] >= r
      pre_clusters$rank[index] <- pre_clusters$rank[index] + 1
      pre_clusters$rank[i] <- r
      if(r == 1){
        c1 = pre_clusters$unique_clusters[i]
      }
    }
  }
  return(pre_clusters)
}
get_rank_post_test2 <- function(tmp, clusters, pre_clusters){
  clusters$rank <- NA
  clusters$rank[1:nrow(pre_clusters)] <- pre_clusters$rank
  
  if(nrow(clusters) > nrow(pre_clusters)){
    for(i in (nrow(pre_clusters)+1):nrow(clusters) ){
      K <- i-1
      
      pr <- 1/sqrt(1:(K+1))
      pr <- pr/sum(pr)
      r = sample(K+1, size = 1, prob = pr)
      if(r == K+1){
        clusters$rank[i] <- r
      } else {
        index <- clusters$rank[1:K] >= r
        clusters$rank[index] <- clusters$rank[index] + 1
        clusters$rank[i] <- r
        if(r == 1){
          c1 = clusters$unique_clusters[i]
        }
      }
    }
  }
  return(clusters)
}

get_unm_mean_dist <- function(n,m, Kn, theta, alpha){
  # the first element in the returned vector is Pr(unm = 0)
  Niter <- 100
  # Km_stars <- sapply(1:Niter, function(x) sample_Kn(m, theta+n, alpha))
  Km_stars <- sample_Kn_formula(Niter, m, theta+n, alpha)
  if(alpha > 0){
    Betas <- rbeta(Niter, theta/alpha + Kn, n/alpha - Kn)
  } else {
    Betas <- rep(theta/(theta+n), Niter)
  }
  upper_bound <- 1+max(Km_stars)
  pr <- numeric(upper_bound)
  # expect_v <- 0
  for(i in 1:Niter){
    pr[1 + 0:Km_stars[i]] <- pr[1 + 0:Km_stars[i]] + 
      dbinom(0:Km_stars[i], size = Km_stars[i], prob =Betas[i])
    # expect_v <- expect_v + Km_stars[i]*Betas[i]
  }
  # expect_v <- expect_v/Niter
  expect_v <- mean(Km_stars*Betas)
  pr <- pr/Niter
  list(mean = expect_v, dist = pr)
}

sample_Kn <- function(n, theta, alpha){
  # the first one is a 1. so we need to compute this for 2:n --> n-1 = 1:(n-1)
  Ki <- 1
  if(n>1){
    for(i in 2:n){
      Ki <- Ki+rbinom(1, size = 1, prob = (theta+Ki*alpha)/(i-1+theta))
    }
  }
  Ki
}

sample_Kn_samples <- function(nsamples, n, theta, alpha){
  # the first one is a 1. so we need to compute this for 2:n --> n-1 = 1:(n-1)
  Ki <- rep(1,nsamples)
  if(n>1){
    for(i in 2:n){
      Ki <- Ki+rbinom(nsamples, size = 1, prob = (theta+Ki*alpha)/(i-1+theta))
    }
  }
  Ki
}

# Rcpp::sourceCpp('scripts/fun_lgfact.cpp')
sample_Kn_formula <- function(n_samples, n, theta, alpha){
  # compute the full distribution
  lpr <- Kn_distribution(n, theta, alpha)
  sample(1:n, size = n_samples, prob = exp(lpr), replace=TRUE)
}

Kn_distribution <- function(n, theta, alpha){
  lpr = sapply(1:n, function(x) {
    lgamma(theta/alpha + x) - lgamma(theta/alpha) - lgamma(theta + n) + lgamma(theta)
  })
  lpr + lgfactorial(n,alpha)
}


Kn_emp_distribution <- function(n, theta, alpha, nsamples = 1000){
  Ks <- sample_Kn_samples(nsamples,n,theta,alpha)
  Ks_t <- table(Ks)
  pr <- numeric(n)
  pr[as.numeric(names(Ks_t))] <- as.numeric(Ks_t)
  pr <- pr / sum(pr)
  log(pr)
}

Kn_DP_distribution <- function(n, theta){
  # \frac{\theta^{k}}{(\theta)_{(n)} |s(n,k)|
  # does not work for large n
  lpr = (1:n)*log(theta) - (lgamma(theta+n)-lgamma(theta)) + log(abs(Stirling1.all(n)))
  lpr
}

get_CI <- function(pi){
  get_CI_general(pi, 1:length(pi) -1)
}

get_CI_Kcurve <- function(pi){
  get_CI_general(pi, 1:length(pi))
}

get_CI_general <- function(pi,support_pi){
  pi_sorted = sort(pi, decreasing = T)
  threshold = min(pi_sorted[1:min(which(cumsum(pi_sorted) > 0.95))])
  index_tmp = which(pi>=threshold)
  index = support_pi[index_tmp]
  if(length(index_tmp)==1){
    CI2 <- index
    CI1 <- index
  } else {
    CI2 <- max(index)
    ## this will break if there are two points (say 1 and 2) over the threshold but
    ## the rest of the interval is after.
    CI1 <- index[ min(which(index[-1]-index[-length(index)] == 1)) ]
  }
  c(CI1, CI2)
}

### Stirling numbers:
### code was taken from the copula package (could not install it as it gave some errors)
### https://github.com/cran/copula/blob/master/R/special-func.R
.nacopEnv <- new.env(parent=emptyenv(), hash=FALSE)
# S1.tab <- list()
assign("S1.tab", list(), envir = .nacopEnv) ## S1.tab[[n]][k] == s(n, k)
assign("S1.full.n", 0  , envir = .nacopEnv)

Stirling1 <- function(n,k)
{
  ## NOTA BENE: There's no "direct" method available here
  stopifnot(length(n) == 1, length(k) == 1)
  if (k < 0 || n < k) stop("'k' must be in 0..n !")
  if(n == 0) return(1)
  if(k == 0) return(0)
  S1 <- function(n,k) {
    if(k == 0 || n < k) return(0)
    if(is.na(S <- St[[n]][k])) {
      ## s(n,k) = s(n-1,k-1) - (n-1) * s(n-1,k) for all n, k >= 0
      St[[n]][k] <<- S <- if(n1 <- n-1L)
        S1(n1, k-1) - n1* S1(n1, k) else 1
    }
    S
  }
  if(compute <- (nt <- length(St <- get("S1.tab", .nacopEnv))) < n) {
    ## extend the "table":
    length(St) <- n
    for(i in (nt+1L):n) St[[i]] <- rep.int(NA_real_, i)
  }
  else compute <- is.na(S <- St[[n]][k])
  if(compute) {
    S <- S1(n,k)
    ## store it back:
    assign("S1.tab", St, envir = .nacopEnv)
  }
  S
}

Stirling1.all <- function(n)
{
  stopifnot(length(n) == 1)
  if(!n) return(numeric(0))
  if(get("S1.full.n", .nacopEnv) < n) {
    assign("S1.full.n", n, envir = .nacopEnv)
    unlist(lapply(seq_len(n), Stirling1, n=n))# which fills "S1.tab"
  }
  else get("S1.tab", .nacopEnv)[[n]]
}



posterior_MH_hyper <- function(ms, niter){
  # ms are the data counts
  rs <- cumsum(ms)
  K <- length(ms)
  fun_ordEPPF_PY <- function(par,ms) -ordEPPF_PY(ms, par[1], par[2]) # par[1] = theta, par[2] = alpha
  # parameters of the proposal (sd of the truncated normal)
  sigma_theta <- 0.3; sigma_alpha <- 0.3
  # parameters of the priors
  a_theta <- 0.1; b_theta <- 0.1
  a_alpha <- 1; b_alpha <- 1
  # starting points:
  tmp_nm <- error_safe(optim(par = c(1,0.5), fn = fun_ordEPPF_PY, ms = ms, method = "Nelder-Mead"))
  theta <- tmp_nm$par[1]
  alpha <- tmp_nm$par[2]
  
  thetas <- numeric(niter+1)
  alphas <- numeric(niter+1)
  acc_theta <- 0
  acc_alpha <- 0
  thetas[1] <- theta
  alphas[1] <- alpha

  for(iter in 1:niter){
    ### theta first
    # the proposal is a truncate (positive) normal with sd = sigma_theta
    theta_star <- rnorm(1, mean = theta, sd = sigma_theta)
    while(theta_star < 0) theta_star <- rnorm(1, mean = theta, sd = sigma_theta)
    # the ratio of proposals 
    loga <- pnorm(theta_star/sigma_theta, log = TRUE)-pnorm(theta/sigma_theta, log = TRUE)
    # the ratio of prior
    loga <- loga + dgamma(theta_star, shape = a_theta, rate = b_theta, log = TRUE) -
                dgamma(theta, shape = a_theta, rate = b_theta, log = TRUE) 
    # the ratio of likelihood
    loga <- loga + sum(log(alpha*rs[-K] + theta_star*ms[-1])) - sum(log(alpha*rs[-K] + theta*ms[-1]))
    loga <- loga - (risfact(theta_star+1,sum(ms)-1) - risfact(theta+1,sum(ms)-1) )
    # accept or reject
    if(log(runif(1,0,1))<loga){
      theta <- theta_star
      acc_theta <- acc_theta + 1
    }
    thetas[1+iter] <- theta
    ### now alpha
    # the proposal is a truncate (positive) normal with sd = sigma_alpha
    alpha_star <- rnorm(1, mean = alpha, sd = sigma_alpha)
    while(alpha_star < 0) alpha_star <- rnorm(1, mean = alpha, sd = sigma_alpha)
    # the ratio of proposals 
    loga <- pnorm(alpha_star/sigma_alpha, log = TRUE)-pnorm(alpha/sigma_alpha, log = TRUE)
    # the ratio of prior
    loga <- loga + dbeta(alpha_star, a_alpha, b_alpha, log = TRUE) -
                dbeta(alpha, a_alpha, b_alpha, log = TRUE) 
    # the ratio of likelihood
    loga <- loga + sum(log(alpha_star*rs[-K] + theta*ms[-1])) - sum(log(alpha*rs[-K] + theta*ms[-1]))
    loga <- loga + (sum(risfact(1-alpha_star, ms-1)) - sum(risfact(1-alpha, ms-1)) )
    # accept or reject
    if(log(runif(1,0,1))<loga){
      alpha <- alpha_star
      acc_alpha <- acc_alpha + 1
    }
    alphas[1+iter] <- alpha
  }
  return(list(thetas = thetas, alphas = alphas, acceptance = c(acc_theta,acc_alpha)/niter))
}
