######## library(dplyr); library(mvtnorm)
Generate_graph <- function(M){
  set.seed(111)
  t <- 5
  # Create a data frame from all combinations of the supplied vectors or factors
  misspat <- data.frame(expand.grid(replicate(t, c(1,0), simplify = FALSE)))
  n_posspat <- nrow(misspat)
  misspat <- misspat %>% arrange(desc(rowSums(misspat)))
  misspat <- as.matrix(misspat)
  prob_pat <- 4^rowSums(misspat) + rowSums(matrix(rep(2^(1:t), n_posspat), nrow = n_posspat, byrow = TRUE) * (1-misspat))
  prob_pat[which(rowSums(misspat) == 0)] <- 0
  ind_pat <- sample(2:n_posspat, M-1, replace = FALSE, prob = prob_pat[-1])
  ind_pat <- c(1, sort(ind_pat))
  misspat <- misspat[ind_pat, ]
  parents <- replicate(M, 1, simplify = FALSE)
  if (M > 2)
    for (j in 3:M){
      for (k in 2:(j-1))
        if (sum(misspat[k, ] >= misspat[j, ]) == t)
          parents[[j]] <- c(parents[[j]], k)
    }
  for (j in 2:M)
    if (length(parents[[j]]) >= 2)
      parents[[j]] <- tail(parents[[j]], 2)
  return(list(misspat = misspat, parents = parents))
}

Generate_data <- function(reps, N, misspat, parents, missmech){
  set.seed(reps)
  t	<- 5
  X <- cbind(rep(1, t), c(0:(t-1)))
  beta <- c(6, 0.15)
  mu <- X %*% beta
  if (missmech == 1)
    sig2_e <- 1
  if (missmech == 2)
    sig2_e <- 3
  sig2_L <- 1
  sig2_S <- 1
  sig_LS <- 0.3
  Sig_LS <- matrix(c(sig2_L, sig_LS, sig_LS, sig2_S), 2, 2)
  e	<-	matrix(rnorm(N*t, 0, sqrt(sig2_e)), N, t)
  u 	<-	rmvnorm(N, c(0,0), Sig_LS) 
  y	<-	matrix(rep(mu, N), N, t, byrow = TRUE) +t(X %*% t(u)) + e
  colnames(y) <- c("y1", "y2", "y3", "y4", "y5")
  y_std <- apply(y, MARGIN = 2, standard)
  
  f_1 <- 80 * (y_std[,1] - 0.15) * (y_std[,1] - 0.3) * (y_std[,1] - 1.05) + 1
  f_2 <- 80 * (y_std[,2] - 0.2) * (y_std[,2] - 0.35) * (y_std[,2] - 1.1) + 1
  f_3 <- 80 * (y_std[,3] - 0.25) * (y_std[,3] - 0.4) * (y_std[,3] - 1.05) + 1
  f_4 <- -80 * (y_std[,4] + 0.1) * (y_std[,4] - 0.65) * (y_std[,4] - 0.75) + 1
  f_5 <- -80 * (y_std[,5] + 0.05) * (y_std[,5] - 0.7) * (y_std[,5] - 0.8) + 1
  f <- cbind(f_1, f_2, f_3, f_4, f_5)
  
  # # type 1 
  # if (missmech == 1){
  #   f_1 <- 80 * (y_std[,1] - 0.15) * (y_std[,1] - 0.3) * (y_std[,1] - 1.05) + 1
  #   f_2 <- 80 * (y_std[,2] - 0.2) * (y_std[,2] - 0.35) * (y_std[,2] - 1.1) + 1
  #   f_3 <- 80 * (y_std[,3] - 0.25) * (y_std[,3] - 0.4) * (y_std[,3] - 1.05) + 1
  #   f_4 <- -80 * (y_std[,4] + 0.1) * (y_std[,4] - 0.65) * (y_std[,4] - 0.75) + 1
  #   f_5 <- -80 * (y_std[,5] + 0.05) * (y_std[,5] - 0.7) * (y_std[,5] - 0.8) + 1
  #   f <- cbind(f_1, f_2, f_3, f_4, f_5)
  # }
  # 
  # # type 2
  # if (missmech == 2){
  #   f_1 <- 0.5 - 300 * (y_std[,1] - 0.1) * (y_std[,1] - 0.25) * (y_std[,1] - 0.75) * (y_std[,1] - 0.9) 
  #   f_2 <- 0.5 - 300 * (y_std[,2] - 0.15) * (y_std[,2] - 0.25) * (y_std[,2] - 0.7) * (y_std[,2] - 0.85) 
  #   f_3 <- 0.5 - 300 * (y_std[,3] - 0.1) * (y_std[,3] - 0.2) * (y_std[,3] - 0.7) * (y_std[,3] - 0.85) 
  #   f_4 <- 0.5 - 300 * (y_std[,4] - 0.15) * (y_std[,4] - 0.25) * (y_std[,4] - 0.7) * (y_std[,4] - 0.9) 
  #   f_5 <- 0.5 - 300 * (y_std[,5] - 0.2) * (y_std[,5] - 0.3) * (y_std[,5] - 0.75) * (y_std[,5] - 0.9)
  #   f <- cbind(f_1, f_2, f_3, f_4, f_5) 
  # }

  ######################  propensity score ######################
  if (missmech == 1 | missmech == 2){
    # O_j = P(R=r_j | y)/P(R=r_PAj | y)
    O <- matrix(NA, nrow = N, ncol = M)
    for (j in 2:M){
      O[, j] <- rowSums(f * matrix(rep(misspat[j, ], N), nrow = N, byrow = TRUE))
      cat(j, range(O[,j]), range(exp(O[,j])), "\n")
    }
    O <- exp(O)
    # Q_j = P(R=r_j | y)/P(R=1_d | y) 
    #     = P(R=r_j | y)/P(R=r_PAj | y) * P(R=r_PAj | y)/P(R=1_d | y)
    #     = O_j * Q_PAj = O_j * sum_{s\in PAj} Q_s
    Q <- matrix(NA, nrow = N, ncol = M)
    Q[, 1] <- 1
    for (j in 2:M){
      Q_PAj <- rep(0, N)
      for (s in parents[[j]])
        Q_PAj <- Q_PAj + Q[, s]
      Q[, j] <- O[, j] * Q_PAj
      cat(j, range(Q[,j]), "\n")
    }
  }

  ######################  missing data  ######################
  prop <- Q/rowSums(Q)
  pattern <- rep(NA, N)
  miss <- matrix(NA, nrow = N, ncol = t)
  for (i in 1:N){
    pattern[i] <- sample(x = 1:M, size = 1, prob = prop[i, ])
    miss[i, ] <- misspat[pattern[i], ]
  }
  miss[which(miss == 0)] <- NA
  data <- y * miss
  data <- data.frame(data)
  data_full <- data.frame(y)
  return(list(data_full = data_full, data = data, pattern = pattern, Q = Q, O = O))
}
  

if (FALSE){
  
  table(pattern)
  
  for (j in 2:M){
    cat(j, range(Q[pattern==1, j]), "\n")
  }
  
  
  
}

# 
# head(data_full)
# head(data)
# head(pattern)
# sum(pattern == 1)
# table(pattern)
# colSums(is.na(data))
# 
# 
# plot(c(1:5), data[1, ], type = "p", ylim = range(data, na.rm = T))
# for (i in 2:N){
#   for (j in 1:5)
#     if (!is.na(data[i,j]))
#       points(j, data[i,j])
# }
# plot(c(1:5), data_full[1, ], type = "p", ylim = range(data, na.rm = T))
# for (i in 2:N){
#   points(1:5, data_full[i, ])
# }
# 




