######## library(caret); library(glmnet); library(nloptr); library(plyr); library(pracma); library(prodlim)
######## standardization function
standard <- function(x){
  a <- range(x, na.rm=T)[1]
  b <- range(x, na.rm=T)[2]
  if (a == b) return(x) else return((x-a)/(b-a))
}

####### negative loglikelihood function of data (given parameter of interest) (to be minimized)
loglik <- function(theta, y, w){
  X <- matrix(c(rep(1, ncol(y)), 0:(ncol(y)-1)), ncol = 2)
  sigma2_L <- theta[1]^2
  sigma2_S <- theta[2]^2
  sigma_LS <- theta[3]
  sigma2_e <- theta[4]^2
  beta <- theta[c(5,6)]
  Sigma <- matrix(c(sigma2_L, sigma_LS, sigma_LS, sigma2_S), 2, 2)
  V <- X %*% Sigma %*% t(X) + sigma2_e * diag(rep(1, 5))
  m <- X %*% beta
  l <- 0
  y <- y[which(w != 0), ]
  w <- w[which(w != 0)]
  l <- mean(log(dmvnorm(x = y, mean = m, sigma = V)) * w)
  return(-l)
  # for (i in 1:nrow(y))
  #   if (w[i] != 0)
  #     l <- l + log(dmvnorm(x = y[i, ], mean = m, sigma = V)) * w[i]
  # return(-l)
}

# loglik_grad <- function(theta, y, w){
#   return(pracma::grad(f = loglik, x0 = theta, y = y, w = w))
# }



####### cross validation for (either tailored or entropy) #########
crossvalidation_miss1 <- function(reps, arg, data, pattern, basis_list, rough_list){
  set.seed(reps)
  if (arg$missmech == 1)
    theta_true <- c(1, 1, 0.3, 1, 6, 0.15)
  if (arg$missmech == 2)
    theta_true <- c(1, 1, 0.3, sqrt(3), 6, 0.15)
  if (is.null(arg$maxit)) maxit <- 1000 else maxit <- arg$maxit
  if (is.null(arg$nfld)) nfld <- 5 else nfld <- arg$nfld
  if (is.null(arg$lams_init)) lams_init <- c(1e3, 1e-8) else lams_init <- arg$lams_init
  if (is.null(arg$nlam)) nlam <- 12 else nlam <- arg$nlam
  flds <- createFolds(y = c(1:arg$N), k = nfld)
  lams <- exp(seq(log(lams_init[1]), log(lams_init[2]), len = nlam))
  loss_cv <- matrix(0, nrow = nlam, ncol = arg$M)
  O <- matrix(0, nrow = arg$N, ncol = arg$M)
  Q <- matrix(0, nrow = arg$N, ncol = arg$M)
  Q[, 1] <- 1 
  for (j in 2:arg$M){
    PAj <- arg$parents[[j]]
    for (ifld in 1:nfld){
      # cross validation : training data split
      ind_valid <- flds[[ifld]]; ind_train <- setdiff(c(1:arg$N), ind_valid)
      pattern_valid <- pattern[ind_valid]; pattern_train <- pattern[ind_train]
      basis_valid <- basis_list[[j]][ind_valid, ]; basis_train <- basis_list[[j]][ind_train, ]
      Q_train <- Q[ind_train, ]; Q_valid <- Q[ind_valid, ]
      N_train <- length(pattern_train); N_valid <- length(pattern_valid)
      
      # problem setup
      lambda <- Parameter(pos=TRUE)
      alpha <- Variable(ncol(basis_list[[j]]))
      if (arg$loss == "sequential"){
        loss1 <- exp(basis_train[pattern_train == 1, ] %*% alpha) * rowSums(Q_train[pattern_train == 1, PAj, drop=FALSE])
        lossj <- -basis_train[pattern_train == j, ] %*% alpha
        loss <- (sum(loss1) + sum(lossj))/N_train
      }
      if (arg$loss == "restricted"){
        lossPAj <- exp(basis_train[pattern_train %in% PAj, ] %*% alpha)
        lossj <- -basis_train[pattern_train == j, ] %*% alpha
        loss <- (sum(lossPAj) + sum(lossj))/N_train
      }
      if (arg$loss == "entropy"){
        # lossPAj <- log(1+exp(basis_train[pattern_train %in% PAj, ] %*% alpha))
        # lossj <- log(1+exp(-basis_train[pattern_train == j, ] %*% alpha))
        lossPAj <- logistic(basis_train[pattern_train %in% PAj, ] %*% alpha)
        lossj <- logistic(-basis_train[pattern_train == j, ] %*% alpha)
        loss <- (sum(lossPAj) + sum(lossj))/N_train
      }
      penalty <- sum(abs(alpha) * rough_list[[j]])
      loss <- loss + lambda * penalty
      obj_train <- Minimize(loss)
      prob_train <- Problem(obj_train)
      # problem solver
      for (ilam in 1:nlam){
        value(lambda) <- lams[ilam]
        result_train <- solve(prob_train, solver="ECOS", max_iters=maxit, feastol=1e-4, reltol_inacc=1e-4) # faster
        if (grepl("error" , result_train$status) | grepl("unbounded" , result_train$status)){
          repeat{if(system2(command = "mkdir", arg = "lockdir",stderr) == 0){break}}
          cat("reps",reps,"j",j,"ilam",ilam,"ifld",ifld,"solver:ECOS:status",result_train$status,"iterations",result_train$num_iters,"\n",file=paste(arg$path,"/",arg$loss,"_failure.txt",sep=""),append=TRUE)
          system2(command = "rmdir", arg = "lockdir")
          result_train <- solve(prob_train, solver="SCS", max_iters=maxit) # slower
        }
        if (grepl("unbounded" , result_train$status)){
          repeat{if(system2(command = "mkdir", arg = "lockdir",stderr) == 0){break}}
          cat("reps",reps,"j",j,"ilam",ilam,"ifld",ifld,"solver:SCS:status",result_train$status,"iterations",result_train$num_iters,"\n",file=paste(arg$path,"/",arg$loss,"_failure.txt",sep=""),append=TRUE)
          system2(command = "rmdir", arg = "lockdir")
          loss_cv[ilam, j] <- Inf
          next
        }
        # validation
        alpha0 <- as.numeric(result_train$getValue(alpha))
        if (arg$loss == "sequential"){
          loss1 <- exp(basis_valid[pattern_valid == 1, ] %*% alpha0) * rowSums(Q_valid[pattern_valid == 1, PAj, drop=FALSE])
          lossj <- -basis_valid[pattern_valid == j, ] %*% alpha0
          loss <- (sum(loss1) + sum(lossj))/N_valid
        }
        if (arg$loss == "restricted"){
          lossPAj <- exp(basis_valid[pattern_valid %in% PAj, ] %*% alpha0)
          lossj <- -basis_valid[pattern_valid == j, ] %*% alpha0
          loss <- (sum(lossPAj) + sum(lossj))/N_valid
        }
        if (arg$loss == "entropy"){
          lossPAj <- log(1+exp(basis_valid[pattern_valid %in% PAj, ] %*% alpha0))
          lossj <- log(1+exp(-basis_valid[pattern_valid == j, ] %*% alpha0))
          loss <- (sum(lossPAj) + sum(lossj))/N_valid
        }
        loss_cv[ilam, j] <- loss + loss_cv[ilam, j]
      }
    }
    
    # find the optimal 
    ilam <- which.min(loss_cv[, j])
    
    # problem setup
    lambda <- Parameter(pos=TRUE)
    alpha <- Variable(ncol(basis_list[[j]]))
    if (arg$loss == "sequential"){
      loss1 <- exp(basis_list[[j]][pattern == 1, ] %*% alpha) * rowSums(Q[pattern == 1, PAj, drop=FALSE])
      lossj <- -basis_list[[j]][pattern == j, ] %*% alpha
      loss <- (sum(loss1) + sum(lossj))/N
    }
    if (arg$loss == "restricted"){
      lossPAj <- exp(basis_list[[j]][pattern %in% PAj, ] %*% alpha)
      lossj <- -basis_list[[j]][pattern == j, ] %*% alpha
      loss <- (sum(lossPAj) + sum(lossj))/N
    }
    if (arg$loss == "entropy"){
      # lossPAj <- log(1+exp(basis_list[[j]][pattern %in% PAj, ] %*% alpha))
      # lossj <- log(1+exp(-basis_list[[j]][pattern == j, ] %*% alpha))
      lossPAj <- logistic(basis_list[[j]][pattern %in% PAj, ] %*% alpha)
      lossj <- logistic(-basis_list[[j]][pattern == j, ] %*% alpha)
      loss <- (sum(lossPAj) + sum(lossj))/N
    }
    penalty <- sum(abs(alpha) * rough_list[[j]])
    loss <- loss + lambda * penalty
    obj <- Minimize(loss)
    prob <- Problem(obj)
    
    # problem solver
    value(lambda) <- lams[ilam]
    result <- solve(prob, solver="ECOS", max_iters=maxit, feastol=1e-4, reltol_inacc=1e-4) # faster
    if (grepl("error" , result$status) | grepl("unbounded" , result$status)){
      repeat{if(system2(command = "mkdir", arg = "lockdir",stderr) == 0){break}}
      cat("reps",reps,"j",j,"ilam",ilam,"ifld",ifld,"solver:ECOS:status",result$status,"iterations",result$num_iters,"\n",file=paste(arg$path,"/",arg$loss,"_failure.txt",sep=""),append=TRUE)
      system2(command = "rmdir", arg = "lockdir")
      result <- solve(prob, solver="SCS", max_iters=maxit) # slower
    }
    if (grepl("unbounded" , result$status)){
      repeat{if(system2(command = "mkdir", arg = "lockdir",stderr) == 0){break}}
      cat("reps",reps,"j",j,"ilam",ilam,"ifld",ifld,"solver:SCS:status",result$status,"iterations",result$num_iters,"\n",file=paste(arg$path,"/",arg$loss,"_failure.txt",sep=""),append=TRUE)
      system2(command = "rmdir", arg = "lockdir")
    }
    alpha0 <- as.numeric(result$getValue(alpha))

    # propensity odds estiamtes
    O[, j] <- exp(basis_list[[j]] %*% alpha0)
    # update Q
    Q[, j] <- rowSums(Q[, PAj, drop = FALSE]) * O[, j]
    if (ilam == 1 | ilam == nlam){
      repeat{if(system2(command = "mkdir", arg = "lockdir",stderr) == 0){break}}
      cat("reps", reps, "j", j, ":", "ilam", ilam, "\n", file = paste(arg$path,"/boundary_lambda_",arg$loss,".txt",sep = ""), append = TRUE)
      system2(command = "rmdir", arg = "lockdir")
    }
  }
  # output
  w <- rep(0, arg$N); w[pattern == 1] <- rowSums(Q)[pattern == 1]
  repeat{if(system2(command = "mkdir", arg = "lockdir",stderr) == 0){break}}
  write.table(O, file = paste(arg$path,"/",reps,"/O_",arg$loss,".txt",sep = ""), row.names = F, col.names = F)
  write.table(Q, file = paste(arg$path,"/",reps,"/Q_",arg$loss,".txt",sep = ""), row.names = F, col.names = F)
  cat(w, file = paste(arg$path,"/",reps,"/w_",arg$loss,".txt",sep = ""), append = TRUE)
  system2(command = "rmdir", arg = "lockdir")
  lm <- nloptr(rep(1, 6), eval_f = loglik, y = data, w = w, opts=list(algorithm="NLOPT_LN_COBYLA", maxeval=1000, xtol_abs=1e-4))
  theta <- as.numeric(lm$solution)
  
  # asymptotic variance
  if (arg$loss == "sequential"){
    sd <- sd_score(N=arg$N, M=arg$M, data=data, basis_list=basis_list, rough_list=rough_list, pattern=pattern, theta=theta, O=O, w=w)
    CI_lower <- theta - qnorm(0.975)*sd; CI_upper <- theta + qnorm(0.975)*sd
    CI_cover <- as.numeric((CI_lower < theta_true) & (theta_true < CI_upper))
    repeat{if(system2(command = "mkdir", arg = "lockdir",stderr) == 0){break}}
    cat(sd, "\n", file=paste(arg$path,"/sd_",arg$loss,".txt",sep=""), append=TRUE)
    cat(sd, "\n", file=paste(arg$path,"/",reps,"/sd_",arg$loss,".txt",sep=""), append=TRUE)
    cat(CI_cover, "\n", file=paste(arg$path,"/CI_cover_",arg$loss,".txt",sep=""), append=TRUE)
    cat(CI_cover, "\n", file=paste(arg$path,"/",reps,"/CI_cover_",arg$loss,".txt",sep=""), append=TRUE)
    system2(command = "rmdir", arg = "lockdir")
  }

  return(theta)
}


sd_score <- function(N, M, data, basis_list, rough_list, pattern, theta, O, w){
  # require package "pracma", "glmnet" (default alpha=1, lasso)
  lams <- exp(seq(log(1), log(1e-5), len = 100))
  # fully-observed pattern
  ind1 <- which(pattern == 1)
  #  psi_theta is derivative of loglikelihood
  psi <- matrix(0, nrow=N, ncol=length(theta))
  for (i in ind1)
    psi[i, ] <- pracma::grad(f = loglik, x0 = theta, y = data[i, ], w = 1)
  # conditional expectation
  u <- NULL
  u[[1]] <- psi
  for (j in 2:M){
    u[[j]] <- matrix(0, nrow=N, ncol=length(theta))
    PAj <- parents[[j]]
    if (length(PAj) == 1){
      indPAj <- which(pattern == PAj)
      for (l in 1:length(theta)){
        fit_l <- cv.glmnet(y=u[[PAj]][indPAj, l], x=basis_list[[j]][indPAj, ], family="gaussian", penalty.factor=rough_list[[j]], intercept = F)
        u[[j]][, l] <- predict(fit_l, newx = basis_list[[j]], s = "lambda.min")
      }
    }else{
      indPAj <- which(pattern %in% PAj)
      ind_yob <- which(colSums(is.na(data[pattern == j, ])) == 0)
      #Xj <- as.matrix(cbind(data[, ind_yob], data[, ind_yob]^2))
      Xj <- as.matrix(cbind(data[, ind_yob]))
      fit_PAj <- cv.glmnet(y = pattern[indPAj], x = Xj[indPAj, ], family = "multinomial", lambda = lams)
      c <- predict(fit_PAj, newx = Xj, s = "lambda.min", type = "response")[, , 1]
      for (k in 1:length(PAj)){
        s_k <- PAj[k]
        ind_sk <- which(pattern == s_k)
        for (l in 1:length(theta)){
          fit_kl <- cv.glmnet(y=u[[s_k]][ind_sk, l], x=basis_list[[j]][ind_sk, ], family="gaussian", penalty.factor=rough_list[[j]], intercept = F)
          u_kl <- predict(fit_kl, newx = basis_list[[j]], s = "lambda.min")
          u[[j]][, l] <- u[[j]][, l] + u_kl * c[, k]
        }
      }
    }
  }
  # Generate paths from j back to s
  paths <- list()
  paths[[1]] <- list(1)
  for (j in 2:M){
    paths[[j]] <- list()
    PAj <- parents[[j]]
    for (k in 1:length(PAj)){
      s_k <- PAj[k]
      for (l in 1:length(paths[[s_k]])){
        pathl <- c(paths[[s_k]][[l]], j)
        paths[[j]] <- append(paths[[j]], list(pathl))
      }
    }
  }
  # Dtheta is the second derivative of loglikelihood. (hessian matrix)
  Dtheta <- hessian(f=loglik, x0=theta, y=data, w=w)
  inv_Dtheta <- solve(Dtheta)
  # Var_psi is the asymptotic variance of estimator (weighted sum of psi) 
  obj <- matrix(0, nrow=N, ncol=length(theta))
  obj[ind1, ] <- obj[ind1, ] + w[ind1]*psi[ind1, ] 
  for (j in 2:M){
    pathj <- paths[[j]]
    for (k in 1:length(pathj)){
      nodes <- rev(pathj[[k]])
      O_prod <- rep(1, N)
      for (l in 1:(length(nodes)-1)){
        node_l <- nodes[l]
        indl <- which(pattern == node_l)
        indPAl <- which(pattern %in% parents[[node_l]])
        obj[indPAl, ] <- obj[indPAl, ] - O[indPAl, node_l] * u[[node_l]][indPAl, ] * O_prod[indPAl]
        obj[indl, ] <- obj[indl, ] + u[[node_l]][indl, ] * O_prod[indl]
        O_prod <- O_prod * O[, node_l] 
      }
    }
  }

  # sqrt(n) * (theta_hat - theta0) -> N(0, invDtheta %*% Var_psi %*% invDtheta)
  Var_psi <- cov(obj)*(N-1)/N #we want E[obj^2],but cov gives var/(n-1)
  Var_psi <- Var_psi + colMeans(obj) %*% t(colMeans(obj))
  Var_theta <- inv_Dtheta %*% Var_psi %*% inv_Dtheta
  sd_theta <- sqrt(diag(Var_theta)/N)
  return(sd_theta = sd_theta)
}









