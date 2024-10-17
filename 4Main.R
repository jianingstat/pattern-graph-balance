
###################### 0.1 Initialization ######################
n_reps <- 1000
N <- 1000; M <- 8; vars <- c(1:5); norder <- 4; nbasis <- 4; cont <- "additive"; cate <- "additive"; missmech <- 1
path1 <- paste(getwd(), "/result", sep = ""); dir.create(path=path1, showWarnings = FALSE)
packagelist <- c("fda","prodlim","caret","nloptr","doParallel","doRNG","parallel","glmnet","pracma","plyr","dplyr","mvtnorm","CVXR"); lapply(packagelist, require, character.only = TRUE)
maxcore <- detectCores()-2; cl <- makeCluster(maxcore, outfile="log.txt"); registerDoParallel(cl)


###################### 0.2 Generate Graph ######################
source("1Data.R"); source("2Basis.R"); source("3Utility-CVXR.R")
set.seed(111)
tmp1 <- Generate_graph(M = M)
misspat <- tmp1$misspat
parents <- tmp1$parents
# equivalently
misspat <- matrix(c(c(1,1,1,1,1), 
                    c(0,1,1,1,1),
                    c(1,0,1,1,1),
                    c(1,1,1,1,0),
                    c(1,1,0,0,1),
                    c(1,0,1,1,0),
                    c(1,1,0,1,0),
                    c(1,1,0,0,0)), nrow = 8, ncol = 5, byrow = TRUE)
parents <- list(c(1), c(1), c(1), c(1), c(1), c(3,4), c(1,4), c(5,7))


foreach (reps = 1:n_reps, .packages = packagelist) %dorng%{
  source("1Data.R", local=TRUE); source("2Basis.R", local=TRUE); source("3Utility-CVXR.R", local=TRUE)
  path2 <- paste(path1, "/", reps, sep=""); dir.create(path = path2, showWarnings = FALSE)

  ###################### 1 Generate Data ######################
  tmp2 <- Generate_data(reps = reps, N = N, misspat = misspat, parents = parents, missmech = missmech)
  data_full <- tmp2$data_full
  data <- tmp2$data
  pattern <- tmp2$pattern
  Q_true <- tmp2$Q
  O_true <- tmp2$O
  
  ###################### 2 Generate Basis ######################
  arg <- list(N=N, M=M, vars=vars, missmech=missmech, parents=parents, norder=norder, nbasis=nbasis, cont=cont, cate=cate)
  tmp3 <- Generate_basis_list(data=data, pattern=pattern, arg=arg)
  basis_list <- tmp3$basis_list; rough_list <- tmp3$rough_list

  ###################### 3.1 Estimate from Full data ######################
  w_full <- rep(1, N)
  lm.1 <- nloptr(rep(1,6), eval_f = loglik, y=data_full, w=w_full, opts=list(algorithm="NLOPT_LN_COBYLA", maxeval=1e4, xtol_abs=1e-4))
  theta_full <- as.numeric(lm.1$solution)
  repeat{if(system2(command="mkdir", arg="lockdir", stderr=NULL)==0){break}}
  cat(theta_full, "\n", file=paste(path1,"/theta_full.txt",sep=""), append=TRUE)
  cat(theta_full, "\n", file=paste(path2,"/theta_full.txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
  
  ###################### 3.2 Estimate from Complete case ######################
  w_cc <- rep(0, N); w_cc[pattern==1] <- 1
  lm.2 <- nloptr(rep(1,6), eval_f = loglik, y=data, w=w_cc, opts=list(algorithm="NLOPT_LN_COBYLA", maxeval=1e4, xtol_abs=1e-4))
  theta_cc <- as.numeric(lm.2$solution)
  repeat{if(system2(command="mkdir", arg="lockdir", stderr=NULL)==0){break}}
  cat(theta_cc, "\n", file=paste(path1,"/theta_cc.txt",sep=""), append=TRUE)
  cat(theta_cc, "\n", file=paste(path2,"/theta_cc.txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
  
  ###################### 3.3 Estimate from IPW with true propensity ######################
  w_true <- rep(0, N)
  w_true[pattern==1] <- rowSums(Q_true[pattern==1, ])
  lm.3 <- nloptr(rep(1,6), eval_f = loglik, y=data, w=w_true, opts=list(algorithm="NLOPT_LN_COBYLA", maxeval=1e4, xtol_abs=1e-4))
  theta_trueweight <- as.numeric(lm.3$solution)
  repeat{if(system2(command="mkdir", arg="lockdir", stderr=NULL)==0){break}}
  cat(theta_trueweight, "\n", file=paste(path1,"/theta_trueweight.txt",sep=""), append=TRUE)
  cat(theta_trueweight, "\n", file=paste(path2,"/theta_trueweight.txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
  
  
  #arg_cv <- append(arg, list(maxit=100, nfld=5, lams_init=c(1e-3, 1e-6), nlam=4, path=path1))
  #arg_cv <- append(arg, list(maxit=1000, nfld=5, lams_init=c(1, 1e-5), nlam=10, path=path1))
  arg_cv <- append(arg, list(maxit=1000, nfld=5, lams_init=c(1, 1e-7), nlam=20, path=path1))
  ###################### 3.4 Estimate from IPW by weights from sequential loss ######################
  arg_cv$loss <- "sequential"
  theta_sequential <- crossvalidation_miss1(reps=reps, arg=arg_cv, data=data, pattern=pattern, basis_list=basis_list, rough_list=rough_list)
  repeat{if(system2(command="mkdir", arg="lockdir", stderr=NULL)==0){break}}
  cat(theta_sequential, "\n", file=paste(path1,"/theta_sequential.txt",sep=""), append=TRUE)
  cat(theta_sequential, "\n", file=paste(path2,"/theta_sequential.txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
  
  ###################### 3.5 Estimate from IPW by weights from restricted loss ######################
  arg_cv$loss <- "restricted"
  theta_restricted <- crossvalidation_miss1(reps=reps, arg=arg_cv, data=data, pattern=pattern, basis_list=basis_list, rough_list=rough_list)
  repeat{if(system2(command="mkdir", arg="lockdir", stderr=NULL)==0){break}}
  cat(theta_restricted, "\n", file=paste(path1,"/theta_restricted.txt",sep=""), append=TRUE)
  cat(theta_restricted, "\n", file=paste(path2,"/theta_restricted.txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
  
  ###################### 3.6 Estimate from IPW by weights from entropy loss ######################
  arg_cv$loss <- "entropy"
  theta_entropy <- crossvalidation_miss1(reps=reps, arg=arg_cv, data=data, pattern=pattern, basis_list=basis_list, rough_list=rough_list)
  repeat{if(system2(command="mkdir", arg="lockdir", stderr=NULL)==0){break}}
  cat(theta_entropy, "\n", file=paste(path1,"/theta_entropy.txt",sep=""), append=TRUE)
  cat(theta_entropy, "\n", file=paste(path2,"/theta_entropy.txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
}



