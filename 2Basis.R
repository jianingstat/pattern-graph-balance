######## library(fda); library(prodlim)
######## generate basis function for all group pairs
Generate_basis_list <- function(data=NULL, pattern=NULL, arg=NULL){
  data_std <- data
  vars_cont <- as.numeric(which(!sapply(data, is.integer)))
  for (k in 1:length(vars_cont))
    data_std[,vars_cont[k]] <- standard(data[,vars_cont[k]])
  basis_list <- NULL; rough_list <- NULL
  for(j in 2:arg$M){
    # j and ancestors of pattern j
    # not just pattern j and Parents of j : which(pattern %in% parents[[j]] | pattern==j)
    # nor pattern1j <- which(pattern==1 | pattern==j)
    varj <- which(!sapply(data[pattern==j, ], anyNA))
    pattern_ACj <- which(rowSums(is.na(data[, varj]))==0)
    ### pattern_ANj <- sort(unique(pattern[which(rowSums(is.na(data[, varj]))==0)]))
    ### pattern_ACj <- which(pattern %in% pattern_ANj)
    tmp <- Generate_basis(data=data_std[pattern_ACj, ], arg=arg)
    basis_list[[j]] <- matrix(NA, nrow=nrow(data), ncol=ncol(tmp$basis))
    basis_list[[j]][pattern_ACj, ] <- tmp$basis
    rough_list[[j]] <- tmp$rough
  }
  return(list(basis_list=basis_list, rough_list=rough_list))
}

######## generate basis functions and their roughness for a pair of patterns: 1_d and r
Generate_basis <- function(data=NULL, arg=NULL){
  ### set tolerance to categorical baisis functions (indicator functions)
  ### if (is.null(arg$tolerance_cate_ratio)) tolerance_cate_ratio <- 1 else tolerance_cate_ratio <- arg$tolerance_cate_ratio
  
  # Only use the commonly observed variables that are of interest
  # Separate continuous and categorical variables
  ### vars <- intersect(arg$vars, which(!sapply(data, anyNA)))
  vars <- as.numeric(which(!sapply(data, anyNA)))
  vars_cont <- intersect(vars, as.numeric(which(!sapply(data, is.integer))))
  vars_cate <- setdiff(vars, vars_cont)
  nvars_cont <- length(vars_cont)
  
  # Generate basis for continuous variables
  if (nvars_cont==0){# if there is no continuous variable, simply create constant function.
    basis_cont <- matrix(1, nrow=nrow(data), ncol=1)
    rough_cont <- 0
  }else{
    # Generate basis for each continuous variable; record roughness (w.r.t each degree) for further calculation.
    basis_cont <- NULL; rough_cont <- NULL
    rough_012 <- matrix(rep(list(NULL), 3*nvars_cont), nrow=3, ncol=nvars_cont)
    for (k in 1:nvars_cont){
      ### dta <- standard(data[, vars_cont[k]])
      dta <- data[, vars_cont[k]]
      if (arg$nbasis==arg$norder) 
        basis_obj <- create.monomial.basis(rangeval=c(0,1), nbasis=arg$nbasis)
      else
        basis_obj <- create.bspline.basis(rangeval=c(0,1), nbasis=arg$nbasis, norder=arg$norder)
      ### rangeval=range(dta)
      basis_cont[[k]] <- eval.basis(evalarg=dta, basisobj=basis_obj)
      for (degree in 0:2)
        rough_012[[(degree+1), k]] <- eval.penalty(basisobj=basis_obj, Lfdobj=degree)
    }
    # orthogonalize the basis functions
    basis_orth <- NULL
    if (arg$cont=="additive" | nvars_cont==1){    # Assemble basis functions (Additive basis functions)
      for (k in 1:nvars_cont){
        rough_orth_k <- rep(0, arg$nbasis)
        coef_orth_k <- matrix(0, nrow=arg$nbasis, ncol=arg$nbasis)
        # find the non-degenerate basis (not linear, nor constant)
        for (i in 1:arg$nbasis)
          if (rough_012[[3, k]][i, i]!=0){
            rough_orth_k[i] <- rough_012[[3, k]][i, i]
            coef_orth_k[i, i] <- 1
            i0 <- i
            break}
        # orthogonalize the rest
        if (i0 < arg$nbasis){
          for (i in (i0+1):arg$nbasis){
            ei <- rep(0, arg$nbasis); ei[i] <- 1; coef_orth_k[i, i] <- 1
            for (j in 1:(i-1))
              if (rough_orth_k[j]!=0){
                proj <- -as.numeric(ei %*% rough_012[[3, k]] %*% coef_orth_k[j, ])/rough_orth_k[j]
                basis_cont[[k]][, i] <- basis_cont[[k]][, i] + proj*basis_cont[[k]][, j] 
                coef_orth_k[i, ] <- coef_orth_k[i, ] + proj*coef_orth_k[j, ]
              }
            rough_orth_k[i] <- as.numeric(coef_orth_k[i, ] %*% rough_012[[3, k]] %*% coef_orth_k[i, ])
          }
        }
        # due to computation error, degenerate functions may has non-zero "rough_orth_k" now. But it's very small.
        basis_orth <- cbind(basis_orth, basis_cont[[k]][, which(rough_orth_k > 1e-5)]) 
        rough_cont <- c(rough_cont, rough_orth_k[which(rough_orth_k > 1e-5)]) 
      }
      # add the degenerate case (constant & linear function)
      for (k in nvars_cont:1){
        basis_orth <- cbind(standard(data[, vars_cont[k]]), basis_orth)
        rough_cont <- c(0, rough_cont) 
      }
      basis_orth <- cbind(rep(1, nrow(basis_orth)), basis_orth)
      rough_cont <- c(0, rough_cont)
    }
    # Tensor product basis functions
    if (arg$cont=="tensor_product" & nvars_cont>1){
      basis_tensor <- matrix(NA, nrow=nrow(data), ncol=arg$nbasis^nvars_cont)
      rough_tensor <- matrix(0, nrow=arg$nbasis^nvars_cont, ncol=arg$nbasis^nvars_cont)
      basis_tensor[, 1:arg$nbasis] <- basis_cont[[1]]
      for (k in 2:nvars_cont)
        for (i in 1:nrow(data))
          basis_tensor[i, 1:(arg$nbasis^k)] <- kronecker(basis_cont[[k]][i, ], basis_tensor[i, 1:arg$nbasis^(k-1)])
      for (k in 1:nvars_cont)
        for (l in 1:nvars_cont){# sum over all possible combinations of v1,v2,...,vp with constraint v1+v2+...+vp=2
          degree <- rep(0, nvars_cont)
          degree[k] <- degree[k] + 1
          degree[l] <- degree[l] + 1
          rough_kl <- 1
          factorials <- 1
          for (i in 1:nvars_cont){# Consider the partial derivative for each variable, use kronecker product to generate the matrix
            # Degree is degree[i](start from 0), but index of list start from 1
            rough_kl <- kronecker(rough_012[[degree[i]+1, i]], rough_kl)
            factorials <- factorials*factorial(degree[i])
          }
          rough_tensor <- rough_tensor + rough_kl*2/factorials
          rough_tensor <- (rough_tensor + t(rough_tensor))/2   # It should be symmetric matrix, but may have computation error. Correct them by hand
        }
      # reordering the basis functions
      basis_tensor <- cbind(basis_tensor[, which(diag(rough_tensor)==0)], basis_tensor[, which(diag(rough_tensor)!=0)])
      rough_tensor <- rough_tensor[c(which(diag(rough_tensor)==0), which(diag(rough_tensor)!=0)), c(which(diag(rough_tensor)==0), which(diag(rough_tensor)!=0))]
      # find the non-degenerate basis (not linear, nor constant)
      rough_orth_tensor <- rep(0, ncol(basis_tensor))
      coef_orth_tensor <- matrix(0, nrow=ncol(basis_tensor), ncol=ncol(basis_tensor))
      for (i in 1:ncol(basis_tensor))
        if (rough_tensor[i, i]!=0){
          rough_orth_tensor[i] <- rough_tensor[i, i]
          coef_orth_tensor[i, i] <- 1
          i0 <- i
          break}
      # orthogonalize the rest
      if (i0 < ncol(basis_tensor)){
        for (i in (i0+1):ncol(basis_tensor)){
          ei <- rep(0, ncol(basis_tensor)); ei[i] <- 1; coef_orth_tensor[i, i] <- 1
          for (j in 1:(i-1)){
            if (rough_orth_tensor[j]!=0){
              proj <- -as.numeric(ei %*% rough_tensor %*% coef_orth_tensor[j, ])/rough_orth_tensor[j]
              basis_tensor[, i] <- basis_tensor[, i] + proj*basis_tensor[, j] 
              coef_orth_tensor[i, ] <- coef_orth_tensor[i, ] + proj*coef_orth_tensor[j, ]
            }
          }
          rough_orth_tensor[i] <- as.numeric(coef_orth_tensor[i, ] %*% rough_tensor %*% coef_orth_tensor[i, ])
        }
      }
      # due to computation error, degenerate functions may has non-zero "rough_orth_tensor" now.
      rough_orth_tensor[which(rough_orth_tensor < 1e-5)] <- 0
      basis_orth <- basis_tensor
      rough_cont <- rough_orth_tensor
    }
    basis_cont <- basis_orth
  }
  ### For stablization, choose the roughness for linear/constant function to be the smallest tolerance or 0.1 (which is smaller).
  # For stablization, choose the roughness penalty for linear/constant function to be the smallest tolerance.
  # Consider categorical variables
  if (length(vars_cate)==0){
    ### rough_cont[which(rough_cont==0)] <- min(rough_cont[which(rough_cont!=0)], 0.1) # for stablization
    rough_cont[which(rough_cont==0)] <- min(rough_cont[which(rough_cont!=0)]) # for stablization
    return(list(basis=basis_cont, rough=rough_cont))
  }
  # set the relative tolerance of categorical varaibles 
  if (arg$cate=="additive"){  # Additive basis functions for each categorical variable separately
    rough_cont[which(rough_cont==0)] <- min(rough_cont[which(rough_cont!=0)]) # for stablization
    data_cate <- as.matrix(data[, vars_cate])
    basis <- cbind(basis_cont, data_cate)
    rough <- c(rough_cont, rep(min(rough_cont), ncol(data_cate)))
    ### rough <- c(rough_cont, rep(tolerance_cate_ratio*min(rough_cont), ncol(data_cate)))
    return(list(basis=basis, rough=rough))
  }
  if (arg$cate=="combination"){# Additive basis functions for each combinations of categorical variables
    ### rough_cont[which(rough_cont==0)] <- min(rough_cont[which(rough_cont!=0)], 0.1) # for stablization
    rough_cont[which(rough_cont==0)] <- min(rough_cont[which(rough_cont!=0)]) # for stablization
    data_cate <- data.frame(data[, vars_cate])
    uni_cate <- unique(data_cate)
    n_cate <- nrow(uni_cate)
    basis_cate <- matrix(0, nrow=nrow(data), ncol=(n_cate-1))# intercept is for individuals where all categorical variables are in their respective reference levels. # take last one as reference.
    for(i in 1:nrow(basis_cont)){
      label <- row.match(data_cate[i, ], uni_cate)
      if (label != 1)# the first combination is set as the baseline.
        basis_cate[i, label] <- 1
    }
    basis <- cbind(basis_cont, basis_cate)
    rough <- c(rough_cont, rep(min(rough_cont), n_cate-1))
    ### rough <- c(rough_cont, rep(tolerance_cate_ratio*min(rough_cont), n_cate-1))
    return(list(basis=basis, rough=rough))
  }
  if (arg$cate=="tensor_product"){  # Basis functions for each combination of categorical variables 
    # Label each row with the unique dummy variable
    data_cate <- data.frame(data[, vars_cate])
    uni_cate <- unique(data_cate)
    n_cate <- nrow(uni_cate)
    basis <- matrix(0, nrow=nrow(data), ncol=n_cate*ncol(basis_cont))
    for(i in 1:nrow(basis)){
      label <- row.match(data_cate[i, ], uni_cate)
      left <- (label-1)*ncol(basis_cont) + 1
      right <- label*ncol(basis_cont)
      basis[i, c(left:right)] <- basis_cont[i, ]
    }
    rough <- rep(rough_cont, n_cate)
    rough[which(rough==0)] <- min(rough_cont[which(rough_cont!=0)]) # for stablization
    return(list(basis=basis, rough=rough))
  }
}