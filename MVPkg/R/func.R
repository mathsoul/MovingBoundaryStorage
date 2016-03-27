# Using Lagrangian to get the optimal solution to the vanilla problem
# which doesn't have a penalty term
# ||...||^2 + lambda(t(mu)w - rho) + eta(t(1)w - 1) 


# get solution for the optimization problems ------------------------------
getVanillaSolution = function(cov, mu_vec = NULL, rho = NULL, lb = NULL, ub = NULL, sum_w = "one"){
  n_stocks = ncol(cov)
  if(sum_w == "one"){
    A = matrix(1,nrow = 1,ncol = n_stocks)
    rhs = 1
    sense = "="
  }else{
    A = matrix(0, nrow = 1,ncol = n_stocks)
    rhs = 0
    sense = "="
  }
  
  if(!is.null(mu_vec) & !is.null(rho)){
    A = rbind(A,mu_vec)
    rhs = c(rhs,rho)
    sense = c(sense,">=")
  }
  
  obj = rep(0,n_stocks)
  
  model = list()
  model$A = A
  model$rhs = rhs
  model$sense = sense
  model$obj = obj
  model$Q = cov
  
  if(!is.null(lb)){
    model$lb = lb
  }else{
    model$lb = rep(-Inf,n_stocks)
  }
  
  if(!is.null(ub)){
    model$ub = ub
  }else{
    model$ub = rep(Inf,n_stocks)
  }
  
  params = list(BarConvTol = 1e-16,OutputFlag = 0)
  result = gurobi(model,params)
  
  return(result$x)
  
}

# get solution by using inverse covariance matrix 
getMultSol = function(inv_cov,mu_vec,rho){
  n_stocks = ncol(inv_cov)
  mu_vec = matrix(mu_vec, ncol = 1)
  one_vec = matrix(1,nrow = n_stocks,ncol = 1)
  
  #   w = (rho * (inv_cov %*% mu_vec %*% t(one_vec) %*% inv_cov %*% one_vec -
  #                inv_cov %*% one_vec %*% t(mu_vec) %*% inv_cov %*% one_vec) +
  #     inv_cov %*% one_vec %*% t(mu_vec) %*% inv_cov %*% mu_vec - 
  #     inv_cov %*% mu_vec %*% t(one_vec) %*% inv_cov %*% mu_vec)/
  #     as.numeric(t(mu_vec) %*% inv_cov %*% mu_vec %*% t(one_vec) %*% inv_cov %*% one_vec -
  #        t(mu_vec) %*% inv_cov %*% one_vec %*% t(one_vec) %*% inv_cov %*% mu_vec)
  
  # another way
  w =  (inv_cov %*% (mu_vec %*% t(one_vec) - one_vec %*% t(mu_vec)) %*% 
          inv_cov %*% (rho * one_vec - mu_vec))/
    as.numeric(t(mu_vec) %*% inv_cov %*% mu_vec %*% t(one_vec) %*% inv_cov %*% one_vec -
                 t(mu_vec) %*% inv_cov %*% one_vec %*% t(one_vec) %*% inv_cov %*% mu_vec)
  
  return(w)
}

getMatSol = function(cov,mu_vec){
  inv_cov = solve(cov)
  n_stocks = ncol(inv_cov)
  mu_vec = matrix(mu_vec, ncol = 1)
  
  B = cbind(mu_vec,rep(1,n_stocks))
  # B = B[,2,drop = FALSE]
  
  C = solve(t(B) %*% inv_cov %*% B)
  D = inv_cov %*% B
  return(D %*% C)
}


getSolution = function(inv_cov, mu_vec, rho){
  n_stocks = ncol(inv_cov)
  mu_vec = matrix(mu_vec, ncol = 1)
  
  B = cbind(mu_vec,rep(1,n_stocks))
  
  C = solve(t(B) %*% inv_cov %*% B)
  D = inv_cov %*% B
  
  w = D %*% C %*% matrix(c(rho,1), ncol = 1)
  
  return(w)
}



getCLassoSolution = function(log_tau_vec, rho, mu_vec, cov_mat, sum_w = "one"){
  n_stocks = ncol(cov_mat)
  
  I_I = cbind(diag(n_stocks), -diag(n_stocks))
  
  if(sum_w == "one"){
    A = matrix(1,nrow = 1,ncol = n_stocks)
    rhs = 1
    sense = "="
  }else{
    A = matrix(0, nrow = 1,ncol = n_stocks)
    rhs = 0
    sense = "="
  }
  
  if(!is.null(mu_vec) & !is.null(rho)){
    A = rbind(A, mu_vec)
    rhs = c(rhs, rho)
    sense = c(sense, ">=")
  }
  
  A = A %*% I_I
  
  Q = t(I_I) %*% cov_mat %*% I_I
  
  if(length(log_tau_vec) == 1){
    log_tau_vec = rep(log_tau_vec, n_stocks)
  }
  
  obj = c(10^log_tau_vec, 10^log_tau_vec)
  # obj = c(rep(0, n_stocks), 2 * 10^log_tau_vec)
  
  model = list()
  model$A = A
  model$rhs = rhs
  model$sense = sense
  model$obj = obj
  model$Q = Q
  model$lb = rep(0, 2*n_stocks)
  
  params = list(BarConvTol = 1e-16,OutputFlag = 0)
  result = gurobi(model,params)
  
  w = I_I %*% result$x
  return(as.numeric(w))
}

# getCLassoSolutionCone = function(log_tau_vec, rho, mu_vec, cov_mat, sum_w = "one"){
#   n_stocks = ncol(cov_mat)
#   
#   I_I = cbind(diag(n_stocks), matrix(0, n_stocks, n_stocks))
#   
#   if(sum_w == "one"){
#     A = matrix(1,nrow = 1,ncol = n_stocks)
#     rhs = 1
#     sense = "="
#   }else{
#     A = matrix(0, nrow = 1,ncol = n_stocks)
#     rhs = 0
#     sense = "="
#   }
#   
#   if(!is.null(mu_vec) & !is.null(rho)){
#     A = rbind(A, mu_vec)
#     rhs = c(rhs, rho)
#     sense = c(sense, "=")
#   }
#   
#   A = A %*% I_I
#   
#   Q = t(I_I) %*% cov_mat %*% I_I
#   
#   if(length(log_tau_vec) == 1){
#     log_tau_vec = rep(log_tau_vec, n_stocks)
#   }
#   
#   obj = c(rep(0, n_stocks), 10^log_tau_vec)
#   # obj = c(rep(0, n_stocks), 2 * 10^log_tau_vec)
#   
#   cones = list()
#   
#   for(i in 1:n_stocks){
#     cones = c(cones, list(list(n_stocks + i, i)))
#   }
#   
#   model = list()
#   model$A = A
#   model$rhs = rhs
#   model$sense = sense
#   model$obj = obj
#   model$Q = Q
#   model$lb = c(rep(-Inf, n_stocks), rep(0, n_stocks))
#   model$cones = cones
#   
#   params = list(BarConvTol = 1e-16, OutputFlag = 0)
#   result = gurobi(model, params)
# 
#   w = I_I %*% result$x
#   
#   print(c(result$objval, getNonZero(w,Para$tol)))
#   return(as.numeric(w))
# }

getCLassoInEqSolutionCone = function(log_tau_vec, rho, mu_vec, cov_mat, sum_w = "one"){
  n_stocks = ncol(cov_mat)
  
  I_I = cbind(diag(n_stocks), matrix(0, n_stocks, n_stocks))
  
  if(sum_w == "one"){
    A = matrix(1,nrow = 1,ncol = n_stocks)
    rhs = 1
    sense = "="
  }else{
    A = matrix(0, nrow = 1,ncol = n_stocks)
    rhs = 0
    sense = "="
  }
  
  if(!is.null(mu_vec) & !is.null(rho)){
    A = rbind(A, mu_vec)
    rhs = c(rhs, rho)
    sense = c(sense, ">=")
  }
  
  A = A %*% I_I
  
  Q = t(I_I) %*% cov_mat %*% I_I
  
  if(length(log_tau_vec) == 1){
    log_tau_vec = rep(log_tau_vec, n_stocks)
  }
  
  obj = c(rep(0, n_stocks), 10^log_tau_vec)
  # obj = c(rep(0, n_stocks), 2 * 10^log_tau_vec)
  
  cones = list()
  
  for(i in 1:n_stocks){
    cones = c(cones, list(list(n_stocks + i, i)))
  }
  
  model = list()
  model$A = A
  model$rhs = rhs
  model$sense = sense
  model$obj = obj
  model$Q = Q
  model$lb = c(rep(-Inf, n_stocks), rep(0, n_stocks))
  model$cones = cones
  
  params = list(BarConvTol = 1e-16, OutputFlag = 0)
  result = gurobi(model, params)
  
  w = I_I %*% result$x
  penalty = sum(10^log_tau_vec * result$x[(n_stocks + 1):(2*n_stocks)])
  n_non_zero = getNonZero(w, Para$tol)
  
#   eigen_est = eigen(cov_mat)
#   projection = t(w) %*% eigen_est$vectors
#   plot(cumsum(abs(projection)), type = 'l', main = log_tau_vec[1])
  
#   if(train_length < Para$n_stocks){
#     sum1 = sum(abs(projection[1:train_length]))
#     sum2 = sum(abs(projection)) - sum1
#     print(round(c(sum1,sum2),2))
#   }
 
  # print(c(result$objval, penalty, result$objval - penalty, round(n_non_zero)))
  return(as.numeric(w))
}

getDualCLassoSolutionCone = function(log_tau_vec, s2, mu_vec, cov_mat, rotate, sum_w = "one"){
  n_stocks = ncol(cov_mat)
  
  w_idx = 1:n_stocks
  t_idx = w_idx + n_stocks
  z_idx = t_idx + n_stocks
  s_idx = max(z_idx) + 1
  
  Gamma = chol(cov_mat)
  
  n_var = max(s_idx)
  n_lconstr = n_stocks + 2
  
  A = matrix(0, nrow = n_lconstr, ncol = n_var)
  # t(w) %*% 1 = 1
  row_idx = 1
  if(sum_w == "one"){
    A[row_idx, w_idx] = 1
  }

  # s = sqrt(s2)
  row_idx = row_idx + 1
  A[row_idx, s_idx] = 1
  # z = Gamma %*% w
  row_idx = row_idx + 1:n_stocks
  A[row_idx, z_idx] = diag(n_stocks)
  A[row_idx, w_idx] = -Gamma
  
  if(sum_w == "one"){
    rhs = c(1, sqrt(s2), rep(0, n_stocks))
  }else{
    rhs = c(0, sqrt(s2), rep(0, n_stocks))
  }
  
  sense = rep("=", n_lconstr)
  
  if(length(log_tau_vec) == 1){
    log_tau_vec = rep(log_tau_vec, n_stocks)
  }
  
  obj = rep(0, n_var)
  obj[w_idx] = - mu_vec
  obj[t_idx] = 10^log_tau_vec
  
  cones = list(as.list(c(s_idx, z_idx)))
  
  for(i in 1:n_stocks){
    if(rotate == FALSE){
      cones = c(cones, list(list(t_idx[i], w_idx[i])))
    }else{
      cones = c(cones, list(list(t_idx[i], z_idx[i])))
    }
    
  }
  
  lb = rep(-Inf, n_var)
  lb[t_idx] = 0
  
  model = list()
  model$A = A
  model$rhs = rhs
  model$sense = sense
  model$obj = obj
  model$lb = lb
  model$cones = cones
  
  params = list(OutputFlag = 0, BarHomogeneous = 1) #,
  result = gurobi(model, params)
  
  if(result$status != "OPTIMAL"){
    # print("result is unbounded or infeasible")
    w = getMinVarSol(inv_cov = solve(cov_mat))
    return(as.numeric(w))
  }
  
  stopifnot(length(result$x) > 0)
  
  
  w = result$x[w_idx]
  return(as.numeric(w))
}

# getDCLassoWEigenSolution = function(log_tau_vec, s2, mu_vec, cov_mat){
#   n_stocks = ncol(cov_mat)
#   
#   w_idx = 1:n_stocks
#   t_idx = w_idx + n_stocks
#   z_idx = t_idx + n_stocks
#   x_idx = z_idx + n_stocks
#   s_idx = max(x_idx) + 1
#   
#   Gamma = chol(cov_mat)
#   eigen_result = eigen(cov_mat/Para$train_length)
#   eigen_vectors = eigen_result$vectors
#   
#   n_var = max(s_idx)
#   n_lconstr = 2 * n_stocks + 2
#   
#   A = matrix(0, nrow = n_lconstr, ncol = n_var)
#   # t(w) %*% 1 = 1
#   row_idx = 1
#   A[row_idx, w_idx] = 1
#   # s = sqrt(s2)
#   row_idx = row_idx + 1
#   A[row_idx, s_idx] = 1
#   # z = Gamma %*% w
#   row_idx = row_idx + 1:n_stocks
#   A[row_idx, z_idx] = diag(n_stocks)
#   A[row_idx, w_idx] = -Gamma
#   # x = eigen_vectors %*% w
#   row_idx = row_idx + n_stocks
#   A[row_idx, x_idx] = diag(n_stocks)
#   A[row_idx, w_idx] = -eigen_vectors
#   
#   rhs = c(1, sqrt(s2), rep(0, 2 * n_stocks))
#   sense = rep("=", n_lconstr)
#   
#   if(length(log_tau_vec) == 1){
#     log_tau_vec = rep(log_tau_vec, n_stocks)
#   }
#   
#   obj = rep(0, n_var)
#   obj[w_idx] = -mu_vec
#   obj[t_idx] = 10^log_tau_vec
#   
#   cones = list(as.list(c(s_idx, z_idx)))
#   
#   for(i in 1:n_stocks){
#       cones = c(cones, list(list(t_idx[i], x_idx[i])))
#   }
#   
#   lb = rep(-Inf, n_var)
#   lb[t_idx] = 0
#   
#   model = list()
#   model$A = A
#   model$rhs = rhs
#   model$sense = sense
#   model$obj = obj
#   model$lb = lb
#   model$cones = cones
#   
#   params = list(OutputFlag = 0, BarHomogeneous = 1) #,
#   result = gurobi(model, params)
#   
#   stopifnot(result$status != "INF_OR_UNBD")
#   
#   w = result$x[w_idx]
#   return(as.numeric(w))
# }

getDCLassoWEigenSolution = function(log_tau_vec, s2, mu_vec, cov_mat, sum_w = "one"){
  n_stocks = ncol(cov_mat)
  
  w_idx = 1:n_stocks
  t_idx = w_idx + n_stocks
  z_idx = t_idx + n_stocks
  x_idx = z_idx + n_stocks
  s_idx = max(x_idx) + 1
  
  Gamma = chol(cov_mat)
  eigen_result = eigen(cov_mat/Para$train_length)
  
  # Y = diag(sqrt(eigen_result$values)) %*% t(eigen_result$vectors)
  Y = t(eigen_result$vectors)
  
  n_var = max(s_idx)
  n_lconstr = 2 * n_stocks + 2
  
  A = matrix(0, nrow = n_lconstr, ncol = n_var)
  # t(w) %*% 1 = 1
  row_idx = 1  
  if(sum_w == "one"){
    A[row_idx, w_idx] = 1
  }
  
  # s = sqrt(s2)
  row_idx = row_idx + 1
  A[row_idx, s_idx] = 1
  # z = Gamma %*% w
  row_idx = row_idx + 1:n_stocks
  A[row_idx, z_idx] = diag(n_stocks)
  A[row_idx, w_idx] = -Gamma
  # x = t(eigen_vectors) %*% w
  row_idx = row_idx + n_stocks
  A[row_idx, x_idx] = diag(n_stocks)
  A[row_idx, w_idx] = -Y
  
  if(sum_w == "one"){
    rhs = c(1, sqrt(s2), rep(0, 2 * n_stocks))
  }else{
    rhs = c(0, sqrt(s2), rep(0, 2 * n_stocks))
  }
  
  sense = rep("=", n_lconstr)
  
  if(length(log_tau_vec) == 1){
    log_tau_vec = rep(log_tau_vec, n_stocks)
  }
  
  obj = rep(0, n_var)
  obj[w_idx] = -mu_vec
  obj[t_idx] = 10^log_tau_vec
  
  cones = list(as.list(c(s_idx, z_idx)))
  
  for(i in 1:n_stocks){
    cones = c(cones, list(list(t_idx[i], x_idx[i])))
  }
  
  lb = rep(-Inf, n_var)
  lb[t_idx] = 0
  
  model = list()
  model$A = A
  model$rhs = rhs
  model$sense = sense
  model$obj = obj
  model$lb = lb
  model$cones = cones
  
  params = list(OutputFlag = 0, BarHomogeneous = 1) #,
  result = gurobi(model, params)
  
  stopifnot(result$status != "INF_OR_UNBD")
  
  w = result$x[w_idx]
  return(as.numeric(w))
}

getDL2Solution = function(log_tau, s2, mu_vec, cov_mat, sum_w = "one"){
  n_stocks = ncol(cov_mat)
  
  w_idx = 1:n_stocks
  z_idx = w_idx + n_stocks
  t_idx = max(z_idx) + 1
  s_idx = max(t_idx) + 1
  
  Gamma = chol(cov_mat)
  
  n_var = max(s_idx)
  n_lconstr = n_stocks + 2
  
  A = matrix(0, nrow = n_lconstr, ncol = n_var)
  # t(w) %*% 1 = 1
  row_idx = 1  
  if(sum_w == "one"){
    A[row_idx, w_idx] = 1
  }
  
  # s = sqrt(s2)
  row_idx = row_idx + 1
  A[row_idx, s_idx] = 1
  # z = Gamma %*% w
  row_idx = row_idx + 1:n_stocks
  A[row_idx, z_idx] = diag(n_stocks)
  A[row_idx, w_idx] = -Gamma
  
  if(sum_w == "one"){
    rhs = c(1, sqrt(s2), rep(0, n_stocks))
  }else{
    rhs = c(0, sqrt(s2), rep(0, n_stocks))
  }
  
  sense = rep("=", n_lconstr)
  
  obj = rep(0, n_var)
  obj[w_idx] = -mu_vec
  obj[t_idx] = 10^log_tau
  
  # s2 >= wSigmaw
  cones = list(as.list(c(s_idx, z_idx)))
  
  # t^2 >= ||w||_2^2
  cones = c(cones, list(as.list(c(t_idx, w_idx))))
  
  lb = rep(-Inf, n_var)
  lb[t_idx] = 0
  
  model = list()
  model$A = A
  model$rhs = rhs
  model$sense = sense
  model$obj = obj
  model$lb = lb
  model$cones = cones
  
  params = list(OutputFlag = 0, BarHomogeneous = 1) #,
  result = gurobi(model, params)
  
  stopifnot(result$status != "INF_OR_UNBD")
  
  w = result$x[w_idx]
  return(as.numeric(w))
}


getSolutionOfRegs = function(Y, factors, sum_alpha_LB){
  n_stocks = ncol(Y)
  
  dim_factors = dim(factors)
  if(is.null(dim_factors)){
    n_factors = 1
  }else{
    n_factors = dim_factors[2]
  }
  
  n_coef = n_factors + 1 # add intercept
  
  Q = matrix(0, nrow = n_coef * n_stocks, ncol = n_coef * n_stocks)
  obj = rep(NA, n_coef * n_stocks)
  A = matrix(0, nrow = 1, ncol = n_coef * n_stocks)
  rhs = sum_alpha_LB
  sense = ">="
  lb = rep(-Inf, n_coef * n_stocks)
  
  X = cbind(1, factors)
  
  for(i in 1:n_stocks){
    idx = 1:n_coef + (i - 1) * n_coef
    Q[idx,idx] = t(X) %*% X
    obj[idx] = -2 * Y[,i] %*% X
    A[1,idx[1]] = 1
  }
  
  model = list(Q = Q, obj = obj, A = A, rhs = rhs, sense = sense, lb = lb)
  
  params = list(OutputFlag = 0)
  result = gurobi(model, params)
  
  sol = matrix(result$x, nrow = n_coef)
  # print(sol)
  Y_pred = X %*% sol
  
  var_unExp = apply(Y_pred - Y, MARGIN = 2, var)
  
  w = (sol[1,]/var_unExp)/sum(sol[1,]/var_unExp)
  
  return(w)
}

# get some statistics from data -------------------------------------------
# Generate Q = 1/n_time \sum_{t=1}^{n_time} t(r_t)%*%r_t 
getQ = function(ret_mat){
  n_time = nrow(ret_mat)
  Q = t(ret_mat)%*% ret_mat/n_time
  
  return(Q)
}

getRhoVec = function(data,Para){
  rho_vec = rep(NA,Para$n_reform)
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    train_idx = (t-Para$train_length+1):t
    ret = as.matrix(data[train_idx,])
    
    rho_vec[i] = mean(ret,na.rm = TRUE)
  }
  return(rho_vec)
}

# Test weights satisfies constraints or not
testWeights = function(w,mu_vec,rho){
  ret_weighted = as.numeric(w %*% mu_vec)
  if(all.equal(sum(w),1)&all.equal(ret_weighted,rho)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

# calculate inverse covariance matrix -------------------------------------
solveInvCovGL2 = function(cov, log_tau){
  # Our problem is:
  # Maximize log det M – tr C M – rho ||M||^2
  # Where C is observed covariance and M is the desired regularized inverse covariance.
  # Taking derivatives, we want
  # M^{-1} – C – 2*rho*M = 0
  # This corresponds to (pre-multiplying by M)
  # I – M C – 2 * rho * M * I * M = 0
  # It also corresponds to (post-multiplying by M)
  # I – C M – 2 * rho * M * I * M = 0
  # Adding, and dividing by 2,
  # (-0.5 C) M + M (-0.5 C) – M * (2 * rho * I) * M + I = 0
  # This is in the form required for the Thm referenced in the paper.
  # http://scholarworks.gsu.edu/cgi/viewcontent.cgi?article=1045&context=math_theses
  # A^t X + X A − X S X + Q = 0.
  n = ncol(cov)
  A = -0.5 * cov
  S = 2 * 10^log_tau * diag(n)
  Q = diag(n)
  
  M = solveRiccati(A, S, Q)
  return(M)
}

solveInvCovCL2 = function(cov, log_tau){
  # objective is to minimize
  # w^t C w + tau ||w||_2^2
  # with some constraints
  # It is the same as w^t (C + tau I) w
  
  n = ncol(cov)
  
  cov_inv = solve(cov + 10^log_tau * diag(n))

  return(cov_inv)
}

solveInvCovByMultBestConst = function(inv_cov_est, inv_cov_true, use = "off diag", cheat_on_diag = TRUE){
  n_stocks = ncol(inv_cov_est)
  diag_idx = which(as.logical(diag(n_stocks)))
  
  if(use == "off diag"){
    off_diag_est = inv_cov_est[-diag_idx]
    off_diag_true = inv_cov_true[-diag_idx]
    
    fit = summary(lm(off_diag_est ~ 0 + off_diag_true))
  }else if(use == "all"){
    fit = summary(lm(c(inv_cov_est) ~ 0 + c(inv_cov_true)))
  }else{
    stop("use is wrong.")
  }
  
  inv_cov = inv_cov_est
  inv_cov = inv_cov_est/fit$coefficients[1,1]
  
  if(cheat_on_diag == TRUE){
    inv_cov[diag_idx] = inv_cov_true[diag_idx]
  }
  
  return(inv_cov)
}


solveInvCovOneByOne = function(data, Para){
  
}


solveInvCovMIQP = function(mat,K){
  n_stocks = ncol(mat)
  n_time = nrow(mat)
  n = n_stocks - 1
  
  inv_cov = diag(n_stocks)
  R2_vec = rep(NA,n_stocks)
  

  for(j in 1:n_stocks){
    y = mat[,j]
    # add intercept to the X
    X = cbind(rep(1,n_time),mat[,-j])
    # use one regress on all the other. Thus - 1.
    
    result = MIQPWIntercept(y,X,K)  
    intercept = result$intercept
    beta = result$beta
    z = result$z
    
    idx = which(abs(z - 1) < 1e-10)
    stopifnot(length(idx) == K)
    
    beta = beta[idx]
    # because we take j stock out, we need to modify idx 
    idx = idx + (idx >= j)
    
    R2 = result$R2
    exp_var = result$obj
    inv_cov[j,idx] = - beta
    inv_cov[j,] = inv_cov[j,]/exp_var
    
    
    R2_vec[j] = result$R2
  }
  
  return(result = list(inv_cov = inv_cov, R2_vec = R2_vec))
}

# get inverse covariance matrix one line by another
solveInvCovPValue = function(data, Para, p){
  inv_cov = diag(Para$n_stocks)
  for(j in 1:Para$n_stocks){
    y = as.matrix(data[, j])
    X = as.matrix(data[, -j])
    result = VarSelectWPvalue(y, X, p)  
    
    inv_cov[j, -j] = - result$beta
    inv_cov[j, ] = inv_cov[j, ]/ result$var_unExp
  }
  return(inv_cov)
}

solveInvCovStep = function(data, Para, direction){
  inv_cov = diag(Para$n_stocks)
  for(j in 1:Para$n_stocks){
    result = VarSelectStep(data, j, direction)  
    
    inv_cov[j, -j] = - result$beta
    inv_cov[j, ] = inv_cov[j, ]/ result$var_unExp
  }
  return(inv_cov)
}

solveInvCovGlasso = function(data, Para, inv_cov_glasso){
  inv_cov = diag(Para$n_stocks)
  for(j in 1:Para$n_stocks){
    result = VarSelectGlasso(data, inv_cov_glasso, j)  
    
    inv_cov[j, -j] = - result$beta
    inv_cov[j, ] = inv_cov[j, ]/ result$var_unExp
  }
  return(inv_cov)
}

solveInvCovIndiReg = function(data, Para, log_tau, type){
  inv_cov = diag(Para$n_stocks)
  
  for(j in 1:Para$n_stocks){
    y = data[, j]
    X = t(data[, -j])
    
    if(type == "L1"){
      result = l1.reg(X, y, 10^log_tau)
    }else if(type == "L2"){
      result = l2.reg(X, y, 10^log_tau)
    }else{
      stop("type is not correct")
    }
    
    inv_cov[j, -j] = - result$estimate
    inv_cov[j, ] = inv_cov[j, ]/ mean(result$residual^2)
  }
  
  return(inv_cov)
}

solveInvCovAverage = function(inv_cov, avg_ratio){
  n_stocks = ncol(inv_cov)
  inv_cov_avg = avg_ratio * inv_cov + (1 - avg_ratio) * diag(n_stocks)
  return(inv_cov_avg)
}

# covariance and inverse covariance ------------------------------------------
modifyCov = function(cov, sim_cov_type, n_group = 2, group_vec = NULL){
  n_stocks = ncol(cov)
  
  if(sim_cov_type == "original"){
    cov_modified = cov
  }else if(sim_cov_type == "diag"){
    cov_modified = diag(diag(cov))
  }else if(sim_cov_type == "mix"){
    cov_modified = matrix(0, n_stocks, n_stocks)
    
    if(!is.null(n_group)){
      group_size = n_stocks / n_group
      stopifnot(is.wholenumber(group_size))
      start_idx = seq(1, n_stocks, by = group_size)
      end_idx = start_idx + group_size - 1
    
      for(i in seq_along(start_idx)){
        idx = start_idx[i]:end_idx[i]
        cov_modified[idx, idx] = cov[idx, idx]
      }
    }else{
      group_names = unique(group_vec)
      
      for(i in seq_along(group_names)){
        name = group_names[i]
        idx = which(group_vec == name)
        cov_modified[idx, idx] = cov[idx, idx]
      }
    }
  }
  
  return(cov_modified)
}



# TS represents time series

getInvCovTSWPValue = function(p, data, Para, use_type, mu_type, cov_type){
  inv_cov_TS = array(NA, dim = c(Para$n_stocks, Para$n_stocks, Para$n_reform))
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    inv_cov_TS[,,i] = solveInvCovPValue(ret, Para, p)
  }
  
  return(inv_cov_TS)
}

getInvCovTSWUBGlasso = function(inv_cov_glasso_TS, data, Para, use_type, mu_type, cov_type){
  inv_cov_TS = array(NA, dim = c(Para$n_stocks, Para$n_stocks, Para$n_reform))
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
  
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    inv_cov_glasso = inv_cov_glasso_TS[,,i]
    inv_cov_TS[,,i] = solveInvCovGlasso(ret, Para, inv_cov_glasso)
  }
  
  return(inv_cov_TS)
}


getInvCovTSWIndiReg = function(log_tau, reg_type, data, Para, use_type, mu_type, cov_type){
  inv_cov_TS = array(NA, dim = c(Para$n_stocks, Para$n_stocks, Para$n_reform))
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    inv_cov_TS[,,i] = solveInvCovIndiReg(ret, Para, log_tau, reg_type)
  }
  
  return(inv_cov_TS)
}

getInvCovTSAverage = function(avg_ratio, data, Para, use_type, mu_type, cov_type){
  inv_cov_TS = array(NA, dim = c(Para$n_stocks, Para$n_stocks, Para$n_reform))
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    inv_cov = solve(cov_mat)
    inv_cov_TS[,,i] = solveInvCovAverage(inv_cov, avg_ratio)
  }
  
  return(inv_cov_TS)
}

getInvCovTSClip = function(tau, data, Para, use_type, mu_type, cov_type, funcky = FALSE, max_eig_value = 3000){
  inv_cov_TS = array(NA, dim = c(Para$n_stocks, Para$n_stocks, Para$n_reform))
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    eig = eigen(cov_mat)
    eig$values[eig$values <= tau] = tau
    
    if(funcky == TRUE){
      eig$values[1] = max_eig_value
    }
    
    
    inv_cov = eig$vectors %*% diag(1/eig$values) %*% t(eig$vectors)
    inv_cov_TS[,,i] = inv_cov
  }
  
  return(inv_cov_TS)
}


getInvCovTSKMeans = function(n_Kmeans, data, Para, use_type, mu_type, cov_type){
  inv_cov_TS = array(NA, dim = c(Para$n_stocks, Para$n_stocks, Para$n_reform))
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    cluster = ccfkms(t(ret), n_Kmeans)
    inv_cov_TS[,,i] = solve(modifyCov(cov_mat, "mix", n_group = NULL, group_vec = cluster$cl))
  }
  
  return(inv_cov_TS)
}

getCovTS = function(data, Para, use_type, mu_type, cov_type){

  cov_TS = array(NA, dim = c(Para$n_stocks, Para$n_stocks, Para$n_reform))
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    cov_TS[,,i] = cov_mat
  }
  
  return(cov_TS)
}

recoverBeta = function(inv_cov){
  var_unExp = diag(1/diag(inv_cov))
  return(var_unExp %*% inv_cov)
}

getAlpha = function(inv_cov, mu_vec){
  n_stocks = ncol(inv_cov)
  mu_vec = matrix(mu_vec, nrow = n_stocks)
  alpha_vec = recoverBeta(inv_cov) %*% mu_vec
  return(alpha_vec)
}


# condition number, eigen values and Lagrangian multiplier ----------------
getEigenValue = function(Para,cov_TS){
  eigen_values = list()
  eigen_values$train = matrix(NA,nrow = Para$n_reform, ncol = Para$n_stocks)
  eigen_values$prev_all = eigen_values$train
  
  for(i in 1:Para$n_reform){
    eigen_result = eigen(cov_TS$train[,,i], symmetric = TRUE)
    eigen_values$train[i,] = eigen_result$values
    
    eigen_result = eigen(cov_TS$prev_all[,,i], symmetric = TRUE)
    eigen_values$prev_all[i,] = eigen_result$values
  }
  return(eigen_values)
}

getLambdaD2C = function(log_tau_vec, w, mu_vec, cov_mat){
  n_stocks = length(w)
  rhs =  mu_vec - 10^log_tau_vec *sign(w)
  A = cbind(2 * cov_mat %*% w, 1)
  
  idx = sample(1:n_stocks, size = 2, replace = FALSE, prob = abs(w))
  
  while(kappa(A[idx,]) > 1e4){
    idx = sample(1:n_stocks, size = 2, replace = FALSE, prob = abs(w))
  }
  
  lambda = solve(A[idx,], rhs[idx])      
  
  log_tau_vec = log_tau_vec - log(lambda[1],10)
  
  if(min(abs(w)[idx])< 1e-3){
    warning("Only one non zero weight!")
  }
  
  non_zero_idx = which(abs(w) > 1e-1)
  
  error = A %*% lambda - rhs
  
  if(max(abs(error[non_zero_idx]) > 1e-1)){
    warning("Error is larger than 1e-1")
    # plot(abs(w)[non_zero_idx], error[non_zero_idx])
  }
  
  return(log_tau_vec)
}

getLambdaC2D = function(log_tau_vec, w, mu_vec, cov_mat){
  n_stocks = length(w)
  
  rhs = -2 * cov_mat %*% w - 10^log_tau_vec *sign(w)
  A = cbind(-mu_vec, 1)
  
  idx = sample(1:n_stocks, size = 2, replace = FALSE, prob = abs(w))
  
  while(kappa(A[idx,]) > 1e4){
    idx = sample(1:n_stocks, size = 2, replace = FALSE, prob = abs(w))
  }
  
  lambda = solve(A[idx,], rhs[idx])      

  log_tau_vec = log_tau_vec - log(lambda[1],10)
  
  if(min(abs(w)[idx])< 1e-3){
    warning("Only one non zero weight!")
  }
  
  non_zero_idx = which(abs(w) > 1e-1)
  
  error = A %*% lambda - rhs
  
  if(max(abs(error[non_zero_idx]) > 1e-1)){
    warning("Error is larger than 1e-1")
    # plot(abs(w)[non_zero_idx], error[non_zero_idx])
  }
  
  return(log_tau_vec)
}

# get return functions --------------------------------------------------------
getRet = function(data, Para, rho_vec, use_type, mu_type, cov_type, idx = NULL, sum_w = "one", daily = FALSE, data_daily = data_daily){
  month_ret = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){
    if(daily == FALSE){
      t = Para$reform_idx[i]
    }else{
      t = Para$reform_idx_daily[i]
    }
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily = daily, data_daily = data_daily)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    if(is.null(idx)){
      mu_vec = mu_vec
      cov_mat = cov_mat
    }else{
      mu_vec = mu_vec[idx]
      cov_mat = cov_mat[idx, idx]
    }
    # for some rho it may be not possible
    rho = getRho(rho_vec, ret)
    
    w = getVanillaSolution(cov_mat, mu_vec, rho, sum_w = sum_w)
    
    t = Para$reform_idx[i]
    
    next_period_idx = (t + 1):(t + Para$period_length)
    
    if(is.null(idx)){
      month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
    }else{
      month_ret[next_period_idx] = getNextPeriodRet(data[, idx], next_period_idx, w)
    }
    
  }
  
  return(month_ret)
}

getMaxSharpeRet = function(data, Para, use_type, mu_type, cov_type, sum_w = "one", daily = FALSE, data_daily){
  month_ret = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){
    if(daily == FALSE){
      t = Para$reform_idx[i]
    }else{
      t = Para$reform_idx_daily[i]
    }
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily = daily, data_daily = data_daily)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    w = solve(cov_mat, mu_vec)/sum(solve(cov_mat, mu_vec))
    
    t = Para$reform_idx[i]
    
    next_period_idx = (t + 1):(t + Para$period_length)
    
    next_period_idx = (t + 1):(t + Para$period_length)
    
    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  return(month_ret)
}


getDemeanRet = function(data, Para, rho_vec, use_type, mu_type, cov_type, cov_demean_type = "cov(ret)", sum_w = "one", daily, data_daily = data_daily){
  month_ret = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){#
    if(daily == FALSE){
      t = Para$reform_idx[i]
    }else{
      t = Para$reform_idx_daily[i]
    }
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily = daily, data_daily = data_daily)
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    mu_mean = rowMeans(ret)
    ret_demean = ret - mu_mean
    
    mu_vec = colMeans(ret_demean) * 252/12
    if(cov_demean_type == "cov(ret_demean)"){
      cov_mat = cov(ret_demean) * 252/12
    }
    
    w = getVanillaSolution(cov_mat, mu_vec, rho = rho_vec[i], sum_w = sum_w)
    
    t = Para$reform_idx[i]
    
    next_period_idx = (t + 1):(t + Para$period_length)
    
    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  
  return(month_ret)
}


getTrueRet = function(data, Para, rho_vec, use_type, mu_type, cov_type, idx = NULL, sum_w = "one"){
  return(getRet(data, Para, rho_vec, use_type,"true", "true", idx, sum_w = sum_w))
}

getVanillaRet = function(data, Para, rho_vec, use_type, mu_type, cov_type, idx = NULL, sum_w = "one"){
  return(getRet(data, Para, rho_vec, use_type, "train", "train", idx, sum_w = sum_w))
}

get1byNRet = function(data, Para){
  month_ret = rep(NA, Para$n_time)
  w = rep(1/Para$n_stocks, Para$n_stocks)
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    next_period_idx = (t + 1):(t + Para$period_length)
    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  
  return(month_ret)
}


getMinVarTrueRet =function(data, Para){
  month_ret = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    w = getVanillaSolution(Para$cov)
    next_period_idx = (t + 1):(t + Para$period_length)
    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  return(month_ret)
}

getNoshortRet = function(data, Para, rho_vec = NULL, use_type, mu_type, cov_type){
  month_ret = rep(NA, Para$n_time)
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    # for some rho it may be not possible
    rho = getRho(rho_vec, ret)
    
    next_period_idx = (t + 1):(t + Para$period_length)
    
    if(rho < min(mu_vec) | rho > max(mu_vec)){
      month_ret[next_period_idx] = NA
    }else{
      lb = rep(0,Para$n_stocks)
      w = getVanillaSolution(cov_mat, mu_vec, rho, lb)
      month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
    }
  }
  
  return(month_ret)
}

getMinVarRet = function(data, Para, use_type, mu_type, cov_type){
  month_ret = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]

    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    mu_vec = NULL
    rho = NULL
    
    w = getVanillaSolution(cov_mat)
    
    next_period_idx = (t + 1):(t + Para$period_length)
    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  
  return(month_ret)
}

getCL2Ret = function(log_tau_TS, data, Para, rho_vec = NULL, use_type, mu_type, cov_type, sum_w = "one", daily = FALSE, data_daily = data_daily, pred = FALSE){
  month_ret = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){
    if(daily == FALSE){
      t = Para$reform_idx[i]
    }else{
      t = Para$reform_idx_daily[i]
    }
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily = daily, data_daily = data_daily, pred = pred)
    
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    rho = getRho(rho_vec, ret)
    
    if(is.null(log_tau_TS)){
      log_tau = 1/2 * log(sum(diag(cov_mat)) * (1/Para$train_length + 1/Para$period_length), 10)
    }else{
      if(length(log_tau_TS) > 1){
        log_tau = log_tau_TS[i]
      }else{
        log_tau = log_tau_TS
      }
    }

    cov_CL2 = cov_mat + 10^log_tau * diag(Para$n_stocks)
    
    w = getVanillaSolution(cov_CL2 , mu_vec, rho, sum_w = sum_w)
    
    
    t = Para$reform_idx[i]
    next_period_idx = (t + 1):(t + Para$period_length)
    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  
  return(month_ret)
}

getModMatRet = function(mod_TS, data, Para, rho_vec = NULL, use_type, mu_type, cov_type, sum_w = "one", daily = FALSE, data_daily = data_daily, type = "MeanVar", s2_vec = NULL, mult_var = 1, AR1_coef = 0){
  month_ret = rep(NA,Para$n_time)
  mu_max_sharpe = rep(NA, Para$n_time)
  s2_max_sharpe = rep(NA, Para$n_time)
  
  for(i in 1:Para$n_reform){
    if(daily == FALSE){
      t = Para$reform_idx[i]
    }else{
      t = Para$reform_idx_daily[i]
    }
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily = daily, data_daily = data_daily)
    
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    rho = getRho(rho_vec, ret)
    
    if(is.null(s2_vec)){
      s2 = sum(cov_mat)/(Para$n_stocks)^2
    }else{
      s2 = s2_vec[i]
    }
    
    cov_mod = cov_mat + mod_TS[i,,]
    if(type == "MeanVar"){
      w = getVanillaSolution(cov_mod , mu_vec, rho, sum_w = sum_w)
    }else if(type == "MaxSharpe"){
      w = solve(cov_mod, mu_vec)/sum(solve(cov_mod, mu_vec))
    }else if(type == "Dlasso"){
#       if(s2 < s2_min_var + 1e-4){ # make sure Dlasso has solution
#         # print(c(s2_min_var + 1e-4, s2))
#         s2 = s2_min_var + 1e-4
#       }
      w = getDualCLassoSolutionCone(log_tau_vec = -Inf, s2, mu_vec, cov_mod, rotate = FALSE, sum_w)
    }else if(type == "Dlasso_DTV"){
      s2_min_var = 1/sum(solve(cov_mod))
      s2 = s2_min_var * mult_var
      w = getDualCLassoSolutionCone(log_tau_vec = -Inf, s2, mu_vec, cov_mod, rotate = FALSE, sum_w)
    }else if(type == "MaxSharpeAR1"){
      w_max_sharpe = solve(cov_mod, mu_vec)/sum(solve(cov_mod, mu_vec))
      mu_max_sharpe[i] = sum(mu_vec * w_max_sharpe)
      # print(mu_max_sharpe[i])
      if(i == 1){
        w = w_max_sharpe
      }else{
        rho = AR1_coef * mu_max_sharpe[i] + (1 - AR1_coef) * mu_max_sharpe[i-1]
        w = getVanillaSolution(cov_mod , mu_vec, rho, sum_w = sum_w)
      }
    }else if(type == "DualMaxSharpeAR1"){
      w_max_sharpe = solve(cov_mod, mu_vec)/sum(solve(cov_mod, mu_vec))
      mu_max_sharpe[i] = sum(mu_vec * w_max_sharpe)
      s2_max_sharpe[i] = t(w_max_sharpe) %*% cov_mod %*% w_max_sharpe
      # print(mu_max_sharpe[i])
      print(s2_max_sharpe[i])
      if(i == 1){
        w = w_max_sharpe
      }else{
        s2 = AR1_coef * s2_max_sharpe[i] + (1 - AR1_coef) * s2_max_sharpe[i-1]
        w = getDualCLassoSolutionCone(log_tau_vec = -Inf, s2, mu_vec, cov_mod, rotate = FALSE, sum_w)
      }
    }
    
    t = Para$reform_idx[i]
    next_period_idx = (t + 1):(t + Para$period_length)
    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  
  return(month_ret)
}


getMinVarClassoRet = function(log_tau, data, Para, use_type, mu_type, cov_type){
  
  month_ret = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    # Because this is minimal variance portfolio, mu_vec and rho should be NULL
    w =  getCLassoInEqSolutionCone(ret, mu_vec = NULL, rho = NULL, log_tau)
    
    next_period_idx = (t + 1):(t + Para$period_length)
    
    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  
  return(month_ret)
}

getClassoRet = function(log_tau_TS, data, Para, rho_vec = NULL, use_type, mu_type, cov_type, sum_w = "one", daily = FALSE, data_daily = data_daily, pred = FALSE){
  month_ret = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){
    if(daily == FALSE){
      t = Para$reform_idx[i]
    }else{
      t = Para$reform_idx_daily[i]
    }
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily = daily, data_daily = data_daily, pred = pred)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    rho = getRho(rho_vec, ret)
    
    if(is.null(log_tau_TS)){
      log_tau_vec = getLogTauVecDlasso(ret, cov_mat, Para, tau_type, prob)
    }else{
      if(length(log_tau_TS) == 1){
        log_tau_vec = rep(log_tau_TS, Para$n_stocks)
      }else{
        log_tau_vec = rep(log_tau_TS[i], Para$n_stocks)
      }
    }
    
    # w = getCLassoSolutionCone(log_tau_vec, rho, mu_vec, cov_mat)
    w = getCLassoInEqSolutionCone(log_tau_vec, rho, mu_vec, cov_mat, sum_w = sum_w)
    
    t = Para$reform_idx[i]
    
    next_period_idx = (t + 1):(t + Para$period_length)
    
    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  
  return(month_ret)
}

getDualClassoRet = function(log_tau_vec = NULL, data, Para, s2_vec = NULL, use_type, mu_type, cov_type, tau_type = NULL, prob = 0.95, mult = NULL, rotate = FALSE, sum_w = "one", daily = daily, data_daily = data_daily){
  
  
  log_tau_vec0 = log_tau_vec # keep a copy of the input log_tau_vec
  
  
  month_ret = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){
    if(daily == FALSE){
      t = Para$reform_idx[i]
    }else{
      t = Para$reform_idx_daily[i]
    }
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily =  daily, data_daily = data_daily)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    if(is.null(log_tau_vec0)){
      log_tau_vec = getLogTauVecDlasso(ret, cov_mat, Para, tau_type, prob, mult)
    }else{
      if(length(log_tau_vec0) == 1){
        log_tau_vec0 = rep(log_tau_vec0, Para$n_reform)
      }
      log_tau_vec = rep(log_tau_vec0[i], Para$n_stocks)
    }
    
    if(is.null(s2_vec)){
      s2 = getInSampleS2(mu_vec, cov_mat, log_tau_vec)
    }else{
      s2 = s2_vec[i]
    }
    
    if(tau_type == "eigen" | tau_type == "eigen unit"){
      w = getDCLassoWEigenSolution(log_tau_vec, s2, mu_vec, cov_mat, sum_w = sum_w)
    }else{
      w = getDualCLassoSolutionCone(log_tau_vec, s2, mu_vec, cov_mat, rotate = rotate, sum_w = sum_w)
    }
    
    t = Para$reform_idx[i]
    
    next_period_idx = (t + 1):(t + Para$period_length)
    
    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  
  return(month_ret)
}

getDL2Ret = function(log_tau = NULL, data, Para, s2_vec = NULL, use_type, mu_type, cov_type, sum_w = "one", mult = mult){
  log_tau0 = log_tau
  month_ret = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    if(is.null(log_tau0)){
      log_tau = log(mult * sqrt(sum(diag(cov_mat))/Para$train_length), 10)
      # print(round(10^log_tau, 2))
    }
    
    if(is.null(s2_vec)){
      s2 = getInSampleS2(mu_vec, cov_mat, log_tau_vec)
    }else{
      s2 = s2_vec[i]
    }
    
    w = getDL2Solution(log_tau, s2, mu_vec, cov_mat, sum_w = sum_w)
    
    next_period_idx = (t + 1):(t + Para$period_length)
    
    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  
  return(month_ret)
}


getRetFromInvCovTS = function(inv_cov_TS, data, Para, rho_vec = NULL, use_type, mu_type, cov_type){
  month_ret = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]

    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    inv_cov = inv_cov_TS[,,i]
    
    # for some rho it may be not possible
    rho = getRho(rho_vec, ret)
    
    w = getSolution(inv_cov, mu_vec, rho)
    
    next_period_idx = (t + 1):(t + Para$period_length)

    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  
  return(month_ret)
}

getRetFromAlphaVarUnExp = function(data, Para, use_type, alpha_type, var_unExp_type){
  month_ret = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    
    train_info = switch(alpha_type,
                        true = getTrainInfo(data, Para, t, use_type, "true", "true"),
                        train = getTrainInfo(data, Para, t, use_type, "train", "train"))
    
    alpha_vec = getAlpha(solve(train_info$cov_mat), train_info$mu_vec)
    
    train_info = switch(var_unExp_type,
                        true = getTrainInfo(data, Para, t, use_type, "true", "true"),
                        train = getTrainInfo(data, Para, t, use_type, "train", "train"))
    
    var_unExp = 1/diag(solve(train_info$cov_mat))
    
    w = (alpha_vec/var_unExp)/sum(alpha_vec/var_unExp)
    
    next_period_idx = (t + 1):(t + Para$period_length)
    
    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  
  return(month_ret)
}

# getRetPCA = function(n_PCA, PCA_type, data, Para, use_type, mu_type, cov_type){
#   month_ret = rep(NA, Para$n_time)
#   
#   alpha_TS = matrix(NA, nrow = Para$n_reform, ncol = Para$n_stocks)
#   varUE_TS = matrix(NA, nrow = Para$n_reform, ncol = Para$n_stocks)
#   R2_TS = matrix(NA, nrow = Para$n_reform, ncol = Para$n_stocks)
#   w_TS = matrix(NA, nrow = Para$n_reform, ncol = Para$n_stocks)
#   
#   for(i in 1:Para$n_reform){
#     t = Para$reform_idx[i]
#     
#     train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
#     mu_vec = train_info$mu_vec
#     cov_mat = train_info$cov_mat
#     ret = train_info$ret
#     use_idx = train_info$use_idx
#     
#     if(PCA_type == "all"){
#       PCA = prcomp(ret)
#       PCA_first = ret %*% PCA$rotation[,1:n_PCA]
#       
#       for(j in 1:Para$n_stocks){
#         y = ret[,j]
#         X = PCA_first
#         fit = summary(lm(y ~ X))
#         alpha_TS[i,j] = fit$coefficients[1,1]
#         varUE_TS[i,j] = mean(fit$residuals^2)
#         R2_TS[i,j] = fit$r.squared
#       }
#     }else if(PCA_type == "other"){
#       for(j in 1:Para$n_stocks){
#         y = ret[,j]
#         X = ret[,-j]
#         PCA = prcomp(X)
#         X = X %*% PCA$rotation[,1:n_PCA]
#         # y = ret[,j] - RF[use_idx,1]
#         # X = Factors[use_idx,]
#         fit = summary(lm(y ~ X))
#         alpha_TS[i,j] = fit$coefficients[1,1]
#         varUE_TS[i,j] = mean(fit$residuals^2)
#         R2_TS[i,j] = fit$r.squared
#       }
#     }
#     
#     w = (alpha_TS[i,]/varUE_TS[i,])/sum(alpha_TS[i,]/varUE_TS[i,])
#     w_TS[i,] = w
#     
#     next_period_idx = (t + 1):(t + Para$period_length)
#     
#     month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
#   }
#   
#   return(month_ret)
# }

getPCARet = function(data, Para, use_type, mu_type, cov_type, sum_w = "one", daily = FALSE, data_daily = data_daily){
  month_mat_MV = matrix(NA, Para$n_time, Para$n_stocks)
  month_mat_MS = matrix(NA, Para$n_time, Para$n_stocks)
  
  for(i in 1:Para$n_reform){
    if(daily == FALSE){
      t = Para$reform_idx[i]
    }else{
      t = Para$reform_idx_daily[i]
    }
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily = daily, data_daily = data_daily)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    eigen_cov = eigen(cov_mat)
    one_vec = rep(1, Para$n_stocks)
    
    for(j in 1:Para$n_stocks){
      sub_idx = 1:j
      
      w = getSubsetPCASolution(mu_vec, eigen_cov, sub_idx)
      
      t = Para$reform_idx[i]
      next_period_idx = (t + 1):(t + Para$period_length)
      
      month_mat_MV[next_period_idx, j] = getNextPeriodRet(data, next_period_idx, w$MV)
      month_mat_MS[next_period_idx, j] = getNextPeriodRet(data, next_period_idx, w$MS)
    }
  }
  
  return(list(MV = month_mat_MV, MS = month_mat_MS))
}
  

#RGL = regularized
getPCAAlphaRGLRet = function(n_PCA, sum_alpha_LB, data, Para, use_type, mu_type, cov_type){
  month_ret = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    PCA = prcomp(ret)
    factors = ret %*% PCA$rotation[,1:n_PCA]
    
    w = getSolutionOfRegs(Y = ret, factors, sum_alpha_LB)
    hist(w)
    next_period_idx = (t + 1):(t + Para$period_length)
    
    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  
  return(month_ret)
}

getFactorsAlphaRGLRet = function(factors, sum_alpha_LB, data, Para, use_type, mu_type, cov_type, cheat_type = NULL){
  month_ret = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    
    if(is.null(cheat_type)){
      train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type)
      mu_vec = train_info$mu_vec
      cov_mat = train_info$cov_mat
      ret = train_info$ret
      use_idx = train_info$use_idx
    }else{
      use_idx = (t + 1):(t + Para$period_length)
      ret = data[use_idx,]
    }
    
    
    sub_factors = as.matrix(factors)[use_idx,]
    
    w = getSolutionOfRegs(Y = ret, sub_factors, sum_alpha_LB)
    # hist(w)
    next_period_idx = (t + 1):(t + Para$period_length)
    
    month_ret[next_period_idx] = getNextPeriodRet(data, next_period_idx, w)
  }
  
  return(month_ret)
}


getOneByNEigenVectorRet = function(data, Para, use_type, mu_type, cov_type, daily = FALSE, data_daily = data_daily){
  ret_mat = matrix(NA, Para$n_time, Para$n_stocks)
  for(i in 1:Para$n_reform){
    if(daily == FALSE){
      t = Para$reform_idx[i]
    }else{
      t = Para$reform_idx_daily[i]
    }
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily = daily, data_daily = data_daily)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    eigen_result = eigen(cov_mat)
    negative_idx = which(mu_vec %*% eigen_result$vectors < 0)
    eigen_result$vectors[,negative_idx] = - eigen_result$vectors[,negative_idx]
    
    w_mat = t(apply(eigen_result$vectors[,(Para$n_stocks:1)], MARGIN = 1, function(x) cumsum(x)/(1:length(x))))
    
    w_mat = apply(w_mat, MARGIN = 2, function(x) x/sum(x))
    
    t = Para$reform_idx[i]
    
    next_period_idx = (t + 1):(t + Para$period_length)
    
    ret_mat[next_period_idx,] = apply(w_mat, MARGIN = 2, getNextPeriodRet, data = data, next_period_idx = next_period_idx)
    
  }
  return(ret_mat)
}

getAllKRet = function(data, Para, rho_vec, use_type, mu_type, cov_type, sum_w = "one", daily = FALSE, data_daily){
  ret_mat = matrix(NA, Para$n_time, Para$n_stocks)
  n_LB_vec = rep(NA, Para$n_reform)
  log_tau_mat = matrix(NA, Para$n_reform, Para$n_stocks)
  
  for(i in 1:Para$n_reform){
    if(daily == FALSE){
      t = Para$reform_idx[i]
    }else{
      t = Para$reform_idx_daily[i]
    }
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily = daily, data_daily = data_daily)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    rho = getRho(rho_vec, ret)
    
    t = Para$reform_idx[i]
    next_period_idx = (t + 1):(t + Para$period_length)
    
    log_tau = 5
    w = getCLassoInEqSolutionCone(log_tau, rho, mu_vec, cov_mat, sum_w = "one")
    n_LB_vec[i] = getNonZero(w, Para$tol)
    
    print(c(i,n_LB_vec[i]))
    log_tau_mat[i, 1:n_LB_vec[i]] = log_tau
    ret_mat[next_period_idx, 1:n_LB_vec[i]] = getNextPeriodRet(data, next_period_idx, w)
    
    if(n_LB_vec[i] < Para$n_stocks){
      for(K in (n_LB_vec[i] + 1):Para$n_stocks){
        # print(K)
        log_tau = getLogTauForFixedK(K, log_tau, rho, mu_vec, cov_mat, sum_w, Para$tol)
        log_tau_mat[i, K] = log_tau
        w = getCLassoInEqSolutionCone(log_tau, rho, mu_vec, cov_mat, sum_w = "one")
        ret_mat[next_period_idx, K] = getNextPeriodRet(data, next_period_idx, w)
      }
    }
  }
  
  result = list(log_tau_mat = log_tau_mat, n_LB_vec = n_LB_vec, ret_mat = ret_mat)
  
  return(result)
}


getInSampleBestKRet = function(ret_mat, Para, type = "prevAll"){
  month_ret = rep(NA, Para$n_time)
  best_K = rep(NA, Para$n_reform)
  
  for(i in 1:Para$n_reform){
    t = Para$reform_idx[i]
    if(type == "prevAll"){
      use_idx = 1:t
    }else{
      use_idx = (t - Para$period_length + 1):t
    }
    
    next_period_idx = (t + 1):(t + Para$period_length)
    
    if(i == 1){#first period log_tau is 5
      month_ret[next_period_idx] = ret_mat[next_period_idx,1]
    }else{
      n_non_NA = sum(!is.na(ret_mat[use_idx,1]))
      if(n_non_NA > 1){ #in order to calculate sd, we need at least 2 observation
        Sharpe_vec = apply(ret_mat[use_idx,], 2, sharpeRatio)
        # plot(Sharpe_vec, ylim = c(0.2, 0.5), main = i)
        best_K[i] = which.max(Sharpe_vec)
        # abline(v = best_K[i], col = 2)
      }else{
        best_K[i] = 1
      }
      
      month_ret[next_period_idx] = ret_mat[next_period_idx, best_K[i]]
    }
  }
  
  result = list(month_ret = month_ret, best_K = best_K)
  return(result)
}

# Sharpe Ratio functions --------------------------------------------------
# Sharpe Ratio of different types of training set
sharpeMarkowitz1 = function(Q_type_seq,mu_type_seq,
                            data,Para,rho_vec = NULL){
  n_Q_type = length(Q_type_seq)
  n_mu_type = length(mu_type_seq)
  
  Sharpe = matrix(NA,n_Q_type,n_mu_type)
  rownames(Sharpe) = paste("Q",Q_type_seq)
  colnames(Sharpe) = paste("mu",mu_type_seq)
  for(i in seq_along(Q_type_seq)){
    for(j in seq_along(mu_type_seq)){
      Sharpe[i,j] = 
        analyseVanillaWQMutype(data,Para,
                               Q_type = Q_type_seq[i],mu_type = mu_type_seq[j],rho_vec)
    }
  }
  return(Sharpe)
}

# Sharpe Ratio of different length of training set
sharpeMarkowitz2 = function(train_cov_from_vec,train_cov_to_vec,
                            train_mu_from_vec,train_mu_to_vec,
                            data,Para,rho_vec = NULL){
  n_Q_from = length(train_cov_from_vec)
  n_mu_from = length(train_mu_from_vec)
  
  Sharpe = matrix(NA,n_Q_from,n_mu_from)
  colnames(Sharpe) = paste("mu",train_mu_from_vec,":",train_mu_to_vec)
  rownames(Sharpe) = paste("cov",train_cov_from_vec,":",train_cov_to_vec)
  
  for(i in 1:n_Q_from){
    for(j in 1:n_mu_from){
      Sharpe[i,j] = analyseTrain(data,Para,
                                train_cov_from_vec[i],train_cov_to_vec[i],
                                train_mu_from_vec[j],train_mu_to_vec[j],rho_vec)
    }
  }
  return(Sharpe)
}

# Sharpe Ratio of different types of training set
sharpeMarkowitz3 = function(cov_type_seq,mu_type_seq,
                            data,Para,rho_vec = NULL){
  n_cov_type = length(cov_type_seq)
  n_mu_type = length(mu_type_seq)
  
  Sharpe = matrix(NA,n_cov_type,n_mu_type)
  rownames(Sharpe) = paste("cov",cov_type_seq)
  colnames(Sharpe) = paste("mu",mu_type_seq)
  for(i in seq_along(cov_type_seq)){
    for(j in seq_along(mu_type_seq)){
      Sharpe[i,j] = 
        analyseVanillaWCovMutype(data,Para,
                               cov_type = cov_type_seq[i],mu_type = mu_type_seq[j],rho_vec)
    }
  }
  return(Sharpe)
}


sharpeTrue = function(rho_vec,data,Para){
  month_ret = getTrueRet(rho_vec,data,Para)  
  
  displayResult(month_ret,"true",Para$train_length)
  Sharpe = sharpeRatio(month_ret)
  return(Sharpe)
}

sharpeMinVarTrue = function(data,Para){
  month_ret = getMinVarTrueRet(data,Para)
  
  displayResult(month_ret,"minimal variance true",Para$train_length)
  Sharpe = sharpeRatio(month_ret)
  return(Sharpe)
}

sharpe1byN = function(data,Para){
  month_ret = get1byNRet(data,Para)
  
  displayResult(month_ret,"1/N",Para$train_length)
  Sharpe = sharpeRatio(month_ret)
  return(Sharpe)
}

sharpeNoShort = function(data,Para,rho_vec = NULL){
  month_ret = getNoShortRet(data,Para,rho_vec)
  
  Sharpe = sharpeRatio(month_ret)
  displayResult(month_ret,"No short",Para$train_length)
  
  return(Sharpe)
}

sharpeMinVar = function(data,Para){
  month_ret = getMinVarRet(data,Para)
  
  Sharpe = sharpeRatio(month_ret)
  displayResult(month_ret,"Minimal Variance",Para$train_length)
  
  return(Sharpe)
}

sharpeMinVarFixTau = function(log_tau_seq,data,Para){
  n_tau = length(log_tau_seq)
  
  Sharpe = rep(NA,n_tau)
  names(Sharpe) = log_tau_seq
  non_zero_mean = rep(NA,n_tau)
  non_zero_min = rep(NA,n_tau)
  non_zero_max = rep(NA,n_tau)
  
  pb = txtProgressBar(min = 1, max = n_tau)
  for(j in 1:n_tau){
    setTxtProgressBar(pb,j)
    log_tau = log_tau_seq[j]
    month_ret = getMinVarClassoRet(log_tau,data,Para)
    
    Sharpe[j] = sharpeRatio(month_ret)
    non_zero_mean[j] = mean(n_non_zero, na.rm = TRUE)
    non_zero_min[j] = min(n_non_zero, na.rm = TRUE)
    non_zero_max[j] = max(n_non_zero, na.rm = TRUE)
    
    
  }
  
  return(Sharpe)
}

sharpeFixTau = function(log_tau_seq,data,Para,rho_vec = NULL){
  n_tau = length(log_tau_seq)
  
  Sharpe = rep(NA,n_tau)
  names(Sharpe) = log_tau_seq
  non_zero_mean = rep(NA,n_tau)
  non_zero_min = rep(NA,n_tau)
  non_zero_max = rep(NA,n_tau)
  
  pb = txtProgressBar(min = 1, max = n_tau)
  for(j in 1:n_tau){
    setTxtProgressBar(pb,j)
    log_tau = log_tau_seq[j]
    month_ret = getClassoRet(log_tau,data,Para,rho_vec)
    
    Sharpe[j] = sharpeRatio(month_ret)
    non_zero_mean[j] = mean(n_non_zero, na.rm = TRUE)
    non_zero_min[j] = min(n_non_zero, na.rm = TRUE)
    non_zero_max[j] = max(n_non_zero, na.rm = TRUE)
  }
  
  return(Sharpe)
}

SharpeRidgeFixTau = function(log_tau_seq,data,Para,rho_vec = NULL){
  n_tau = length(log_tau_seq)
  
  Sharpe = rep(NA,n_tau)
  names(Sharpe) = log_tau_seq
  non_zero_mean = rep(NA,n_tau)
  non_zero_min = rep(NA,n_tau)
  non_zero_max = rep(NA,n_tau)
  
  pb = txtProgressBar(min = 1, max = n_tau)
  for(j in 1:n_tau){
    setTxtProgressBar(pb,j)
    log_tau = log_tau_seq[j]
    month_ret = getRidgeRet(log_tau,data,Para,rho_vec)
      
    Sharpe[j] = sharpeRatio(month_ret)
    non_zero_mean[j] = mean(n_non_zero, na.rm = TRUE)
    non_zero_min[j] = min(n_non_zero, na.rm = TRUE)
    non_zero_max[j] = max(n_non_zero, na.rm = TRUE)
  }
  
  return(Sharpe)
}






# cos(mv_vec,1) -----------------------------------------------------------
simCos = function(mu_type_seq,data,Para){
  n_mu_type = length(mu_type_seq)
  
  cos_mat = matrix(NA, Para$n_reform - 1,n_mu_type)
  rownames(cos_mat) = 1:(Para$n_reform - 1)
  colnames(cos_mat) = paste("mu",mu_type_seq)
  for(j in seq_along(mu_type_seq)){
    cos_mat[,j] = analyseCosMutype(data,Para,mu_type = mu_type_seq[j])
  }

  return(cos_mat)
}


# Plot --------------------------------------------------------------------
frontierSharpePlot = function(title,Para,Sharpe,Sharpe_mat,rho_seq,log_tau_seq,Sim = NULL){
  F_col = 1:20
  p_col = 1:2
  n_tau = length(log_tau_seq)
  
  plot(x = NA, y = NA,xlab = "rho",
       ylab = "Sharpe Ratio", main = paste(Para$data_name,title),
       xlim = c(0,7),ylim = c(0.1,0.6))
  
  for(i in 1:n_tau){
    points(x = rho_seq, y = Sharpe_mat[,i], col = F_col[i], pch = 17)
  }
  
  name = paste0("Sharpe$F_mu_true_train_",Para$train_length)
  
  points(x = rho_seq, y = eval(parse(text = name)), col = F_col[1], pch = 19)
  # points(x = rho_seq, y = Sharpe$F_no_short, col = F_col[2], pch = 19)
  # points(x = rho_seq, y = Sharpe$F_prev_all, col = F_col[3], pch = 19)
  
  if(is.null(Sim)){
    n = 1
    
    lines = c(Sharpe$one_by_N,Sharpe$min_var)
    n_lines = length(lines)
    abline(h = lines,col = 1:n_lines)
    abline(v = mean(Para$rho_vec),col = 1 + n_lines)
    
    legend("right", legend = c("Markowitz",log_tau_seq,"1/N","Min Var"), col = c(1:n,F_col[1:n_tau],1:n_lines), lty = c(rep(NA,n_tau + n),rep(1,n_lines)), pch = c(rep(19,n),rep(17,n_tau),rep(NA,n_lines)))
    
    # ,"No Short","Previsou All"
  }else{
    name = paste0("Sharpe$F_mu_true_true_",Para$train_length)
    points(x = rho_seq, y = eval(parse(text = name)), col = F_col[4], pch = 19)
    n = 2
    
    lines = c(Sharpe$one_by_N,Sharpe$min_var,Sharpe$true_min_var)
    n_lines = length(lines)
    abline(h = lines,col = 1:n_lines)
    abline(v = mean(Para$rho_vec),col = 1 + n_lines)
    
    legend("right", legend = c("Markowitz","True",log_tau_seq,"1/N","Min Var","True Min Var"), col = c(1:n,F_col[1:n_tau],1:n_lines), lty = c(rep(NA,n_tau + n),rep(1,n_lines)), pch = c(rep(19,n),rep(17,n_tau),rep(NA,n_lines)))
    
    # ,"No Short","Previsou All"
  }
}

frontierSDMuPlot = function(title,Para,mu,mu_mat,SD,SD_mat,Sim = NULL){
  rho_seq = as.numeric(rownames(mu$F_Classo))
  log_tau_seq = colnames(mu_mat)
  
  F_col = 1:20
  p_col = 1:3
  n_tau = length(log_tau_seq)
  
  plot(x = NA, y = NA, xlab = "standard deviation",
       ylab = "expected return",
       main = paste(Para$data_name,title),
       xlim = c(0,10),ylim = range(rho_seq))
  
  for(i in 1:n_tau){
    lines(x = SD_mat[,i], y = mu_mat[,i], col = F_col[i], type = "b", pch = 17)
  }
  
  name1 = paste0("SD$F_mu_train_train_",Para$train_length)
  name2 = paste0("mu$F_mu_train_train_",Para$train_length)
  
  lines(x = eval(parse(text = name1)), y = eval(parse(text = name2)), col = F_col[1], type = "b", pch = 19)
  lines(x = SD$F_no_short, y = mu$F_no_short, col = F_col[2], type = "b", pch = 19)
  lines(x = SD$F_prev_all, y = mu$F_prev_all, col = F_col[3], type = "b", pch = 19)
  
  if(is.null(Sim)){
    n = 3
    points(x = c(SD$one_by_N,SD$min_var),
           y = c(mu$one_by_N,mu$min_var), col = p_col, pch = 3)
    abline(h = mean(Para$rho_vec), col = 1)
    
    legend("topleft", legend = c("Markowitz","No Short","Previous All",log_tau_seq,"1/N","Minimal Variance"), col = c(1:n,F_col[1:n_tau],p_col), pch = c(rep(19,n),rep(17,n_tau),3,3))
    
  }else{
    name1 = paste0("SD$F_mu_true_true_",Para$train_length)
    name2 = paste0("mu$F_mu_true_true_",Para$train_length)
    
    
    points(x = eval(parse(text = name1)), y = eval(parse(text = name2)), col = F_col[4], pch = 19)
    n = 4
    points(x = c(SD$one_by_N,SD$min_var),
           y = c(mu$one_by_N,mu$min_var), col = p_col, pch = 3)
    abline(h = mean(Para$rho_vec), col = 1)
    
    legend("right", legend = c("Markowitz","No Short","Previous All","True",log_tau_seq,"1/N","Min Var","True Min Var"), col = c(1:n,F_col[1:n_tau],p_col), pch = c(rep(19,n),rep(17,n_tau),rep(3,3)))
  }
}


matrixHeatPlot = function(mat,Para,breaks,title){
  # make n_breaks the smallest even number that is larger than length(breaks)
  # I want the palette to be symmetric
  n_breaks = length(breaks) - length(breaks)%%2 
  palette = divPalette(n_breaks, name = "PRGn")
  
  # data for heat plot
  data_HP = data.frame(
    Row = rep(1:Para$n_stocks,Para$n_stocks),
    Col = rep(1:Para$n_stocks,each = Para$n_stocks),
    Z =  cut(mat,breaks = breaks)
  )
  
  idx = which(table(data_HP$Z)>0)
  palette = palette[idx]
  
  ggplot(data =  data_HP, aes(x = Row, y = Col)) + 
    geom_tile(aes(fill = Z), colour = "white") +
    scale_fill_manual(values = palette) + ggtitle(title)
}

histCoefAlongTime = function(coef_mat,reform_idx,tol){
  x_min = floor(min(coef_mat))
  x_max = ceiling(max(coef_mat))
  breaks = seq(x_min,x_max,by = 0.2)
  
  n_reform = length(reform_idx)
  for(i in 1:n_reform){
    idx = which(abs(coef_mat[i,])>tol)
    hist(coef_mat[i,idx],breaks = breaks,
         main = paste0("Hist of weights at month ",reform_idx[i]),
         xlab = "weights")
    # readline(prompt = "Press Enter to continue.")
  }
}


scatterPlot = function(x, y, title, xlab = "true", ylab = "L2", limit = NULL, diag_idx = NULL){
  if(is.null(limit)){
    limit = range(c(x, y))
  }
  
  if(is.null(diag_idx)){
    plot(x, y, xlim = limit, ylim = limit, main = title, xlab = xlab, ylab = ylab)
  }else{
    plot(x[diag_idx], y[diag_idx], xlim = limit, ylim = limit, main = title, xlab = xlab, ylab = ylab)
    points(x[-diag_idx], y[-diag_idx], col = 2)
  }
  
  abline(a = 0, b = 1, col = 2 )
}

bestSharpePlot = function(mu, SD, methods, title, best_Sharpe){
  rho_seq = rownames(mu$F_Classo)
  if(rho_seq[1] == "oneByN"){
    ylimit = c(0.5, 2.5)
  }else{
    ylimit = range(as.numeric(rho_seq))
  }
  
  n_methods = length(methods)

  mu_mat = NULL
  SD_mat = NULL
  for(i in 1:n_methods){
    m = methods[i]
    best_col = as.numeric(best_Sharpe["Index",m])
    mu_mat = cbind(mu_mat, eval(parse(text = paste0("mu$F_",m)))[,best_col])
    SD_mat = cbind(SD_mat, eval(parse(text = paste0("SD$F_",m)))[,best_col])
  }
  
  matplot(x = SD_mat, y = mu_mat, type = 'p', lty = 1, pch = 1:n_methods,lwd = 2, xlab = "standard deviation", ylab = "expected return", main = paste(Para$data_name,title), xlim = c(2,10), ylim = ylimit)
  
  idx = n_methods + 1
  
  name1 = paste0("SD$F_mu_train_train_",Para$train_length)
  name2 = paste0("mu$F_mu_train_train_",Para$train_length)
  
  lines(x = eval(parse(text = name1)), y = eval(parse(text = name2)), type = "b", col = idx, pch = idx, lwd = 2)
  
  idx = idx + 1 
  
  name1 = paste0("SD$F_mu_true_train_",Para$train_length)
  name2 = paste0("mu$F_mu_true_train_",Para$train_length)
  
  lines(x = eval(parse(text = name1)), y = eval(parse(text = name2)), type = "b", col = idx, pch = idx, lwd = 2)
  idx = idx + 1
  
  name1 = paste0("SD$F_mu_true_true_",Para$train_length)
  name2 = paste0("mu$F_mu_true_true_",Para$train_length)
  
  lines(x = eval(parse(text = name1)), y = eval(parse(text = name2)), type = "b", col = idx, pch = idx, lwd = 2)
  
#   idx = idx + 1
#   lines(x = SD$F_Demean[,"cov(ret)"], y = mu$F_Demean[,"cov(ret)"], type = "b", col = idx, pch = idx, lwd = 2)
#   
#   idx = idx + 1
#   lines(x = SD$F_Demean[,"cov(ret_demean)"], y = mu$F_Demean[,"cov(ret_demean)"], type = "b", col = idx, pch = idx, lwd = 2)
  
#   idx = idx + 1
#   lines(x = SD$F_mu_train_train_daily, y = mu$F_mu_train_train_daily, type = "b", col = idx, pch = idx, lwd = 2)
  
  legend_names = c(paste(methods, best_Sharpe["Parameter",]), "train train", "true train", "true true") #, "demean cov(ret)", "demean cov(ret_demean)", , "daily Markowitz"
  legend("topleft",legend = legend_names, col = 1:idx, lty = 1, lwd = 2,pch = 1:idx)
}

corPlot = function(cor_ret, methods, rho_seq, best_Sharpe){
  n_methods = length(methods)
  n_comb = choose(n_methods,2)
  
  cor_best = matrix(NA, nrow = length(rho_seq), ncol = n_comb)
  rownames(cor_best) = rho_seq
  colnames(cor_best) = rep(NA, n_comb)
  col_idx = 1
  
  for(i in 1:(n_methods-1)){
    for(j in (i+1): n_methods){
      method1 = methods[i]
      method2 = methods[j]
      
      idx1 = as.numeric(best_Sharpe["Index", method1])
      idx2 = as.numeric(best_Sharpe["Index", method2])
      para1 = best_Sharpe["Parameter", method1]
      para2 = best_Sharpe["Parameter", method2]
      
      cor_array = eval(parse(text = paste0("cor_ret$", method1, "_", method2)))
      if(is.null(cor_array)){
        cor_array = eval(parse(text = paste0("cor_ret$", method2, "_", method1)))
        cor_best[,col_idx] = cor_array[idx2, idx1, ]
      }else{
        cor_best[,col_idx] = cor_array[idx1, idx2, ]
      }
      
      colnames(cor_best)[col_idx] = paste(method1, para1, method2, para2, sep = "_")
      col_idx = col_idx + 1
    }
  }
  if(rho_seq == "oneByN"){
    matplot(x = 1, y = cor_best, ylim = c(0.6, 1), ylab = "correlation", xlab = "rho", main = "Correlation between methods", type = "b",lty = 1, pch = 1:n_comb, lwd = 2)
  }else{
    matplot(x = rho_seq, y = cor_best, ylim = c(0.6, 1), ylab = "correlation", xlab = "rho", main = "Correlation between methods", type = "b",lty = 1, pch = 1:n_comb, lwd = 2)
  }
  
  legend("topleft", legend = colnames(cor_best), col = 1:n_comb, lty = 1, pch = 1:n_comb, lwd = 2)
}

muIntPlot = function(mu_vec, SD_vec, mult, title){
  plot(mu_vec, main = title, ylab = "mu")
  for(i in 1:Para$n_stocks){
    mu = mu_vec[i]
    SD = SD_vec[i]
    lb = mu - mult * SD
    ub = mu + mult * SD
    lines(x = c(i,i), y = c(lb, ub), type = "l", col = 2)
  }
}


# Variable Selection ------------------------------------------------------
VarSelectWPvalue = function(y, X, p){
  beta = rep(0, ncol(X))
  
  fit_all = summary(lm(y ~ X))
  p_values = fit_all$coefficients[,4]
  idx = which(p_values[-1] <= p)
  if(length(idx) > 0){
    fit_select = summary(lm(y ~ X[,idx]))
    # the first element of coefficients is intercept
    beta[idx] = fit_select$coefficients[-1,1]
    var_unExp = sum(fit_select$residuals^2)/nrow(X)
  }else{
    # beta are all 0 which is the same as the initialization
    var_unExp = mean((y - mean(y))^2)
  }
  
  
  return(result = list(beta = beta, var_unExp = var_unExp))
}

VarSelectStep = function(data, i, direction){
  n_stocks = ncol(data)
  colnames(data) = paste0("x", 1:n_stocks)
  
  full = lm(paste0("x", i, " ~ ."), data = data)
  null = lm(paste0("x", i, " ~ 1"), data = data)
  
  
  if(direction == "both"){
    result = step(full, scope=list(upper = full, lower = null), direction = "both",
                  trace = FALSE) 
  }else if(direction == "forward"){
    result = step(null, scope=list(upper=full, lower=null), direction = "forward",
                  trace = FALSE)
  }else if(direction == "backward"){
    result = step(full, scope=list(upper=full, lower=null), direction = "backward",
                  trace = FALSE)
  }else{
    stop("direction is wrong")
  }
  
  names = names(result$coefficients)
  idx = as.numeric(gsub(pattern = "x", replacement = "", names[-1]))
  
  beta = rep(0, n_stocks)
  beta[idx] = result$coefficients[-1]
  
  beta = beta[-i]
  
  var_unExp = sum(result$residuals^2)/nrow(data)
  
  return(result = list(beta = beta, var_unExp = var_unExp))
}

VarSelectGlasso = function(data, inv_cov_glasso, i){
  n_stocks = ncol(data)
  y = as.matrix(data[,i])
  X = as.matrix(data[,-i])
  
  non_zero_idx = which(abs(inv_cov_glasso[i,-i]) > 1e-16)

  if(length(non_zero_idx) == 0){
    beta = rep(0, n_stocks - 1)
    var_unExp = var(y)
  }else{
    X = X[, non_zero_idx]
    result = summary(lm(y ~ X))
    
    beta = rep(0, n_stocks - 1)
    beta[non_zero_idx] = result$coefficients[-1,1]
    
    var_unExp = sum(result$residuals^2)/nrow(data) 
  }
  
  return(result = list(beta = beta, var_unExp = var_unExp))
}




# get in sample s2 estimation ---------------------------------------------
getGamma = function(cov_mat){
  # t(Gamma) %*% Gamma = as.numeric(t(one_vec) %*% inv_cov %*% one_vec) * inv_cov  - inv_cov %*% one_vec %*% t(one_vec) %*% inv_cov
  n_stocks = ncol(cov_mat)
  
  one_vec = matrix(1, nrow = n_stocks, 1)
  
  inv_cov = solve(cov_mat)
  
  M = as.numeric(t(one_vec) %*% inv_cov %*% one_vec) * inv_cov  - inv_cov %*% one_vec %*% t(one_vec) %*% inv_cov
  
  eigen_result = eigen(M)
  eigen_vec = eigen_result$values
  eigen_vec[n_stocks] = 0
  P = eigen_result$vectors
  
  Gamma = diag(sqrt(eigen_vec)) %*% t(P)
  Gamma = Gamma[-n_stocks,]
  
  return(Gamma)
}


getMaxVarWRTMu = function(cov_mat, mu_vec, log_tau_vec, tau_type = "LInf indi t", prob = 0.99, proposed){
  n_stocks = ncol(cov_mat)
  
  one_vec = matrix(1, nrow = n_stocks, 1)
  Gamma = getGamma(cov_mat)
  
  w_0 = rep(1, n_stocks)/n_stocks
  w_min_var = as.numeric(getMinVarSol(solve(cov_mat)))
  
  log_tau_vec = getLogTauVecDlasso(ret, cov_mat, Para, tau_type, prob = 0.99)
  
  w_0 = rep(1, n_stocks)/n_stocks
  w_min_var = as.numeric(getMinVarSol(solve(cov_mat)))
  
  
  # Conic programming begins
  mu_idx = 1:n_stocks
  y_idx = 1:(n_stocks - 1) + max(mu_idx)
  z_idx = max(y_idx) + 1
  
  n_var = max(z_idx)
  n_lconstr = n_stocks
  
  A = matrix(0, nrow = n_lconstr, ncol = n_var)
  # z = sqrt(t) (w_0 - w_min_var)^t mu
  row_idx = 1
  A[row_idx, z_idx] = 1
  A[row_idx, mu_idx] = - (w_0 - w_min_var)
  # y = Gamma %*% mu
  row_idx = max(row_idx) + 1:(n_stocks-1)
  A[row_idx, y_idx] = diag(n_stocks-1)
  A[row_idx, mu_idx] = -Gamma*sqrt(proposed)
  
  rhs = rep(0, n_lconstr)
  sense = rep("=", n_lconstr)
  
  if(length(log_tau_vec) == 1){
    log_tau_vec = rep(log_tau_vec, n_stocks)
  }
  
  obj = rep(0, n_var)
  
  cones = list(as.list(c(z_idx, y_idx)))
  
  lb = rep(-Inf, n_var)
  lb[mu_idx] = mu_vec - 10^log_tau_vec
  lb[z_idx] = 0
  
  ub = rep(+Inf, n_var)
  ub[mu_idx] = mu_vec + 10^log_tau_vec
  
  model = list()
  model$A = A
  model$rhs = rhs
  model$sense = sense
  model$obj = obj
  model$ub = ub
  model$lb = lb
  model$cones = cones
  
  params = list(OutputFlag = 0, BarHomogeneous = 1) #,
  result = gurobi(model, params)
  
  mu = result$x[mu_idx]
  
  #plug mu back into mean variance to see what's the optimal variance
  w_mean_var = getVanillaSolution(cov_mat, mu, mean(mu))
  if(is.numeric(w_mean_var)){
    s2 = t(w_mean_var) %*% cov_mat %*% w_mean_var
  }else{
    s2 = NA
  }
  return(c(s2))
}

getInSampleS2 = function(mu_vec, cov_mat, log_tau_vec){
  n_stocks = length(mu_vec)
  if(length(log_tau_vec) == 1){
    log_tau_vec = rep(log_tau_vec, n_stocks)
  }
  
#   # type 1: 1/N portfolio's variance
  s2 = sum(cov_mat)/n_stocks^2
  
#   # type 2: mean variance portfolio's in sample variance with rho = mean(mu_vec)
#   rho = mean(mu_vec)
#   w_mean_var = getVanillaSolution(cov_mat, mu_vec, rho)
#   s2 = w_mean_var %*% cov_mat %*% w_mean_var
  
#   # type 3: in-sample max Sharpe portfolio's variance
#   w_max_Sharpe = solve(cov_mat, mu_vec)/sum(solve(cov_mat, mu_vec))
#   s2 = w_max_Sharpe %*% cov_mat %*% w_max_Sharpe
  
#   # type 4: mean variance portfolio's in sample variance with rho = mean(mu_vec) + mean(10^log_tau_vec) and mu = mu_vec
#   mu_max = mu_vec + 10^log_tau_vec
#   rho = mean(mu_max)
#   w_mean_var = getVanillaSolution(cov_mat, mu_vec, rho)
#   s2 = w_mean_var %*% cov_mat %*% w_mean_var
#   print(paste("mu_vec, rho_max", round(s2,2)))
    
#   # type 5: mean variance portfolio's in sample variance with rho = mean(mu_vec) + mean(10^log_tau_vec) and mu = mu_vec + 10^log_tau_vec
#   mu_max = mu_vec + 10^log_tau_vec
#   rho = mean(mu_max)
#   w_mean_var = getVanillaSolution(cov_mat, mu_max, rho)
#   s2 = w_mean_var %*% cov_mat %*% w_mean_var
#   print(paste("mu_max, rho_max", round(s2,2))) 
  
  # type 6: the maximum of mean variance portfolio's in sample variance with mu_vec in the uncertainty set
#   set.seed(070298)
#   N = 1000
#   s2_Sim = rep(NA, N)
#   cl <- makeCluster(8)
#   registerDoParallel(cl)
#   
#   result = foreach(k = 1:N,.options.RNG = 070298,.packages = pkg_names)%dorng%{
#     source("code/funcPar.R")
#     mu_Sim= runif(n_stocks, min = mu_vec - 10^log_tau_vec, max = mu_vec + 10^log_tau_vec)
#     rho = mean(mu_Sim)
#     w_mean_var = getVanillaSolution(cov_mat, mu_Sim, rho)
#     s2 = w_mean_var %*% cov_mat %*% w_mean_var
#   }
#   
#   stopCluster(cl)
#   
#   s2_Sim = unlist(result)
#   
# #   for(i in 1:N){
# #     mu_Sim= runif(n_stocks, min = mu_vec - 10^log_tau_vec, max = mu_vec + 10^log_tau_vec)
# #     rho = mean(mu_Sim)
# #     w_mean_var = getVanillaSolution(cov_mat, mu_Sim, rho)
# #     s2_Sim[i] = w_mean_var %*% cov_mat %*% w_mean_var
# #   }
#   
#   s2 = max(s2_Sim)
#   print(paste("simulation", s2))
  
  # type 7: use optim to solve for the maximum
  # optim thinks this function is discontinuous
  
  # type 8: the largest difference between ret_min_var and 1/N
  # mu_sign = mu_0 + sign(w_min_var - 1/N) * tau
#   w0 = rep(1/n_stocks, n_stocks)
#   w_min_var = rowMeans(solve(cov_mat))/sum(solve(cov_mat))
#   mu_sign = mu_vec + sign(w_min_var - w0) * 10^log_tau_vec
#   rho = mean(mu_sign)
#   w_mean_var = getVanillaSolution(cov_mat, mu_sign, rho)
#   s2 = w_mean_var %*% cov_mat %*% w_mean_var
#   print(paste("mu_sign, rho", round(s2,2))) 
  
  return(as.numeric(s2 + 1e-4))
}


getS2seq = function(rho_type, data, Para, use_type, mu_type, cov_type){
  cov_TS = getCovTS(data, Para, use_type, mu_type, cov_type)
  s2_TS_min_var = rep(NA, Para$n_reform)
  
  for(i in 1:Para$n_reform){
    s2_TS_min_var[i] = 1/sum(solve(cov_TS[,,i]))
  }
  
  s2_max_Sharpe_all = (Para$mu_vec %*% solve(Para$cov) %*% t(Para$mu_vec))/(sum(solve(Para$cov) %*% t(Para$mu_vec)))^2
  
  if(rho_type == "oneByN"){
    return(s2_seq = "oneByN")
  }else if(rho_type == "constant"){
    if(Para$n_stocks == 10){
      s2_seq = (seq(4, 5, by = 0.25))^2
    }else{
      s2_seq = (seq(2, 5, by = 0.25))^2
    }
    
  }

  return(s2_seq)
}




# get tau ---------------------------------------------------------
getLogTauVecDlasso = function(ret, cov_mat, Para, tau_type, prob, mult = NULL){
  log_tau_vec = NULL
  
  df_train = Para$train_length - 1
  df_next = Para$period_length - 1
  
  prob_one = prob^(1/Para$period_length)
  t_value = abs(qt( (1 - prob_one)/2 , df_train))
  
  mean_vec = apply(ret, 2, mean)
  var_vec = diag(cov_mat)
  max_vec = apply(ret, 2, max)
  min_vec = apply(ret, 2, min)
  
  deviation_mat = rbind(abs(max_vec - mean_vec), abs(min_vec - mean_vec))
  max_deviation = max(deviation_mat)
  max_deviation_vec = apply(deviation_mat, 2, max)
  
  if(tau_type == "LInf all t"){
    log_tau_vec = log(mult * sqrt(max(var_vec)/(df_train + 1)), 10)
  }else if(tau_type == "LInf indi t"){
    log_tau_vec = log( mult * sqrt(var_vec/(df_train + 1)), 10)
  }else if(tau_type == "LInf all NP"){
    log_tau_vec = log(max_deviation, 10)
  }else if(tau_type == "LInf indi NP"){
    log_tau_vec = log(max_deviation_vec, 10)
  }else if(tau_type == "sd mu_vec"){
    log_tau_vec = log(mult * sd(mean_vec), 10)
  }else if(tau_type == "sd demean"){
    demean_ret = ret - rowMeans(ret)
    log_tau_vec = log(mult/sqrt(df_train + 1) * apply(demean_ret, 2, sd), 10)
  }else if(tau_type == "rotate"){
    log_tau_vec = log(mult/sqrt(Para$train_length), 10)
  }else if(tau_type == "eigen"){
    eigen_result = eigen(cov_mat/Para$train_length)
    log_tau_vec = log(mult * sqrt(eigen_result$values), 10)
  }else if(tau_type == "eigen unit"){
    log_tau_vec = log(mult, 10)
  }
  
  return(log_tau_vec)
}

getTauWellAndPoor = function(df, N, poor_space, draw_type){
  if(draw_type == "boot"){
    draw_type = "Boot"
  }else if(draw_type == "wish"){
    draw_type = "Wish"
  }
  
  getDiffEigenValVec = match.fun(paste0("getDiffEigenValVec", draw_type))
  
  diff_val_mod =  getDiffEigenValVec(df, N, poor_space)
  
  tau = max(diff_val_mod)
  
  return(tau)
}


getDiffEigenValVecWish = function(df, N, poor_space){
  cov_true = cov(df)
  cov_true = t(poor_space) %*% cov_true %*% poor_space
  
  p = nrow(df)
  
  n_stocks = ncol(cov_true)
  
  type = getWishartType(cov_true, p, N)
  cov_true = cov_true + 10^(-10) * diag(n_stocks)
  
  result = foreach(k = 1:N,.options.RNG = 070298,.packages = pkg_names)%dorng%{
    cov_est = rWishartAllCases(p, cov_true, type)/p
    
    diff = cov_true - cov_est
    eigen_diff = eigen(diff)
    eigen_values = eigen_diff$values
  }
  
  result = array(unlist(result), dim = c(n_stocks,N))
  
  diff = apply(result, 1, mean)
  
  return(diff)
}



getDiffEigenValIWishart = function(N, cov_true, cov_est, pct_vec, Para, daily){
  if(daily == FALSE){
    p = Para$train_length - 1 
  }else{
    if(is.null(Para$train_length_daily)){
      p = Para$train_length * 21 - 1 
    }else{
      p = Para$train_length_daily - 1
    }
  }
  
  diff_mat = matrix(NA, nrow = Para$n_stocks, ncol = N)
  
  cov_true = cov_true + 10^(-10) * diag(n_stocks)
  
  cov_est = rWishartAllCases(p, cov_true, type)/p
  for(i in 1:N){
    cov_IW = solve(rwish(p,solve(cov_est*p)))
    
    diff = cov_IW - cov_est
    eigen_diff = eigen(diff)
    diff_mat[,i] = eigen_diff$values
  }
  
  diff = apply(diff_mat, 1, mean)
  
  diff_qtl = quantile(diff, pct_vec, na.rm = TRUE)
  
  return(diff_qtl)
}

getCorTrueAndEstEigenVector = function(cov_true, p, N){
  result = getCorTrueAndEstEigenVectorAllSims(cov_true, p, N)
  
  cor_mean = rowMeans(result)
  
  return(cor_mean)
}

getCorTrueAndEstEigenVectorAllSims = function(cov_true, p, N){
  n_stocks = ncol(cov_true)
  
  type = getWishartType(cov_true, p, N)
  
  cov_true = cov_true + 10^(-10) * diag(n_stocks)
  eigen_true = eigen(cov_true)
  
  result = foreach(k = 1:N,.options.RNG = 070298,.packages = pkg_names)%dorng%{
    cov_est = rWishartAllCases(p, cov_true, type)/p
    
    eigen_est = eigen(cov_est)
    
    cor_eigen_vectors = cor(eigen_true$vectors, eigen_est$vectors)
    
    cor_vec = abs(diag(cor_eigen_vectors))
  }
  
  result = array(unlist(result),dim = c(n_stocks,N))
  
  return(result)
}

getP = function(Para, daily){
  if(daily == FALSE){
    p = Para$train_length - 1 
  }else{
    if(is.null(Para$train_length_daily)){
      p = Para$train_length * 21 - 1 
    }else{
      p = Para$train_length_daily - 1
    }
  }
  return(p)
}

getPoorSpace = function(cov_true, cor_true, alpha){
  eigen_true = eigen(cov_true)
  
  poor_idx = which(cor_true < alpha)
  
  if(length(poor_idx) > 0){
    poor_space = eigen_true$vectors[,poor_idx]
  }else{
    poor_space = diag(ncol(cov_true))
  }
  
  return(poor_space)
}

getLogTauForFixedK = function(K, log_tau_u, rho, mu_vec, cov_mat, sum_w, tol){
  log_tau_l = -5
  log_tau_diff = log_tau_u - log_tau_l
  tau_diff = 10^log_tau_u - 10^log_tau_l
  n_non_zero_l = getNonZeroClasso(log_tau_l, rho, mu_vec, cov_mat, tol)
  n_non_zero_u = getNonZeroClasso(log_tau_u, rho, mu_vec, cov_mat, tol)
  
  while((K > n_non_zero_u | K < n_non_zero_l) & (log_tau_diff > 10^-2 | tau_diff > 0.01)){
    log_tau_m = mean(c(log_tau_l,log_tau_u))
    n_non_zero_m = getNonZeroClasso(log_tau_m, rho, mu_vec, cov_mat, tol)
#     display = matrix(c(log_tau_l, log_tau_m, log_tau_u,
#                        n_non_zero_l, n_non_zero_m, n_non_zero_u),
#                      nrow = 2,byrow = TRUE)
#     print(display)
    if(n_non_zero_m > K){
      log_tau_l = log_tau_m
      n_non_zero_l = n_non_zero_m
    }else if(n_non_zero_m < K){
      log_tau_u = log_tau_m
      n_non_zero_u = n_non_zero_m
    }else{
      log_tau_l = mean(c(log_tau_l,log_tau_m))
      n_non_zero_l = getNonZeroClasso(log_tau_l, rho, mu_vec, cov_mat, tol)
      log_tau_u = mean(c(log_tau_u,log_tau_m))
      n_non_zero_u = getNonZeroClasso(log_tau_u, rho, mu_vec, cov_mat, tol)
    }
    log_tau_diff = log_tau_u - log_tau_l
    tau_diff = 10^log_tau_u - 10^log_tau_l
  }
  
  return(log_tau_u)
}

# additional functions ---------------------------------------------------

getCI = function(Mean, SD, mult){ #CI confidence interval
  return(c(Mean + mult * SD, Mean - mult *SD))
}


retSummary = function(ret, Para){
  ret_annual = matrix(ret[-(1:Para$train_length)], nrow = 12)
  
  mu_annual = apply(ret_annual, MARGIN = 2, FUN = mean)
  SD_annual = apply(ret_annual, MARGIN = 2, FUN = sd)
  Sharpe_annual = apply(ret_annual, MARGIN = 2, FUN = sharpeRatio)
  
  return(cbind(mu_annual, SD_annual, Sharpe_annual))
}

getMonthDataFromDaily = function(data_daily){
  dates_daily = rownames(data_daily)
  dates = unique(substr(dates_daily, 1, 6))
  n_time = length(dates)
  n_stocks = ncol(data_daily)
  data = as.data.frame(matrix(NA, n_time, n_stocks))
  rownames(data) = dates
  for(i in 1:n_time){
    start_idx = min(grep(dates[i], dates_daily))
    end_idx = max(grep(dates[i], dates_daily))
    raw_ret = 1 + data_daily[start_idx:end_idx,, drop = FALSE]/100
    data[i,] = (apply(raw_ret, MARGIN = 2, prod) - 1) * 100
    # data[i,] = colSums(data_daily[start_idx:end_idx,, drop = FALSE])
  }
  return(data)
}

getAnualData = function(data){
  dates = rownames(data)
  dates_annual = unique(substr(dates, 1, 4))
  n_time = length(dates_annual)
  n_stocks = ncol(data)
  data_annual = as.data.frame(matrix(NA, n_time, n_stocks))
  rownames(data_annual) = dates_annual
  for(i in 1:n_time){
    start_idx = min(grep(dates_annual[i], dates))
    end_idx = max(grep(dates_annual[i], dates))
    raw_ret = 1 + data[startge_idx:end_idx,, drop = FALSE]/100
    data_annual[i,] = (apply(raw_ret, MARGIN = 2, prod) - 1) * 100
    # data_annual[i,] = colSums(data[start_idx:end_idx,, drop = FALSE])
  }
  return(data_annual)
}

getSubsetData = function(data, start_date, end_date, train_length){
  dates = rownames(data)
  start_idx = min(grep(pattern = start_date, x = dates)) - train_length
  end_idx = max(grep(pattern = end_date, x = dates)) 
  if(is.infinite(start_idx)){
    start_idx = 1
  }
  
  if(is.infinite(end_idx)){
    end_idx = nrow(data)
  }
  return(data[start_idx:end_idx, , drop = FALSE])
}

getTrainInfo = function(data, Para, t, use_type, mu_type, cov_type, daily = FALSE, data_daily = data_daily, pred = FALSE){
  if(pred == FALSE){
    if(daily == FALSE){
      use_idx = switch(use_type,
                       prev_one = (t - Para$train_length + 1):t,
                       prev_all = 1:t)
      
      ret = as.matrix(data[use_idx,])
      
      mu_vec = switch(mu_type,
                      true = Para$mu_vec,
                      train = colMeans(ret,na.rm = TRUE))
      
      cov_mat = switch(cov_type,
                       true = Para$cov,
                       train = cov(ret, use = "complete") + 10^(-10) * diag(Para$n_stocks))
    }else{
      use_idx = switch(use_type,
                       prev_one = (t - Para$train_length_daily + 1):t,
                       prev_all = 1:t)
      
      ret = as.matrix(data_daily[use_idx,])
      
      mu_vec = switch(mu_type,
                      true = Para$mu_vec,
                      train = colMeans(ret,na.rm = TRUE) * 252/12)
      
      cov_mat = switch(cov_type,
                       true = Para$cov,
                       train = (cov(ret, use = "complete") +10^(-10) * diag(Para$n_stocks)) * 252/12)
      
    }
    
    result = list(ret = ret, mu_vec = mu_vec, cov_mat = cov_mat, use_idx = use_idx)
    
  }else{
    data = data - pred_month_ret
    result = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily = daily, data_daily = data_daily, pred = FALSE)
  }
 
  return(result)
}


getBestSharpeInfo = function(methods, Sharpe){
  n_methods = length(methods)
  
  best_Sharpe = matrix(NA, nrow = 2, ncol =  n_methods)
  rownames(best_Sharpe) = c("Parameter", "Index")
  colnames(best_Sharpe) = methods

  for(i in 1:n_methods){
    method_i = methods[i]
    Sharpe_i = eval(parse(text = paste0("Sharpe$F_", method_i)))
    
    idx = which.max(apply(Sharpe_i, MARGIN = 2, FUN = max, na.rm = TRUE))
    best_Sharpe[1,i] = names(idx)
    best_Sharpe[2,i] = idx
  }
  
  return(best_Sharpe)
}

corTwoSim = function(ret_mat_1, ret_mat_2){
  rho_seq = dimnames(ret_mat_1)$rho
  n_rho = length(rho_seq)
  log_tau_1 = dimnames(ret_mat_1)$seq
  log_tau_2 = dimnames(ret_mat_2)$seq
  # cov_ret = array(NA, dim = c(length(log_tau_1), length(log_tau_2), n_rho), dimnames = list(log_tau_1, log_tau_2, rho_seq))
  cor_ret = array(NA, dim = c(length(log_tau_1), length(log_tau_2), n_rho), dimnames = list(log_tau_1, log_tau_2, rho_seq))
  
  n_dim = length(dim(ret_mat_1))
  if(n_dim == 4){
    N = dim(ret_mat_1)[4]
  }else{
    # this will make dimnames disappear
    dim(ret_mat_1) = c(dim(ret_mat_1),1)
    dim(ret_mat_2) = c(dim(ret_mat_2),1)
    N = 1
  }
  
  for(i in seq_along(log_tau_1)){
    for(j in seq_along(log_tau_2)){
      for(k in seq_along(rho_seq)){
        # cov_vec = rep(NA, N)
        cor_vec = rep(NA, N)
        for(l in 1:N){
          # cov_vec[l] = cov(ret_mat_1[,k,i,l], ret_mat_2[,k,j,l], use = "complete")
          cor_vec[l] = cor(ret_mat_1[,k,i,l], ret_mat_2[,k,j,l], use = "complete")
        }
        # cov_ret[i,j,k] = mean(cov_vec)
        cor_ret[i,j,k] = mean(cor_vec)
      }
    }
  }
  # return(list(cov_mat = cov_ret, cor = cor_ret))
  return(cor_ret)
}


solveRiccati = function(A, S, Q){
  # We need equation of the form
  # A^t X + X A − X S X + Q = 0.
    
  n = ncol(Q)
  H = matrix(0, nrow = 2*n, ncol = 2*n)
  H[1:n, 1:n] = A
  H[1:n, (n + 1):(2*n)] = -S
  H[(n + 1):(2*n), 1:n] = -Q
  H[(n + 1):(2*n), (n + 1):(2*n)] = -t(A)
  
  eigen_result = eigen(H)
  
  idx = which(Re(eigen_result$values) < 0)
  
  P = eigen_result$vectors[,idx]
  
  P11 = P[1:n,]
  P21 = P[(n + 1):(2*n),]
  
  P = P21 %*% solve(P11)
  
  return(P)
}

getMinVarSol = function(inv_cov){
  w_MinVar = rowSums(inv_cov)/sum(rowSums(inv_cov))
  return(w_MinVar)
}

getMeanVarSol = function(inv_cov, mu_vec){
  mu_vec = matrix(mu_vec, ncol = 1)
  w_MeanVar = inv_cov %*% mu_vec/ sum(inv_cov %*% mu_vec)
  return(w_MeanVar)
}

vecNorm = function(vec,p){
  return(sum(abs(vec)^p,na.rm = TRUE)^(1/p))
}

sharpeRatio = function(vec){
  vec = vec #- RF[,1]
  sharpe_ratio = mean(vec,na.rm = TRUE)/sd(vec,na.rm = TRUE)
  return(sharpe_ratio)
}

annualySharpe = function(vec){
  annual_sharpe = apply(matrix(vec, 12), MARGIN = 2, FUN = sharpeRatio)
  return(annual_sharpe)
}

getNextPeriodRet = function(data,next_period_idx,weights){
  next_period_ret_mat = data[next_period_idx,,drop = FALSE]
  cum_ret = apply(next_period_ret_mat,MARGIN = 2,cumsum)
  next_period_cum_ret = cum_ret %*% weights
  next_period_ret = diff(c(0,next_period_cum_ret))
  return(next_period_ret)
  
  # It turns out Bordie's paper use log return
  # The following commented out code use raw return which is useless here
  # cum_ret = apply(1+next_period_ret_mat/100,MARGIN = 2,cumprod)
  # oneByN_cum_ret = rowMeans(cum_ret)
  # oneByN_cum_ret_lag = c(1,oneByN_cum_ret[1:(period_length-1)])
  # oneByN_month_ret[next_period_idx] = (oneByN_cum_ret/oneByN_cum_ret_lag - 1)*100
  
}


getNonZero = function(vec,tol){
  return(sum(abs(vec)>tol))
}

cosXandY = function(x,y){
  txy = t(x) %*% y
  norm_x = vecNorm(x,2)
  norm_y = vecNorm(y,2)
  cos_xy = txy/norm_x/norm_y
  return(cos_xy)
}

getDir = function(data_name, name, file_type,Sim = NULL, cov_type = NULL, rho_type = NULL, load_daily  = FALSE, train_length = NULL, period_length = NULL, use_type = NULL, author = "Brodie", draw_type = "boot"){
  if(load_daily){
    file_name = paste0("Result/", author, "/", data_name, Sim, cov_type, rho_type,
                       " ", name, " ", use_type, " ",train_length, " ", period_length, draw_type,
                       "_daily", ".",file_type)
  }else{
    file_name = paste0("Result/", author, "/",data_name, Sim, cov_type, rho_type,
                       " ", name, " ",use_type, " ",train_length, " ", period_length, draw_type, 
                       ".", file_type)
  }
  
  return(file_name)
}



is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
 return(abs(x - round(x)) < tol) 
}

getNonZeroClasso = function(log_tau_vec, rho, mu_vec, cov_mat, tol){
  w = getCLassoInEqSolutionCone(log_tau_vec, rho, mu_vec, cov_mat, sum_w = "one")
  n_non_zero = getNonZero(w, tol)
  return(n_non_zero)
}


# ProbUnExpVar ------------------------------------------------------------
# getGammaRatioWish = function(cov_true, N, p){
#   n_stocks = ncol(cov_true)
#   cov_true = cov_true + 10^(-10) * diag(n_stocks)
#   
#   type = getWishartType(cov_true, p, N)
#   eigen_true = eigen(cov_true)
#   
#   result = foreach(k = 1:N,.options.RNG = 070298,.packages = pkg_names)%dorng%{
#     cov_est = rWishartAllCases(p, cov_true, type)/p
#     
#     eigen_est = eigen(cov_est)
#     
#     var_real = diag(t(eigen_est$vectors) %*% cov_true %*% eigen_est$vectors)
#     
#     ratio_vec = var_real/eigen_est$values
#   }
#   
#   result = array(unlist(result), dim = c(n_stocks,N))
#   
#   return(result)
# }

getGammaRatioWish = function(cov_true, N, p){
  n_stocks = ncol(cov_true)
  type = getWishartType(cov_true, p, N)
  
  cov_true = cov_true + 10^(-10) * diag(n_stocks)
  eigen_true = eigen(cov_true)
  
  result = foreach(k = 1:N,.options.RNG = 070298,.packages = pkg_names)%dorng%{
    cov_est = rWishartAllCases(p, cov_true, type)/p
    
    eigen_est = eigen(cov_est)
    
    var_real = diag(t(eigen_est$vectors) %*% cov_true %*% eigen_est$vectors)
    
    if(type == "degenerate"){
      ratio_vec = rep(Inf, n_stocks)
      ratio_vec[1:p] = var_real[1:p]/eigen_true$values[1:p]
    }else if(type == "invertible"){
      ratio_vec = var_real/eigen_true$values
    }
    
  }
  
  result = array(unlist(result), dim = c(n_stocks,N))
  
  return(result)
}

getTopKProbFromIndicatorMat = function(indicator_mat){
  N = ncol(indicator_mat)
  
  top_k = apply(indicator_mat, MARGIN = 2, FUN = cumprod)
  freq = rowSums(top_k)
  prob = (c(N, freq) - c(freq, 0))/N
  
  return(prob)
}

# get mod_matrix for poor_space from all eigenvectors to none
# getAllModMat = function(df, N, draw_type){
#   cov_true = cov(df)
#   n_stocks = ncol(cov_true)
#   eigen_true = eigen(cov_true)
#   
#   tau_vec = rep(0, n_stocks + 1)
#   mod_array = array(0, dim = c(n_stocks, n_stocks, n_stocks + 1))
#   
#   for(i in 1:(n_stocks - 1)){
#     poor_idx = i : n_stocks
#     poor_space = eigen_true$vectors[,poor_idx]
#     tau_vec[i]  = getTauWellAndPoor(df, N, poor_space, draw_type)
#   }
#   
#   for(i in 1:(n_stocks - 1)){
#     poor_idx = i : n_stocks
#     poor_space = eigen_true$vectors[,poor_idx]
#     mod_array[,,i] = tau_vec[i] * poor_space %*% t(poor_space)
#   }
#   
#   return(result = list(tau_vec = tau_vec, mod_array = mod_array))
# }

getAllModMat = function(df, N, draw_type){
  cov_true = cov(df)
  n_stocks = ncol(cov_true)
  eigen_true = eigen(cov_true)
  
  tau_vec = rep(0, n_stocks + 1)
  mod_array = array(0, dim = c(n_stocks, n_stocks, n_stocks + 1))
  
  result = foreach(i = 1:(n_stocks - 1),.options.RNG = 070298,.packages = pkg_names)%dorng%{
    poor_idx = i : n_stocks
    poor_space = eigen_true$vectors[,poor_idx]
    getTauWellAndPoor(df, N, poor_space, draw_type)
  }
  
  tau_vec[1:(n_stocks - 1)] = c(unlist(result))
  
  for(i in 1:(n_stocks - 1)){
    poor_idx = i : n_stocks
    poor_space = eigen_true$vectors[,poor_idx]
    mod_array[,,i] = tau_vec[i] * poor_space %*% t(poor_space)
  }

  return(result = list(tau_vec = tau_vec, mod_array = mod_array))
}

getMeanModMat = function(mod_array, prob){
  n_stocks = length(prob) - 1
  
  mod_mat = matrix(0, n_stocks, n_stocks)
  for(i in 1:(n_stocks + 1)){
    mod_mat = mod_mat + mod_array[,,i] * prob[i]
  }
  
  return(mod_mat)
}

#Probabilistic unexpected variance 
getTauVecGR = function(N, eps_vec, all_mod_mat_TS, all_tau_TS, Para, data, data_daily, daily = FALSE, top_K = TRUE){
  n_eps = length(eps_vec)
  
  tau_mat = matrix(NA, Para$n_reform, n_eps)
  colnames(tau_mat) = eps_vec
  
  mod_array = array(NA, dim = c(Para$n_reform, n_eps, Para$n_stocks, Para$n_stocks), dimnames = c(NULL, eps_vec = eps_vec, NULL, NULL))
  
  for(i in 1:Para$n_reform){
    if(daily == FALSE){
      t = Para$reform_idx[i]
    }else{
      t = Para$reform_idx_daily[i]
    }
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily = daily, data_daily = data_daily)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    eigen_cov = eigen(cov_mat)
    
    p = nrow(ret) - 1
  
    ratio_est = getGammaRatioWish(cov_mat, N, p)
    
    for(j in seq_along(eps_vec)){
      eps = eps_vec[j]
      
      indicator_mat = (abs(ratio_est - 1) < eps)
      
      if(top_K == TRUE){
        prob_est = getTopKProbFromIndicatorMat(indicator_mat)
        
        tau_mat[i,j] = sum(all_tau_TS[i,] * prob_est)
        mod_array[i,j,,] = getMeanModMat(all_mod_mat_TS[i,,,], prob_est)
        
        test1 = getMeanModMat(all_mod_mat_TS[i,,,], prob_est)
      }else{
        prob_mat = getProbMatFromIndicatorMat(indicator_mat)
        
        mod_values = c(all_tau_TS[i,] %*% prob_mat)
        
        mod_array[i,j,,] = eigen_cov$vectors %*% diag(mod_values) %*% t(eigen_cov$vectors)
        
        test2 = eigen_cov$vectors %*% diag(mod_values) %*% t(eigen_cov$vectors)
      }
      
    }
  }
  
  result = list(tau_mat = tau_mat, mod_array = mod_array)
  
  return(result)
}

getAllModMatTS = function(data, Para, N, data_daily, daily, draw_type = "boot"){
  all_mod_mat_TS = array(NA, dim = c(Para$n_reform, Para$n_stocks, Para$n_stocks, Para$n_stocks + 1))
  all_tau_TS = matrix(NA, nrow = Para$n_reform, ncol = Para$n_stocks + 1)
  
  for(i in 1:Para$n_reform){
    if(daily == FALSE){
      t = Para$reform_idx[i]
    }else{
      t = Para$reform_idx_daily[i]
    }
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily = daily, data_daily = data_daily)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    p = nrow(ret) - 1
    
    result = getAllModMat(ret, N, draw_type)
    all_tau_TS[i,] = result$tau_vec
    all_mod_mat_TS[i,,,] = result$mod_array
  }
  
  result = list(all_tau_TS = all_tau_TS, all_mod_mat_TS = all_mod_mat_TS)
  
  return(result)
}

getAvgTargetRet = function(data, Para, use_type, mu_type, cov_type, daily = FALSE, data_daily = data_daily){
  avg_target_ret = matrix(NA, Para$n_reform, 2)
  
  for(i in 1:Para$n_reform){
    if(daily == FALSE){
      t = Para$reform_idx[i]
    }else{
      t = Para$reform_idx_daily[i]
    }
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily = daily, data_daily = data_daily)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    avg_target_ret[i,1] = mean(mu_vec, na.rm = TRUE)
    avg_target_ret[i,2] = max(mu_vec, na.rm = TRUE)
  }
  return(avg_target_ret)
  
}

getMonthResult = function(ret, FUN){
  FUN = match.fun(FUN)
  ret = matrix(ret, 12)
  vec = apply(ret, 2, FUN)
  return(vec)
}

getPosNegCov = function(data, type){
  if(type == "all"){
    cov_mat = cov(as.matrix(data), use = "complete")
  }else{
    n_time = nrow(data)
    mean_ret = rowMeans(data, na.rm = TRUE)
    rank_ret = rank(mean_ret)/n_time
    idx = switch(type, 
                 pos = which(mean_ret >= 0),
                 neg = which(mean_ret <= 0),
                 bottom20 = which(rank_ret <= 0.2),
                 bottom30 = which(rank_ret <= 0.3),
                 top20 = which(rank_ret >= 0.8),
                 top30 = which(rank_ret >= 0.3))
    cov_mat = cov(as.matrix(data[idx,]), use = "complete")
  }
  return(cov_mat)
}

getCumSharpe = function(ret){
  n_time = length(ret)
  first_non_NA = min(which(!is.na(ret)))
  
  cumSharpe = rep(NA, n_time)
  for(i in (first_non_NA + 11):n_time){
    idx = 1:i
    cumSharpe[i] = sharpeRatio(ret[idx])
  }
  return(cumSharpe)
}

printRound2 = function(vec){
  print(round(vec,2))
}

getIndexMat = function(mat){
  index_mat = matrix(NA, nrow = 2, ncol = length(mat))
  rownames(index_mat) = c("row", "col")
  index_mat[1,] = rep(1:nrow(mat), ncol(mat))
  index_mat[2,] = rep(1:ncol(mat), each = nrow(mat))
  return(index_mat)
}

getRho = function(rho_vec, ret){
  if(is.null(rho_vec)){
    rho = mean(ret,na.rm = TRUE)
  }else{
    rho = rho_vec[i]
  }
  
  return(rho)
}


rWishartAllCases = function(p, cov_mat, type = "invertible"){
  if(type == "degenerate"){
    rMat = rWishartDegenerate(p, cov_mat)
  }else if(type == "invertible"){
    rMat = rwish(p, cov_mat)
  }
  return(rMat)
}

rWishartDegenerate = function(p, cov_mat){
  n = ncol(cov_mat)
  X = rmvnorm(p, mean = rep(0,n) ,sigma = cov_mat)
  rMat = t(X) %*% X
  return(rMat)
}

getWishartType = function(cov_mat, p, N){
  n = ncol(cov_mat)
  rank_mat = rankMatrix(cov_mat)
  if(rank_mat < n | p < n){
    type = "degenerate"
  }else{
    type = "invertible"
  }
  return(type)
}


getProbMatFromIndicatorMat = function(indicator_mat){
  N = ncol(indicator_mat)
  n_stocks = nrow(indicator_mat)
  
  prob_mat = matrix(0, nrow = n_stocks + 1, ncol = n_stocks)
  
  for(i in 1:N){
    FALSE_idx = which(!indicator_mat[,i])
    min_FALSE_idx = min(FALSE_idx)
    prob_mat[min_FALSE_idx + 1, FALSE_idx] = prob_mat[min_FALSE_idx + 1, FALSE_idx] + 1/N
  }
  
  return(prob_mat)
}

getSubsetPCASolution = function(mu_vec, eigen_cov, sub_idx){
  one_vec = rep(1, ncol(eigen_cov$vectors))
  
  inner_prod_MV = c(one_vec %*% eigen_cov$vectors)
  inner_prod_MS = c(mu_vec %*% eigen_cov$vectors)
  eigen_values = eigen_cov$values
  eigen_values[-sub_idx] = Inf
  
  w_MV = eigen_cov$vectors %*% (inner_prod_MV/eigen_values)
  w_MV = w_MV/sum(w_MV)
  
  w_MS = eigen_cov$vectors %*% (inner_prod_MS/eigen_values)
  w_MS = w_MS/sum(w_MS)
  
  return(list(MV = w_MV, MS = w_MS))
}

gammaRatioStat = function(df, idx, cov_true, eigen_true){
  cov_est = cov(df[idx,], use = "complete")
  eigen_est = eigen(cov_est)
  gamma_ratio = diag(t(eigen_est$vectors) %*% cov_true %*% eigen_est$vectors) / eigen_true$values
  if(length(idx) < ncol(df)){
    gamma_ratio[(length(idx) + 1) : ncol(df)] = Inf
  }
  return(gamma_ratio)
}

getGammaRatioBoot = function(df, N){
  cov_true = cov(df, use = "complete")
  eigen_true = eigen(cov_true)
  
  # it is not trivial to do parallel in boot function, so I comment it out
  boot_result = boot(df, statistic = gammaRatioStat, R = N, # parallel = "snow"
                     cov_true = cov_true, eigen_true = eigen_true)
  
  gamma_ratio = t(boot_result$t)
  
  return(gamma_ratio)
}



wSubsetPCAStat = function(df, idx, mu_true, cov_true, eigen_true, side){
  gamma_ratio = gammaRatioStat(df, idx, cov_true, eigen_true)
  
  if(side == 1){
    I_vec = (gamma_ratio <= 1 + eps)
  }else{
    I_vec = (abs(gamma_ratio - 1) <= eps)
  }
  
  sub_idx = which(I_vec)
  
  if(length(sub_idx) == 0){
    sub_idx = 1
  }
  
  w = getSubsetPCASolution(mu_true, eigen_true, sub_idx)
  
  return(c(w$MV, w$MS))
}


getGammaRatioSubsetPCARet = function(eps, N, data, Para, use_type, mu_type, cov_type, sum_w = "one", daily = FALSE, data_daily = data_daily, side = 2){
  month_ret_MV = rep(NA,Para$n_time)
  month_ret_MS = rep(NA,Para$n_time)
  
  for(i in 1:Para$n_reform){
    if(daily == FALSE){
      t = Para$reform_idx[i]
    }else{
      t = Para$reform_idx_daily[i]
    }
    
    train_info = getTrainInfo(data, Para, t, use_type, mu_type, cov_type, daily = daily, data_daily = data_daily)
    mu_vec = train_info$mu_vec
    cov_mat = train_info$cov_mat
    ret = train_info$ret
    
    # it is not trivial to do parallel in boot function, so I comment it out
    boot_result = boot(ret, wSubsetPCAStat, R = N, mu_true = mu_vec, cov_true = cov_mat, eigen_true = eigen(cov_mat), side = side) #, parallel = "snow"
    
    w = colMeans(boot_result$t)
    w_MV = w[1:Para$n_stocks]
    w_MS = w[-(1:Para$n_stocks)]
  
    t = Para$reform_idx[i]
    
    next_period_idx = (t + 1):(t + Para$period_length)
    
    month_ret_MV[next_period_idx] = getNextPeriodRet(data, next_period_idx, w_MV)
    month_ret_MS[next_period_idx] = getNextPeriodRet(data, next_period_idx, w_MS)
  }
  
  return(list(MV = month_ret_MV, MS = month_ret_MS))
}


getDiffEigenValVecBoot = function(df, N, poor_space){
  cov_true = cov(df)
  # it is not trivial to do parallel in boot function, so I comment it out
  boot_result = boot(df, diffEigenValVecStat, R = N, cov_true = cov_true, poor_space = poor_space) #, parallel = "snow"
  diff_sim_result = boot_result$t
  diff_vals = colMeans(diff_sim_result)
  
  return(diff_vals)
}  
  
diffEigenValVecStat =  function(df, idx, cov_true, poor_space){
  cov_est = cov(df[idx,])
  diff_mat = t(poor_space) %*% (cov_est - cov_true) %*% poor_space
  diff_mat = diff_mat + 10^(-10) * diag(ncol(diff_mat))
  eigen_diff = eigen(diff_mat)
  return(Re(eigen_diff$values))
}

