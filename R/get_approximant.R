#library(tidyverse)
#library(data.table)
#library(Matrix)
#library(lpSolve)

# example use 
# # load the example data
# example <- readRDS("data/example1.R")
# get_approximant_p_infinity(P = example$P, M = example$M)

get_approximant <- function(P, M, norm = "infinity") {
  norm <- as.character(norm)
  stopifnot(norm %in% c("infinity", "1", "2"))
  if(norm == "infinity") return(get_approximant_p_1orinfinity(P = P, M = M, norm = norm))
  if(norm == "1") return(get_approximant_p_1orinfinity(P = P, M = M, norm = norm))
  if(norm == "2") stop()
}


get_approximant_p_1orinfinity <- function(P, M, norm = "infinity") {
  norm <- as.character(norm)
  stopifnot(norm %in% c("infinity", "1", "2"))
  p_is_infinity <- (norm == "infinity")
  m <- nrow(P)
  
  # create the linear program
  
  # variables arrange all in a long vector: 
  # Q_ij: the transition probability from state i to state j
  # xi_ij: slacks for the absolute difference between P_ij and Q_ij
  # l: the maximum column sum of Q
  n_vars <- 2*m^2+1
  Q_idx <- 1:m^2              # index assumes byrow = TRUE
  xi_idx <- (m^2+1):(2*m^2)   # index assumes byrow = TRUE
  l_idx <- 2*m^2+1
  
  
  # helper function identify idx from Q_ij indices in each row with ones 
  idx_from_ij <- function(i, j, cols)  { (i-1)*cols + j }
  i_from_idx <- function(idx, cols) { floor((idx-1) / cols) + 1  }
  j_from_idx <- function(idx, cols) { (idx-1)%%cols + 1}
  unroll_mat <- function(mat) as.vector(t(mat))
  
  
  # objective function: involves only l 
  objective.fn <- Matrix::sparseMatrix(x=1, i = 1, j = l_idx, dims=c(1,n_vars))
  
  # constraints:
  # 0.   abs(P_ij - Q_ij) <= xi_ij
  # 0a.  P_ij - Q_ij <= xi_ij  -->  Q_ij + xi_ij >= P_ij
  coeffs.0a <- Matrix::sparseMatrix(
    x = rep(1, length(Q_idx) + length(xi_idx)),
    i = c(1:length(Q_idx), 1:length(Q_idx)), 
    j = c(idx_from_ij(i=rep(1:m, each = m), j = rep(1:m, m), cols = m),
          length(Q_idx) + idx_from_ij(i=rep(1:m, each = m), j = rep(1:m, m), cols = m)),
    dims = c(length(Q_idx), n_vars)
  )
  rhs.0a = unroll_mat(P) # row major unroll 
  sense.0a = rep(">=", length(Q_idx))
  
  # 0b. -P_ij + Q_ij <= xi_ij --> -Q_ij + x_ij  >= -P_ij
  coeffs.0b <- Matrix::sparseMatrix(
    x = rep(c(-1, 1), each = length(Q_idx)),
    i = c(1:length(Q_idx), 1:length(Q_idx)), 
    j = c(idx_from_ij(i=rep(1:m, each = m), j = rep(1:m, m), cols = m),
          length(Q_idx) + idx_from_ij(i=rep(1:m, each = m), j = rep(1:m, m), cols = m)),
    dims = c(length(Q_idx), n_vars)
  )
  rhs.0b = unroll_mat(-P) # row major unroll 
  sense.0b = rep(">=", length(Q_idx))
  
  # 1. Q is a rows sum to 1 
  coeffs.1 <- Matrix::sparseMatrix(
    x = rep(1, length(Q_idx)),
    i = rep(1:m, each = m), 
    j = idx_from_ij(i=rep(1:m, each = m), j = rep(1:m, m), cols = m), 
    dims = c(m, n_vars)
  )
  rhs.1 = rep(1, m)
  sense.1 = rep("=", m)
  
  # 2. Q non negative 
  # also xi non negative 
  # also l non negative 
  coeffs.2 <- Matrix::sparseMatrix(
    x = rep(1, n_vars),
    i = 1:n_vars, 
    j = 1:n_vars, 
    dims = c(n_vars, n_vars)
  )
  rhs.2 <- rep(0, n_vars)
  sense.2 = rep(">=", n_vars)
  
  # 3. Q is a matrix with structural zeros Q <= M
  coeffs.3 <- Matrix::sparseMatrix(
      x = rep(1, length(Q_idx)),
      i = Q_idx, 
      j = Q_idx, 
      dims = c(length(Q_idx), n_vars)
    )
  rhs.3 <- unroll_mat(M)
  sense.3 = rep("<=", length(Q_idx))
  #browser()
  if(p_is_infinity) {
    # 4. all row sums for xi's (including their max) is less than l
    # OR -(row sum of xi's) +  l >= 0 
    coeffs.4 <- Matrix::sparseMatrix(
      x = rep(-1, length(xi_idx)),
      i = rep(1:m, each = m), 
      j = length(Q_idx) + idx_from_ij(i=rep(1:m, each = m), j = rep(1:m, m), cols = m), 
      dims = c(m, n_vars)
    )
    coeffs.4[1:m, l_idx] <- 1 
    rhs.4 = rep(0, m)
    sense.4 = rep(">=", m)  
  } else {
    # 4. all column sums for xi's (including their max) is less than l
    # OR -(col sum of xi's) +  l >= 0 
    coeffs.4 <- Matrix::sparseMatrix(
      x = rep(-1, length(xi_idx)),
      i = rep(1:m, each = m), 
      j = length(Q_idx) + idx_from_ij(i=rep(1:m, m), j = rep(1:m, each = m), cols = m), 
      dims = c(m, n_vars)
    )
    coeffs.4[1:m, l_idx] <- 1 
    rhs.4 = rep(0, m)
    sense.4 = rep(">=", m)  
  }
  
  
  # Bundle all coeffs, sense, and rhs together 
  coeffs <- rbind(coeffs.0a, coeffs.0b, coeffs.1, coeffs.2, coeffs.3, coeffs.4)
  sense <- c(sense.0a, sense.0b, sense.1, sense.2, sense.3, sense.4)
  rhs <- c(rhs.0a, rhs.0b, rhs.1, rhs.2, rhs.3, rhs.4)
  
  ## Solve: 
  sol <- lpSolve::lp(
    direction = "min", 
    objective.in = as.vector(objective.fn), 
    const.mat = as.matrix(coeffs), 
    const.dir = sense, 
    const.rhs = rhs
  )
  
  Q <- matrix(data = sol$solution[Q_idx], ncol = m, byrow = TRUE)
  colnames(Q) <- colnames(P)
  rownames(Q) <- rownames(P)
  return(list(LP= sol, Q = Q))
}
