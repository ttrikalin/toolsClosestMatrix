library(tidyverse)
library(data.table)
library(Matrix)
library(gurobi)


l1 <- 0.20
l2 <- 0.15

# System 
# S1 -> S2 -> S3

# generator matrix
A <- matrix(
  c(-l1 , l1, 0, 
    0, -l2, l2, 
    0, 0, 0), nrow = 3, byrow = TRUE)

# transition matrix
P <- expm(A)
N <- nrow(P)
unrolled_P <- as.vector(P)  # in gurobi all variables in a vector

unrolled_index <- function(i, j) (i - 1) * N + j # this is the index of the element in the vector

# a column of ones
ones <- matrix(rep(1, N), nrow = N)

# a column of zeros
zeros <- matrix(rep(0, N), nrow = N)





## Gurobi model



n <- nrow(P)
N <- n * n

# Flatten Q and P into vectors (row-major)
Pvec <- as.vector(t(P))
Wvec <- as.vector(t(W))

# Indices of variables that are allowed (i.e., W == 1)
active_idx <- which(Wvec == 1)

# Variable count
num_vars <- length(active_idx)

# Build the quadratic objective matrix Q: identity matrix for (Q_ij - P_ij)^2
Qobj <- diag(num_vars)
cvec <- -2 * Pvec[active_idx]

# Equality constraints: rows must sum to 1
# We'll build an A matrix with a row per row sum constraint
A <- matrix(0, nrow = n, ncol = num_vars)

for (i in 1:n) {
  for (j in 1:n) {
    idx <- (i - 1) * n + j
    if (Wvec[idx] == 1) {
      var_idx <- which(active_idx == idx)
      A[i, var_idx] <- 1
    }
  }
}

model <- list()
model$Q <- Qobj
model$obj <- cvec
model$A <- A
model$rhs <- rep(1, n)
model$sense <- rep("=", n)
model$lb <- rep(0, num_vars)
model$ub <- rep(1, num_vars)
model$vtype <- rep("C", num_vars)
model$modelsense <- "min"

# Solve
result <- gurobi(model)

# Reconstruct full Q matrix
Qvec <- rep(0, N)
Qvec[active_idx] <- result$x
Qmat <- matrix(Qvec, nrow = n, byrow = TRUE)

Qmat
