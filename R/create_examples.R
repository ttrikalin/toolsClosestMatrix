library(Matrix)

rm(list = ls())

#  this is a simple example of a system with a matrix P that we want to approximate 
# by a matrix Q that has structural zeros.

# A-> B -> C

l1 <- 0.03
l2 <- 0.015

# generator matrix
A <- matrix(
  c(-l1 , l1, 0, 
    0, -l2, l2, 
    0, 0, 0), nrow = 3, byrow = TRUE)

# The matrix that we wish to approximate has an non-0 A->C transition 
P <- expm(A)
colnames(P) <- c("toA", "toB", "toC")
rownames(P) <- c("fromA", "fromB", "fromC")
P
P["fromA","toC"]!=0

# The matrix with the structural constraint to enforce a A->C transition
M <- matrix(rep(1, length(P)), ncol = ncol(P))
colnames(M) <- colnames(P)
rownames(M) <- rownames(P)
M["fromA","toC"] <- 0
M

saveRDS(list(A=A, P=P, M=M), "data/example1.R")
