# read functions
source("R/functions.R")
source("R/get_approximant.R")

# example A -> B -> C, embeddable 
ex1 <- readRDS("data/example1.R")

# infinity norm
Q <- get_approximant(P = ex1$P, M = ex1$M, norm = "infinity")
ex1$P; Q
