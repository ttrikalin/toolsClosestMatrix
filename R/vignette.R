# read functions
source("R/functions.R")
source("R/get_approximant.R")

# example A -> B -> C, embeddable 
ex1 <- readRDS("data/example1.rds")

# infinity norm -- sol1 
result1 <- get_approximant(P = ex1$P, M = ex1$M, norm = "infinity")
result1$LP
ex1$P; result1$Q

# infinity norm -- sol2
M2 <- ex1$M
M2[1,1] <- ex1$P[1,1]  # do not let it increase Q[1,1] above P[1,1]
result2 <- get_approximant(P = ex1$P, M = M2, norm = "infinity")
result2$LP  # same objective value!
ex1$P; result2$Q

# the above makes sense -- simply move mass. 
# you only see which rows need some change 
# and you simply redistribute in that row's sum  
# so the infinity norm is not doing anything useful -- any redistribution works. 

# 1 norm -- sol3  -- again, simply moves mass. 
result3 <- get_approximant(P = ex1$P, M = ex1$M, norm = "1")
result3$LP
ex1$P; result3$Q


## No difference between the matrices when there is no constraints 
Msame <- ex1$M
Msame[1,3] <- 1  
result4 <- get_approximant(P = ex1$P, M = Msame, norm = "infinity")
result4$LP  # same objective value!
ex1$P; result4$Q



## This is the spectral norm approximant 
result5 <- get_approximant(P = ex1$P, M = ex1$M, norm = "2")
result5$result$value  
ex1$P; result5$Q
