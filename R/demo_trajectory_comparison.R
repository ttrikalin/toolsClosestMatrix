#' Demonstration script for comparing Markov chain trajectories (mean only)
#' 
#' This script demonstrates how to use the trajectory plotting functions
#' to compare the behavior of transition matrices P and Q.

# Load required libraries
library(Matrix)

# Source the functions
source("R/functions.R")

# Load the example data
  example <- readRDS("data/example1.rds")
  P <- as.matrix(example$P)
  M <- example$M
  
  cat("Loaded example data:\n")
  cat("Matrix P (original):\n")
  print(P)
  cat("\nMatrix M (structural constraints):\n")
  print(M)
  
  # Get the approximant matrix Q using the existing function
  cat("\nCalculating approximant matrix Q...\n")
  result <- get_approximant(P, M, norm = "frobenius")
  Q <- result$Q
  
  cat("Matrix Q (approximant):\n")
  print(Q)
  
  # Check if matrices are stochastic
  cat("\nChecking stochasticity:\n")
  cat("P is row stochastic:", is_row_stochastic(P), "\n")
  cat("Q is row stochastic:", is_row_stochastic(Q), "\n")
  
  # Set up initial state and parameters
  initial_state <- c(1, 0, 0)  # Start in state A
  state_names <- c("A", "B", "C")
  n_steps <- 200
  
  cat("\nCalculating mean trajectories...\n")
  cat("Initial state:", state_names[which(initial_state == 1)], "\n")
  cat("Number of time steps:", n_steps, "\n")
  
  # Create comprehensive comparison plots
  cat("\nCreating plots...\n")
  result_plots <- plot_markov_comparison(
    P = P, 
    Q = Q, 
    initial_state = initial_state,
    n_steps = n_steps,
    state_names = state_names,
    plot_type = "both",
    title = "Example 1: A->B->C System"
  )
  
  # Display summary statistics
  cat("\nSummary of results:\n")
  cat("Final probabilities using P:\n")
  final_P <- result_plots$mean_trajectory_P[, ncol(result_plots$mean_trajectory_P)]
  for (i in 1:length(state_names)) {
    cat(sprintf("  %s: %.4f\n", state_names[i], final_P[i]))
  }
  
  cat("\nFinal probabilities using Q:\n")
  final_Q <- result_plots$mean_trajectory_Q[, ncol(result_plots$mean_trajectory_Q)]
  for (i in 1:length(state_names)) {
    cat(sprintf("  %s: %.4f\n", state_names[i], final_Q[i]))
  }
  
  # Calculate differences
  cat("\nDifferences (P - Q):\n")
  differences <- final_P - final_Q
  for (i in 1:length(state_names)) {
    cat(sprintf("  %s: %.4f\n", state_names[i], differences[i]))
  }
  
  cat("\nMaximum absolute difference:", max(abs(differences)), "\n")
  
