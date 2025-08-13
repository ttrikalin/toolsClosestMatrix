#' Check if a matrix is stochastic
#' 
#' This function checks whether a matrix is stochastic (row stochastic, column stochastic, 
#' or doubly stochastic) based on the specified type.
#' 
#' @param matrix A numeric matrix to check
#' @param type Character string specifying the type of stochastic matrix to check:
#'   - "row": row stochastic (each row sums to 1)
#'   - "column": column stochastic (each column sums to 1) 
#'   - "doubly": doubly stochastic (both row and column stochastic)
#'   - "any": checks if matrix is either row or column stochastic
#' @param tolerance Numeric tolerance for floating point comparisons (default: 1e-10)
#' @return Logical value indicating whether the matrix is stochastic of the specified type
#' 
#' @examples
#' # Create a row stochastic matrix
#' P <- matrix(c(0.5, 0.3, 0.2,
#'               0.2, 0.6, 0.2,
#'               0.1, 0.4, 0.5), nrow = 3, byrow = TRUE)
#' 
#' is_stochastic(P, "row")     # TRUE
#' is_stochastic(P, "column")  # FALSE
#' is_stochastic(P, "doubly")  # FALSE
#' is_stochastic(P, "any")     # TRUE
#' 
#' @export
is_stochastic <- function(matrix, type = "row", tolerance = 1e-10) {
  # Input validation
  if (!is.matrix(matrix) || !is.numeric(matrix)) {
    stop("Input must be a numeric matrix")
  }
  
  if (!type %in% c("row", "column", "doubly", "any")) {
    stop("Type must be one of: 'row', 'column', 'doubly', 'any'")
  }
  
  # Check if all elements are non-negative
  if (any(matrix < -tolerance)) {
    return(FALSE)
  }
  
  # Check row stochastic (each row sums to 1)
  is_row_stochastic <- all(abs(rowSums(matrix) - 1) < tolerance)
  
  # Check column stochastic (each column sums to 1)
  is_column_stochastic <- all(abs(colSums(matrix) - 1) < tolerance)
  
  # Return result based on type
  switch(type,
         "row" = is_row_stochastic,
         "column" = is_column_stochastic,
         "doubly" = is_row_stochastic && is_column_stochastic,
         "any" = is_row_stochastic || is_column_stochastic)
}

#' Check if a matrix is row stochastic
#' 
#' Convenience function to check if a matrix is row stochastic.
#' 
#' @param matrix A numeric matrix to check
#' @param tolerance Numeric tolerance for floating point comparisons (default: 1e-10)
#' @return Logical value indicating whether the matrix is row stochastic
#' 
#' @examples
#' P <- matrix(c(0.5, 0.3, 0.2,
#'               0.2, 0.6, 0.2,
#'               0.1, 0.4, 0.5), nrow = 3, byrow = TRUE)
#' is_row_stochastic(P)  # TRUE
#' 
#' @export
is_row_stochastic <- function(matrix, tolerance = 1e-10) {
  is_stochastic(matrix, "row", tolerance)
}

#' Check if a matrix is column stochastic
#' 
#' Convenience function to check if a matrix is column stochastic.
#' 
#' @param matrix A numeric matrix to check
#' @param tolerance Numeric tolerance for floating point comparisons (default: 1e-10)
#' @return Logical value indicating whether the matrix is column stochastic
#' 
#' @examples
#' P <- matrix(c(0.5, 0.2, 0.1,
#'               0.3, 0.6, 0.4,
#'               0.2, 0.2, 0.5), nrow = 3, byrow = TRUE)
#' is_column_stochastic(P)  # TRUE
#' 
#' @export
is_column_stochastic <- function(matrix, tolerance = 1e-10) {
  is_stochastic(matrix, "column", tolerance)
}

#' Check if a matrix is doubly stochastic
#' 
#' Convenience function to check if a matrix is doubly stochastic.
#' 
#' @param matrix A numeric matrix to check
#' @param tolerance Numeric tolerance for floating point comparisons (default: 1e-10)
#' @return Logical value indicating whether the matrix is doubly stochastic
#' 
#' @examples
#' P <- matrix(c(0.5, 0.3, 0.2,
#'               0.3, 0.4, 0.3,
#'               0.2, 0.3, 0.5), nrow = 3, byrow = TRUE)
#' is_doubly_stochastic(P)  # TRUE
#' 
#' @export
is_doubly_stochastic <- function(matrix, tolerance = 1e-10) {
  is_stochastic(matrix, "doubly", tolerance)
}


row_normalize <- function(A){
  for(r in 1:nrow(A)) {
    rs <- sum(A[r,])
    A[r,] <- A[r,] / rs
  }
  return(A)
}


#' Calculate theoretical trajectory using matrix powers
#' 
#' Calculates the theoretical probability distribution over time using P^t * initial_state.
#' 
#' @param P Transition matrix
#' @param initial_state Initial state distribution
#' @param n_steps Number of time steps
#' @param state_names Names for the states (optional)
#' @return Matrix where rows are states and columns are time points
#' 
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.1, 0.9), nrow = 2, byrow = TRUE)
#' initial <- c(1, 0)
#' theoretical <- calculate_theoretical_trajectory(P, initial, 100)
#' 
#' @export
calculate_theoretical_trajectory <- function(P, initial_state, n_steps, state_names = NULL) {
  # Input validation
  if (!is.matrix(P) || !is.numeric(P)) {
    stop("P must be a numeric matrix")
  }
  if (length(initial_state) != nrow(P)) {
    stop("Length of initial_state must match number of rows in P")
  }
  
  # Set state names if not provided
  if (is.null(state_names)) {
    state_names <- paste0("State_", 1:nrow(P))
  }
  
  n_states <- nrow(P)
  time_points <- 0:n_steps
  
  # Initialize probability matrix
  prob_matrix <- matrix(0, nrow = n_states, ncol = n_steps + 1)
  prob_matrix[, 1] <- initial_state
  
  # Calculate P^t * initial_state for each time step
  P_power <- diag(n_states)  # P^0 = I
  for (t in 1:n_steps) {
    P_power <- P_power %*% P
    prob_matrix[, t + 1] <- initial_state %*% P_power
  }
  
  rownames(prob_matrix) <- state_names
  colnames(prob_matrix) <- time_points
  
  return(prob_matrix)
}

#' Plot and compare Markov chain trajectories (mean only)
#' 
#' Creates comprehensive plots comparing mean trajectories from two transition matrices.
#' Uses theoretical calculations (no simulation) for efficiency.
#' 
#' @param P Original transition matrix
#' @param Q Approximant transition matrix
#' @param initial_state Initial state distribution
#' @param n_steps Number of time steps to calculate
#' @param state_names Names for the states (optional)
#' @param plot_type Type of plot: "individual", "comparison", or "both"
#' @param colors Color palette for states
#' @param title Plot title
#' @return List containing theoretical calculations
#' 
#' @examples
#' P <- matrix(c(0.8, 0.2, 0.1, 0.9), nrow = 2, byrow = TRUE)
#' Q <- matrix(c(0.9, 0.1, 0.05, 0.95), nrow = 2, byrow = TRUE)
#' initial <- c(1, 0)
#' result <- plot_markov_comparison(P, Q, initial, n_steps = 100)
#' 
#' @export
plot_markov_comparison <- function(P, Q, initial_state, n_steps = 100, 
                                   state_names = NULL, plot_type = "both", 
                                   colors = NULL, title = NULL) {
  
  # Set state names if not provided
  if (is.null(state_names)) {
    state_names <- paste0("State_", 1:nrow(P))
  }
  
  # Set colors if not provided
  if (is.null(colors)) {
    colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
  }
  
  # Calculate theoretical trajectories (mean)
  cat("Calculating mean trajectories for matrix P...\n")
  mean_P <- calculate_theoretical_trajectory(P, initial_state, n_steps, state_names)
  
  cat("Calculating mean trajectories for matrix Q...\n")
  mean_Q <- calculate_theoretical_trajectory(Q, initial_state, n_steps, state_names)
  
  # Create plots
  if (plot_type %in% c("individual", "both")) {
    # Individual plots for P and Q
    par(mfrow = c(2, 1))
    
    # Plot for matrix P
    plot(0:n_steps, mean_P[1, ], type = "l", col = colors[1], lwd = 2,
         ylim = c(0, 1), xlab = "Time Step", ylab = "Probability",
         main = paste("Mean Trajectories using Matrix P", ifelse(is.null(title), "", paste(" -", title))),
         cex.main = 0.9)
    
    for (i in 2:nrow(mean_P)) {
      lines(0:n_steps, mean_P[i, ], col = colors[i], lwd = 2)
    }
    
    legend("topright", legend = state_names, col = colors[1:nrow(mean_P)], 
           lwd = 2, cex = 0.8, title = "States")
     
    # Plot for matrix Q
    plot(0:n_steps, mean_Q[1, ], type = "l", col = colors[1], lwd = 2,
         ylim = c(0, 1), xlab = "Time Step", ylab = "Probability",
         main = paste("Mean Trajectories using Matrix Q", ifelse(is.null(title), "", paste(" -", title))),
         cex.main = 0.9)
    
    for (i in 2:nrow(mean_Q)) {
      lines(0:n_steps, mean_Q[i, ], col = colors[i], lwd = 2)
    }
    
    legend("topright", legend = state_names, col = colors[1:nrow(mean_Q)], 
           lwd = 2, cex = 0.8, title = "States")
  }
  
  if (plot_type %in% c("comparison", "both")) {
    # Comparison plot
    if (plot_type == "both") {
      par(mfrow = c(1, 1))
    }
    
    plot(0:n_steps, mean_P[1, ], type = "l", col = colors[1], lwd = 2,
         ylim = c(0, 1), xlab = "Time Step", ylab = "Probability",
         main = paste("Comparison: P vs Q Mean Trajectories", ifelse(is.null(title), "", paste(" -", title))),
         cex.main = 0.9)
    
    # Add all other states for P
    for (i in 2:nrow(mean_P)) {
      lines(0:n_steps, mean_P[i, ], col = colors[i], lwd = 2)
    }
    
    # Add Q trajectories with dashed lines
    for (i in 1:nrow(mean_Q)) {
      lines(0:n_steps, mean_Q[i, ], col = colors[i], lwd = 2, lty = 2)
    }
    
    legend("topright", legend = c(paste(state_names, "(P)"), paste(state_names, "(Q)")),
           col = rep(colors[1:nrow(mean_P)], 2), 
           lwd = 2, lty = rep(c(1, 2), each = nrow(mean_P)), cex = 0.8)
  }
  
  # Reset plot parameters
  par(mfrow = c(1, 1))
  
  # Return results
  return(list(
    mean_trajectory_P = mean_P,
    mean_trajectory_Q = mean_Q,
    state_names = state_names,
    initial_state = initial_state
  ))
}
