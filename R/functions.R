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
