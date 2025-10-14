#' @title Radial Basis Function (Gaussian) kernel
#'
#' @description Compute the RBF kernel matrix between two sets of row-vectors.
#'
#' @param X1 Numeric matrix of size n1 x p. Each row is a p-dimensional point.
#' @param X2 Numeric matrix of size n2 x p. Each row is a p-dimensional point.
#' @param sigma Positive scalar, kernel bandwidth. 
#'
#' @return Numeric matrix K of size n1 x n2 with entries \eqn{K_{ij} = \exp(-\|X_{1,i} - X_{2,j}\|^2 / (2\sigma^2))}.
#' @export
#'
#' @examples
#' X1 <- matrix(rnorm(6), nrow = 3)
#' X2 <- matrix(rnorm(8), nrow = 4)
#' K  <- rbf_kernel(X1, X2, sigma = 1)
#' 
rbf_kernel <- function(X1, X2, sigma) {
  inner_product <- X1 %*% t(X2)         
  X1_sq <- rowSums(X1^2)               
  X2_sq <- rowSums(X2^2)              
  sq_dist <- outer(X1_sq, X2_sq, "+") - 2 * inner_product  
  sq_dist[sq_dist < 0] <- 0          
  K <- exp(-sq_dist / (2 * sigma^2))  
  return(K)
}
