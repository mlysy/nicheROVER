# single multivariate normal draw
.rmvn <- function(mean, sigma) {
  (rnorm(length(mean)) %*% chol(sigma))[1,] + mean
}
