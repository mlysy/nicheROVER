# 2d example
d <- 2 # number of dimensions
V <- crossprod(matrix(rnorm(4),d,d))
mu <- rnorm(d)
plot(ellipse(mu, V), type = "l")
points(niche.runif(1e4, mu, V), col = "brown", pch = ".")
