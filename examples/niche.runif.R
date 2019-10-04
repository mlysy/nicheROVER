# 2d example
V <- crossprod(matrix(rnorm(4),2,2))
mu <- rnorm(2)
plot(ellipse(mu, V), type = "l")
points(niche.runif(1e4, mu, V), col = "brown", pch = ".")
