# spherical case: compare Monte Carlo method to analytic formula

d <- 2 # 2D example
mA <- rnorm(d)
mB <- rnorm(d)
sigA <- rexp(1)
SigA <- sigA^2 * diag(d)
sigB <- rexp(1)
SigB <- sigB^2 * diag(d)

# plot circles
ellA <- ellipse(mA, SigA)
ellB <- ellipse(mB, SigB)
plot(0, type = "n",
     xlim = range(ellA[,1], ellB[,1]),
     ylim = range(ellA[,2], ellB[,2]), xlab = "x", ylab = "y")
lines(ellA, col = "red")
lines(ellB, col = "blue")
legend("topright", legend = c("niche A", "niche B"),
       fill = c("red", "blue"), bg = "white")

# compare niche calculations
overlap.sphere(mA, sigA, mB, sigB)
overlap.unif(mA, SigA, mB, SigB, nprob = 1e5)
