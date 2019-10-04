# for each species, size of 95% niche region using sample variance
tapply(1:nrow(fish), fish$species, function(ind) {
  X <- fish[ind,2:4] # all measurements for given species
  Sighat <- var(X) # sample variance
  niche.size(Sigma = Sighat, alpha = .95)
})
