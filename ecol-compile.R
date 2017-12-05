##############################################3

# script for compiling an R package
# martin lysy, 2014

require(devtools)
require(knitr)

# do this on a clean workspace
rm(list = ls())

###############################################

# ok.  new package building instructions.

pkg.name <- "nicheROVER"
pkg.path <- "C:/Users/Jerome/Documents/R/Rpackages/nicheROVER/"
# source file containing all functions, documentation, examples
# written in roxygen2 format
pkg.src <- "ecol-package2.R"
# file containing all datasets
pkg.data <- "ecol-data.RData"
# description file
pkg.descr <- "ecol-DESCRIPTION"
# package vignette to be generated with knitr
pkg.vign <- "ecol-vignette.Rmd"

###########################################33

# create package skeleton
# loading workspace objects was the only way i could figure out to include
# both functions and datasets into the package...
tmp.ls <- c(ls(), "tmp.ls")
load(pkg.data)
source(pkg.src)
package.skeleton(list = ls()[!ls() %in% tmp.ls],
                 name = pkg.name, force = TRUE)

###############################################

# create the documentation

# separate the package file from the main doc source file,
# and place each of these into the package folder
con <- file(pkg.src, "r")
zz <- readLines(con)
close(con)
ind <- which(sapply(zz, grepl, pattern = "NULL"))[1]
# source file
cat(zz[(ind+1):length(zz)], sep = "\n",
    file = paste0(pkg.path, pkg.name, "/R/", pkg.src))
# package file
cat(zz[1:ind], sep = "\n",
    file = paste0(pkg.path, pkg.name, "/R/",
      pkg.name, "-package.R"))

# add the description file to package folder
file.copy(from = paste0(pkg.path, pkg.descr),
          to = paste0(pkg.path, pkg.name, "/DESCRIPTION"),
          overwrite = TRUE)

# get the version number from description file.
# good to have around for install command coming up.
con <- file(pkg.descr, "r")
zz <- readLines(con)
close(con)
ind <- which(sapply(zz, grepl, pattern = "Version:"))[[1]]
pkg.vers <- strsplit(zz[ind], split = " ")[[1]][2]
pkg.vers

# add vignette to package folder
dir.create(path = paste0(pkg.path, pkg.name, "/vignettes"))
file.copy(from = paste0(pkg.path, pkg.vign),
          to = paste0(pkg.path, pkg.name, "/vignettes/", pkg.vign),
          overwrite = TRUE)

# roxygenize all of the documentation together
document(pkg = paste0(pkg.path, pkg.name), clean = TRUE)

##########################################################3

# build, install, check package (for windows)

# build package
cmd <- paste0(R.home(component = "bin"), "/R CMD build ",
              pkg.name)
compiled <- system(cmd, intern = FALSE)

# install package
cmd <- paste0(R.home(component = "bin"), "/R CMD INSTALL ",
              pkg.name, "_", pkg.vers, ".tar.gz")
compiled <- system(cmd, intern = FALSE)

# check package
run.ex <- TRUE # set to true for full CRAN check, false to skip examples (faster)
cmd <- paste0(R.home(component = "bin"), "/R CMD check ",
              ifelse(run.ex, "--as-cran ", "--no-examples "),
              pkg.name, "_", pkg.vers, ".tar.gz")
compiled <- system(cmd, intern = FALSE)

##########################################################

# tests

# first, quit R...
q(save = "no")

require(nicheROVER)

browseVignettes(package = "nicheROVER")

example(niche.plot, ask = FALSE)

example(overlap.plot, ask = FALSE)

example(nicheROVER)

help(package = "nicheROVER")

