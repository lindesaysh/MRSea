detach(package:MRSea)
library(devtools)

# setwd('C:\Users\...\Dropbox\Marine Scotland\package\MRSea')

#setwd('../Dropbox/Marine Scotland/package_postWorkshop/MRSea')
#setwd('/Users/Lindesay/Dropbox/Marine Scotland/package_postWorkshop/MRSea/')
setwd('/Users/Lindesay/Dropbox/Research/MarineScotland_LSH/package/MRSea/')
#setwd('Dropbox/MarineScotland_LSH/package/MRSea/')

document(pkg="MRSea/.")  #updates documentation modified via roxygen2

# change namespace file 
#  this should execute operating system commands to handle the namespace change
#shell("del NAMESPACE")
#shell("copy NAMESPACE.txt NAMESPACE")

build(pkg="MRSea/.", binary=TRUE)

# re-install package
# restart workspace
# load package
require(MRSea)
?makesplineParams
?runSALSA1D
?runSALSA2D


# terminal/command prompt code for checking package and making manual
# note the check function creates the manual too
# R CMD Check MRSea

# the no clean means a hidden folder with the tex file is created.
# R CMD Rd2pdf --no-clean