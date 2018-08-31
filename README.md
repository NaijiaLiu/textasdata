# textasdata
## For the package:

library(MCMCpack)
library(topicmodels)
library(devtools)

devtools::build()
devtools::load_all()
devtools::document()

## There are two functions in the package and this is the way to get manual
?Marginalikelihood 
?visualizeTMM