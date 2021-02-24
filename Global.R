#' Global Definition
#'

#' @docType package
#' @name NomiX
#' @importFrom dplyr %>%

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("."))

pack = c("data.table",
         "magrittr",
         "plyr",
         "parallel",
         "dplyr",
         "OptimalCutpoints",
         "stringr",
         "utils")

for ( i in 1:length(pack)){
  pkg = pack[i]

  if(!require(pkg, character.only = TRUE)){
    install.packages(pack[i], dependencies = TRUE)
  }
}

