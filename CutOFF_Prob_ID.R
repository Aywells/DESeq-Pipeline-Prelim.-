#' Cutoff for probability of having the disease
#'
#' @param Prob vector of probability of the Individual to be affected by the disease
#' @param Pheno vector (binary) of the individuals affected by the disease
#' @param method a method for the cutoff selection
#' @param conf.level the confidence level for the CI
#' @importFrom stringr str_remove_all str_replace_all str_split
#' @importFrom OptimalCutpoints optimal.cutpoints control.cutpoints summary.optimal.cutpoints

#' @return the cutoff prob for individual to be affected by the disease
#' @export


CutOff_Prob_ID <- function(Pheno, Prob,
                           method = "MaxSp",
                           conf.level = 0.90){

  data = data.frame(Pheno = Pheno, Prob = Prob)
  optimal.cutpoint <- OptimalCutpoints::optimal.cutpoints(X = "Prob", status = "Pheno",
                                                          tag.healthy = 0,
                                                          methods = method,
                                                          pop.prev = NULL,
                                                          control = OptimalCutpoints::control.cutpoints(),
                                                          ci.fit = TRUE,
                                                          conf.level = conf.level,
                                                          trace = FALSE,
                                                          data = data)
  # OptimalCutpoints::plot.optimal.cutpoints(optimal.cutpoint)
  CF = optimal.cutpoint[[1]]$Global$optimal.cutoff$cutoff

  return(CF)
}
