##### Cutoff definition for each gene


##########################################################################
# Sensitivity equal to Specificity Method ("SpEqualSe"): Covariate gene
##########################################################################

#' Cutoff for each gene
#'
#' @param probs vector of probability of the Individual to be affected by the disease
#' @param Disease vector (binary) of the individuals affected by the disease
#' @param method for the cutoff definition.
#' @importFrom stringr str_remove_all str_replace_all str_split
#' @importFrom OptimalCutpoints optimal.cutpoints control.cutpoints summary.optimal.cutpoints

#' @return the cutoff prob for each gene
#' @export

cutoff = function( probs, Disease, method = "MaxDOR" ){
  probs %<>% as.numeric()
  data = data.frame(Gene = probs, Pheno = Disease)
  optimal.cutpoint.SpEqualSe <- OptimalCutpoints::optimal.cutpoints(X = "Gene", status = "Pheno", tag.healthy = 0,
                                                                    methods = method, pop.prev = NULL,
                                                                    control = OptimalCutpoints::control.cutpoints(),
                                                                    ci.fit = TRUE, conf.level = 0.99,
                                                                    trace = FALSE, data = data)
  Cutoff_tmp = OptimalCutpoints::summary.optimal.cutpoints(optimal.cutpoint.SpEqualSe)
  Cutoff = Cutoff_tmp[[1]]$Global$optimal.cutoff$cutoff
  AUC = Cutoff_tmp$p.table$Global$AUC_CI %>%
    stringr::str_remove_all(., " ") %>%
    stringr::str_replace_all(., "\\(", ",") %>%
    stringr::str_remove_all(., "\\)") %>%
    stringr::str_split(., ",", simplify = TRUE) %>%
    as.numeric()
  names(AUC) = c("AUC", "AUC_l", "AUC_u")
  return(list(Cutoff, AUC))
}
