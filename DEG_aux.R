####################################################################
### Define GDAs for a disease.
### Internal function.
### Function based ONLY in the DEGs.
### In the future, the idea is to implement it based on other GDAs.
####################################################################
### Author: Deisy Morselli Gysi
### Date: 12/01/2019
### Version: 1
####################################################################

#' Auxiliar function for DEGs
#'
#' @param i internal paramenter for the parallel function.
#' @param Disease_vector Is a binary (0,1) vector if the patient has or not the disease under study.
#' @param Disease_name The name of the disease under study.
#' @param Data data.frame containing all the normalized (z-scored) gene expression of all patients.
#' @param method for the cutoff definition.

#' @export
#' @importFrom stats glm
#' @importFrom magrittr "%>%" "%<>%"
#'
DEG_aux = function(i,
                   Disease_vector = Disease_vector,
                   Disease_name = Disease_name,
                   Data = Data, method = "MaxDOR"){
  `%>%` <- magrittr::`%>%`
  G_ID = names(Data)[i]
  Gene = Data[,i] %>% as.vector()
  Disease_vector =  as.vector(Disease_vector)
  fit = stats::glm(Disease_vector ~ Gene, family = "binomial") %>%
    summary()
  b0 = fit$coefficients[1,1]
  fc = fit$coefficients[2,1]
  p = fit$coefficients[2,4]

  probs = 1/(1+exp(-(b0 + fc*Gene)))
  cutoffs = NomiX::cutoff(probs, Disease_vector, method = method)

  Out = cbind.data.frame(Gene = G_ID,
                         b0 = b0,
                         FC = fc,
                         pval = p,
                         Disease = Disease_name,
                         cutoff = cutoffs[[1]],
                         AUC = cutoffs[[2]][1],
                         AUC_l = cutoffs[[2]][2],
                         AUC_u = cutoffs[[2]][3])
  return(Out)
}
