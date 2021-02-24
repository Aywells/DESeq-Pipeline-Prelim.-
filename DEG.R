####################################################################
### Define GDAs for a disease
### Function based ONLY in the DEGs.
### In the future, the idea is to implement it based on other GDAs.
####################################################################
### Author: Deisy Morselli Gysi
### Date: 12/01/2019
### Version: 1
####################################################################

#' DEGs: Define Differentially Expressed Genes
#' @description This function defines the set of DEGs for each disease. Should be run every 96 new patients. So we can update our model.
#'
#' @param Disease_vector Is a binary (0,1) vector if the patient has or not the disease under study.
#' @param Disease_name The name of the disease under study.
#' @param Data data.frame containing all the normalized (z-scored) gene expression of all patients.
#' @param k Number of threads to use for parallel function.
#' @param method for the cutoff definition.

#' @return a data.frame containing the Gene name, Fold Change, the p-value and the p value FDR adjusted for the selected disease.
#' @export
#' @importFrom stats glm p.adjust
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom parallel makeCluster clusterExport clusterApplyLB stopCluster
#' @importFrom data.table rbindlist
#'
#' @examples
#'
#' #load("~/Desktop/PostDoc/00_Projects/NomiX/Method_dev/01_Codes/99_Package/01_Data/example.RData")
#' #load("~/Desktop/PostDoc/00_Projects/NomiX/Method_dev/01_Codes/99_Package/01_Data/Pheno.RData")
#' #DEGs = DEG(Data = example[,1:100],
#' #Disease_vector = Pheno$Asthma,
#' #Disease_name = "Asthma",
#' #k = 5)
#'
#'
DEG = function(Data, Disease_vector, k, Disease_name, method = "MaxDOR"){
  `%>%` <- magrittr::`%>%`
  cl <- parallel::makeCluster(k)
  # parallel::clusterExport(cl, "Disease_vector")
  # parallel::clusterExport(cl, "Data")
  parallel::clusterExport(cl, "DEG_aux")
  parallel::clusterExport(cl, "cutoff")
  parallel::clusterExport(cl, "%>%")

  len = ncol(Data)

  Output = parallel::clusterApplyLB(cl, 1:len, DEG_aux,  Disease_vector = Disease_vector,
                                    Disease_name = Disease_name,
                                    Data = Data, method = method)
  Output = data.table::rbindlist(Output)
  Output$padj = stats::p.adjust(Output$p, "fdr")
  parallel::stopCluster(cl)

  return(Output)
}

