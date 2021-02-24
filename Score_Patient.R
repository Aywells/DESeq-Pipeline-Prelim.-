

#' Scoring the Patients
#'
#' @param Disease is a string with the disease to be tested.
#' @param RNAseq the RNAseq (raw counts) from the patient.
#' @param dir the directory of the most recent DEGs to compare against.
#' @param normalize TRUE or FALSE if the RNAseq needs to be Z-normalized.
#'
#' @importFrom data.table fread
#' @importFrom dplyr select filter
#' @importFrom utils globalVariables
#' @return a score.
#' @export

Score_Patient = function(Disease,
                         RNAseq,
                         dir = "~/Desktop/PostDoc/00_Projects/NomiX/Method_dev/99_Test/01_Update/",
                         normalize = F){
  `%>%` <- magrittr::`%>%`
  FC = Gene = padj = NULL
  if(normalize == T){
    RNAseq = scale(RNAseq)
  }

  DEGS = system(paste0("ls ", dir,"/00_FC/FC_", Disease, ".csv"), intern = T) %>%
    data.table::fread()
  DEGS_Name = paste0(dir, "/01_Genes/DEG_", Disease, ".csv") %>%
    data.table::fread()

  DEGS = subset(DEGS, DEGS$Gene %in% DEGS_Name$Gene)

  DEGS = subset(DEGS, DEGS$Gene %in% names(RNAseq))

  RNAseq = subset(RNAseq, select =   names(RNAseq)[names(RNAseq) %in% DEGS$Gene])

  RNAseq %<>% t() %>% as.data.frame()
  names(RNAseq) = "ID"
  RNAseq$Gene = row.names(RNAseq)
  DEG_ID = suppressMessages( suppressWarnings(  plyr::join(DEGS, RNAseq, type = "inner")) )
  Probs = data.frame(Gene = DEG_ID$Gene,
                     Prob = 1/(1+exp(-(DEG_ID$b0 + DEG_ID$FC*DEG_ID$ID))),
                     CO = DEG_ID$cutoff)
  SCR_PR = sum(Probs$Prob > Probs$CO) / nrow(Probs)

  return( list(Prob_Affected = data.frame(Prob = SCR_PR),
               Genes = Probs))
}
