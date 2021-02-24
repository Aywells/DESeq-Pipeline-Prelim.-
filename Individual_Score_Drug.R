

#' Calculates the prob of a patient be affected by the disease and also its drugs.
#'
#' @param Disease Name of the disease to be tested
#' @param RNAseq Normalized RNAseq for each patient
#' @param dir Where the data is stored
#' @importFrom utils read.table
#'
#' @export

Score_Individual = function(Disease = "Asthma",
                            RNAseq,
                            dir = "~/Desktop/PostDoc/00_Projects/NomiX/Method_dev/99_Test/01_Update/"){

  DB = Gene=cutoff_disease = NULL
  cutoff_disease = utils::read.table(paste0(dir, "/02_CutOff/", Disease, "_CutOff.txt"))
  AUX = Score_Patient(Disease, RNAseq = RNAseq,
                      normalize = F,
                      dir = dir)
  Affected = ifelse(AUX[[1]]$Prob>cutoff_disease, "Affected", "Healthy")

  Out = data.frame(Prob = AUX[[1]]$Prob, Affected)

  Your_Genes = subset(AUX[[2]], AUX[[2]]$Prob > AUX[[2]]$CO, select = Gene)
  message("Scoring Drugs\n")
  load("~/Desktop/PostDoc/00_Projects/NomiX/Method_dev/01_Codes/99_Package/01_Data/DrugBank.RData")
  Genes2Drug = subset(DB, DB$Gene_Target != "", select = c("Gene_Target", "ID"))

  Drugs_4U = Enrichment(Background = AUX[[2]]$Gene, Genes2Drug, Genes = Your_Genes)
  Drugs_4U

  DB_map = DB %>% select(ID, Name, Approved) %>% unique()
  names(DB_map)[1] = "Drug"
  Drugs_4U[[2]]<- plyr::join(DB_map, Drugs_4U[[2]]) %>% na.exclude()

  return(list(Affected = Out, Drugs_4U = Drugs_4U))
}


