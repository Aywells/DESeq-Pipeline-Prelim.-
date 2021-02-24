#' Drug-Gene Enrichment test
#'
#' @param Background The complete list of tested genes (DEGs from the diasease).
#' @param Genes2Drug a list (from drug bank) of "Gene_Target" and drug "ID".
#' @param Genes the list of significative DEGs of each patient
#' @importFrom magrittr `%>%` `%<>%`
#' @importFrom dplyr group_by summarise
#' @importFrom plyr join
#' @importFrom stats prop.test
#'
#' @export

Enrichment = function (Background, Genes2Drug, Genes)
{
  `%>%` <- magrittr::`%>%`
  `%<>%` <- magrittr::`%<>%`
  Drug = n = NULL
  Genes2Drug = data.frame(Gene.symbol = Genes2Drug$Gene_Target,
                          Drug = Genes2Drug$ID, stringsAsFactors = F) %>% unique()
  Genes2Drug = subset(Genes2Drug, Genes2Drug$Gene.symbol %in%
                        Background) %>% unique()

  Exp = Genes2Drug %>%
    dplyr::group_by(Drug) %>%
    dplyr::summarise( n = n())

  Obs_aux = subset(Genes2Drug, Genes2Drug$Gene.symbol %in% Genes$Gene)

  Obs = Obs_aux %>%
    dplyr::group_by(Drug) %>%
    dplyr::summarise( n = n())


  if(nrow(Obs) > 0){
    names(Exp)[2] = "Exp"
    names(Obs)[2] = "Obs"

    Exp_Obs = suppressMessages( plyr::join(Exp, Obs))
    Exp_Obs[is.na(Exp_Obs)]<-0
    Exp_Obs$ExpTotal = sum(Exp_Obs$Exp)- Exp_Obs$Exp
    Exp_Obs$ObsTotal = sum(Exp_Obs$Obs)- Exp_Obs$Obs

    ## Loop for fisher's
    Exp_Obs$p = apply(Exp_Obs,1, FUN = function(x) {
      e = c(x[3],x[2]) %>% matrix(., ncol= 2) %>% as.numeric()
      o = c(x[5],x[4]) %>% matrix(., ncol= 2) %>% as.numeric()
      p = suppressWarnings(stats::fisher.test(rbind(e,o), alternative = "greater")$p.value)
      return(p)
    })
  } else if(nrow(Obs) == 0){
    Exp_Obs$Obs = 0
    Exp_Obs$p = 0.9999
  }
  Exp_Obs$Obs = 0
  Exp_Obs$p = 0.9999

  Exp_Obs$Pad = p.adjust(Exp_Obs$p)
  Your_Drugs = subset(Exp_Obs, Exp_Obs$Pad<0.05)

  Exp_Obs = Exp_Obs[order(c(Exp_Obs$Pad, Exp_Obs$p), decreasing = TRUE),]
  Exp_Obs$Pad = ifelse(Exp_Obs$Pad == 1, 0.9999, Exp_Obs$Pad)

  return(list (Your_Drugs, Exp_Obs))
}
