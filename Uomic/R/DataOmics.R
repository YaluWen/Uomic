#' Simulated omics data for 500 subjects. It has gene expression, methylation levels, and SNPs data for gene \emph{RB1}. 
#' 
#'
#' @format A data list with three elements
#' \describe{
#'   \item{Expr}{gene expression level for gene \emph{RB1}, and it is simulated using the InterSIM software}
#'   \item{Met}{methylation levels for all CpG sites within gene \emph{RB1}, and they are simulated using the InterSIM software}
#'   \item{SNPs}{genomic data for all genetic variants within gene \emph{RB1}, and they were simulated using the hapgen2 software with CEU in HapMap project as reference}
#' }
"Omics"