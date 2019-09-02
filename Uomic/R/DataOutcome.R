#' Simulated outcomes under different models.
#'
#'
#' @format A data frame with 500 rows and 8 variables:
#' \describe{
#'   \item{all}{all layers of omics data contributed to disease risk (i.e., 5 SNPs, methylation levels at 5 CpG sites and gene expression levels).}
#'   \item{gen}{5 SNPs were selected to be causal.}
#'   \item{met}{methylation levels at 5 CpG sites were selected to be causal.}
#'   \item{expr}{gene expression level was associated with the outcome.}
#'   \item{gen_met}{5 SNPs and methylation levels at 5 CpG sites were associated with the outcome.}
#'   \item{gen_expr}{5 SNPs and gene expression level were associated with the outcome.}
#'   \item{met_expr}{methylation levels at 5 CpG sites and gene expression level were associated with the outcome.}
#'   \item{null}{none of the omics data were associated with the outcome}
#' }
"outcome"