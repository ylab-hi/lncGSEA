#' \code{lncGSEA} package
#'
#' lncRNA associated GSEA
#'
#' @docType package
#' @name lncGSEA
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("NES", "logP.sign", "name",
                                                        "pathway", "logP", "cohort",
                                                        "pathway", "Status", "NES",
                                                        "FDR"))
