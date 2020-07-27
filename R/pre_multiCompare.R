#' Enriched pathways results of multiple lncRNAs per cohort or multiple cohort per lncRNA
#'
#' Combine the enriched pathways' results of multiple lncRNAs in one cohort
#' or one lncRNAs in multiple cohort
#'
#' @param files A list of data frames which stores enriched pathway results for multi-lncRNA per cohort
#'              or multi-cohort per lncRNA
#' @param compare A character string to specify comparison among lncRNAs or among cohorts
#'
#' @import data.table
#' @return A data frame ready for heatmap by plot_compareCohort
#' @example
#' prad <- pre_multiCompare(files = list("ENST00000625256_PRAD_cor.txt",
#'                                   "ENST00000561978_PRAD_cor.txt",
#'                                   "ENST00000561519_PRAD_cor.txt"),
#'                      compare = "lncRNA")
#' @export


pre_multiCompare <- function(files, compare = c("lncRNA", "cohort")){

    res <- vector(mode="list", length = length(files))
    res <- lapply(files, function(x) read.delim(x))
    compare <- match.arg(compare)
    if (compare == "lncRNA") {
        for (i in 1:length(res)){
        res[[i]]$cohort <- strsplit(files[[i]],"_")[[1]][1]
    }
}
    if (compare == "cohort") {
        for (i in 1:length(res)){
            res[[i]]$cohort <- strsplit(files[[i]],"_")[[1]][2]
        }
}


    res.df <- data.table::rbindlist(res)
    res.df <- res.df[, c("pathway","padj","NES","cohort")]
    res.df$Status <- ifelse(res.df$NES > 0, "Upregulated", ifelse(res.df$NES < 0, "Downregulated", ""))
    res.df$Status <- factor(res.df$Status, levels = c("Upregulated", "Downregulated", ""))

    res.df$FDR <- ifelse(res.df$padj > 0.05, "> 0.05", ifelse(res.df$padj < 0.001, "< 0.001", "< 0.05"))
    res.df$FDR <- factor(res.df$FDR, levels = c("> 0.05", "< 0.05", "< 0.001"))

    return(res.df)

}


