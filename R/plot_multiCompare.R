#' Create data frame for plot_compareCohort
#'
#' Combine outputs of multiple cohort for a lncRNA from lnc_GSEA for heatmap plot
#'
#' @param lncRNA A character string which is the name of a transcript id, e.g.,"ENST00000561519"
#' @param cohorts A upper case string showing one of TCGA cohorts, e.g, "PRAD"
#' @param geneset A character string which gives the path to a gene set gmt file
#' @import data.table
#' @import ggplot2
#' @return A data frame contains GSEA results of one lncRNA in multiple studies
#' @example
#' cohorts <- c("PRAD","LUAD","BRCA")
#' arlnc1.df <- pre_compareCohort(lncRNA="ENST00000561519", cohorts = cohorts)
#' @export

pre_compareCohort <- function(lncRNA, cohorts, geneset = NULL){
    # apply pre_gsea to multiple cohorts ---
    files <- sapply(cohorts, pre_gsea, lncRNA)
    print(length(files))

    res <- vector(mode = "list", length = length(files))
    for (i in 1:length(files)){
        # name the first column of each data frame by its cohort name
        names(files[[i]])[1] <- cohorts[i]
        res[[i]] <- lnc_gsea(tid_cohort = files[[i]],geneset)
        res[[i]]$cohort <- cohorts[i]
    }
   print(length(res))

   res.df <- data.table::rbindlist(res)
   res.df <- res.df[, c("pathway","padj","NES","cohort")]

   res.df$Status <- ifelse(res.df$NES > 0, "Upregulated", ifelse(res.df$NES < 0, "Downregulated", ""))
   res.df$Status <- factor(res.df$Status, levels = c("Upregulated", "Downregulated", ""))

   res.df$FDR <- ifelse(res.df$padj > 0.05, "> 0.05", ifelse(res.df$padj < 0.001, "< 0.001", "< 0.05"))
   res.df$FDR <- factor(res.df$FDR, levels = c("> 0.05", "< 0.05", "< 0.001"))

   return(res.df)
}




#' plot enriched pathways for one lncRNA in multiple cohorts
#'
#' Visualize the enriched pathways' results in multi-cohorts for interested lncRNA
#'
#' @param gseares A list of data frames which stores interested lncRNA's enriched pathway results for multiple cohorts
#' @import ggplot2
#'
#' @return A heatmap plot
#' @example
#' cohorts <- c("PRAD","LUAD","BRCA")
#' arlnc1.df <- pre_compareCohort(lncRNA="ENST00000561519", cohorts = cohorts)
#' plot_multiCompare(arlnc1.df)
#' @export


plot_multiCompare <- function(gseares) {
     # create plot
    p <- ggplot2::ggplot(data = gseares, mapping = ggplot2::aes(x = cohort,
                                        y =  pathway,
                                        colour = Status,
                                        shape = Status,
                                        fill = NES,
                                        size = FDR)) +
        geom_point() +
        scale_shape_manual(values = c(24,25,1)) +
        scale_fill_gradient2(high = "coral3", low = "deepskyblue3", mid = "white") +
        scale_color_manual(values = c("coral3","deepskyblue3","white")) +
        scale_size_manual(values=c(0,3,4)) +
        theme_bw() +
        guides(size = guide_legend(override.aes = list(shape=17))) +
        theme(plot.title = element_text(hjust = 0.5),
              axis.title.y =  element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size=8, angle = 45, hjust = 1, face = "plain", color='black'),
              axis.text.y = element_text(face = "plain",size=8, color='black'))

    p

}
