#' process expression data for sg_GSEA function
#'
#' Find TCGA samples appear both in RNAseq expression (FPKM) data frame and lncRNAs
#' data frame, if the sample ids are different, use tcga.meta.info to match.
#'
#' @param cohort A character string which is the short name for cohorts in TCGA
#' @param t_id A character string of transcript id, either start with ENST or T
#'
#' @return A data frame of lncRNA's transcript and other genes' expression in interested cohort
#'
#' @import data.table
#' @import stringr
#' @import dplyr
#'
#' @details cohort should be one of studies in TCGA, for example, "PRAD" can be used as
#'          the function's input to stand for prostate cancer; t_id can be from mitranscriptome
#'          database or reflnc database.
#'
#' @example
#' pre_gsea("PRAD", "T280760")
#' pre_gsea("BRCA","ENST00000430998")
#'
#' @export
pre_gsea <- function(cohort, t_id){

    #dataDir <- path.expand("./data")

    # interested cohort ---
    cohort <- toupper(cohort)
    cohort.file <- paste0("./data/", cohort, ".FPKM.txt")

    # for test -----
    #cohort.file <- system.file("extdata", paste0(cohort,".FPKM.txt.gz"), package= "lncGSEA")

    if (file.exists(cohort.file)) {
    cancer <- as.data.frame(data.table::fread(cohort.file))
    } else {stop(cohort.file, " does not exist!")}

    # filter out low expressed genes -----
    cancer <- cancer[rowMeans(cancer[,-1])>=1, ]
    rownames(cancer) <- cancer$GeneID
    cancer$GeneID <- NULL


    # transpose ---
    cancerT <- as.data.frame(t(cancer))
    cancerT <- data.frame(patient_id = row.names(cancerT), cancerT)
    # interested t_id ---
    reflnc <- data.table::fread(system.file("extdata", "RefLnc_lncRNA_tumor_sample_FPKM_tid.gz",package = "lncGSEA"))
    mich <- data.table::fread(system.file("extdata", "mitranscriptome.expr.fpkm_tid.tsv.gz",package = "lncGSEA"))

    if (any(grepl(t_id, reflnc[[1]]))) {
        dataFiles <- "./data/RefLnc_lncRNA_tumor_sample_FPKM.gz"
        # for test ----
        # dataFiles <- system.file("extdata", "RefLnc_lncRNA_tumor_sample_FPKM.gz", package = "lncGSEA")
        meta <- data.table::fread(system.file("extdata", "tcga.meta.file.txt.gz", package = "lncGSEA"))

    } else if (any(grepl(t_id, mich[[1]]))) {
        # for test ---
        # dataFiles <- system.file("extdata", "mitranscriptome.expr.fpkm.tsv.gz", package = "lncGSEA")
        dataFiles <- "./data/mitranscriptome.expr.fpkm.tsv.gz"
        meta <- data.table::fread(system.file("extdata", "library_info.txt.gz", package = "lncGSEA"))

    } else { message("Can not find the input t_id in the package database.
                     Users should provide their lncRNAs expression matrix for
                     TCGA samples.")}

    # create data.table with header ---
    dt <- data.table::fread(dataFiles)
    names(dt)[1] <- "transcript_id"
    #dt <- dt[dt$transcript_id==t_id, ] # change enst00098899.1 to enst00098899
    dt <- dt[grep(t_id, dt$transcript_id), ]
    if (nrow(dt) > 1){
        dt <- unique(dt)
    }

    # transpose ---
    dt <- as.data.frame(dt)
    rownames(dt) <- dt[,1]
    dt[,1] <- NULL
    t_dt <- as.data.frame(t(dt))

    # test if sample name end with "gdc_realn_rehead"
    if (any(grepl("gdc_realn_rehead", rownames(t_dt)))){
        t_dt$library_id <- stringr::str_split(rownames(t_dt),"_",simplify = TRUE)[,1]
        t_dt.join <- dplyr::left_join(t_dt, meta, by = "library_id")
        t_dt.join.sub <- t_dt.join[,c(3,1)]
        tid_cohort <- na.omit(dplyr::full_join(t_dt.join.sub, cancerT,
                                               by = c("ALIQUOT_BARCODE" = "patient_id")))
    } else {
        t_dt$library_id <- rownames(t_dt)
        t_dt.join <- dplyr::left_join(t_dt, meta, by = "library_id")
        #t_dt.join.sub <- t_dt.join[,c("tcga_legacy_sample_id",t_id)]
        t_dt.join.sub <- t_dt.join[,c(13,1)]
        tid_cohort <- na.omit(dplyr::full_join(t_dt.join.sub, cancerT,
                                               by = c("tcga_legacy_sample_id" = "patient_id")))
    }

    names(tid_cohort)[1] <- cohort
    if (nrow(tid_cohort) > 0){
        return(tid_cohort)
    } else {message("There are 0 rows in the dataframe, please try another cohort.")}


}
