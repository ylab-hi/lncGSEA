#' Filter out low expression genes
#'
#' Filter out low expression genes by mean of FPKM smaller than an cutoff
#'
#'@param df A data frame, rows are genes, columns are samples
#'@param col A integer showing column number, df[[col]] is not considered for filtering
#'@param cutoff A numeric value as expression threshold
#'@return A data frame with genes' expression larger than cutoff
#'@example
#'cohort <- filter_mean(cohort,1,1)
#'@export

filter_mean <- function(df, col, cutoff){
    df <- df[rowMeans(df[,-col]) >= cutoff, ]
    return(df)
}


#'Set a column as row names
#'
#'Set a column as rownames and remove that column in the data frame
#'
#'@param df A data frame
#'@param col A integer, which column chosen to be as row names
#'@return A data frame
#'@example
#'cohort <- col2rownames(cohort,1)
#'@export

col2rownames <- function(df, col){
    rownames(df) <- df[[col]]
    df[[col]] <- NULL
    return(df)
}



#'Transpose an expression data frame or matrix
#'
#'Transpose an expression data frame or matrix and add rownames as its first column 'patient_id'
#'
#'@param df A data frame
#'@return A transposed data frame
#'
#'@example
#'cohortT <- transdf(cohort)
#'
#'@export

transdf <- function(df){
    datT <- as.data.frame(t(df))
    datT <- data.frame(patient_id = row.names(datT), datT)
}



#' Function to join customized lncRNA expression with other genes' expression
#'
#' Join lncRNA and other genes' expression, prepare a data frame ready for lnc_gsea
#'
#'@param lncRNA A character string showing the file which stores a lncRNA expression provided by users
#'@param cohort.file A character string, path to the file which stores other genes' expression
#'@import dplyr
#'@return A data frame ready for lnc_gsea
#'@example
#'lnctest <- custom_lnc("igf1ras1.custom.txt","../PRAD.FPKM.txt")
#'
#'@export

custom_lnc <- function(lncRNA, cohort.file){
    lnc <- read.delim(lncRNA, sep = "\t", stringsAsFactors = FALSE,
                      header = TRUE, check.names = FALSE)
    cohort <- as.data.frame(data.table::fread(cohort.file))

    # filter genes whose mean FPKM <1
    cohort <- filter_mean(cohort, 1, 1)
    # set first column as rownames
    cohort <- col2rownames(cohort, 1)
    # transpose
    cohortT <- transdf(cohort)
    novelT <- transdf(novel)

    # join lncRNA and other genes expression matrix
    novel.cohort <- dplyr::left_join(novelT, cohortT, by = "patient_id")
    novel.cohort <- na.omit(novel.cohort)

    return(novel.cohort)
}



