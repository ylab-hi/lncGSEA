context("test multiple lncRNAs/cohort results")

library(lncGSEA)

test_that("First column name is pathway and ncol is 6", {
    df.test1 <- pre_gsea("PRAD", "T136792", "/Users/yren/Projects/lncGSEA/data/")
    lnc.test1 <- lnc_gsea(df.test1, metric = "cor", cor.method = "pearson")

    df.test2 <- pre_gsea("PRAD", "T136804", "/Users/yren/Projects/lncGSEA/data/")
    lnc.test2 <- lnc_gsea(df.test2, metric = "cor", cor.method = "pearson")

    prad <- pre_multiCompare(files = list("T136792_PRAD_cor.txt",
                                          "T136804_PRAD_cor.txt"),
                             compare = "lncRNA")
    expect_equal(names(prad)[1], "pathway")
    expect_equal(ncol(prad), 6)
    expect_equal(prad[1,4], "T136792")
})



test_that("plot_multiCompare works", {
    df.test1 <- pre_gsea("PRAD", "T136792", "/Users/yren/Projects/lncGSEA/data/")
    lnc.test1 <- lnc_gsea(df.test1, metric = "cor",cor.method = "pearson")

    df.test2 <- pre_gsea("PRAD", "T136804", "/Users/yren/Projects/lncGSEA/data/")
    lnc.test2 <- lnc_gsea(df.test2, metric = "cor",cor.method = "pearson")

    prad <- pre_multiCompare(files = list("T136792_PRAD_cor.txt",
                                          "T136804_PRAD_cor.txt"),
                             compare = "lncRNA")

    pp <- plot_multiCompare(prad)
    tmp <- tempfile("plot", fileext = ".png")
    ggsave(tmp, plot=pp)
    expect_true(TRUE)
})
