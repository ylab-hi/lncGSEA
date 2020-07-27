context("test lnc_gsea results")
library(lncGSEA)

test_that("First column name is pathway and ncol is 8", {
    df.test <- pre_gsea("PRAD", "T136792")
    lnc.test <- lnc_gsea(df.test, metric = "cor", cor.method = "pearson")

    expect_equal(names(lnc.test)[1], "pathway")
    expect_equal(ncol(lnc.test), 8)
})
