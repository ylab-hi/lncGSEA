context("test pre_gsea output")

library(lncGSEA)

test_that("First column name is PRAD and Second column is T136792", {
    df.test <- pre_gsea("PRAD", "T136792")

    expect_equal(names(df.test)[1], "PRAD")
    expect_equal(names(df.test)[2], "T136792")
})