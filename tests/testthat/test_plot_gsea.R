context("test plot_gsea results")

library(lncGSEA)

test_that("plot_gsea works", {
    df.test  <- pre_gsea("PRAD", "T136792", "/Users/yren/Projects/lncGSEA/inst/extdata/")
    lnc.test <- lnc_gsea(df.test, metric = "cor", cor.method = "pearson")
    pp <- plot_gsea("T136792_PRAD_cor.txt", direction = "both")
    tmp <- tempfile("plot", fileext = ".png")
    ggsave(tmp, plot=pp)
    expect_true(TRUE)
})
