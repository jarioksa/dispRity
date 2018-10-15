context("lda.test")

## Test
test_that("lda.test works", {

    load("lda_test_data.Rda")

    ## Sanitizing
    expect_error(lda.test(data = "no_matrix", subset = 10))
    wrong_matrix <- lda_test_data$procrustes
    expect_error(lda.test(data = wrong_matrix, subset = 10))
    wrong_disparity <- dispRity(lda_test_data$procrustes, metric = mean)
    expect_error(lda.test(data = wrong_disparity, subset = 10))

    ## Right output
    expect_is(
        lda.test()
        , "class")
    expect_equal(
        dim(lda.test())
        , dim)
})
