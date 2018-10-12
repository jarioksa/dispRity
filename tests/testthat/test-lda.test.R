context("lda.test")

## Test
test_that("lda.test works", {

    ## Sanitizing
    expect_error(lda.test())

    ## Right output
    expect_is(
        lda.test()
        , "class")
    expect_equal(
        dim(lda.test())
        , dim)
})
