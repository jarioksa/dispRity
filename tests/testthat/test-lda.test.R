context("lda.test")

## Factorise subsets
test_that("factorise.subsets works", {

    ## Simple factorisation return
    ordinated_matrix <- matrix(data = rnorm(90), nrow = 10)
    expect_warning(data <- custom.subsets(ordinated_matrix, list(c(1:4), c(5:10))))
    test <- factorise.subsets(data)
    expect_is(test, "list")
    expect_equal(unlist(lapply(test, length)), 10)
    expect_equal(unname(unlist(lapply(test, unique))), c("1", "2"))

    ## Double elements factorisation return
    require(geomorph)
    data(plethodon)
    procrustes <- geomorph::gpagen(plethodon$land)
    geomorph_df <- geomorph::geomorph.data.frame(procrustes, species = as.factor(lda_test_data$species), morpho = as.factor(lda_test_data$morpho))
    data <- geomorph.ordination(geomorph_df, ordinate = FALSE)

    test <- factorise.subsets(data)
    expect_is(test, "list")
    expect_equal(names(test), c("species", "morpho"))
    expect_equal(unlist(lapply(test, length)), c("species" = 40, "morpho" = 40))
    expect_equal(unname(unlist(lapply(test, unique))), c("Jord", "Teyah", "group1", "group2", "group3"))


    ## Complex factorisation (with NAs)
    ordinated_matrix <- matrix(data = rnorm(90), nrow = 10,
         dimnames = list(letters[1:10]))
    data <- custom.subsets(ordinated_matrix,
         list("A" = c("a", "b", "c", "d"), "B" = c("e", "f", "g", "h", "i", "j"),
              "C" = c("a", "c", "d", "f", "h")))

    test <- factorise.subsets(data)
    expect_is(test, "list")
    expect_equal(names(test), c("A", "B", "C"))
    expect_equal(unlist(lapply(test, length)), c("A" = 10, "B" = 10, "C" = 10))
    expect_equal(unname(unlist(lapply(test, unique))), c(1, 0, 0, 1, 1, 0))
})



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
