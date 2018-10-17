context("lda.test")

load("lda_test_data.Rda")

## Factorise subsets
test_that("factorise.subsets works", {

    ## Simple factorisation return
    ordinated_matrix <- matrix(data = rnorm(90), nrow = 10)
    expect_warning(data <- custom.subsets(ordinated_matrix, list(c(1:4), c(5:10))))
    test <- factorise.subsets(data)
    expect_is(test, "list")
    expect_equal(unlist(lapply(test, length)), 10)
    expect_equal(unname(unlist(lapply(test, unique))), as.factor(c("1", "2")))

    ## Double elements factorisation return
    require(geomorph)
    data(plethodon)
    procrustes <- geomorph::gpagen(plethodon$land, print.progress = FALSE)
    geomorph_df <- geomorph::geomorph.data.frame(procrustes, species = as.factor(lda_test_data$species), morpho = as.factor(lda_test_data$morpho))
    data <- geomorph.ordination(geomorph_df, ordinate = FALSE)

    test <- factorise.subsets(data)
    expect_is(test, "list")
    expect_equal(names(test), c("species", "morpho"))
    expect_equal(unlist(lapply(test, length)), c("species" = 40, "morpho" = 40))
    expect_equal(as.vector(unname(unlist(lapply(test, unique)))), c("Jord", "Teyah", "group1", "group2", "group3"))


    ## Complex factorisation (with NAs)
    ordinated_matrix <- matrix(data = rnorm(90), nrow = 10, dimnames = list(letters[1:10]))
    data <- custom.subsets(ordinated_matrix,
         list("A" = c("a", "b", "c", "d"), "B" = c("e", "f", "g", "h", "i", "j"),
              "C" = c("a", "c", "d", "f", "h")))

    test <- factorise.subsets(data)
    expect_is(test, "list")
    expect_equal(names(test), c("A", "B", "C"))
    expect_equal(unlist(lapply(test, length)), c("A" = 10, "B" = 10, "C" = 10))
    expect_equal(unname(unlist(lapply(test, unique))), as.factor(c(1, 0, 0, 1, 1, 0)))
})

## Run one lda
test_that("run.one.lda works", {
    data <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]), species = rep(c("s","c","v"), rep(50,3)))
    data_matrix <- data[, -ncol(data)]
    factor <- as.factor(data[, ncol(data)])

    set.seed(1)
    test <- run.one.lda(factor, data_matrix, prior = FALSE, train = 50, CV = FALSE, fun.type = MASS::lda)
    expect_is(test, "list")
    expect_equal(names(test), c("fit", "predict", "training"))
    expect_is(test$fit, "lda")

})

## Test
test_that("lda.test works", {

    ## Plethodon test
    require(geomorph)
    data(plethodon)
    procrustes <- geomorph::gpagen(plethodon$land, print.progress = FALSE)
    geomorph_df <- geomorph.data.frame(procrustes, species = as.factor(lda_test_data$species), morpho = as.factor(lda_test_data$morpho))
    data_disparity <- geomorph.ordination(geomorph_df, ordinate = FALSE)

    ## Iris test
    data_df <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]), species = rep(c("s","c","v"), rep(50,3)))
    prior_iris <- rep(1/3, 3)


    ## Sanitizing
    expect_error(lda.test(data = "no_matrix", train = 10))
    wrong_matrix <- lda_test_data$procrustes
    expect_error(lda.test(data = wrong_matrix, train = 10))
    wrong_disparity <- dispRity(lda_test_data$procrustes, metric = mean)
    expect_error(lda.test(data = wrong_disparity, train = 10))
    expect_error(lda.test(data = data_df, train = "151"))
    expect_error(lda.test(data = data_df, train = 151))
    expect_error(lda.test(data = data_df, train = 10, prior = "a"))
    expect_error(lda.test(data = data_df, train = 10, prior = c(0.5, 0.5)))
    expect_error(lda.test(data = data_disparity, train = 10, prior = c(0.5, 0.5)))
    expect_error(lda.test(data = data_disparity, train = 10, prior = list(c(0.5, 0.5), c(1))))
    expect_error(lda.test(data = data_df, train = 10, type = "a"))
    expect_warning(expect_error(lda.test(data = data_df, train = 10, type = mean)))
    expect_error(lda.test(data = data_df, train = 10, bootstrap = "a"))
    expect_error(lda.test(data = data_df, train = 10, bootstrap = c(1,2,3)))
    expect_error(lda.test(data = data_df, train = 10, CV = "x"))
    expect_error(lda.test(data = data_df, train = 10, LASSO = "x"))
    expect_error(lda.test(data = data_df, train = 10, PGLS = "x"))

    ## Running the two examples
    test <- lda.test(data_df, train = 50, bootstraps = 7)
    expect_is(test, c("dispRity", "lda-test"))
    expect_equal(names(test), c("species", "call"))
    expect_equal(length(test[[1]]), 7)
    expect_equal(names(test[[1]][[1]]), c("fit", "predict", "training"))

    expect_warning(test <- lda.test(data_disparity, train = 10, bootstraps = 5))
    expect_is(test, c("dispRity", "lda-test"))
    expect_equal(names(test), c("species", "morpho", "call"))
    expect_equal(length(test[[1]]), 5)
    expect_equal(names(test[[1]][[1]]), c("fit", "predict", "training"))    

})
