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


test_that("extract.lda.test works" ,{

    lda_test <- lda_test_data$lda_test

    ## Works with sd
    test <- extract.lda.test(lda_test, what = "prior", where = "fit", deviation = sd, cent.tend = mean)

    expect_is(test, "list")
    expect_equal(names(test), c("cent.tend", "deviation"))
    expect_equal(names(test[[1]]), c("species", "morpho"))

    expect_equal(test$cent.tend$species, c("Jord" = 0.5, "Teyah" = 0.5))
    expect_equal(round(test$deviation$species, digits = 7), c("Jord" = 0.2645751, "Teyah" = 0.2645751))
    expect_equal(round(test$cent.tend$morpho, digits = 7), c("group1" = 0.2666667, "group2" = 0.3333333, "group3" = 0.4000000))
    expect_equal(round(test$deviation$morpho, digits = 8),c("group1" = 0.20816660, "group2" = 0.05773503, "group3" = 0.17320508))

    ## Works with means
    test <- extract.lda.test(lda_test, what = "x", where = "predict", deviation = c(0.05, 0.95), cent.tend = mean)

    expect_is(test, "list")
    expect_equal(names(test), c("cent.tend", "deviation"))
    expect_equal(names(test[[1]]), c("species", "morpho"))

    expect_equal(length(test$cent.tend$species), 30)
    expect_equal(dim(test$deviation$species), c(2,30))
    expect_equal(length(test$cent.tend$morpho), 30)
    expect_equal(dim(test$deviation$morpho), c(2,30))

    ## Works without arguments
    test <- extract.lda.test(lda_test, what = "class", where = "predict")
    expect_is(test, "list")
    expect_equal(names(test), c("species", "morpho"))
    expect_equal(levels(unlist(c(test[[1]]))), c("Jord", "Teyah"))
    expect_equal(levels(unlist(c(test[[2]]))), c("group1", "group2", "group3"))

    ## Works without what or where
    test1 <- extract.lda.test(lda_test, where = "training")
    test2 <- extract.lda.test(lda_test, what = "training")

    expect_equal(names(test1), names(test2))
    expect_equal(test1[[1]], test2[[1]])
    expect_equal(test1[[2]], test2[[2]])

    ## Works with just the cent.tend
    test1 <- extract.lda.test(lda_test, what = "means", where = "fit", cent.tend = sd)
    expect_equal(round(test1$cent.tend$species, digits = 7), c("Jord" = 0.2056696, "Teyah" = 0.2056992))
    expect_equal(round(test1$cent.tend$morpho, digits = 7), c("group1" = 0.2060365, "group2" = 0.2060382, "group3" = 0.2057907))

    ## Works just with the deviation
    test2 <- extract.lda.test(lda_test, what = "means", where = "fit", deviation = sd)
    expect_equal(test1$cent.tend$species, test2$deviation$species)
    expect_equal(test1$cent.tend$morpho, test2$deviation$morpho)
})


test_that("accuracy.score works" ,{

    lda_test <- lda_test_data$lda_test
    prediction <- lda_test[[1]][[1]]$predict$class
    un_trained <- lda_test$support$factors[[1]][-lda_test[[1]][[1]]$training]

    ## Works by default
    expect_equal(round(accuracy.score(prediction, un_trained, return.table = FALSE), 7), 0.5666667)
    expect_is(accuracy.score(prediction, un_trained, return.table = TRUE), "table")
    expect_equal(accuracy.score(prediction, un_trained, return.table = TRUE), table(prediction, un_trained))
})



test_that("apply.accuracy.score works" ,{

    lda_test <- lda_test_data$lda_test

    ## Works by default
    test1 <- apply.accuracy.score(lda_test)
    expect_equal(names(test1), names(lda_test)[1:2])
    expect_equal(unique(unlist(lapply(test1, length))), 3)
    expect_equal(round(test1[[1]], 7), c("1" = 0.5666667, "2" = 0.6666667, "3" =0.7000000))
    expect_equal(round(test1[[2]], 7), c("1" = 0.6666667, "2" = 0.7666667, "3" =0.6666667))

    test1 <- apply.accuracy.score(lda_test, return.table = TRUE)
    expect_equal(names(test1), names(lda_test)[1:2])
    expect_equal(lapply(test1, length), list("species" = 3, "morpho" = 3))
    expect_equal(unique(unlist(lapply(test1, lapply, class))), "table")
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
    expect_error(lda.test(data = data_df, train = c(1,2)))
    expect_error(lda.test(data = data_disparity, train = c(10, 1000)))
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
    expect_equal(names(test), c("species", "support"))
    expect_equal(length(test[[1]]), 7)
    expect_equal(names(test[[1]][[1]]), c("fit", "predict", "training"))

    expect_warning(test <- lda.test(data_disparity, train = 10, bootstraps = 5))
    expect_is(test, c("dispRity", "lda-test"))
    expect_equal(names(test), c("species", "morpho", "support"))
    expect_equal(length(test[[1]]), 5)
    expect_equal(names(test[[1]][[1]]), c("fit", "predict", "training"))    

})
