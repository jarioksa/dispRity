#' @title Run a LDA or CBA test
#'
#' @description Running a linear discriminant analysis or canonical variable analysis test on a \code{dispRity} object
#'
#' @param data A \code{dispRity} object with attributed subsets or a \code{matrix} with factors as the last column.
#' @param train The size of the training dataset (by default the size is @@@)
#TODO: get a default subset size selector? Maybe something like the Silverman's rule? 
#' @param prior A \code{vector} of prior expectations (if missing the prior is set to the observed proportions in the training set).
#' @param bootstrap Optional, the number of bootstrap replicates to run (is missing, no bootstraps are run).
#' @param CV Logical, whether to perform cross-validation or not (default is \code{FALSE}).
#' @param LASSO Logical, whether to perform a LASSO analysis.
#' @param PGLS Logical, whether to perform a PGLS analysis.
#' @param ... Optional arguments to be passed to \code{\link[MASS]{lda}}.
#' 
#' @examples
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

# MASS, lda

# stop("DEBUG lda.test")
# source("lda.test_fun.R")
# library(geomorph)
# load("../tests/testthat/lda_test_data.Rda")
# source("sanitizing.R")
# data(plethodon)
# procrustes <- geomorph::gpagen(plethodon$land)
# geomorph_df <- geomorph.data.frame(procrustes, species = as.factor(lda_test_data$species), morpho = as.factor(lda_test_data$morpho))
# data <- geomorph.ordination(geomorph_df, ordinate = FALSE)

## Iris training
# data <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]), species = rep(c("s","c","v"), rep(50,3)))
# prior <- 


lda.test <- function(data, train, prior, bootstrap, CV = FALSE, LASSO = FALSE, PGLS = FALSE, ...) {

    match_call <- match.call()

    ## Sanitizing
    ## Data
    data_class <- check.class(data, c("dispRity", "matrix", "data.frame"))
    if(data_class != "dispRity") {
        ## Check if the last column is a factor
        last_col <- ncol(data)
        test_factor <- as.factor(data[, last_col])
        if(length(levels(test_factor)) == nrow(data)) {
            stop.call(msg.pre = "The last column of ", call = match_call$data, msg = " does not contain factors.")
        }

        ## Convert into the table and factors
        data_matrix <- data[, -last_col]
        factors <- test_factor

    } else {
        ## Check if the dispRity object has subsets
        if(is.null(data$subsets) || length(data$subsets) == 1) {
            stop.call(match_call$data, " must contain subsets. Use the custom.subsets() function.")
        }

        ## Convert into the data table + factors
        factors <- factorise.subsets(data)
        data_matrix <- data$matrix
    }

    ## Subsets
    check.class(train, c("numeric", "integer"))
    check.length(train, 1, " must a single numeric value the number of elements to use for training.")
    if(train > nrow(data_matrix)) {
        stop.call(msg.pre = paste0("The training set size (", train ,") cannot be bigger than the number of elements in "), call = match_call$data, msg = paste0(" (", nrow(data_matrix), ")."))
    }

    ## Prior
    if(missing(prior)) {
        ## Calculate the prior later one
        prior <- FALSE
    } else {
        ## If only one factor is test
        if(class(factors) != "list") {
            ## Right class
            check.class(prior, c("numeric", "integer"))
            ## Right length
            check.length(prior, length(levels(factors)), " must be equal to the number of levels to test.")
        } else {
            ## Right list
            check.class(prior, "list", " must be a list of the same length of the number of categories to test.")
            check.length(prior, length(factors)," must be a list of the same length of the number of categories to test.")

            ## Right elements
            test_class <- unlist(lapply(prior, class))
            test_class <- test_class == "numeric" | test_class == "integer"
            if(any(!test_class)) {
                stop.call("prior", paste0(" ", which(!test_class), " must be numeric."))
            }

            ## Right lengths
            test_length <- mapply(function(x,y) return(unname(x) == unname(y)), lapply(prior, length), lapply(factors, function(x) length(levels(x))))
            if(any(!test_length)) {
                stop.call("prior", paste0(" ", which(!test_length), " is not the right length."))
            }
        }
    }

    ## Bootstrap
    if(!missing(bootstrap)) {
        ## Do bootstrap
        do_bootstrap <- TRUE
        check.class(bootstrap, c("integer", "numeric"))
        check.length(bootstrap, 1, " must be a single numeric value.")
    } else {
        ## Don't bootstrap
        do_bootstrap <- FALSE
    }

    ## CV, LASSO and PGLS
    check.class(CV, "logical")
    check.class(LASSO, "logical")
    check.class(PGLS, "logical")

    ## Run the LDA
    lda_out <- lapply(factors, run.one.LDA, data_matrix, prior = prior, train = train, CV = CV)

    run.one.LDA(factors[[1]], data_matrix, prior = prior, train = train, CV = CV)

    ## Handle colinear warnings



    # #TODO: Add the option to update models



    # ## If estimate prior
    # prior <- as.numeric((table(training)/sum(table(training))))

    # ## First we select a subset of the dataset for training
    # subset <- sample(1:nrow(data), subset)
    # training <- data[subset, ncol(data)]

    # ## Get the number classes from the data
    # classes <- unique(data[, ncol(data)])

    # ## Second we set the prior


    # ## Fitting the LDA
    # lda_fit <- MASS::lda(x = data[, -ncol(data)], grouping = data[, ncol(data)],
    #                      prior = prior, subset = subset, CV = CV, ...)

    # ## Predicting the fit
    # lda_predict <- predict(lda_fit, data[-subset, -ncol(data)])

    # ## Return the predictions results and the fit
    # return(list("fit" = lda_fit, "predict" = lda_predict, "training" = subset, "data" = data))

}