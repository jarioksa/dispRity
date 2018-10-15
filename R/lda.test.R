#' @title Run a LDA or CBA test
#'
#' @description Running a linear discriminant analysis or canonical variable analysis test on a \code{dispRity} object
#'
#' @param data A \code{dispRity} object with attributed subsets or a \code{matrix} with factors as the last column.
#' @param subset The size of the training dataset (by default the size is @@@)
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

# stop("DEBUG test.lda")
# load("../tests/testthat/lda_test_data.Rda")
# source("sanitizing.R")

test.lda <- function(data, subset, prior, bootstrap, CV = FALSE, LASSO = FALSE, PGLS = FALSE, ...) {

    match_call <- match.call()

    ## Sanitizing
    ## Data
    data_class <- check.class(data, c("dispRity", "matrix", "data.frame"))
    if(data_class == "matrix") {
        ## Check if the last column is a factor
        last_col <- ncol(data)
        test_factor <- as.factor(data[, last_col])
        if(length(levels(test_factor)) == nrow(data)) {
            stop("")
        }
    } else {
        ## Check if the dispRity object has 

    }


    #TODO: Add the option to update models

    ## First we select a subset of the dataset for training
    subset <- sample(1:nrow(data), subset)
    training <- data[subset, ncol(data)]

    ## Get the number classes from the data
    classes <- unique(data[, ncol(data)])

    ## Second we set the prior
    if(missing(prior)) {
        prior <- as.numeric((table(training)/sum(table(training))))
    }

    ## Fitting the LDA
    lda_fit <- MASS::lda(x = data[, -ncol(data)], grouping = data[, ncol(data)],
                         prior = prior, subset = subset, CV = CV, ...)

    ## Predicting the fit
    lda_predict <- predict(lda_fit, data[-subset, -ncol(data)])

    ## Return the predictions results and the fit
    return(list("fit" = lda_fit, "predict" = lda_predict, "training" = subset, "data" = data))

}