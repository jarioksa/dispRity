#' @title Run a discriminant analysis type
#'
#' @description Running a discriminant analysis or canonical variable analysis test on a \code{dispRity} object
#'
#' @param data A \code{dispRity} object with attributed subsets or a \code{matrix} with factors as the last column.
#' @param train The size of the training dataset. This can be either a single value (to apply to all factors) or a vector of values equal to the number of factors.

#TODO: get a default subset size selector? Maybe something like the Silverman's rule? 

#' @param prior A \code{vector} of prior expectations (if missing the prior is set to the observed proportions in the training set).
#' @param type The type of discriminant function to use: either \code{"linear"} (default), \code{"quadratic"} or any function of similar to the generic type \code{\link[MASS]{lda}}.
#' @param bootstraps Optional, the number of bootstraps replicates to run (is missing, no bootstrapss are run).
#' @param CV Logical, whether to perform cross-validation or not (default is \code{FALSE}).
#' @param LASSO Logical, whether to perform a LASSO analysis.
#' @param PGLS Logical, whether to perform a PGLS analysis.
#' @param all.levels Logical, whether to always include all levels in the training dataset (\code{TRUE} - default) or not (\code{FALSE}). Not including speeds up the calculation but can generate errors when apply the discriminant function.
#' @param ... Optional arguments to be passed to the discriminant function.

# TODO: add a parallel option (to be applied at the wrapper level of run.multi.lda)]]

#' @return
#' A list of class \code{"dispRity"} and \code{"lda.test"} containing the results for each tested factor(s) and each bootstraps along with a \code{$support} element containing meta-data.
#' 
#' For example, for two tested factors with three bootstrap replicates each:
#' 
#' \code{lda_out$factor1$bootstraps1}
#' \code{               $bootstraps2}
#' \code{               $bootstraps3}
#' \code{       $factor2$bootstraps1}
#' \code{               $bootstraps2}
#' \code{               $bootstraps3}
#' \code{       $support}
#' 
#' @examples
#' ## The classic iris example:
#' 
#' ## A data.frame with the factor as the last column
#' iris_data <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),
#'                         species = rep(c("s","c","v"), rep(50,3)))
#' 
#' ## The lda test with a training sample size of 50 specimen
#' iris_test <- lda.test(iris_data, train = 50)
#' 
#' ## Printing the results
#' iris_test
#' 
#' ## Summarising and plotting the results
#' summary(iris_test)
#' plot(iris_test)
#' 
#' \dontrun{
#' ## A geomorph example
#' set.seed(1)
#' require(geomorph)
#'
#' ## Loading the plethodon dataset
#' data(plethodon)
#' 
#' ## Running a Procrustes superimposition
#' procrustes <- geomorph::gpagen(plethodon$land)
#' 
#' ## Creating two factors: the morpho groups and a random group
#' morpho_group <- as.factor(c(rep("group1", 10), rep("group2", 10), rep("group3", 20)))
#' random_group <- as.factor(sample(c("random1", "random2", "random3"), 40, replace = TRUE))
#' 
#' ## Combining them into a geomorph.data.frame object
#' geomorph_df <- geomorph.data.frame(procrustes, random = random_group, morpho = morpho_group)
#' 
#' ## Converting it into a dispRity object without ordination
#' disparity_data <- geomorph.ordination(geomorph_df, ordinate = FALSE)
#' 
#' ## Applying a simple lda test with a training sample size of 10 specimen bootstrapped 100 times
#' plethodon_test <- lda.test(disparity_data, train = 10, bootstraps = 100)
#' 
#' ## Printing the results
#' plethodon_test
#' 
#' ## Summarising and plotting the results
#' summary(plethodon_test)
#' plot(plethodon_test)
#' }
#'
#' @seealso \code{\link{dispRity}}, \code{\link{test.dispRity}}, \code{\link{summary.dispRity}}, \code{\link{plot.dispRity}}.
#' 
#' @author Thomas Guillerme
#' @export

# stop("DEBUG lda.test")
# source("lda.test_fun.R")
# library(geomorph)
# library(dispRity)
# load("../tests/testthat/lda_test_data.Rda")
# source("sanitizing.R")
# data(plethodon)
# procrustes <- geomorph::gpagen(plethodon$land)
# geomorph_df <- geomorph.data.frame(procrustes, species = as.factor(lda_test_data$species), morpho = as.factor(lda_test_data$morpho))
# data <- geomorph.ordination(geomorph_df, ordinate = FALSE)
# lda_test <- lda.test(data, train = 10, bootstraps = 3)

## Iris training
# data <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]), species = rep(c("s","c","v"), rep(50,3))) 

# fun.type = "linear"
# CV = FALSE
# LASSO = FALSE
# PGLS = FALSE
# bootstraps = 4

lda.test <- function(data, train, prior, type = "linear", bootstraps, CV = FALSE, LASSO = FALSE, PGLS = FALSE, all.levels = TRUE, ...) {

    match_call <- match.call()

    ## Sanitizing
    ## Data
    data_class <- check.class(data, c("dispRity", "matrix", "data.frame"))
    if(data_class != "dispRity") {
        ## Check if the last column is a factor
        last_col <- ncol(data)
        test_factor <- list(as.factor(data[, last_col]))
        if(length(levels(test_factor[[1]])) == nrow(data)) {
            stop.call(msg.pre = "The last column of ", call = match_call$data, msg = " does not contain factors.")
        }

        ## Convert into the table and factors
        data_matrix <- data[, -last_col]
        factors <- test_factor
        names(factors) <- colnames(data)[last_col]

    } else {
        ## Check if the dispRity object has subsets
        if(is.null(data$subsets) || length(data$subsets) == 1) {
            stop.call(match_call$data, " must contain subsets. Use the custom.subsets() function.")
        }

        ## Convert into the data table + factors
        factors <- factorise.subsets(data)
        if(class(factors) != "list") {
            factors <- list(factors)
        }
        data_matrix <- data$matrix
    }

    ## Train
    check.class(train, c("numeric", "integer"))

    if(length(train) > 1) {
        length_required <- length(factors)
        check.length(train, length_required, paste0(" must be either a single numeric value to apply to all factors or a vector of numeric values equal to the number of factors (", length_required, ")."), errorif = FALSE)
    } 
    if(any(train > nrow(data_matrix))) {
        stop.call(msg.pre = paste0("The training set size cannot be bigger than the number of elements in "), call = match_call$data, msg = paste0(" (", nrow(data_matrix), ")."))
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

    ## Type
    type_class <- check.class(type, c("character", "function"))
    if(type_class == "character") {
        allowed_classes <- c("linear", "quadratic")
        check.method(type, allowed_classes, "type argument must be a function or")

        if(type == "linear") {
            fun.type <- MASS::lda
        }
        if(type == "quadratic") {
            fun.type <- MASS::qda
        }
    } else {
        ## User defined function
        fun.type <- type
    }

    ## Bootstrap
    if(!missing(bootstraps)) {
        ## Do bootstraps
        check.class(bootstraps, c("integer", "numeric"))
        check.length(bootstraps, 1, " must be a single numeric value.")
    } else {
        ## Don't bootstraps
        bootstraps <- 1
    }

    ## CV, LASSO and PGLS
    check.class(CV, "logical")
    check.class(LASSO, "logical")
    check.class(PGLS, "logical")
    check.class(all.levels, "logical")

    ## Run the ldas
    lda_out <- lapply(factors, run.multi.lda, data_matrix, prior = prior, train = train, CV = CV, fun.type = fun.type, bootstraps = bootstraps, all.levels = all.levels, ...)
    # lda_out <- lapply(factors, run.multi.lda, data_matrix, prior = prior, train = train, CV = CV, fun.type = fun.type, bootstraps = bootstraps, all.levels = all.levels)
    
    ## Add lda names
    names(lda_out) <- names(factors)
    

    #TODO: Handle lda warnings! Capture them only once and print them only once.

    # lda_out structure =
    # lda_out$factor1$bootstraps1
    #                $bootstraps2
    #                $bootstraps3
    #        $factor2$bootstraps1
    #                $bootstraps2
    #                $bootstraps3


    #TODO: Add the option to update models
    #TODO: Test Cross Validation
    #TODO: Add PGLS
    #TODO: Add LASSO

    ## Add the support element
    lda_out$support <- list()

    ## Get the call
    lda_out$support$call <- match_call
    ## Get the bootstraps
    lda_out$support$bootstraps <- bootstraps
    ## Get the factors
    lda_out$support$factors <- factors
    ## Get the overall accuracy
    lda_out$support$accuracy <- apply.accuracy.score(lda_out)
    ## Get the overall fit
    lda_out$support$prop.trace <- apply.prop.trace(lda_out)

    ## Output the results
    class(lda_out) <- c("dispRity", "lda.test")
    return(lda_out)
}