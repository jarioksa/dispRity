#' @title Prints a \code{dispRity} object.
#'
#' @description Summarises the content of a \code{dispRity} object.
#'
#' @param x A \code{dispRity} object.
#' @param all \code{logical}; whether to display the entire object (\code{TRUE}) or just summarise its contents (\code{FALSE} - default).
#' @param ... further arguments to be passed to \code{print} or to \code{print.dispRity}.
#' 
#' @examples
#' ## Load the disparity data based on Beck & Lee 2014
#' data(disparity)
#' 
#' ## Displaying the summary of the object content
#' disparity
#' print(disparity) # the same
#' print.dispRity(disparity) # the same
#'
#' ## Displaying the full object
#' print.dispRity(disparity, all = TRUE)
#'
#' @seealso \code{\link{custom.subsets}}, \code{\link{chrono.subsets}}, \code{\link{boot.matrix}}, \code{\link{dispRity}}.
#'
#' @author Thomas Guillerme

## DEBUG
# warning("DEBUG dispRity.R")
# library(dispRity)
# source("sanitizing.R")
# source("dispRity.R")
# source("dispRity_fun.R")
# source("dispRity.metric.R")
# source("make.dispRity.R")
# source("fetch.dispRity.R")
# source("boot.matrix.R") ; source("boot.matrix_fun.R")
# source("chrono.subsets.R") ; source("chrono.subsets_fun.R")
# source("custom.subsets.R") ; source("custom.subsets_fun.R")
# data(BeckLee_mat50)
# data(BeckLee_tree)
# data_simple <- BeckLee_mat50
# data_boot <- boot.matrix(BeckLee_mat50, bootstraps = 11, rarefaction = c(5,6))
# data_subsets_simple <- chrono.subsets(BeckLee_mat50, tree = BeckLee_tree,  method = "discrete", time = c(120,80,40,20))
# data_subsets_boot <- boot.matrix(data_subsets_simple, bootstraps = 11, rarefaction = c(5,6))
# data <- dispRity(data_subsets_boot, metric = c(variances))

 
print.dispRity <- function(x, all = FALSE, ...) {

    match_call <- match.call()
    x_name <- match_call$x

    if(all) {
        ## Print everything
        tmp <- x
        class(tmp) <- "list"
        print(tmp)

    } else {

        ## ~~~~~~~
        ## Composite dispRity objects (subclasses)
        ## ~~~~~~~
        if(length(class(x)) > 1) {
            ## randtest
            if(class(x)[2] == "randtest") {

                ## Remove the call (messy)
                remove.call <- function(element) {
                    element$call <- "dispRity::null.test"
                    return(element)
                }
                x <- lapply(x, remove.call)

                if(length(x) == 1) {
                    print(x[[1]])
                } else {
                    tmp <- x
                    class(tmp) <- "list"
                    print(tmp) 
                }
                return()
            } 

            if(class(x)[2] == "model.test") {
                cat("Disparity evolution model fitting:\n")
                cat(paste0("Call: ", as.expression(x$call), " \n\n"))
                
                print(x$aic.models)

                cat(paste0("\nUse x$full.details for displaying the models details\n"))
                cat(paste0("or summary(x) for summarising them.\n"))

                return()
            }

            if(class(x)[2] == "model.sim") {

                cat("Disparity evolution model simulation:\n")
                cat(paste0("Call: ", as.expression(x$call), " \n\n"))
                cat(paste0("Model simulated (", x$nsim, " times):\n"))

                print(x$model)

                cat("\n")

                if(!is.null(x$p.value)) {
                    print(x$p.value)
                }

                return()
            }

            if(class(x)[2] == "dtt") {
                if(length(x) != 2){
                    ## Tested dtt
                    cat("Disparity-through-time test (modified from geiger:dtt)\n")
                    cat(paste0("Call: ", as.expression(x$call), " \n\n"))

                    cat(paste0("Observation: ", x$MDI , "\n\n"))

                    cat(paste0("Model: ", x$call$model , "\n"))
                    cat(paste0("Based on ", length(x$sim_MDI) , " replicates\n"))
                    cat(paste0("Simulated p-value: ", x$p_value , "\n"))
                    cat(paste0("Alternative hypothesis: ", x$call$alternative , "\n\n"))

                    print(c("Mean.dtt" = mean(x$dtt), "Mean.sim_MDI" = mean(x$sim_MDI), "var.sim_MDI" = var(x$sim_MDI)))

                    cat(paste0("\nUse plot.dispRity() to visualise."))
                    return()
                } else {
                    ## raw dtt
                    ## Fake an object with no attributes
                    x_tmp <- x
                    class(x_tmp) <- "list"
                    print(x_tmp)
                    cat(paste0("- attr(*, \"class\") = \"dispRity\" \"dtt\"\n"))
                    cat(paste0("Use plot.dispRity to visualise."))
                }
            }

            if(class(x)[2] == "lda.test") {

                ## Tested lda
                cat("Discriminant test:\n")
                cat(paste0("Call: ", as.expression(x$support$call), " \n\n"))

                ## Check the bootstraps
                is_bootstrapped <- ifelse(x$support$bootstraps > 1, TRUE, FALSE)

                ## Rounding digits
                round_digit <- 3

                ## Accuracy
                cat("Overall accuracy:\n")
                if(is_bootstrapped) {
                    median <- summarise.extract.list(x$support$accuracy$score, fun = stats::median, rounding = round_digit)
                    sd <- summarise.extract.list(x$support$accuracy$score, fun = stats::sd, rounding = round_digit)
                    print_df <- cbind(median, sd)
                    print_df <- as.data.frame(print_df)
                } else {
                    print_df <- cbind(lapply(x$support$accuracy$score, round, digits = round_digit))
                    print_df <- as.data.frame(print_df)
                    colnames(print_df) <- ""
                }

                if(!is.null(x$support$accuracy$test)) {
                    ## Getting the number of columns before the test
                    n_col_prior <- ncol(print_df)

                    token.converter <- function(p) {
                        p <- as.numeric(format(p, digits = 4))

                        if(p < 2.2e-16) {
                            return(paste0("< 2.2e-16", " ***"))
                        } else {
                            if(p < 0.001) {
                                return(paste0(p, " ***"))
                            } else {
                                if(p < 0.01) {
                                   return(paste0(p, " **")) 
                                } else {
                                    if(p < 0.05) {
                                        return(paste0(p, " *")) 
                                    } else {
                                        if(p < 0.1) {
                                            return(paste0(p, " *"))
                                        } else {
                                            return(paste0(p, ""))
                                        }
                                    }
                                }
                            }
                        }
                    }

                    ## Convert p_values
                    p_values <- lapply(x$support$accuracy$p_value, token.converter)
                    print_df$test <- p_values
                    ## Add the test name
                    colnames(print_df)[n_col_prior + 1] <- x$support$accuracy$test

                }

                print(print_df)
                cat("\n")

                ## Trace
                cat("Proportion of trace:\n")
                if(is_bootstrapped) {
                    median <- summarise.extract.list(x$support$prop.trace, fun = stats::median, rounding = round_digit)
                    sd <- summarise.extract.list(x$support$prop.trace, fun = stats::sd, rounding = round_digit)
                    prop_trace <- mapply(cbind, median, sd, SIMPLIFY = FALSE)
                    print_df <- lapply(prop_trace, function(x) {colnames(x) <- c("median", "sd") ; return(x)})
                } else {
                    clean.table <- function(table, round_digit) {
                        colnames(table) <- ""
                        return(round(table, digits = round_digit))
                    }
                    print_df <- lapply(x$support$prop.trace, clean.table, round_digit)
                }
                print(print_df)

                cat("Use summary.dispRity() or plot.dispRity() for displaying the full results.\n")
                return()
            }
        }

        
        ## ~~~~~~~
        ## Simple dispRity objects
        ## ~~~~~~~

        if(length(x$call) == 0) {
            if(!is.null(x$matrix) && class(x$matrix) == "matrix") {
                cat(" ---- dispRity object ---- \n")
                dims <- dim(x$matrix)
                cat(paste0("Contains only a matrix ", dims[1], "x", dims[2], "."))
            } else {
                cat("Empty dispRity object.\n")
            }
            return()
        }

        cat(" ---- dispRity object ---- \n")

        ## Print the matrix information
        if(any(names(x$call) == "subsets") && length(x$subsets) != 1) {
            ## Get the number of subsets (minus origin)
            subsets <- names(x$subsets)

            ## Check if there is more than one subset
            if(length(subsets) != 1) {

                ## Get the method
                method <- x$call$subsets
                if(length(method) != 1) {
                    method <- paste(method[1], " (", method[2],")", sep = "")
                }
                if(method == "customised") {
                    cat(paste(length(subsets), method, "subsets for", nrow(x$matrix), "elements"))    
                } else {
                    cat(paste(length(subsets), method, "time subsets for", nrow(x$matrix), "elements"))
                }
                if(length(x$call$dimensions) != 0) cat(paste(" with", x$call$dimensions, "dimensions"), sep = "")
                cat(":\n")
                if(length(subsets) > 5) {
                    cat("    ",paste(subsets[1:5], collapse=", "),"...\n")
                } else {
                    cat("    ",paste(subsets, collapse=", "), ".\n", sep="")
                }
            }
        } else {
            cat(paste(nrow(x$matrix), "elements"))
            if(length(x$call$dimensions) != 0) cat(paste(" with", x$call$dimensions, "dimensions"), sep = "")
            cat(".\n")
        }
        
        ## Print the bootstrap information
        if(any(names(x$call) == "bootstrap")) {
            if(x$call$bootstrap[[1]] != 0) {
                cat(paste("Data was bootstrapped ", x$call$bootstrap[[1]], " times (method:\"", x$call$bootstrap[[2]], "\")", sep = ""))
            }
            if(!is.null(x$call$bootstrap[[3]])) {
                if(x$call$bootstrap[[3]][[1]] == "full") {
                    cat(" and fully rarefied")
                } else {
                    cat(paste(" and rarefied to ", paste(x$call$bootstrap[[3]], collapse = ", "), " elements", sep = ""))
                }
            }
            cat(".\n")
        }

        ## Print the disparity information
        if(any(names(x$call) == "disparity")) {
            
            #metrics <- as.character(x$call$disparity$metrics)
            #strsplit(strsplit(metrics, split = "c(", fixed = TRUE)[[1]], split = ")", fixed = TRUE)[[2]][1]

            cat(paste("Disparity was calculated as:", paste(as.character(x$call$disparity$metrics$name), collapse = ", ")))
            cat(".\n")
        }
    }
}