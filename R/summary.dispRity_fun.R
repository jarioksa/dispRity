## Converts one or more CI into quantile probabilities
CI.converter <- function(CI) {
    sort(c(50-CI/2, 50+CI/2)/100)
}

# Wrapper function for summarising a single rarefaction
get.summary <- function(disparity_subsets_rare, cent.tend, quantiles) {

    output <- list()

    ## Summarising NA
    if(is.na(disparity_subsets_rare[[1]])) {
        if(!missing(cent.tend)) {
            output$cent.tend <- NA
        }
        if(!missing(quantiles)) {
            output$quantiles <- rep(NA, length(quantiles) * 2)
        }
        return(output)
    }

    ## Summarising normal data
    if(!missing(cent.tend)) {
        output$cent_tend <- cent.tend(as.vector(disparity_subsets_rare))
    }
    if(!missing(quantiles)) {
        output$quantiles <- quantile(as.vector(disparity_subsets_rare), probs = CI.converter(quantiles))
    }
    return(output)
}

## lapply wrapper function for summarising a single subset
lapply.summary <- function(disparity_subsets, cent.tend, quantiles) {
    return(lapply(disparity_subsets[-1], get.summary, cent.tend, quantiles))
}

## lapply wrapper for getting elements
lapply.get.elements <- function(subsets, bootstrapped = TRUE) {
    if(bootstrapped){
        return(unlist(lapply(subsets[-1], nrow)))
    } else {
        return(unlist(lapply(subsets, nrow)))
    }
}

## lapply wrapper for getting the disparity observed values
lapply.observed <- function(disparity) {
    output <- if(any(is.na(disparity$elements))) {
        return(NA)
    } else {
        return(c(disparity$elements))
    }
}

## mapply wrapper for getting the disparity observed values
mapply.observed <- function(disparity, elements) {
    return(c(disparity, rep(NA, (length(elements)-1))))
}

## Get digits for table (shifts the decimal point to contain a maximum of four characters)
get.digit <- function(column) {
    if(max(nchar(round(column)), na.rm = TRUE) <= 4) {
        return(4 - max(nchar(round(column)), na.rm = TRUE))
    } else {
        return(0)
    }
}


round.column <- function(column, digits) {
    ## Get the digits value
    digits <- ifelse(digits != "default", digits, get.digit(as.numeric(column)))
    ## Round the table
    if(class(column) != "group") {
        return(round(as.numeric(column), digits = digits))
    } else {
        return(round(as.numeric(as.character(column)), digits = digits))
    }
}


## Function for digits the results
digits.fun <- function(results_table, digits, model.test = FALSE) {

    ## model.test
    start_column <- ifelse(model.test, 1, 3)

    ## Apply the digits
    rounded_table <- as.matrix(results_table[,c(start_column:ncol(results_table))])
    rounded_table <- apply(rounded_table, 2, round.column, digits)
    results_table[,c(start_column:ncol(results_table))] <- rounded_table

    return(results_table)
}


## match parameters lapply (for model.test summaries)
match.parameters <- function(one_param, full_param) {
    ## Which parameters are present:
    available_param <- full_param %in% names(one_param)
    if(any(!available_param)) {
        ## Making the missing parameters NAs
        output_param <- rep(NA, length(full_param))
        names(output_param) <- full_param
        output_param[available_param] <- one_param[order(match(names(one_param), full_param))]
        return(output_param)
    } else {
        return(one_param)
    }
}

## Summary for model.test
summary.model.test <- function(data, quantiles, cent.tend, recall, digits, match_call){
    ## Extracting the AICs and the log likelihoods
    base_results <- cbind(data$aic.models, "log.lik" = sapply(data$full.details, function(x) x$value))

    ## Extracting the additional parameters
    parameters <- lapply(data$full.details, function(x) x$par)

    # MP: allow summaries to work on a single model
    if(length(parameters) == 1)  {
        param.tmp <- c(parameters[[1]])
        #names(param.tmp) <- rownames(parameters[[1]])
        parameters <- list(param.tmp)
        }
    base_results <- cbind(base_results, "param" = unlist(lapply(parameters, length)))
    
    ## Get the full list of parameters
    
       names_list <- lapply(parameters, names)
       full_param <- unique(unlist(names_list))

    output_table <- cbind(base_results, do.call(rbind, lapply(parameters, match.parameters, full_param)))
   
    ## Rounding
    summary_results <- digits.fun(output_table, digits, model.test = TRUE)

    return(summary_results)
}

## Summary for model.sim
summary.model.sim <- function(data, quantiles, cent.tend, recall, digits, match_call) {
    ## Extract the central tendencies
    simulation_data_matrix <- sapply(data$simulation.data$sim, function(x) x$central_tendency)

    ## Get the quantiles
    simulation_results <- apply(simulation_data_matrix, 1, get.summary, cent.tend = cent.tend, quantiles = quantiles)
    simulation_results <- cbind(do.call(rbind, lapply(simulation_results, function(X) rbind(X$cent_tend[[1]]))),
                                do.call(rbind, lapply(simulation_results, function(X) rbind(X$quantiles))))
    colnames(simulation_results)[1] <- as.character(match_call$cent.tend)

    ## Output table
    output_table <- cbind("subsets" = rev(data$simulation.data$fix$subsets),
                          "n" = data$simulation.data$fix$sample_size,
                          "var" = unname(data$simulation.data$fix$variance),
                          simulation_results)
    rownames(output_table) <- seq(1:nrow(output_table))
    return(output_table)
}


summary.lda.test <- function(data, quantiles, cent.tend, recall, match_call, digits, ...) {

    ## Prediction:

    ## Default digits handling
    if(digits == "default") {
        digits = 3
    }

    ## Checking whether the data is bootstrapped
    is_bootstrapped <- ifelse(data$support$bootstraps > 1, TRUE, FALSE)

    ## Get the observed priors
    prior_obs <- lapply(data$support$factors, table)
    prior_rel <- unlist(lapply(prior_obs, function(counts) return(counts/sum(counts))))
    prior_obs <- unlist(prior_obs)

    ## Get the groups names
    group_names <- names(prior_obs)

    ## Get the training prior
    priors <- extract.lda.test(data, what = "prior", where = "fit")

    ## Get the counts
    counts <- extract.lda.test(data, what = "counts", where = "fit")

    ## Get the posteriors
    posteriors <- extract.lda.test(data, what = "class", where = "predict")

    # Get the posterior tables
    # warning("The posteriors are based on the classes only, not on the posterior -> to change")
    posterior_tables <- get.posterior.tables(posteriors, convert.prior.table(data$support$factors, data$support$bootstraps), extract.lda.test(data, what = "training"))

    ## Get posterior sums
    posteriors_counts <- lapply(posterior_tables, lapply, function(table) apply(table, 1, sum))
    posteriors_ratio <- lapply(posteriors_counts, lapply, function(value) return(value/sum(value)))
    posteriors_ratio <- lapply(lapply(posteriors_ratio, lapply, round, digits = digits), function(factor) do.call(cbind, factor))

    ## Check quantile class
    quantile_fun <- ifelse(class(quantiles) == "function", TRUE, FALSE) 

    if(is_bootstrapped) {
        ## Priors
        priors_cent_tend <- summarise.extract.list(priors, fun = cent.tend, rounding = digits)
        if(quantile_fun) {
            priors_spread <- summarise.extract.list(priors, fun = quantile, rounding = digits)
        } else {
            priors_spread <- summarise.extract.list(priors, fun = quantile, rounding = digits, probs = CI.converter(quantiles))
        }
        
        
        ## Prior chunk
        prior_chunk <- rbind(unlist(priors_cent_tend), matrix(unlist(priors_spread), ncol = length(group_names)))
        ## Prior chunk names
        prior_chunk_names <- paste0("prior.", as.expression(match_call$cent.tend))
        col_names <- colnames(priors_spread[[1]])
        if(is.null(col_names)){
            prior_chunk_names <- c(prior_chunk_names, paste0("prior.", as.expression(match_call$quantiles)))
        } else {
            prior_chunk_names <- c(prior_chunk_names, paste0("prior.", col_names))
        }

        ## Posteriors
        post_cent_tend <- summarise.extract.list(posteriors_ratio, fun = cent.tend, rounding = digits)
        if(quantile_fun) {
            post_spread <- summarise.extract.list(posteriors_ratio, fun = quantile, rounding = digits)
        } else {
            post_spread <- summarise.extract.list(priors, fun = quantile, rounding = digits, probs = CI.converter(quantiles))
        }

        ## Posterior chunk
        post_chunk <- rbind(unlist(post_cent_tend), matrix(unlist(post_spread), ncol = length(group_names)))
        ## Posterior chunk names
        post_chunk_names <- paste0("post.", as.expression(match_call$cent.tend))
        col_names <- colnames(post_spread[[1]])
        if(is.null(col_names)){
            post_chunk_names <- c(post_chunk_names, paste0("post.", as.expression(match_call$quantiles)))
        } else {
            post_chunk_names <- c(post_chunk_names, paste0("post.", col_names))
        }

    } else {
        ## Get the priors
        prior_chunk <- unlist(priors)
        prior_chunk_names <- "prior"
        ## Get the posteriors
        post_chunk <- unlist(posteriors_ratio)
        post_chunk_names <- "posterior"
    }

    ## Make the output table
    posterior_results <- matrix(c(prior_obs, prior_rel), nrow = 2, dimnames = list(c("obs.prior", "proportion"), group_names), byrow = TRUE)
    posterior_results <- rbind(posterior_results, prior_chunk, post_chunk)
    rownames(posterior_results)[-c(1,2)] <- c(prior_chunk_names, post_chunk_names)


    ## Group means:

    ## Get the group means
    group_means <- extract.lda.test(data, what = "means", where = "fit")

    if(is_bootstrapped) {
        ## Splitting the matrices
        split.matrix<-function(one_group_mean, bootstraps) {
            start <- seq(from = 1, to = ncol(one_group_mean), by = ncol(one_group_mean)/bootstraps)
            end <- (start + ncol(one_group_mean)/bootstraps) - 1
            return(lapply(as.list(data.frame(rbind(start, end))), function(margin, matrix) return(matrix[, c(margin[1]:margin[2])]), matrix = one_group_mean))
        } 

        ## list per bootstraps
        group_means <- lapply(group_means, split.matrix, bootstraps = data$support$bootstraps)

        ## Merging the matrices (and extracting the means)
        merge.matrices <- function(one_group_mean, fun) {
            return(apply(simplify2array(one_group_mean), 1:2, fun))
        }

        ## Apply the central tendency on all the groups
        group_means <- lapply(group_means, merge.matrices, fun = cent.tend)
        group_means_out <- do.call(rbind, group_means)
        rownames(group_means_out) <- paste(rep(names(group_means), unlist(lapply(group_means, nrow))), rownames(group_means_out), sep = ".")
    } else {
        ## Get the output in the right format (non-bootstrapped)
        group_means_out <- do.call(rbind, group_means)
        rownames(group_means_out) <- paste(rep(names(group_means), unlist(lapply(group_means, nrow))), rownames(group_means_out), sep = ".")
    }


    ## coefficients:
    coefficients <- NULL


    return(list("prediction" = posterior_results, "group_means" = group_means_out))


    #$Prediction
    #          total mean.prior sd.prior mean.factor.1.post factor.2.post
    # factor.1 0.33  0.3333333 0.3333333 
    # factor.2 0.3333333 0.3333333 0.3333333 
    #TG: check if transposing the matrix is not nicer
    #       factor.1 factor.2 factor.3
    # total
    # mean.prior
    # sd.prior

    #$Group_means
    #          var1     var2    var3
    # factor.1 0.33  0.3333333 0.3333333 
    # factor.2 0.3333333 0.3333333 0.3333333 

    #$Coefficient discriminant
    #          LD1     LD2    L3
    # var1 0.33  0.3333333 0.3333333 
    # var2 0.3333333 0.3333333 0.3333333 

}









# ## Saving results function
# save.results.seq.test <- function(model, results) {
#     save_out <- match(results, names(summary(model)))
#     return(summary(model)[save_out])
# }

# ## Transform a matrix (usually the coefficient results) into a list
# matrix.to.list <- function(matrix, no.intercept = TRUE) {
#     ## Remove the intercept column (if needed)
#     if(no.intercept != FALSE) {
#         if(any(rownames(matrix) == "(Intercept)")) {
#             matrix <- matrix[-which(rownames(matrix) == "(Intercept)"),]
#         }
#     }

#     ## Transform the matrix into a list
#     output <- as.list(matrix)
#     ## Adding names (if necessary)
#     if(is.null(names(output))) {
#         names(output) <- colnames(matrix)
#     }
#     return(output)
# }

# ## Relists a list (recursive) with element names
# relist.names <- function(element, elements_names) {
#     output <- list(element)
#     names(output) = elements_names
#     return(output)
# }

# ## Gets results table from results elements
# get.results.table <- function(element, results_elements, cent.tend, quantiles, comparisons, match_call, is.distribution) {
#     ## Get the data
#     if(is.null(element)) {
#         element = 1
#     }
#     data <- lapply(lapply(results_elements, lapply, `[[`, element), unlist)
    
#     ## Central tendency
#     results_table <- matrix(data = lapply(data, cent.tend), ncol = 1, dimnames = list(comparisons))
#     if(!is.null(match_call$cent.tend)) {
#         colnames(results_table) <- paste(match_call$cent.tend)
#     } else {
#         colnames(results_table) <- "mean"
#     }

#     if(is.distribution != FALSE) {
#         ## Quantiles
#         results_quantiles <- lapply(data, quantile, probs = CI.converter(quantiles))
#         ## Create table
#         results_table <- cbind(results_table, matrix(data = unlist(results_quantiles), nrow = length(comparisons), byrow = TRUE, dimnames = list(c(), c(names(results_quantiles[[1]])))))
#     }
#     return(results_table)
# }

# summary.seq.test <- function(data, quantiles, cent.tend, recall, digits, results, match_call) {
    
#     ## SANITIZING
#     ## quantiles
#     check.class(quantiles, "numeric", " must be any value between 1 and 100.")
#     ## remove warnings
#     if(any(quantiles < 1)) {
#         stop("quantiles(s) must be any value between 1 and 100.")
#     }
#     if(any(quantiles > 100)) {
#         stop("quantiles(s) must be any value between 1 and 100.")
#     }

#     ## Check if is distribution
#     is.distribution <- ifelse(length(data$models[[1]]) == 1, FALSE, TRUE)

#     ## SAVING THE RESULTS
#     models_results <- lapply(data$models, lapply, save.results.seq.test, results)
#     intercepts_results <- lapply(data$intercepts, unlist)
#     comparisons <- unique(unlist(names(data$models)))

#     ## Getting the slopes
#     if(is.distribution != FALSE) {
#         ## Getting the coefficients
#         ## Transforming the coefficients into a list
#         results_coefficients <- lapply(lapply(models_results, lapply, `[[`, 1), lapply, matrix.to.list, no.intercept = TRUE)

#         ## Creating the tables for each element in the matrices
#         elements_list_matrix <- as.list(names(results_coefficients[[1]][[1]]))
#         table_matrix <- lapply(elements_list_matrix, get.results.table, results_coefficients, cent.tend = cent.tend, quantiles = quantiles, comparisons = comparisons, match_call = match_call, is.distribution = is.distribution)
#         names(table_matrix) <- elements_list_matrix

#         ## Check if there are any other results to output
#         if(length(results) > 1) {
#             coefficients_matrix <- which(unlist(lapply(models_results[[1]][[1]], class)) == "matrix")
#             other_results <- results[-which(results == "coefficients")]
#             ## Extracting the other elements
#             results_list <- lapply(models_results, lapply, `[[`, -coefficients_matrix)
#             ## Rename and relist the elements
#             elements_names <- names(models_results[[1]][[1]][-coefficients_matrix])
#             results_list <- lapply(results_list, lapply, relist.names, elements_names)

#             ## Creating the tables for each element in the list
#             elements_list_list <- as.list(names(results_list[[1]][[1]]))
#             table_list <- lapply(elements_list_list, get.results.table, results_list, cent.tend = cent.tend, quantiles = quantiles, comparisons = comparisons, match_call = match_call, is.distribution = is.distribution)
#             names(table_list) <- elements_list_list

#             table_matrix <- append(table_matrix, table_list)
#         }

#         ## Apply p-adjust (if necessary)
#         if(!is.null(data$correction)) {
#             table_matrix$`Pr(>|t|)` <- apply(table_matrix$`Pr(>|t|)`, 2, p.adjust, method = data$correction)
#         }

#     } else {
#         ## Creating the table for the first model
#         table_matrix <- models_results[[1]][[1]][[1]]
#         ## Creating the table for the other models
#         table_tmp <- matrix(unlist(lapply(models_results[-1], lapply, `[[`, "coefficients"), use.names = FALSE), ncol = ncol(table_matrix), byrow = TRUE)
#         ## Combining the tables
#         table_matrix <- rbind(table_matrix, table_tmp)
#         ## Removing the intercept
#         table_matrix <- table_matrix[-1,]
#         ## Adding the rownames
#         rownames(table_matrix) <- comparisons
#         ## Check if there are any other results to output
#         if(length(results) > 1) {
#             coefficients_matrix <- which(unlist(lapply(models_results[[1]][[1]], class)) == "matrix")
#             ## Extracting the other elements
#             results_list <- as.matrix(unlist(lapply(models_results, lapply, `[[`, -coefficients_matrix)))
#             ## Adding the elements names
#             colnames(results_list) <- results[-which(results == "coefficients")]
#             ## Combing it to the table
#             table_matrix <- cbind(table_matrix, results_list)
#         }

#         ## Apply p-adjust (if necessary)
#         if(!is.null(data$correction)) {
#             table_matrix[,which(colnames(table_matrix) == "Pr(>|t|)")] <- p.adjust(table_matrix[,which(colnames(table_matrix) == "Pr(>|t|)")], method = data$correction)
#         }
#     }

#     ## Getting the initial intercepts
#     if(is.distribution == TRUE) {
#         ## Get the first intercept
#         initial_intercept <- lapply(lapply(models_results, lapply, `[[`, 1), lapply, matrix.to.list, no.intercept = FALSE)[1]
#         ## Creating the tables for each element in the matrices
#         elements_list_matrix <- as.list(names(initial_intercept[[1]][[1]]))
#         ## Remove the NAs (i.e. slopes)
#         elements_list_matrix <- elements_list_matrix[-which(is.na(elements_list_matrix))]
#         intercept_list <- lapply(elements_list_matrix, get.results.table, initial_intercept, cent.tend = cent.tend, quantiles = quantiles, comparisons = comparisons[1], match_call = match_call, is.distribution = is.distribution)
#         ## Convert list into a matrix
#         initial_intercept <- matrix(unlist(intercept_list), ncol = ncol(intercept_list[[1]]), nrow = length(elements_list_matrix), byrow = TRUE,
#             dimnames = list(c(elements_list_matrix), c(colnames(intercept_list[[1]]))))

#         ## Add the predicted
#         predicted_intercept <- get.results.table(NULL, intercepts_results[-1], cent.tend = cent.tend, quantiles = quantiles, comparisons = comparisons[-1], match_call = match_call, is.distribution = is.distribution)
#         intercept_matrix <- list("Initial" = initial_intercept, "Predicted" = predicted_intercept)
#     } else {
#         ## Get the first intercept
#         initial_intercept <- models_results[[1]][[1]][[1]]
#         ## Add the other intercepts (estimated)
#         intercept_matrix <- matrix(NA, nrow = (length(intercepts_results)-1), ncol = ncol(initial_intercept))
#         intercept_matrix[,1] <- unlist(intercepts_results[-1])
#         ## Bind the two tables
#         intercept_matrix <- rbind(initial_intercept, intercept_matrix)
#         ## Remove the slope
#         intercept_matrix <- intercept_matrix[-2,]
#         ## Add rownmaes
#         rownames(intercept_matrix) <- comparisons
#     }

#     ## Combining the tables
#     results_out <- list("Slopes" = table_matrix, "Intercepts" = intercept_matrix)

#     return(results_out)
# }