## mapply function for applying the tests to distribution data
test.mapply <- function(pair_comparison, data, test, ...) {
    return(mapply(test, data[[pair_comparison[[1]]]], data[[pair_comparison[[2]]]], MoreArgs = ..., SIMPLIFY = FALSE))
}

## transforming test into a lapply from a given list of comp for multiple distributions
test.list.lapply.distributions <- function(list_of_comp, data, test, ...) {
    output <- lapply(list_of_comp, test.mapply, data, test, ...)
    return(output)
}

## Creating a list of sequences
set.sequence <- function(length) {
    ## Sequence of only two
    if(length == 2) {
        output <- matrix(data = c(1, 2), nrow = 2, byrow = TRUE)
    } else {
    ## sequence of more
        output <- matrix(data = c(1:(length - 1), 2:length), nrow = 2, byrow = TRUE)
    }
    return(output)
}

## convert a list from character to numeric
convert.to.numeric <- function(list, object) {
    return(lapply(list, match, names(object)))
}

## Getting the names (character) (convert.to.character internal)
names.fun <- function(list, object) {
    return(names(object[c(list)]))
}

## convert a list from numeric to character
convert.to.character <- function(list, object) {
    ## Applying to the list
    return(lapply(list, names.fun, object))
}


## function for repeating the extracted_data names (list to table internal)
rep.names <- function(name, subsets) {
    return(rep(name, subsets))
}


## Convert a list into a table (for aov)
list.to.table <- function(extracted_data, style = "group") {

    ## Get the list of names
    names_list <- as.list(names(extracted_data))
    ## If no names, just get list numbers
    if(length(names_list) == 0) {
        names_list <- as.list(seq(from = 1, to = length(extracted_data)))
    }
    subsets_length <- unlist(lapply(extracted_data, length), recursive = FALSE)

    ## Create the data.frame
    output <- data.frame("data" = unlist(extracted_data), row.names = NULL, "subsets" = unlist(mapply(rep.names, names_list, subsets_length, SIMPLIFY = FALSE)))

    ## Transform groups to numeric
    if(style == "binomial") {
        output$group <- as.numeric(output$group)-1
    }

    return(output)
}

## lapply function for htest.to.vector
get.element <- function(print, htest) {
    return(htest[grep(print, names(htest))][[1]])
}

## Converts an htest into a vector
htest.to.vector <- function(htest, print) {
    ## print is a vector of htest elements to print
    return(unlist(lapply(print, get.element, htest)))
}

## Set the list of comparisons
set.comparisons.list <- function(comp, extracted_data, comparisons) {
    options(warn = -1)

    switch(comp,
        custom = {
            ## get the list of subsets to compare
            comp_subsets <- comparisons            
        },
        pairwise = {
            ## Get the pairs of subsets
            comp_subsets <- combn(1:length(extracted_data), 2)
            ## convert pair subsets table into a list of pairs
            comp_subsets <- unlist(apply(comp_subsets, 2, list), recursive = FALSE)
        },
        sequential = {
            ## Set the list of sequences
            comp_subsets <- set.sequence(length(extracted_data))
            ## convert seq subsets into a list of sequences
            comp_subsets <- unlist(apply(comp_subsets, 2, list), recursive = FALSE)            
        },
        referential = {
            ## Set the list of comparisons as a matrix
            matrix_data <- c(rep(1, length(extracted_data) - 1), seq(from = 2, to = length(extracted_data)))
            comp_subsets <- matrix(matrix_data, ncol = (length(extracted_data) - 1), byrow = TRUE)
            ## convert pair subsets table into a list of pairs
            comp_subsets <- unlist(apply(comp_subsets, 2, list), recursive = FALSE)
        }
    )

    options(warn = 0)

    return(comp_subsets)
}

## Save the comparisons list
save.comparison.list <- function(comp_subsets, extracted_data) {
    ## Saving the list of comparisons
    comparisons_list <- convert.to.character(comp_subsets, extracted_data)
    comparisons_list <- unlist(lapply(comparisons_list, paste, collapse = " : "))
    return(comparisons_list)
}

## Function for lapplying aov type functions
lapply.lm.type <- function(data, test, ...) {
    return(test(data ~ subsets, data = data, ...))
}


## Calculate the central tendency and the quantiles from a table of results
get.quantiles.from.table <- function(table, con.cen.tend, conc.quantiles, ...) {
    return(t(rbind(apply(table, 1, con.cen.tend, ...), apply(table, 1, quantile, probs = conc.quantiles, na.rm = TRUE, ...))))
}

## Returning a table of numeric values
output.numeric.results <- function(details_out, name, comparisons_list, conc.quantiles, con.cen.tend) {
    ## Transforming list to table
    table_temp <- do.call(rbind.data.frame, details_out)

    ## Getting the eventual parameter name
    param_name <- unique(as.character(lapply(details_out, names)))

    ## Calculate the quantiles and the central tendency
    if(!missing(conc.quantiles) && !missing(con.cen.tend)) {
        table_out <- get.quantiles.from.table(table_temp, con.cen.tend, conc.quantiles)
    } else {
        table_out <- table_temp
    }

    ## Getting column names
    if(param_name != "NULL") {
        colnames(table_out)[1] <- paste(name, param_name, sep = ": ")
    } else {
        colnames(table_out)[1] <- name
    }
    ## Getting row names (the comparisons)
    row.names(table_out) <- comparisons_list

    return(table_out)
}



## lapply function for getting the test elements (output.htest.results internal)
lapply.output.test.elements <- function(test_element, details_out, comparisons_list, conc.quantiles, con.cen.tend) {
    if(!missing(conc.quantiles) && !missing(con.cen.tend)) {
        return(output.numeric.results(lapply(lapply(details_out, lapply, htest.to.vector, print = test_element), unlist), test_element, comparisons_list, conc.quantiles, con.cen.tend))
    } else {
        return(output.numeric.results(lapply(lapply(details_out, lapply, htest.to.vector, print = test_element), unlist), test_element, comparisons_list))
    }
}

## Returning a table for htests
output.htest.results <- function(details_out, comparisons_list, conc.quantiles, con.cen.tend, correction) {
    ## Getting the test elements
    test_elements <- unique(unlist(lapply(details_out, lapply, names)))
    ## Selecting the numeric (or) integer elements only
    test_elements <- test_elements[grep("numeric|integer", unlist(lapply(as.list(details_out[[1]][[1]]), class)))]
    ## Remove null.value and the estimates
    remove <- match(c("null.value", "conf.int", "estimate"), test_elements)
    if(any(is.na(remove))) {
        remove <- remove[-which(is.na(remove))]
    }
    if(length(remove) > 0) {
        test_elements <- test_elements[-remove]
    }
    
    ## Get the results
    if(!missing(conc.quantiles) && !missing(con.cen.tend)) {
        table_out <- lapply(as.list(test_elements), lapply.output.test.elements, details_out, comparisons_list, conc.quantiles, con.cen.tend)
    } else {
        table_out <- lapply(as.list(test_elements), lapply.output.test.elements, details_out, comparisons_list)
    }

    ## Applying the correction
    if(correction != "none") {
        ## Find which element contains the p-values
        pvalues <- unlist(lapply(lapply(lapply(table_out, colnames), function(X) X == "p.value"), any))
        if(any(pvalues)) {
            table_out[pvalues][[1]] <- apply(table_out[pvalues][[1]], 2, p.adjust, method = correction, n = nrow(table_out[pvalues][[1]]))
        }
    }

    return(table_out)
}

# ## Handling output for lm multiple tests
# output.lm.results <- function(details_out, conc.quantiles, con.cen.tend) {

#     ## Getting the summaries
#     summaries <- lapply(details_out, summary)
    
#     ## Transforming the list 
#     list_of_results <- list()
#     for(element in 1:length(summaries[[1]][[1]])) {
#         list_of_results[[element]] <- matrix(unlist(lapply(lapply(summaries, `[[`, 1), `[[`, element)), nrow = length(summaries[[1]][[1]][[element]]),
#             dimnames = list(c("subsets", "Residuals")))
#     }

#     ## Get the quantiles
#     list_of_results <- lapply(list_of_results, get.quantiles.from.table, con.cen.tend, conc.quantiles, na.rm = TRUE)

#     ## Name the elements
#     for(element in 1:length(summaries[[1]][[1]])) {
#         colnames(list_of_results[[element]])[[1]] <- names(summaries[[1]][[1]])[[element]]
#     }    

#     return(list_of_results)
# }