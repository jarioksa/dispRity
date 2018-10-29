## Factorising subsets from dispRity objects
factorise.subsets <- function(data) {

    ## Function for splitting the vector by categories
    split.category <- function(category, split_names) {
        output <- lapply(split_names, function(x, category) ifelse(x[[1]] == category, x[[2]], NA), category)
        return(as.vector(na.omit(unlist(output))))
    }

    ## Function for translating one category (unique)
    translate.category.all <- function(category, level, n_elements, subsets_list) {
        ## Empty vector out
        vector_out <- rep(NA, n_elements)
        
        ## Get the sublist of interest
        if(!is.null(category)) {
            sub_list <- unlist(subsets_list[grep(paste0(category, "."), names(subsets_list))], recursive = FALSE)
        } else {
            sub_list <- unlist(subsets_list, recursive = FALSE)
        }

        ## Loop through the replacement elements
        for(one_level in 1:length(sub_list)) {
            vector_out[sub_list[[one_level]]] <- level[one_level]
        }

        return(as.factor(vector_out))
    }

    ## Function for translating one category (not unique)
    translate.category.binary <- function(category, n_elements, subsets_list) {
        ## Empty vector out
        vector_out <- rep(NA, n_elements)

        ## Level name
        level <- names(subsets_list)[category]

        ## Translate
        vector_out[subsets_list[[category]]$elements] <- 1
        vector_out[-subsets_list[[category]]$elements] <- 0
        return(as.factor(vector_out))
    }

    ## Get the number of elements to subset
    n_elements <- dim(data$matrix)[1]

    ## Get the subsets_list
    subsets_list <- lapply(data$subsets, lapply, as.vector)

    ## Get the number of actual elements
    element_names <- names(subsets_list)

    ## Check whether it is a dot named list
    is_dot_named_list <- ifelse(length(grep("\\.", element_names)) == length(element_names), TRUE, FALSE)

    ## Check if the number of categories is unique
    length_unique <- length(unlist(subsets_list))
    are_unique <- ifelse(length_unique > n_elements, FALSE, TRUE)

    ## If not unique, check whether it is because of multiple categories
    if(!are_unique && is_dot_named_list) {
        are_unique <- ifelse(length_unique/n_elements == as.integer(length_unique/n_elements), TRUE, FALSE)
    }

    ## Each category as the right number of elements
    if(are_unique) {
        ## Named categories
        if(is_dot_named_list) {
            ## Split the names
            split_names <- strsplit(element_names, split = "\\.")

            ## Get the number of categories
            categories <- unique(unlist(lapply(split_names, function(x) return(x[[1]]))))

            ## Splitting by categories
            levels <- lapply(as.list(categories), split.category, split_names)

            ## Translate all categories
            subsets_out <- mapply(translate.category.all, as.list(categories), levels, MoreArgs = list(n_elements, subsets_list), SIMPLIFY = FALSE)
            names(subsets_out) <- categories

            return(subsets_out)
        } else {
            ## Translate a unique non-named list
            return(list(translate.category.all(category = NULL, level = element_names, n_elements, subsets_list)))
        }
    } else {
        ## Create binary categories
        output <- lapply(as.list(1:length(subsets_list)), translate.category.binary, n_elements, subsets_list)
        names(output) <- names(subsets_list)
        return(output)
    }
}

## Run a single LDA
run.one.lda <- function(factor, data_matrix, prior, train, CV, fun.type, ...) {

    ## First we select a subset of the dataset for training
    subset <- sample(1:nrow(data_matrix), train)
    training <- factor[subset]

    ## Set the prior (if missing)
    if(!prior[[1]][[1]]){
        ## If estimate prior
        prior <- as.numeric((table(training)/sum(table(training))))
    }

    ## Fitting the LDA
    lda_fit <- fun.type(x = data_matrix, grouping = factor, prior = prior, subset = subset, CV = CV, ...)

    ## Predicting the fit
    lda_predict <- predict(lda_fit, data_matrix[-subset, ])

    ## Return the predictions results and the fit
    return(list("fit" = lda_fit, "predict" = lda_predict, "training" = subset))
}

## Wrap lda runs
run.multi.lda <- function(factor, data_matrix, prior, train, CV, fun.type, bootstraps, ...) {
    ## Replicate the LDAs
    return(replicate(bootstraps, run.one.lda(factor, data_matrix, prior = prior, train = train, CV = CV, fun.type = fun.type, ...), simplify = FALSE))
}

## Getting a specific variable from a lda.test list
extract.lda.test <- function(lda_test, what, where, deviation, cent.tend) {

    ## Remove the calls and the factors
    data_tmp <- lda_test
    data_tmp$support <- NULL

    get.from.bootstrap <- function(one_bootstrap, where, what) {
        if(missing(what)) {
            return(one_bootstrap[[where]])
        } 
        if(missing(where)) {
            return(one_bootstrap[[what]])
        }

        return(one_bootstrap[[where]][[what]])
    }

    get.from.factor <- function(one_factor, where, what) {
        return(lapply(one_factor, get.from.bootstrap, where, what))
    }

    ## Get the elements
    elements <- lapply(data_tmp, get.from.factor, where, what)

    ## Make matrices
    if(all(unlist(lapply(elements, lapply, class)) == "factor")) {
        ## The matrix is a dataframe
        data.frame.no.names <- function(X) {
            data_frame <- data.frame(X)
            colnames(data_frame) <- seq(1:ncol(data_frame))
            return(data_frame)
        }
        matrices <- lapply(elements, data.frame.no.names)
    } else {
        ## The matrix is a matrix...
        matrices <- lapply(elements, function(X) do.call(cbind, X))
    }

    ## Get central tendency and deviation
    if(!missing(cent.tend)){
        central <- lapply(matrices, function(X) apply(X, 1, cent.tend))
    }
    
    if(!missing(deviation)){
        if(class(deviation) == "function") {
            ## Deviation is a deviation
            dev <- lapply(matrices, function(X) apply(X, 1, deviation))
        } else {
            ## Deviation are quantiles
            dev <- lapply(matrices, function(X) apply(X, 1, quantile, probs = deviation))
        }
    }

    if(missing(cent.tend) && missing(deviation)) {
        return(matrices)
    }

    if(missing(cent.tend) && !missing(deviation)) {
        return(list("deviation" = dev))
    }

    if(!missing(cent.tend) && missing(deviation)) {
        return(list("cent.tend" = central))
    }

    return(list("cent.tend" = central, "deviation" = dev))

}

## Applying a function to an extracted lda list
summarise.extract.list <- function(list_lda, fun, rounding = 10, ...) {

    ## Applying the function on one factor (i.e. element of the list)
    apply.to.matrix <- function(one_factor, fun, rounding, ...){
        output <- apply(one_factor, 1, function(one_factor, fun, rounding, ...) round(fun(one_factor, ...), digits = rounding), fun, rounding, ...)
        # output <- apply(one_factor, 1, function(one_factor, fun, rounding) round(fun(one_factor), digits = rounding), fun, rounding) ; warning("DEBUG lda.test_fun")

        ## Getting the output in the right format
        if(class(output) == "numeric") {
            ## If it's a numeric vector, make it into a matrix with levels as rownames
            output <- matrix(output, dimnames = list(names(output)))
        } else {
            ## If it's a matrix, check if the rownames match the input
            if(!all(rownames(output) %in% rownames(one_factor))) {
                output <- t(output)
            }
        }
        return(output)
    }

    apply.to.vector <- function(one_factor, fun, rounding, ...){
        output <- round(fun(one_factor, ...), digits = rounding)
        # output <- round(fun(one_factor), digits = rounding) ; warning("DEBUG lda.test_fun")

        ## Get the output in the right format (a 1xn matrix)
        return(matrix(output, nrow = 1, dimnames = list(NULL, names(output))))
    }

    ## Check the list elements class
    input_class <- unique(unlist(lapply(list_lda, class)))

    if(input_class == "matrix") {
        ## Input is a matrix, output is going to be a list of matrix
        return(lapply(list_lda, apply.to.matrix, fun = fun, rounding = rounding, ...))
    } else {
        ## Input is a vector, output is going to be a list of vector
        return(lapply(list_lda, apply.to.vector, fun = fun, rounding = rounding, ...))

    }
} 

## Accuracy score
accuracy.score <- function(posterior_table) {
    #TODO: Allow scaling
    return(sum(diag(posterior_table))/sum(posterior_table))
}

## Counting the posteriors
posterior.table <- function(posterior, prior, training) {
    return(table(posterior, prior[-training]))

    # priors
    # [1] 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0
    # posteriors
    # [1] 1 1 1 1 1 0 0 0 0 0 0
    # training
    # [1] 1 2 3 4 5 6 7 8 9

    #     1 2 3 4 5 6 7 8 9  11  13  15  17  19
    # [1] 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0
    #     t t t t t t t t t 1 1 1 1 1 0 0 0 0 0 0

    # Correct 1: 1
    # Correct 0: 6

    # 1 that are 0: 0
    # 0 that are 1: 4

    #           0 1 <- priors
    # post -> 0 6 0
    #         1 4 1
}

## Convert the priors into a bootstrap table
convert.prior.table <- function(priors, bootstraps) {
    return(lapply(priors, function(prior) as.data.frame(replicate(bootstraps, prior))))
}

## Getting the list of posterior tables
get.posterior.tables <- function(posteriors, priors, trainings) {
    ## Get the tables for each bootstraps
    mapply.posterior.table <- function(one_posteriors, one_priors, one_trainings) {
        mapply(posterior.table, as.list(one_posteriors), as.list(one_priors), as.list(data.frame(one_trainings)), SIMPLIFY = FALSE)
    }
    return(mapply(mapply.posterior.table, posteriors, priors, trainings, SIMPLIFY = FALSE))
}

## Applying the accuracy score to a whole lda-test object
apply.accuracy.score <- function(lda_test) {

    ## Get the posteriors
    posteriors <- extract.lda.test(lda_test, what = "class", where = "predict")

    ## Get the trainings
    trainings <- extract.lda.test(lda_test, what = "training")

    ## Get the prior in a table format
    priors <- convert.prior.table(lda_test$support$factors, lda_test$support$bootstraps) 

    ## Get the tables for each bootstraps
    posterior_tables <- get.posterior.tables(posteriors, priors, trainings)

    ## Get the accuracy scores
    accuracy_scores <- lapply(posterior_tables, lapply, accuracy.score)

    ## Get the scores in the right format
    output <- lapply(accuracy_scores, unlist)

    return(output)
}

## Applying the proportion of trace to a whole lda-test object
apply.prop.trace <- function(lda_test) {

    ## Get the svds
    svds <- extract.lda.test(lda_test, what = "svd", where = "fit")

    ## Get prop.trace
    get.svds <- function(one_svd) {
        out <- apply(one_svd, 2, function(svd) return(svd^2/sum(svd^2)))
        rownames(out) <- paste0("LD", 1:nrow(out))
        return(out)
    } 
    prop.trace <- lapply(svds, get.svds)
    
    return(prop.trace)
}