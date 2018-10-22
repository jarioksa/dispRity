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

## Accuracy score
accuracy.score <- function(lda.test) {
        ## Get the attribution table
    attribution_table <- table(data$predict$class, data$data[-data$training, ncol(data$data)])

    # ## Get the prediction accuracy
    # get.accuracy <- function(data, attribution_table, scale = scale.accuracy) {
    #     if(!scale) {
    #         ## Return the average number of correct predictions
    #         return(sum(diag(attribution_table))/sum(attribution_table))
    #     } else {
    #         ## Return the number of scaled correct predictions
    #         obs_class_proportion <- table(data$data[, ncol(data$data)])
    #         obs_class_proportion <- obs_class_proportion/sum(obs_class_proportion)

    #         obs_class_proportion <- c(0.333, 0.333, 0.333)            

    #         sum(diag(attribution_table)/obs_class_proportion)

    #         /(sum(attribution_table))

    #         sum(attribution_table*obs_class_proportion)

    #         test <- table(data$predict$class, data$data[-data$training, ncol(data$data)])

    #         sum(diag(attribution_table))/sum(attribution_table)

    #         predict_rand1 <- which(data$predict$class == "random1")
    #         (data$predict$class[predict_rand1] == data$data[-data$training, ncol(data$data)][predict_rand1])
    #     }
    # }

    prediction_accuracy <- mean(data$predict$class == data$data[-data$training, ncol(data$data)])
}


