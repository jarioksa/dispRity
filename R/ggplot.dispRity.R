#' @title Converts dispRity data for \code{\link[ggplot2]{ggplot}}
#'
#' @description Converts the dispRity into data that can be passed to \code{\link[ggplot2]{ggplot}}
#'
#' @param data the dispRity data to convert for \code{\link[ggplot2]{ggplot}}
#' @param raw logical, whether to pass the raw disparity data with no additional arguments (\code{TRUE}) or with the default arguments from \code{\link{plot.dispRity}} (\code{FALSE} - default). See details.
#' @param ... any arguments to be passed from \code{\link{plot.dispRity}} to \code{\link[ggplot2]{ggplot}}. See details.
#' 
#' @return
#' A \code{"gg"}, \code{"ggplot"} object.
#' 
#' @details
#' When using \code{raw = TRUE} and no optional arguments, the functions returns a table of at least two columns containing the groups and disparity values.
#' When using \code{raw = FALSE}, default optional arguments from \code{\link{plot.dispRity}} are also passed to \code{\link[ggplot2]{ggplot}}.
#' Optional arguments (\code{...}) will be passed to the \code{\link[ggplot2]{ggplot}} object.
#' 
#' @examples
#' ## Loading built in data
#' data(disparity)
#' 
#' ## Generating the gg object
#' gg <- ggplot.dispRity(disparity)
#' 
#' ## Plotting the ggplot object
#' plot(gg)
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export
#' @importFrom ggplot2 ggplot

ggplot.dispRity <- function(data, raw = FALSE, ...) {

    ## data
    check.class(data, "dispRity")

    ## raw
    check.class(raw, "logical")

    ## optional arguments
    dots <- list(...)
    optional_args <- ifelse(length(dots) == 0, FALSE, TRUE)

    ## Get the disparity data for plotting
    disparity_data <- do.call(rbind.data.frame, extract.dispRity(data, observed = FALSE))

    ## Convert into a proper data.frame
    split_data <- strsplit(rownames(disparity_data), "\\.")
    ## Get the column IDs
    split_ID <- unlist(lapply(split_data, function(X) return(X[2])))
    ## Get the group values
    split_group <- unlist(lapply(split_data, function(X) return(X[1])))
    if(!any(is.na(as.numeric(split_group)))) {
        split_group <- as.numeric(split_group)
    }

    ## Combining into a proper data.frame
    disparity_data <- cbind(split_group, disparity_data[,1])
    rownames(disparity_data) <- split_ID
    colnames(disparity_data)[1] <- ifelse(data$call$subsets[1] == "continuous", "time", "group")
    colnames(disparity_data)[2] <- "disparity"
    class(disparity_data) <- "data.frame"

    x <- ggplot2::ggplot(disparity_data)

    return(x)
}