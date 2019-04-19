#' Partition a dataset into independent subsets 
#'
#' To get a trustworthy estimate of how a developed treatment rule will perform in independent samples drawn from the same population, it is critical that rule development be performed independently of rule evaluation. Further, it is common to perform model selection to settle on the form of the developed treatment rule and, in this case, it is essential that the ultimately chosen treatment rule is also evaluated on data that did not inform any stage of the model-building. The \code{SplitData()} function partitions a dataset so rule development/validation/evaluation (or development/evaluation if there is no model selection) can quickly be performed on independent datasets. This function is only appropriate for the simple setting where the rows in a given dataset are independent of one another (e.g. the same individuals are not represented with multiple rows).

#' @param data A data frame representing the *development*  dataset used for building a treatment rule
#' @param n.sets A numeric/integer equal to either 3 (if a development/validation/evaluation partition is desired) or 2 (if there is no model-selection and only a development/evaluation partition is desired).
#' @param split.proportions A numeric vector with length equal to \code{n.sets}, providing the proportion of observations in \code{data} that should be assigned to the development/evaluation partitions (if \code{n.sets=2}) or to the development/validation/evaluation partitions (if \code{n.sets=3}). The entries must sum to 1.
#' @return A data.frame equal to \code{data} with an additional column named `partition', which is a factor variable with levels equal to `development' and `evaluation' (if \code{n.sets=2})  or to  `development', `validation', and `evaluation' (if \code{n.sets=3}).
#' @examples
#' set.seed(123)
#' example.split <- SplitData(data=obsStudyGeneExpressions,
#'                                      n.sets=3, split.proportions=c(0.5, 0.25, 0.25))
#' table(example.split$partition)
#' @export

SplitData <- function(data,
                              n.sets=c(3, 2),
                              split.proportions=NULL) {
    if (!is.data.frame(data)) {
        stop("the data must be a data frame")
    }
    if (is.null(split.proportions) == FALSE) {
        if (sum(split.proportions) != 1) {
            stop("split.proportions must be a numeric vector whose entries add up to 1")
        }
    }
    if (length(split.proportions) != n.sets) {
        stop("the length of split.proportions must be equal to n.sets")
    }
    n <- nrow(data)
    if (n.sets == 2) {
        if (is.null(split.proportions)) {
            split.proportions <- c(0.7, 0.3)
        }
        n.development <- ceiling(split.proportions[1] * n)
        idx.development <- sort(sample(x=1:n, size=n.development, replace=FALSE))
        idx.evaluation <- (1:n)[-idx.development]
        partition <- rep(NA, n)
        partition[idx.development] <- "development"
        partition[idx.evaluation] <- "evaluation"
        partition <- factor(partition, levels=c("development", "evaluation"))
    }
    if (n.sets == 3) {
        if (is.null(split.proportions)) {
            split.proportions <- c(0.5, 0.25, 0.25)
        }
        n.development <- ceiling(split.proportions[1] * n)
        n.validation <- floor(split.proportions[2] * n)
        n.evaluation <- n - (n.development + n.validation)
        idx.development <- sort(sample(x=1:n, size=n.development, replace=FALSE))
        idx.validation <- sort(sample((1:n)[-idx.development], size=n.validation, replace=FALSE))
        idx.evaluation <- (1:n)[-sort(c(idx.development, idx.validation))]
        partition <- rep(NA, n)
        partition[idx.development] <- "development"
        partition[idx.validation] <- "validation"
        partition[idx.evaluation] <- "evaluation"
        partition <- factor(partition, levels=c("development", "validation", "evaluation"))
    }
    modified.data <- data.frame(data, partition)
    return(modified.data)
}
