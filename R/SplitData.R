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
