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
        n.training <- ceiling(split.proportions[1] * n)
        idx.training <- sort(sample(x=1:n, size=n.training, replace=FALSE))
        idx.test <- (1:n)[-idx.training]
        partition <- rep(NA, n)
        partition[idx.training] <- "training"
        partition[idx.test] <- "test"
        partition <- factor(partition, levels=c("training", "test"))
    }
    if (n.sets == 3) {
        if (is.null(split.proportions)) {
            split.proportions <- c(0.5, 0.25, 0.25)
        }
        n.training <- ceiling(split.proportions[1] * n)
        n.validation <- floor(split.proportions[2] * n)
        n.test <- n - (n.training + n.validation)
        idx.training <- sort(sample(x=1:n, size=n.training, replace=FALSE))
        idx.validation <- sort(sample((1:n)[-idx.training], size=n.validation, replace=FALSE))
        idx.test <- (1:n)[-sort(c(idx.training, idx.validation))]
        partition <- rep(NA, n)
        partition[idx.training] <- "training"
        partition[idx.validation] <- "validation"
        partition[idx.test] <- "test"
        partition <- factor(partition, levels=c("training", "validation", "test"))
    }
    modified.data <- data.frame(data, partition)
    return(modified.data)
}
