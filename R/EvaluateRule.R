EvaluateRule <- function(data,
                                  BuildRule.object=NULL,
                                  B=NULL,
                                  study.design, #=c("RCT", "observational"),
                                  name.outcome,
                                  type.outcome,
                                  desirable.outcome,
                                  separate.propensity.estimation=TRUE,
                                  clinical.threshold=0,
                                  name.treatment,
                                  names.influencing.treatment,
                                  names.influencing.rule,
                                  propensity.method="logistic.regression",
                                  truncate.propensity.score=TRUE,
                                  truncate.propensity.score.threshold=0.05, 
                                  observation.weights=NULL,
                                  additional.weights=rep(1, nrow(data)),
                                  lambda.choice=c("min", "1se"),
                                  propensity.k.cv.folds=10,
                                  propensity.b.cv.repeats=1,
                                  bootstrap.CI=FALSE,
                                  bootstrap.CI.replications=1000,
                                  bootstrap.type="basic") {
    lambda.choice <- match.arg(lambda.choice)
    if (is.null(propensity.method) & is.null(observation.weights)) {
        stop("need to specify either logistic regression for estimating the propensity score, or a vector of observation weights to use instead")
    }
    if (!is.null(propensity.method) & !is.null(observation.weights)) {
        stop("cannot specify both a method for estimating the propenstiy score (implying it needs to be computed) and a vector of observation weights to use instead")
    }
    if (!is.null(propensity.method)) {
        stopifnot(propensity.method %in% c("logistic.regression"))
    }
    if (!is.null(observation.weights) & !is.numeric(observation.weights))
        stop("if observation weights are provided, they need to be numeric")
    if (!is.data.frame(data)) {
        stop("dataset must be a data frame")
    }
    if (all(data[, name.treatment] %in% c(0, 1)) == FALSE) {
        stop("treatment needs to be coded as a binary indicator")
    }
    if (is.null(BuildRule.object) & is.null(B)) {
        stop("either BuildRule.object or B has to be specified")
    }
    if (is.null(BuildRule.object) == FALSE) {
        prediction.approach <- BuildRule.object$prediction.approach
        if (desirable.outcome != TRUE & (BuildRule.object$prediction.approach %in% c("OWL", "OWL.framework"))) {
            stop("OWL-based procedures assumes that higher values of the outcome variable are better")
        }
        if (prediction.approach %in% c("OWL", "OWL.framework")) {
            warning(paste(prediction.approach, "assumes that larger values of the outcome variable are better"))
        }
    }
    stopifnot(is.logical(truncate.propensity.score))
    stopifnot(is.numeric(truncate.propensity.score.threshold))
    if (truncate.propensity.score == FALSE) {
        truncate.propensity.score.threshold <- 0
    }
    if (truncate.propensity.score.threshold < 0 | truncate.propensity.score.threshold > 0.25) {
        stop("truncate.propensity.scorethreshold should be between 0 and 0.25")
    }
    stopifnot(length(additional.weights) == nrow(data))
    stopifnot(is.numeric(additional.weights))
    stopifnot(is.logical(separate.propensity.estimation))
    n <- nrow(data)
    my.formatted.data <- FormatData(data=data,
                                                   name.outcome=name.outcome,
                                                   name.treatment=name.treatment,
                                                   type.outcome=type.outcome,
                                                   names.influencing.rule=names.influencing.rule,
                                                   names.influencing.treatment=names.influencing.treatment)
    do.one.EvaluateRule <- EvaluateRuleOnce(data=data,
                                               my.formatted.data=my.formatted.data,
                                               BuildRule.object=BuildRule.object,
                                               B=B,
                                               study.design=study.design, 
                                               type.outcome=type.outcome,
                                               desirable.outcome=desirable.outcome,
                                               separate.propensity.estimation=separate.propensity.estimation,
                                               clinical.threshold=clinical.threshold,
                                               names.influencing.treatment=names.influencing.treatment,
                                               names.influencing.rule=names.influencing.rule,
                                               propensity.method=propensity.method,
                                               truncate.propensity.score=truncate.propensity.score,
                                               truncate.propensity.score.threshold=truncate.propensity.score.threshold,
                                               observation.weights=observation.weights,
                                               additional.weights=additional.weights,
                                               lambda.choice=lambda.choice,
                                               propensity.k.cv.folds=propensity.k.cv.folds,
                                               propensity.b.cv.repeats=propensity.b.cv.repeats)
    
    observed.n.test.positives <- do.one.EvaluateRule$n.test.positives
    observed.ATE.test.positives <- do.one.EvaluateRule$ATE.test.positives
    observed.n.test.negatives <- do.one.EvaluateRule$n.test.negatives
    observed.ATE.test.negatives <- do.one.EvaluateRule$ATE.test.negatives
    observed.ABR <- do.one.EvaluateRule$ABR
    if (bootstrap.CI == TRUE) {
        vec.n.test.positives <- rep(NA, bootstrap.CI.replications)
        vec.mean.ATE.test.positives <- rep(NA, bootstrap.CI.replications)
        vec.n.test.negatives <- rep(NA, bootstrap.CI.replications)
        vec.mean.ATE.test.negatives <- rep(NA, bootstrap.CI.replications)
        vec.mean.ABR <- rep(NA, bootstrap.CI.replications)
        for (b in 1:bootstrap.CI.replications) {
            idx.boot <- sample(1:n, size=n, replace=TRUE)
            my.formatted.data.boot <- FormatData(data=data[idx.boot, ],
                                                  name.outcome=name.outcome,
                                                  name.treatment=name.treatment,
                                                  type.outcome=type.outcome,
                                                  names.influencing.rule=names.influencing.rule,
                                                  names.influencing.treatment=names.influencing.treatment)
            if (!is.null(observation.weights)) {
                boot.observation.weights <- observation.weights[idx.boot]
            } else {
                boot.observation.weights <- NULL
            }
            do.one.EvaluateRule.boot <- EvaluateRuleOnce(data=data[idx.boot, ],
                                                            my.formatted.data=my.formatted.data.boot,
                                                            BuildRule.object=BuildRule.object,
                                                            B=B[idx.boot], 
                                                            study.design=study.design, 
                                                            type.outcome=type.outcome,
                                                            desirable.outcome=desirable.outcome,
                                                            additional.weights=additional.weights,
                                                            separate.propensity.estimation=separate.propensity.estimation,
                                                            lambda.choice=lambda.choice,
                                                            clinical.threshold=clinical.threshold,
                                                            names.influencing.treatment=names.influencing.treatment,
                                                            names.influencing.rule=names.influencing.rule,
                                                            propensity.method=propensity.method,
                                                            observation.weights=boot.observation.weights,
                                                            propensity.k.cv.folds=propensity.k.cv.folds,
                                                            propensity.b.cv.repeats=propensity.b.cv.repeats,
                                                            truncate.propensity.score=truncate.propensity.score,
                                                            truncate.propensity.score.threshold=truncate.propensity.score.threshold)
            vec.n.test.positives[b] <- do.one.EvaluateRule.boot$n.test.positives
            vec.mean.ATE.test.positives[b] <- do.one.EvaluateRule.boot$ATE.test.positives
            vec.n.test.negatives[b] <- do.one.EvaluateRule.boot$n.test.negatives
            vec.mean.ATE.test.negatives[b] <- do.one.EvaluateRule.boot$ATE.test.negatives
            vec.mean.ABR[b] <- do.one.EvaluateRule.boot$ABR
        }
        if (bootstrap.type == "basic") {
            if (observed.n.test.positives > 0) {
                basic.bootstrap.CI.ATE.test.positives <- c(2 * observed.ATE.test.positives - quantile(vec.mean.ATE.test.positives, probs=0.975),
                                                                          2 * observed.ATE.test.positives - quantile(vec.mean.ATE.test.positives, probs=0.025))
                names(basic.bootstrap.CI.ATE.test.positives) <- c("2.5%", "97.5%")
            } else {
                basic.bootstrap.CI.ATE.test.positives <- c(NA, NA)
            }
            if (observed.n.test.negatives > 0) {
                basic.bootstrap.CI.ATE.test.negatives <- c(2 * observed.ATE.test.negatives - quantile(vec.mean.ATE.test.negatives, probs=0.975),
                                                                          2 * observed.ATE.test.negatives - quantile(vec.mean.ATE.test.negatives, probs=0.025))
                names(basic.bootstrap.CI.ATE.test.negatives) <- c("2.5%", "97.5%")
            } else {
                basic.bootstrap.CI.ATE.test.negatives <- c(NA, NA)
            }
            basic.bootstrap.CI.ABR <- c(2 * observed.ABR - quantile(vec.mean.ABR, probs=0.975),
                                                                          2 * observed.ABR - quantile(vec.mean.ABR, probs=0.025))
            names(basic.bootstrap.CI.ABR) <- c("2.5%", "97.5%")
            vec.mean.ABR
        } else {
            stop("only basic bootstrap CI is supported for now")
        }
        return(list("recommended.treatment"=do.one.EvaluateRule$recommended.treatment,
                      "vec.mean.ATE.test.positives"=vec.mean.ATE.test.positives,
                      "bootstrap.CI.ATE.test.positives"=basic.bootstrap.CI.ATE.test.positives,
                      "vec.mean.ATE.test.negatives"=vec.mean.ATE.test.negatives,
                      "bootstrap.CI.ATE.test.negatives"=basic.bootstrap.CI.ATE.test.negatives,
                      "bootstrap.CI.ABR"=basic.bootstrap.CI.ABR,
                      "n.test.positives"=observed.n.test.positives,
                      "ATE.test.positives"=observed.ATE.test.positives,
                      "n.test.negatives"=observed.n.test.negatives,
                      "ATE.test.negatives"=observed.ATE.test.negatives,
                      "ABR"=observed.ABR))
    } else {
        return(list("recommended.treatment"=do.one.EvaluateRule$recommended.treatment,
                      "fit.object"=do.one.EvaluateRule$fit.object,
                       "n.test.positives"=observed.n.test.positives,
                      "ATE.test.positives"=observed.ATE.test.positives,
                      "n.test.negatives"=observed.n.test.negatives,
                      "ATE.test.negatives"=observed.ATE.test.negatives,
                      "ABR"=observed.ABR))
    }
}
