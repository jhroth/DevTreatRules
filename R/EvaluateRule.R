#' Evaluate a Treatment Rule 
#'
#' Perform principled evaluation of a treatment rule (using the IPW approach to account for potential confounding) on a dataset that is independent of the development dataset on which the rule was developed, either to perform model selection (with a validation dataset) or to obtain trustworthy estimates of performance for a pre-specified treatment rule (with an evaluation dataset).
#'
#' @param evaluation.data A data frame representing the *validation* or *evaluation* dataset used to estimate the performance of a rule that was developed on an independent development dataset.
#' @param BuildRule.object The object returned by the \code{BuildRule()} function. Defaults to NULL but is required if a treatment rule is not provided in the \code{B} argument. Only one of \code{BuildRule.object} and \code{B} should be specified. 
#' @param B A numeric vector representing a pre-specified treatment rule, which must have length equal to the number of rows in \code{evaluation.data} and elements equal to \code{0/FALSE} indicating no treatment and \code{1/TRUE} indicating treatment. Defaults to \code{NULL} but is required if \code{BuildRule.object} is not specified.  Only one of \code{BuildRule.object} and \code{B} should be specified. 
#' @param study.design Either `observational', `RCT', or `naive'. For the \code{observational} design, the function will use inverse-probability-of-treatment observation weights (IPW) based on estimated propensity scores with predictors \code{names.influencing.treatment}; for the \code{RCT} design, the function will use IPW based on propensity scores equal to the observed sample proportions; for the \code{naive} design, all observation weights will be uniformly equal to 1.
#' @param name.outcome A character indicating the name of the outcome variable in \code{evaluation.data}.
#' @param type.outcome Either `binary' or `continuous', the form of \code{name.outcome}.
#' @param desirable.outcome A logical equal to \code{TRUE} if higher values of the outcome are considered desirable (e.g. for a binary outcome, 1 is more desirable than 0). The \code{OWL.framework} and \code{OWL} approaches to treatment rule estimation require a desirable outcome.
#' @param separate.propensity.estimation A logical equal to \code{TRUE} if propensity scores should be estimated separately in the test-positives and test-negatives subpopulations and equal to \code{FALSE} if propensity scores should be estimated in the combined sample. Default is \code{TRUE}.
#' @param clinical.threshold A numeric equal to a positive number above which the predicted outcome under treatment must be superior to the predicted outcome under control for treatment to be recommended. Only used when \code{BuildRuleObject} was specified and derived from the split-regression or direct-interactions approach. Default is 0.
#' @param name.treatment A character indicating the name of the treatment variable in \code{evaluation.data}.
#' @param names.influencing.treatment A character vector (or element) indicating the names of the variables in \code{evaluation.data} that are expected to influence treatment assignment in the current dataset. Required for \code{study.design=}`observational'. 
#' @param names.influencing.rule A character vector (or element) indicating the names of the variables in \code{evaluation.data} that may influence response to treatment and are expected to be observed in future clinical settings.
#' @param propensity.method One of`logistic.regression', `lasso', or `ridge'. This is the underlying regression model used to estimate propensity scores (for \code{study.design=}`observational'. If \code{bootstrap.CI=TRUE}, then \code{propensity.method} must be `logistic.regression'. Defaults to NULL.
#' @param show.treat.all A logical variable dictating whether summaries for the naive rule that assigns treatment to all observations are reported, which help put the performance of the estimated treatment rule in context. Default is TRUE
#' @param show.treat.none A logical variable dictating whether summaries for the naive rule that assigns treatment to no observations are reported, which help put the performance of the estimated treatment rule in context. Default is TRUE
#' @param truncate.propensity.score A logical variable dictating whether estimated propensity scores less than \code{truncate.propensity.score.threshold} away from 0 or 1 should be truncated to be \code{truncate.propensity.score.threshold} away from 0 or 1.
#' @param truncate.propensity.score.threshold A numeric value between 0 and 0.25.
#' @param observation.weights A numeric vector equal to the number of rows in \code{evaluation.data} that provides observation weights to be used in place of the IPW weights estimated with \code{propensity.method}. Defaults to NULL. Only one of the \code{propensity.method} and \code{observation.weights} should be specified.
#' @param additional.weights A numeric vector of observation weights that will be multiplied by IPW weights in the rule evaluation stage, with length equal to the number of rows in \code{evaluation.data}.. This can be used, for example, to account for a non-representative sampling design or to apply an IPW adjustment for missingness. The default is a vector of 1s.
#' @param lambda.choice Either `min' or `1se', corresponding to the \code{s} argument in \code{predict.cv.glmnet()} from the \code{glmnet} package; only used when \code{propensity.method} or \code{rule.method} is `lasso' or `ridge'. Default is `min'.
#' @param propensity.k.cv.folds An integer dictating how many folds to use for K-fold cross-validation that chooses the tuning parameter when \code{propensity.method} is `lasso' or `ridge'. Default is 10.
#' @param bootstrap.CI Logical indicating whether the ATE/ABR estimates returned by \code{EvaluateRule()} should be accompanied by 95\% confidence intervals based on the bootstrap. Default is \code{FALSE}
#' @param bootstrap.CI.replications An integer specifying how many bootstrap replications should underlie the computed CIs. Default is 1000.
#' @param bootstrap.type One character element specifying the type of bootstrap CI that should be computed. Currently the only supported option is \code{bootstrap.type=}`basic', but this may be expanded in the future.
#' @return A list with the following components
#' \itemize{
#'   \item \code{recommended.treatment}: A numeric vector of 0s and 1s, with length equal to the number of rows in \code{evaluation.data}, where a 0 indicates treatment is not recommended and a 1 indicates treatment is recommended for the corresponding observation in \code{evaluation.data}.
#'   \item \code{fit.object}: A list consisting of one of the following: the propensity scores estimated in the test-positives and in the test-negatives (if \code{separate.propensity.estimation=TRUE}, \code{study.design=}`observational', and \code{observation.weights=NULL}); the propensity scores estimated in the combined sample (if \code{separate.propensity.estimation=FALSE}, \code{study.design=}`observational', and \code{observation.weights=NULL}); and simply is simply null if \code{study.design=}`RCT' (in which case propensity score would just be the inverse of the sample proportion receiving treatment)
#' \item \code{summaries}: a matrix with columns reporting the following summaries of treatment rule performance: the number of observations in \code{evaluation.data} recommended to receive treatment. (\code{n.positives}); the estimated average treatment effect among those recommended to receive treatment (\code{ATE.positives}); the number of observations in \code{evaluation.data} recommended to not receive treatment (\code{n.negatives}); the estimated average treatment effect among those recommended to not receive treatment (\code{ATE.negatives}); the estimated average benefit of using the rule, with the weighted average of ATE.positives and -1 * ATE.negatives where weights are the proportions of test-positives and test-negatives (\code{ABR}). If \code{bootstrap.CI=TRUE}, then 4 additional columns are included, showing the lower bound (LB) and upper bound (UB) of the 95\% CIs for \code{ATE.positives} and \code{ATE.negatives}.
#' }
#' @examples
#' set.seed(123)
#' example.split <- SplitData(data=obsStudyGeneExpressions,
#'                                      n.sets=3, split.proportions=c(0.5, 0.25, 0.25))
#' development.data <- example.split[example.split$partition == "development",]
#' validation.data <- example.split[example.split$partition == "validation",]
#' one.rule <- BuildRule(development.data=development.data,
#'                      study.design="observational",
#'                      prediction.approach="split.regression",
#'                      name.outcome="no_relapse",
#'                      type.outcome="binary",
#'                      desirable.outcome=TRUE,
#'                      name.treatment="intervention",
#'                      names.influencing.treatment=c("prognosis", "clinic", "age"),
#'                      names.influencing.rule=c("age", paste0("gene_", 1:10)),
#'                      propensity.method="logistic.regression",
#'                      rule.method="glm.regression")
#' split.validation <- EvaluateRule(evaluation.data=validation.data,
#'                           BuildRule.object=one.rule,
#'                           study.design="observational",
#'                           name.outcome="no_relapse",
#'                           type.outcome="binary",
#'                           desirable.outcome=TRUE,
#'                           name.treatment="intervention",
#'                           names.influencing.treatment=c("prognosis", "clinic", "age"),
#'                           names.influencing.rule=c("age", paste0("gene_", 1:10)),
#'                           propensity.method="logistic.regression",
#'                           bootstrap.CI=FALSE)
#' split.validation[c("n.positives", "n.negatives",
#'                        "ATE.positives", "ATE.negatives", "ABR")]
#' @export

EvaluateRule <- function(evaluation.data,
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
                                  propensity.method=NULL, #"logistic.regression",
                                  show.treat.all=TRUE,
                                  show.treat.none=TRUE,
                                  truncate.propensity.score=TRUE,
                                  truncate.propensity.score.threshold=0.05, 
                                  observation.weights=NULL,
                                  additional.weights=rep(1, nrow(evaluation.data)),
                                  lambda.choice=c("min", "1se"),
                                  propensity.k.cv.folds=10,
                                  bootstrap.CI=FALSE,
                                  bootstrap.CI.replications=1000,
                                  bootstrap.type="basic") {
    lambda.choice <- match.arg(lambda.choice)
    stopifnot(is.logical(show.treat.all) & is.logical(show.treat.all))
    if (is.null(propensity.method) & is.null(observation.weights)) {
        stop("need to specify either logistic regression for estimating the propensity score, or a vector of observation weights to use instead")
    }
    if (!is.null(propensity.method) & !is.null(observation.weights)) {
        stop("cannot specify both a method for estimating the propensity score (implying it needs to be computed) and a vector of observation weights to use instead")
    }
    if (!is.null(propensity.method) & bootstrap.CI == TRUE) {
        if (!(propensity.method %in% c("logistic.regression"))) {
            stop("if bootstrap.CI=TRUE, then propensity.method must be logistic.regression")
        }
    }
    if (!is.null(observation.weights) & !is.numeric(observation.weights))
        stop("if observation weights are provided, they need to be numeric")
    if (!is.data.frame(evaluation.data)) {
        stop("dataset must be a data frame")
    }
    if (all(evaluation.data[, name.treatment] %in% c(0, 1)) == FALSE) {
        stop("treatment needs to be coded as a binary indicator")
    }
    if (is.null(BuildRule.object) & is.null(B)) {
        stop("either BuildRule.object or B has to be specified")
    }
    if (!is.null(B)) {
        stopifnot(is.numeric(B))
        stopifnot(length(B) == nrow(evaluation.data))
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
    stopifnot(length(additional.weights) == nrow(evaluation.data))
    stopifnot(is.numeric(additional.weights))
    stopifnot(is.logical(separate.propensity.estimation))
    if (show.treat.all == TRUE & show.treat.none == TRUE) {
        nrow.summaries <- 3
        rownames.summaries <- c("estimated.rule", "treat.all", "treat.none")
    } else if (show.treat.all == TRUE & show.treat.none == FALSE) {
        nrow.summaries <- 2
        rownames.summaries <- c("estimated.rule", "treat.all")
    } else if (show.treat.all == FALSE & show.treat.none == TRUE) {
        nrow.summaries <- 2
        rownames.summaries <- c("estimated.rule", "treat.none")
    } else if (show.treat.all == FALSE & show.treat.none == FALSE) {
        nrow.summaries <- 1
        rownames.summaries <- c("estimated.rule")
    }
    n <- nrow(evaluation.data)
    my.formatted.data <- FormatData(data=evaluation.data,
                                                   name.outcome=name.outcome,
                                                   name.treatment=name.treatment,
                                                   type.outcome=type.outcome,
                                                   names.influencing.rule=names.influencing.rule,
                                                   names.influencing.treatment=names.influencing.treatment)
    do.one.EvaluateRule <- EvaluateRuleOnce(data=evaluation.data,
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
                                               propensity.k.cv.folds=propensity.k.cv.folds)
    if (show.treat.all == TRUE) {
        do.one.EvaluateRule.treat.all <- EvaluateRuleOnce(data=evaluation.data,
                                                          my.formatted.data=my.formatted.data,
                                                          BuildRule.object=NULL,
                                                          B=rep(1, nrow(evaluation.data)),
                                                          study.design=study.design, 
                                                          type.outcome=type.outcome,
                                                          desirable.outcome=desirable.outcome,
                                                          separate.propensity.estimation=separate.propensity.estimation,
                                                          clinical.threshold=clinical.threshold,
                                                          names.influencing.treatment=names.influencing.treatment,
                                                          names.influencing.rule=names.influencing.rule,
                                                          propensity.method="logistic.regression",
                                                          truncate.propensity.score=truncate.propensity.score,
                                                          truncate.propensity.score.threshold=truncate.propensity.score.threshold,
                                                          observation.weights=observation.weights,
                                                          additional.weights=additional.weights,
                                                          lambda.choice=lambda.choice,
                                                          propensity.k.cv.folds=propensity.k.cv.folds)
    }
    if (show.treat.none == TRUE) {
        do.one.EvaluateRule.treat.none <- EvaluateRuleOnce(data=evaluation.data,
                                                           my.formatted.data=my.formatted.data,
                                                           BuildRule.object=NULL,
                                                           B=rep(0, nrow(evaluation.data)),
                                                           study.design=study.design, 
                                                           type.outcome=type.outcome,
                                                           desirable.outcome=desirable.outcome,
                                                           separate.propensity.estimation=separate.propensity.estimation,
                                                           clinical.threshold=clinical.threshold,
                                                           names.influencing.treatment=names.influencing.treatment,
                                                           names.influencing.rule=names.influencing.rule,
                                                           propensity.method="logistic.regression",
                                                           truncate.propensity.score=truncate.propensity.score,
                                                           truncate.propensity.score.threshold=truncate.propensity.score.threshold,
                                                           observation.weights=observation.weights,
                                                           additional.weights=additional.weights,
                                                           lambda.choice=lambda.choice,
                                                           propensity.k.cv.folds=propensity.k.cv.folds)
    }
    observed.n.positives <- do.one.EvaluateRule$n.positives
    observed.ATE.positives <- do.one.EvaluateRule$ATE.positives
    observed.n.negatives <- do.one.EvaluateRule$n.negatives
    observed.ATE.negatives <- do.one.EvaluateRule$ATE.negatives
    observed.ABR <- do.one.EvaluateRule$ABR
    if (bootstrap.CI == TRUE) {
        vec.n.positives <- rep(NA, bootstrap.CI.replications)
        vec.mean.ATE.positives <- rep(NA, bootstrap.CI.replications)
        vec.n.negatives <- rep(NA, bootstrap.CI.replications)
        vec.mean.ATE.negatives <- rep(NA, bootstrap.CI.replications)
        vec.mean.ABR <- rep(NA, bootstrap.CI.replications)
        for (b in 1:bootstrap.CI.replications) {
            idx.boot <- sample(1:n, size=n, replace=TRUE)
            my.formatted.data.boot <- FormatData(data=evaluation.data[idx.boot, ],
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
            do.one.EvaluateRule.boot <- EvaluateRuleOnce(data=evaluation.data[idx.boot, ],
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
                                                            truncate.propensity.score=truncate.propensity.score,
                                                            truncate.propensity.score.threshold=truncate.propensity.score.threshold)
            vec.n.positives[b] <- do.one.EvaluateRule.boot$n.positives
            vec.mean.ATE.positives[b] <- do.one.EvaluateRule.boot$ATE.positives
            vec.n.negatives[b] <- do.one.EvaluateRule.boot$n.negatives
            vec.mean.ATE.negatives[b] <- do.one.EvaluateRule.boot$ATE.negatives
            vec.mean.ABR[b] <- do.one.EvaluateRule.boot$ABR
        }
        if (bootstrap.type == "basic") {
            if (observed.n.positives > 0) {
                basic.bootstrap.CI.ATE.positives <- c(2 * observed.ATE.positives - stats::quantile(vec.mean.ATE.positives, probs=0.975),
                                                                          2 * observed.ATE.positives - stats::quantile(vec.mean.ATE.positives, probs=0.025))
                names(basic.bootstrap.CI.ATE.positives) <- c("2.5%", "97.5%")
            } else {
                basic.bootstrap.CI.ATE.positives <- c(NA, NA)
            }
            if (observed.n.negatives > 0) {
                basic.bootstrap.CI.ATE.negatives <- c(2 * observed.ATE.negatives - stats::quantile(vec.mean.ATE.negatives, probs=0.975),
                                                                          2 * observed.ATE.negatives - stats::quantile(vec.mean.ATE.negatives, probs=0.025))
                names(basic.bootstrap.CI.ATE.negatives) <- c("2.5%", "97.5%")
            } else {
                basic.bootstrap.CI.ATE.negatives <- c(NA, NA)
            }
            basic.bootstrap.CI.ABR <- c(2 * observed.ABR - stats::quantile(vec.mean.ABR, probs=0.975),
                                                                          2 * observed.ABR - stats::quantile(vec.mean.ABR, probs=0.025))
            names(basic.bootstrap.CI.ABR) <- c("2.5%", "97.5%")
        } else {
            stop("only basic bootstrap CI is supported for now")
        }
        mat.summaries <- matrix(NA, nrow=nrow.summaries, ncol=9)
        rownames(mat.summaries) <- rownames.summaries
        colnames(mat.summaries) <- c("n.positives", "ATE.positives", "ATE.positives.CI.LB", "ATE.positives.CI.UB", 
                                                    "n.negatives", "ATE.negatives", "ATE.negatives.CI.LB", "ATE.negatives.CI.UB", 
                                                    "ABR")
        mat.summaries["estimated.rule", "n.positives"] <- do.one.EvaluateRule$n.positives
        mat.summaries["estimated.rule", "ATE.positives"] <- do.one.EvaluateRule$ATE.positives
        mat.summaries["estimated.rule", "ATE.positives.CI.LB"] <- basic.bootstrap.CI.ATE.positives[1]
        mat.summaries["estimated.rule", "ATE.positives.CI.UB"] <- basic.bootstrap.CI.ATE.positives[2]
        mat.summaries["estimated.rule", "n.negatives"] <- do.one.EvaluateRule$n.negatives
        mat.summaries["estimated.rule", "ATE.negatives"] <- do.one.EvaluateRule$ATE.negatives
        mat.summaries["estimated.rule", "ATE.negatives.CI.LB"] <- basic.bootstrap.CI.ATE.negatives[1]
        mat.summaries["estimated.rule", "ATE.negatives.CI.UB"] <- basic.bootstrap.CI.ATE.negatives[2]
        mat.summaries["estimated.rule", "ABR"] <- do.one.EvaluateRule$ABR
        
        if (show.treat.all == TRUE) {
            mat.summaries["treat.all", "n.positives"] <- do.one.EvaluateRule.treat.all$n.positives
            mat.summaries["treat.all", "ATE.positives"] <- do.one.EvaluateRule.treat.all$ATE.positives
            mat.summaries["treat.all", "ATE.positives.CI.LB"] <- NA
            mat.summaries["treat.all", "ATE.positives.CI.UB"] <- NA
            mat.summaries["treat.all", "n.negatives"] <- do.one.EvaluateRule.treat.all$n.negatives
            mat.summaries["treat.all", "ATE.negatives"] <- do.one.EvaluateRule.treat.all$ATE.negatives
            mat.summaries["treat.all", "ATE.negatives.CI.LB"] <- NA
            mat.summaries["treat.all", "ATE.negatives.CI.UB"] <- NA
            mat.summaries["treat.all", "ABR"] <- do.one.EvaluateRule.treat.all$ABR
        }
        if (show.treat.none == TRUE) {
            mat.summaries["treat.none", "n.positives"] <- do.one.EvaluateRule.treat.none$n.positives
            mat.summaries["treat.none", "ATE.positives"] <- do.one.EvaluateRule.treat.none$ATE.positives
            mat.summaries["treat.none", "ATE.positives.CI.LB"] <- NA
            mat.summaries["treat.none", "ATE.positives.CI.UB"] <- NA
            mat.summaries["treat.none", "n.negatives"] <- do.one.EvaluateRule.treat.none$n.negatives
            mat.summaries["treat.none", "ATE.negatives"] <- do.one.EvaluateRule.treat.none$ATE.negatives
            mat.summaries["treat.none", "ATE.negatives.CI.LB"] <- NA
            mat.summaries["treat.none", "ATE.negatives.CI.UB"] <- NA
            mat.summaries["treat.none", "ABR"] <- do.one.EvaluateRule.treat.none$ABR
        }
        return(list("recommended.treatment"=do.one.EvaluateRule$recommended.treatment,
                      "fit.object"=do.one.EvaluateRule$fit.object,
                       "summaries"=mat.summaries))
                       ## "bootstrap.CI.ATE.positives"=basic.bootstrap.CI.ATE.positives,
                       ## "bootstrap.CI.ATE.negatives"=basic.bootstrap.CI.ATE.negatives,
                       ## "bootstrap.CI.ABR"=basic.bootstrap.CI.ABR))
                      ## "vec.mean.ATE.positives"=vec.mean.ATE.positives,
                      ## "vec.mean.ATE.negatives"=vec.mean.ATE.negatives,
                      ##                       "n.positives"=observed.n.positives,
                      ## "ATE.positives"=observed.ATE.positives,
                      ## "n.negatives"=observed.n.negatives,
                      ## "ATE.negatives"=observed.ATE.negatives,
                      ## "ABR"=observed.ABR))
    } else {
        mat.summaries <- matrix(NA, nrow=nrow.summaries, ncol=5)
        colnames(mat.summaries) <- c("n.positives", "ATE.positives", "n.negatives", "ATE.negatives", "ABR")
        rownames(mat.summaries) <- rownames.summaries
        mat.summaries["estimated.rule", "n.positives"] <- do.one.EvaluateRule$n.positives
        mat.summaries["estimated.rule", "ATE.positives"] <- do.one.EvaluateRule$ATE.positives
        mat.summaries["estimated.rule", "n.negatives"] <- do.one.EvaluateRule$n.negatives
        mat.summaries["estimated.rule", "ATE.negatives"] <- do.one.EvaluateRule$ATE.negatives
        mat.summaries["estimated.rule", "ABR"] <- do.one.EvaluateRule$ABR
        if (show.treat.all == TRUE) {
            mat.summaries["treat.all", "n.positives"] <- do.one.EvaluateRule.treat.all$n.positives
            mat.summaries["treat.all", "ATE.positives"] <- do.one.EvaluateRule.treat.all$ATE.positives
            mat.summaries["treat.all", "n.negatives"] <- do.one.EvaluateRule.treat.all$n.negatives
            mat.summaries["treat.all", "ATE.negatives"] <- do.one.EvaluateRule.treat.all$ATE.negatives
            mat.summaries["treat.all", "ABR"] <- do.one.EvaluateRule.treat.all$ABR
        }
        if (show.treat.none == TRUE) {
            mat.summaries["treat.none", "n.positives"] <- do.one.EvaluateRule.treat.none$n.positives
            mat.summaries["treat.none", "ATE.positives"] <- do.one.EvaluateRule.treat.none$ATE.positives
            mat.summaries["treat.none", "n.negatives"] <- do.one.EvaluateRule.treat.none$n.negatives
            mat.summaries["treat.none", "ATE.negatives"] <- do.one.EvaluateRule.treat.none$ATE.negatives
            mat.summaries["treat.none", "ABR"] <- do.one.EvaluateRule.treat.none$ABR
        }
        return(list("recommended.treatment"=do.one.EvaluateRule$recommended.treatment,
                       "fit.object"=do.one.EvaluateRule$fit.object,
                       "summaries"=mat.summaries))
                      ##  "n.positives"=observed.n.positives,
                      ## "ATE.positives"=observed.ATE.positives,
                      ## "n.negatives"=observed.n.negatives,
                      ## "ATE.negatives"=observed.ATE.negatives,
                      ## "ABR"=observed.ABR))
    }
}

