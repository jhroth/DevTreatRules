#' Build a Treatment Rule
#'
#' Perform principled development of a treatment rule (using the IPW approach to account for potential confounding) on a development dataset (i.e. training set) that is independent of datasets used for model selection (i.e. validation set) and rule evaluation (i.e. test set).
#'
#' @param development.data A data frame representing the *development* dataset (i.e. training set) used for building a treatment rule.
#' @param study.design Either `observational', `RCT', or `naive'. For the \code{observational} design, the function uses inverse-probability-of-treatment observation weights (IPW) based on estimated propensity scores with predictors \code{names.influencing.treatment}; for the \code{RCT} design, the function uses IPW based on propensity scores equal to the observed sample proportions; for the \code{naive} design, all observation weights will be uniformly equal to 1.
#' @param prediction.approach One of `split.regression', `direct.interactions', `OWL', or `OWL.framework'.
#' @param name.outcome A character indicating the name of the outcome variable in \code{development.data}.
#' @param type.outcome Either `binary' or `continuous', the form of \code{name.outcome}.
#' @param name.treatment A character indicating the name of the treatment variable in \code{development.data}.
#' @param names.influencing.treatment A character vector (or single element) indicating the names of the variables in \code{development.data} that are expected to influence treatment assignment in the current dataset. Required for \code{study.design=}`observational'.
#' @param names.influencing.rule A character vector (or single element) indicating the names of the variables in \code{development.data} that may influence response to treatment and are expected to be observed in future clinical settings.
#' @param desirable.outcome A logical equal to \code{TRUE} if higher values of the outcome are considered desirable (e.g. for a binary outcome, a 1 is more desirable than a 0). The \code{OWL.framework} and \code{OWL} prediction approaches require a desirable outcome.
#' @param rule.method One of `glm.regression', `lasso', or `ridge'. For \code{type.outcome=}`binary', `glm.regression' leads to logistic regression; for a \code{type.outcome=}`continuous', `glm.regression' specifies linear regression. This is the underlying regression model used to develop the treatment rule.
#' @param propensity.method One of `logistic.regression', `lasso', or `ridge'. This is the underlying regression model used to estimate propensity scores for \code{study.design=}`observational'.
#' @param truncate.propensity.score A logical variable dictating whether estimated propensity scores less than \code{truncate.propensity.score.threshold} away from 0 or 1 should be truncated to be no more than \code{truncate.propensity.score.threshold} away from 0 or 1.
#' @param truncate.propensity.score.threshold A numeric value between 0 and 0.25.
#' @param additional.weights A numeric vector of observation weights that will be multiplied by IPW weights in the rule development stage, with length equal to the number of rows in \code{development.data}. This can be used, for example, to account for a non-representative sampling design or to apply an IPW adjustment for missingness. The default is a vector of 1s.
#' @param type.observation.weights Default is NULL, but other choices are `IPW.L', `IPW.L.and.X', and `IPW.ratio', where L indicates \code{names.influencing.treatment}, X indicates \code{names.influencing.rule}. The default behavior is to use the `IPW.ratio' observation weights (propensity score based on X divided by propensity score based on L and X) for \code{prediction.approach=}`split.regression' and to use `IPW.L' observation weights (inverse of propensity score based on L) for the `direct.interactions', `OWL', and `OWL.framework' prediction approaches.
#' @param propensity.k.cv.folds An integer specifying how many folds to use for K-fold cross-validation that chooses the tuning parameters when \code{propensity.method} is `lasso' or `ridge'. Default is 10.
#' @param rule.k.cv.folds An integer specifying how many folds to use for K-fold cross-validation that chooses the tuning parameter when \code{rule.method} is \code{lasso} or `ridge'. Default is 10.
#' @param lambda.choice Either `min' or `1se', corresponding to the \code{s} argument in \code{predict.cv.glmnet()} from the \code{glmnet} package. Only used when \code{propensity.method} or \code{rule.method} is `lasso' or `ridge'. Default is `min'.
#' @param OWL.lambda.seq Used when \code{prediction.approach=}`OWL', a numeric vector that corresponds to the \code{lambdas} argument in the \code{owl()} function from the \code{DynTxRegime} package. Defaults to \code{2^seq(-5, 5, 1)}.
#' @param OWL.kernel Used when \code{prediction.approach=}`OWL', a character equal to either `linear' or `radial'. Corresponds to the \code{kernel} argument in the \code{owl()} function from the \code{DynTxRegime} package. Default is `linear'.
#' @param OWL.kparam.seq Used when \code{prediction.approach=}`OWL' and \code{OWL.kernel=}`radial'.  Corresponds to the \code{kparam} argument in the \code{owl()} function from the \code{DynTxRegime} package. Defaults to \code{2^seq(-10, 10, 1)}.
#' @param OWL.cvFolds Used when \code{prediction.approach=}`OWL', an integer corresponding to the \code{cvFolds} argument in the \code{owl()} function from the \code{DynTxRegime} package. Defaults to 10.
#' @param OWL.verbose Used when \code{prediction.approach=}`OWL', a logical corresponding to the \code{verbose} argument in the \code{owl()} function from the \code{DynTxRegime} package. Defaults to \code{TRUE}.
#' @param OWL.framework.shift.by.min Logical, set to \code{TRUE} by default in recognition of our empirical observation that, with a continuous outcome, OWL framework performs far better in simulation studies when the outcome was shifted to have a minimum of just above 0.
#' @param direct.interactions.center.continuous.Y Logical, set to \code{TRUE} by default in recognition of our empirical observation that, with a continuous outcome, direct-interactions performed far better in simulation studies when the outcome was mean-centered.
#' @param direct.interactions.exclude.A.from.penalty Logical, set to \code{TRUE} by default in recognition of our empirical observation that, with a continuous outcome and lasso/ridge used specified as the \code{rule.method}, direct-interactions performed far better in simulation studies when the coefficient corresponding to the treatment variable was excluded from the penalty function.
#' @return A list with some combination of the following components (depending on specified \code{prediction.approach})
#' \itemize{
#'   \item \code{type.outcome}: The \code{type.outcome} specified above (used by other functions that are based on \code{BuildRule()})
#'   \item \code{prediction.approach}: The \code{prediction.approach} specified above (used by other functions that are based on \code{BuildRule()})
#'   \item \code{rule.method}: The \code{rule.method} specified above (used by other functions that are based on \code{BuildRule()})
#'   \item \code{lambda.choice}: The \code{lambda.choice} specified above (used by other functions that are based on \code{BuildRule()})
#'   \item \code{propensity.score.object}: A list containing the relevant regression object from propensity score estimation. The list has two elements for \code{type.observation.weights=}`IPW.ratio' (the default for \code{prediction.approach=}`split.regression'), has one element for \code{type.observation.weights=}`IPW.L' (the default for `OWL', `OWL.framework' and `direct.interactions'), has one element when \code{type.observation.weights=}`IPW.L.and.X', and is simply equal to NA if \code{study.design=}`RCT' (in which case propensity score would just be the inverse of sample proportion receiving treatment).
#'   \item \code{owl.object}: For \code{prediction.approach=}`OWL' only, the object returned by the \code{owl()} function in the \code{DynTxRegime} package.
#'   \item \code{observation.weights}: The observation weights used for estimating the treatment rule
#'   \item \code{rule.object}: For \code{prediction.approach=}`OWL.framework' or \code{prediction.approach=}`direct.interactions', the regression object returned from treatment rule estimation (to which the \code{coef()} function could be applied, for example)
#'   \item \code{rule.object.control}: For \code{prediction.approach=}`split.regression' the regression object returned from treatment rule estimation (to which the \code{coef()} function could be applied, for example) that estimates the outcome variable for individuals who do not receive treatment.
#'   \item \code{rule.object.treatment}: For \code{prediction.approach=}`split.regression' the regression object returned from treatment rule estimation (to which the \code{coef()} function could be applied, for example) that estimates the outcome variable for individuals who do receive treatment.
#' }

#' @references
#' \itemize{
#' \item Yingqi Zhao, Donglin Zeng, A. John Rush & Michael R. Kosorok  (2012)
#' Estimating individualized treatment rules using outcome weighted learning.
#' Journal of the American Statistical Association,
#' 107:499 1106--1118.

#' \item Shuai Chen, Lu Tian, Tianxi Cai, Menggang Yu (2017)
#' A general statistical framework for subgroup identification and comparative treatment scoring.
#' Biometrics,
#' 73:4: 1199--1209.

#' \item Lu Tian, Ash A. Alizadeh, Andrew J. Gentles, Robert Tibshirani (2014)
#' A simple method for estimating interactions between a treatment and a large number of covariates.
#' Journal of the American Statistical Association,
#' 109:508: 1517--1532.
#' \item Jeremy Roth and Noah Simon (2019).
#' Using propensity scores to develop and evaluate treatment rules with observational data
#' (Manuscript in progress)
#' \item Jeremy Roth and Noah Simon (2019).
#' Elucidating outcome-weighted learning and its comparison to split-regression: direct vs. indirect methods in practice.
#' (Manuscript in progress)
#' }

#' @examples
#' set.seed(123)
#' example.split <- SplitData(data=obsStudyGeneExpressions,
#'                                      n.sets=3, split.proportions=c(0.5, 0.25, 0.25))
#' development.data <- example.split[example.split$partition == "development",]
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
#' coef(one.rule$rule.object.control)
#' coef(one.rule$rule.object.treatment)
#' @import glmnet
#' @import DynTxRegime
#' @export

BuildRule <- function(development.data,
                      study.design, #=c("RCT", "observational"),
                      prediction.approach, #=c("OWL", "OWL.framework", "split.regression", "direct.interactions"),
                      name.outcome,
                      type.outcome, #=c("binary", "continuous"),
                      name.treatment,
                      names.influencing.treatment=NULL,
                      names.influencing.rule,
                      desirable.outcome,
                      rule.method=NULL,
                      propensity.method,
                      additional.weights=rep(1, nrow(development.data)),
                      truncate.propensity.score=TRUE,
                      truncate.propensity.score.threshold=0.05,
                      type.observation.weights=NULL,
                      propensity.k.cv.folds=10,
                      rule.k.cv.folds=10,
                      lambda.choice=c("min", "1se"),
                      OWL.lambda.seq=NULL,
                      OWL.kernel="linear",
                      OWL.kparam.seq=NULL,
                      OWL.cvFolds=10,
                      OWL.verbose=TRUE,
                      OWL.framework.shift.by.min=TRUE,
                      direct.interactions.center.continuous.Y=TRUE,
                      direct.interactions.exclude.A.from.penalty=TRUE) {
    # LOTS of checks
    if (is.data.frame(development.data) == FALSE) {
        stop("dataset must be a data frame")
    }
    if (all(development.data[, name.treatment] %in% c(0, 1)) == FALSE) {
        stop("treatment needs to be coded as a binary indicator")
    }
    if (!(study.design %in% c("RCT", "observational", "naive"))) {
        stop("study design needs to be RCT, observational, or naive")
    }
    if (study.design == "observational") {
        if (is.null(names.influencing.treatment)) {
            stop("names.influencing.treatment needs to be provided for an observational study design")
        }
    }
    if (is.null(type.observation.weights) & study.design == "observational" & (prediction.approach %in% c("OWL", "OWL.framework"))) {
        #type.observation.weights <- "IPW.L"
        type.observation.weights <- "IPW.L.and.X"
    }
    if (is.null(type.observation.weights) & study.design == "observational" & (prediction.approach %in% c("split.regression"))) {
        type.observation.weights <- "IPW.ratio"
    }
    if (is.null(type.observation.weights) & prediction.approach %in% c("direct.interactions")) {
        type.observation.weights <- "IPW.L.and.X"
    }
    lambda.choice <- match.arg(lambda.choice)
    if (min(development.data[, name.outcome], na.rm=TRUE) < 0 & prediction.approach %in% c("OWL", "OWL.framework")) {
        stop("negative values of the outcome are not allowed when fitting OWL or  OWL framework can accommodate this")
    }
    if (is.logical(desirable.outcome) == FALSE) {
        stop("desirable.outcome has to be TRUE or FALSE")
    }
    if (desirable.outcome == FALSE & prediction.approach %in% c("OWL", "OWL.framework")) {
        stop("the OWL approach assumes that larger values of the outcome variable are better")
    }
    stopifnot((type.outcome %in% c("binary", "continuous")) == TRUE)
    if (prediction.approach == "OWL") {
        if (is.null(rule.method) == FALSE) {
            stop("rule.method argument should be NULL when prediction.approach='OWL'")
        }
    } else {
        if (is.null(rule.method) == TRUE || !(rule.method %in% c("glm.regression", "lasso", "ridge"))) {
            stop("rule.method must be one of glm.regression, lasso, ridge")
        }
        if (type.outcome == "binary" & rule.method == "glm.regression") {
            rule.method <- "logistic.regression"
        }
        if (type.outcome == "continuous" & rule.method == "glm.regression" & (prediction.approach %in% c("split.regression", "direct.interactions"))) {
            rule.method <- "linear.regression"
        }
        if (type.outcome == "continuous" & rule.method == "glm.regression" & (prediction.approach %in% c("OWL", "OWL.framework"))) {
            rule.method <- "logistic.regression"
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
    stopifnot(length(additional.weights) == nrow(development.data))
    stopifnot(is.numeric(additional.weights))
    # format the data frame and return different model.matrix objects (as data frames and matrices)
    my.formatted.data <- FormatData(data=development.data,
                                                    name.outcome=name.outcome,
                                                    name.treatment=name.treatment,
                                                    type.outcome=type.outcome,
                                                    names.influencing.rule=names.influencing.rule,
                                                    names.influencing.treatment=names.influencing.treatment)
    n <- nrow(development.data)
    idx.control <- which(my.formatted.data$df.model.matrix.all[, "treatment"] == 0)
    idx.treatment <- which(my.formatted.data$df.model.matrix.all[, "treatment"] == 1)
    if (study.design == "observational") {
        if (prediction.approach %in% c("direct.interactions")) {
            ## response.for.propensity.score <- "fac.treatment.neg.pos"
            ## data.matrix.for.propensity.score <- my.formatted.data$model.matrix.X.times.A.plus.AY
            ## data.df.for.propensity.score <- my.formatted.data$df.model.matrix.X.times.A.plus.AY
            response.for.propensity.score <- "fac.treatment"
            data.matrix.for.propensity.score <- my.formatted.data$model.matrix.all
            data.df.for.propensity.score <- my.formatted.data$df.model.matrix.all
        } else if (prediction.approach %in% c("OWL", "OWL.framework", "split.regression")) {
            response.for.propensity.score <- "fac.treatment"
            data.matrix.for.propensity.score <- my.formatted.data$model.matrix.all
            data.df.for.propensity.score <- my.formatted.data$df.model.matrix.all
        } else {
            stop("invalid prediction.approach")
        }
        # Predict P(T=1 | L)
        propensity.score.L.object <- DoPrediction(data.matrix=data.matrix.for.propensity.score,
                                                         data.df=data.df.for.propensity.score,
                                                         name.response=response.for.propensity.score,
                                                         type.response="binary",
                                                         names.features=names(my.formatted.data$df.model.matrix.L),
                                                         observation.weights=rep(1, n),
                                                         method=propensity.method,
                                                         lambda.choice=lambda.choice,
                                                         k.cv.folds=propensity.k.cv.folds,
                                                         include.intercept=TRUE)
        propensity.score.L.probability <- TruncateProbability(probability=propensity.score.L.object$one.fit.predicted.probability,
                                                                               threshold=truncate.propensity.score.threshold)
        if (type.observation.weights == "IPW.ratio") {
       # Predict P(T=1 | L, X)
        propensity.score.L.and.X.object <- DoPrediction(data.matrix=data.matrix.for.propensity.score,
                                                         data.df=data.df.for.propensity.score,
                                                         name.response=response.for.propensity.score,
                                                         type.response="binary",
                                                         names.features=names(my.formatted.data$df.model.matrix.L.and.X),
                                                         observation.weights=rep(1, n),
                                                         method=propensity.method,
                                                         lambda.choice=lambda.choice,
                                                         k.cv.folds=propensity.k.cv.folds,
                                                         include.intercept=TRUE)
        propensity.score.L.and.X.probability <- TruncateProbability(probability=propensity.score.L.and.X.object$one.fit.predicted.probability,
                                                                                        threshold=truncate.propensity.score.threshold)
            # Predict P(T=1 | X)
            propensity.score.X.object <- DoPrediction(data.matrix=data.matrix.for.propensity.score,
                                                       data.df=data.df.for.propensity.score,
                                                       name.response=response.for.propensity.score,
                                                       type.response="binary",
                                                       names.features=names(my.formatted.data$df.model.matrix.X),
                                                       observation.weights=rep(1, n),
                                                       method=propensity.method,
                                                       lambda.choice=lambda.choice,
                                                       k.cv.folds=propensity.k.cv.folds,
                                                       include.intercept=TRUE)
            propensity.score.X.probability <- TruncateProbability(probability=propensity.score.X.object$one.fit.predicted.probability,
                                                                                   threshold=truncate.propensity.score.threshold)
            obs.weights <- rep(NA, n)
            obs.weights[idx.treatment] <- (propensity.score.X.probability / propensity.score.L.and.X.probability)[idx.treatment]
            obs.weights[idx.control] <- ((1 - propensity.score.X.probability) / (1 - propensity.score.L.and.X.probability))[idx.control]
            obs.weights <- additional.weights * obs.weights
            propensity.score.object <- list(propensity.score.X.object, propensity.score.L.and.X.object)
        } else if (type.observation.weights == "IPW.L.and.X") {
        # Predict P(T=1 | L, X)
            propensity.score.L.and.X.object <- DoPrediction(data.matrix=data.matrix.for.propensity.score,
                                                         data.df=data.df.for.propensity.score,
                                                         name.response=response.for.propensity.score,
                                                         type.response="binary",
                                                         names.features=names(my.formatted.data$df.model.matrix.L.and.X),
                                                         observation.weights=rep(1, n),
                                                         method=propensity.method,
                                                         lambda.choice=lambda.choice,
                                                         k.cv.folds=propensity.k.cv.folds,
                                                         include.intercept=TRUE)
        ## propensity.score.L.and.X.object <- DoPrediction(data.matrix=data.matrix.for.propensity.score,
        ##                                                  data.df=data.df.for.propensity.score,
        ##                                                  name.response=response.for.propensity.score,
        ##                                                  type.response="binary",
        ##                                                  names.features=names(my.formatted.data$df.model.matrix.L.and.X),
        ##                                                  observation.weights=rep(1, n),
        ##                                                  method=propensity.method,
        ##                                                  lambda.choice=lambda.choice,
        ##                                                  k.cv.folds=propensity.k.cv.folds,
        ##                                                  include.intercept=TRUE)
            propensity.score.L.and.X.probability <- TruncateProbability(probability=propensity.score.L.and.X.object$one.fit.predicted.probability,
                                                                                        threshold=truncate.propensity.score.threshold)
            obs.weights <- rep(NA, n)
            obs.weights[idx.treatment] <- (1 / propensity.score.L.and.X.probability)[idx.treatment]
            obs.weights[idx.control] <- (1 / (1 - propensity.score.L.and.X.probability))[idx.control]
            obs.weights <- additional.weights * obs.weights
            propensity.score.object <- list(propensity.score.L.and.X.object)
        } else if (type.observation.weights == "IPW.L") {
            obs.weights <- rep(NA, n)
            obs.weights[idx.treatment] <- (1 / propensity.score.L.probability)[idx.treatment]
            obs.weights[idx.control] <- (1 / (1 - propensity.score.L.probability))[idx.control]
            obs.weights <- additional.weights * obs.weights
            propensity.score.object <- list(propensity.score.L.object)
        } else {
            stop("invalid type of observation weighting specified")
        }
    } else if (study.design == "RCT") {
        propensity.score.empirical <- mean(my.formatted.data$df.model.matrix.all[, "treatment"] == 1)
        obs.weights <- rep(NA, n)
        obs.weights[idx.treatment] <- 1 / propensity.score.empirical
        obs.weights[idx.control] <- 1 / (1 - propensity.score.empirical)
        obs.weights <- additional.weights * obs.weights
        propensity.score.object <- NA
    } else if (study.design == "naive") {
        obs.weights <- rep(1, n)
        obs.weights <- additional.weights * obs.weights
        propensity.score.object <- NA
    }
     # do weighted regression or classification problem using IPW weights
    if (prediction.approach == "OWL") {
        if (study.design == "observational" & type.observation.weights == "IPW.L") {
            predicted.propensity.score.larger.model <- propensity.score.L.probability
        } else if (study.design == "observational" & type.observation.weights == "IPW.L.and.X") {
            predicted.propensity.score.larger.model <- propensity.score.L.and.X.probability
        } else if (study.design == "observational" & type.observation.weights == "IPW.ratio") {
            print("might not want to use IPW ratio here")
        } else if (study.design == "RCT") {
            predicted.propensity.score.larger.model <- NULL
        } else {
            stop("invalid type.observation.weights chosen for observational study design")
        }
        predict.T.with.X.object <- OneDynTxRegime(data.matrix=my.formatted.data$model.matrix.all,
                                                                     data.df=my.formatted.data$df.model.matrix.all,
                                                                     study.design=study.design,
                                                                     names.X=names(my.formatted.data$df.model.matrix.X),
                                                                     name.outcome="outcome", # changed by FormatData()
                                                                     name.treatment="treatment", # changed by FormatData()
                                                                     DynTxRegime.method="OWL",
                                                                     predicted.propensity.score.larger.model=predicted.propensity.score.larger.model,
                                                                     lambda.seq=OWL.lambda.seq,
                                                                     kernel=OWL.kernel,
                                                                     kparam.seq=OWL.kparam.seq,
                                                                     cvFolds=OWL.cvFolds,
                                                                     OWL.verbose=OWL.verbose)
        warning("the OWL approach assumes that larger values of the outcome variable are better")
        return(list("type.outcome"=type.outcome,
                       "propensity.score.object"=propensity.score.object,
                       "prediction.approach"=prediction.approach,
                      "owl.object"=predict.T.with.X.object))
    } else if (prediction.approach == "OWL.framework") {
        if (is.null(rule.method)) {
            stop("need to specify rule method for the OWL framework approach")
        }
        response.for.OWL.framework <- "fac.treatment"
        if (OWL.framework.shift.by.min == TRUE) {
            OWL.weights <- (my.formatted.data$df.model.matrix.all[, "outcome"] - abs(min(my.formatted.data$df.model.matrix.all[, "outcome"])) * 0.999) * obs.weights
        } else {
            OWL.weights <- my.formatted.data$df.model.matrix.all[, "outcome"] * obs.weights
        }
        if (rule.method %in% c("linear.regression")) {
            print("you specified linear.regression as the rule.method, but OWL.framework is predicting the binary treatment indicator. as a result, we changed the rule.method to logistic.regression just for fitting the OWL.framework rule")
            rule.method.for.OWL.framework <- "logistic.regression"
        } else {
            rule.method.for.OWL.framework <- rule.method
        }
        predict.T.with.X.object <- DoPrediction(data.matrix=my.formatted.data$model.matrix.all,
                                                 data.df=my.formatted.data$df.model.matrix.all,
                                                 name.response=response.for.OWL.framework,
                                                 type.response="binary",
                                                 names.features=names(my.formatted.data$df.model.matrix.X),
                                                 observation.weights=OWL.weights,
                                                 method=rule.method.for.OWL.framework,
                                                 lambda.choice=lambda.choice,
                                                 k.cv.folds=rule.k.cv.folds,
                                                 include.intercept=TRUE)
        warning("the OWL framework approach assumes that larger values of the outcome variable are better")
        return(list("type.outcome"=type.outcome,
                       "propensity.score.object"=propensity.score.object,
                       "observation.weights"=OWL.weights,
                       "prediction.approach"=prediction.approach,
                       "lambda.choice"=lambda.choice,  #03/25 added for consistency
                       "rule.method"=rule.method,
                      "rule.object"=predict.T.with.X.object$one.fit))
    } else if (prediction.approach == "split.regression") {
        if (is.null(rule.method)) {
            stop("need to specify rule method for the split regression approach")
        }
        if (type.outcome == "binary") {
            name.outcome.for.split.regression <-"outcome.fac"
        } else {
            name.outcome.for.split.regression <-"outcome"
        }
        predict.Y.with.X.object.control <- DoPrediction(data.matrix=my.formatted.data$model.matrix.all[idx.control, ],
                                                         data.df=my.formatted.data$df.model.matrix.all[idx.control, ],
                                                         name.response=name.outcome.for.split.regression,
                                                         type.response=type.outcome,
                                                         names.features=names(my.formatted.data$df.model.matrix.X),
                                                         observation.weights=obs.weights[idx.control],
                                                         method=rule.method,
                                                         lambda.choice=lambda.choice,
                                                         k.cv.folds=rule.k.cv.folds,
                                                         include.intercept=TRUE)
        predict.Y.with.X.object.treatment <- DoPrediction(data.matrix=my.formatted.data$model.matrix.all[idx.treatment, ],
                                                           data.df=my.formatted.data$df.model.matrix.all[idx.treatment, ],
                                                           name.response=name.outcome.for.split.regression,
                                                           type.response=type.outcome,
                                                           names.features=names(my.formatted.data$df.model.matrix.X),
                                                           observation.weights=obs.weights[idx.treatment],
                                                           method=rule.method,
                                                           lambda.choice=lambda.choice,
                                                           k.cv.folds=rule.k.cv.folds,
                                                           include.intercept=TRUE)
        return(list("type.outcome"=type.outcome,
                      "propensity.score.object"=propensity.score.object,
                      "observation.weights"=obs.weights,
                      "prediction.approach"=prediction.approach,
                      "rule.method"=rule.method,
                      "lambda.choice"=lambda.choice,
                      "rule.object.control"=predict.Y.with.X.object.control$one.fit,
                      "rule.object.treatment"=predict.Y.with.X.object.treatment$one.fit))
        } else if (prediction.approach == "direct.interactions") {
            if (is.null(rule.method)) {
                stop("need to specify rule method for the interactions approach")
            }
            if (type.outcome == "binary") {
                name.outcome.for.interactions.approach <-"outcome.fac"
            } else {
                if (direct.interactions.center.continuous.Y == TRUE) {
                    name.outcome.for.interactions.approach <-"outcome.centered"
                } else {
                    name.outcome.for.interactions.approach <- "outcome"
                }
            }
            predict.Y.with.X.times.A <- DoPrediction(data.matrix=my.formatted.data$model.matrix.X.times.A,
                                                      data.df=my.formatted.data$df.model.matrix.X.times.A.plus.AY,
                                                      name.response=name.outcome.for.interactions.approach,
                                                      type.response=type.outcome,
                                                      names.features=colnames(my.formatted.data$model.matrix.X.times.A), # already includes "treatment.neg.pos, which was handled in FormatData()
                                                      observation.weights=obs.weights,
                                                      method=rule.method,
                                                      lambda.choice=lambda.choice,
                                                      k.cv.folds=rule.k.cv.folds,
                                                      include.intercept=FALSE,
                                                      exclude.A.from.penalty=direct.interactions.exclude.A.from.penalty)
            #print(coef(predict.Y.with.X.times.A$one.fit))
            return(list("type.outcome"=type.outcome,
                        "propensity.score.object"=propensity.score.object,
                        "observation.weights"=obs.weights,
                        "prediction.approach"=prediction.approach,
                        "rule.method"=rule.method,
                        "lambda.choice"=lambda.choice,
                        "rule.object"=predict.Y.with.X.times.A$one.fit,
                        "predict.Y.with.X.times.A"=predict.Y.with.X.times.A))
        }
}
