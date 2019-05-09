#' @import DynTxRegime
#' @import modelObj


Logit <- function(x) log(x / (1 - x))
Expit <- function(x) exp(x) / (1 + exp(x))
NamesToGlmFormula <- function(name.response, names.features, include.intercept) {
    if (include.intercept == TRUE) {
        result <- stats::as.formula(paste0(name.response, " ~ ", paste0(names.features, collapse=" + ")))
    } else if (include.intercept == FALSE) {
        result <- stats::as.formula(paste0(name.response, " ~ 0 +", paste0(names.features, collapse=" + ")))
    }
    return(result)
}

ComputeABR <- function(n.positives, ATE.positives,
                                     n.negatives, ATE.negatives) {
    if (n.positives == 0) {
        result <- -1 * ATE.negatives
    } else if (n.negatives == 0) {
        result <- ATE.positives
    } else {
        result <- (n.positives / (n.positives + n.negatives)) * ATE.positives +
                    (n.negatives / (n.positives + n.negatives)) * (-1 * ATE.negatives)
    }
    return(result)
}

OneCvGlmnet <- function(method, features, response, observation.weights, my.family, my.metric, lambda.seq, k.cv.folds, include.intercept, my.penalty.factor) {
    tryCatch({
        lambda.seq <- NULL
        if (method == "ridge") {
            my.alpha <- 0
        }
        else if (method == "lasso") {
            my.alpha <- 1
        }
        result <- cv.glmnet(x=features,
                            y=response,
                            weights=observation.weights,
                            family=my.family,
                            type.measure=my.metric,
                            lambda=lambda.seq,
                            nfolds=k.cv.folds,
                            alpha=my.alpha,
                            intercept=include.intercept,
                            penalty.factor=my.penalty.factor)
        return(result)
    }, error = function(e) {
        if (method == "ridge") {
            ridge.lambda.max <- max(svd(scale(features))$d)^2
            lambda.seq <- exp(seq(from=log(0.001*ridge.lambda.max), to=log(ridge.lambda.max), length.out=100))
            my.alpha <- 0
        } else if (method == "lasso") {
            if (is.factor(response) == TRUE) {
                stopifnot(length(levels(response)) == 2)
                response.mat <- matrix(NA, nrow=length(response), ncol=1)
                response.mat[response == levels(response)[1], 1] <-  0
                response.mat[response == levels(response)[2], 1] <-  1
                response.mat <- as.matrix(response.mat)
                lasso.lambda.max <- max(abs(t(as.matrix(features)) %*% response.mat))
                lambda.seq <- exp(seq(from=log(0.001*lasso.lambda.max), to=log(lasso.lambda.max), length.out=100))
                my.alpha <- 1
            } else {
                lasso.lambda.max <- max(abs(t(as.matrix(scale(features))) %*% as.vector(scale(response))))
                lambda.seq <- exp(seq(from=log(0.001*lasso.lambda.max), to=log(lasso.lambda.max), length.out=100))
            }
        }
        result <- cv.glmnet(x=features,
                            y=response,
                            weights=observation.weights,
                            family=my.family,
                            type.measure=my.metric,
                            lambda=lambda.seq,
                            nfolds=k.cv.folds,
                            alpha=my.alpha,
                            intercept=include.intercept,
                            penalty.factor=my.penalty.factor)
        warning("there had been an error performing using cv.glmnet, but it was fixed by pre-specifying a lambda sequence")
        return(result)
    })
}

OneDynTxRegime <- function(data.matrix,
                                        data.df,
                                        study.design, #c("RCT", "observational")
                                        names.X,
                                        name.outcome,
                                        name.treatment,
                                        predicted.propensity.score.larger.model=NULL,
                                        DynTxRegime.method="OWL", # c("OWL", "earl")
                                        lambda.seq=NULL,
                                        kernel=c("linear", "radial"),
                                        kparam.seq=NULL,
                                        cvFolds=10,
                                        OWL.verbose=TRUE) {
    kernel <- match.arg(kernel)
    if (is.null(lambda.seq)) {
        lambda.seq <- 2 ^ seq(from=-5, to=5, by=1) 
    }
    if (kernel == "linear") {
        kparam.seq <- NULL
    } else if (is.null(kparam.seq) & kernel == "radial") {
        kparam.seq <- 2 ^ seq(from=-10, to=10, by=1)
    } 
    if (DynTxRegime.method == "OWL") {
        if (study.design == "observational") {
            # re-define outcome for observational study
            get.propensity.model.for.owl <- RedefineOutcomeForOWL(data.df=data.df,
                                                                                     name.outcome=name.outcome,
                                                                                     name.treatment=name.treatment,
                                                                                     predicted.propensity.score.larger.model=predicted.propensity.score.larger.model)
            data.df.for.owl <- data.frame(data.df,
                                                    "outcome.redefined"=get.propensity.model.for.owl$outcome.redefined)
            one.fit <- owl(moPropen=get.propensity.model.for.owl$one.ModelObj.glm.constant.only,
                           kernel=kernel,
                           data=data.df.for.owl,
                           reward=data.df.for.owl[, "outcome.redefined"],
                           txName="treatment",
                           regime=NamesToGlmFormula(name.response="", names.features=names.X, include.intercept=FALSE),
                           lambdas=lambda.seq,
                           kparam=kparam.seq,
                           cvFolds=cvFolds,
                           verbose=OWL.verbose)
        } else if (study.design == "RCT") {
            n <- length(data.df[, name.outcome])
            one.ModelObj.glm.constant.only <- buildModelObj(model= ~ 1,
                                                                                 solver.method="glm",
                                                                                 solver.args=list("weights"=rep(1, n),
                                                                                                      "family"="binomial"),
                                                                                 predict.method="predict.glm",
                                                                                 predict.args=list(type="response"))
            data.df.for.owl <- data.df
            one.fit <- owl(moPropen=one.ModelObj.glm.constant.only,
                           kernel=kernel,
                           data=data.df.for.owl,
                           reward=data.df.for.owl[, "outcome"],
                           txName="treatment",
                           regime=NamesToGlmFormula(name.response="", names.features=names.X, include.intercept=FALSE),
                           lambdas=lambda.seq,
                           kparam.seq=kparam.seq,
                           cvFolds=cvFolds,
                           verbose=OWL.verbose)
        }
        return(list("one.fit"=one.fit))
    } else if (DynTxRegime.method == "earl") {
        stop("Yingqi said that earl is also an appropriate alternative to OWL that is available in DynTxRegime")
    }
}

RedefineOutcomeForOWL <- function(data.df,
                                                  name.outcome,
                                                  name.treatment,
                                                  predicted.propensity.score.larger.model) {
    n <- length(data.df[, name.outcome])
    # fit constant-only GLM model
    one.ModelObj.glm.constant.only <- buildModelObj(model= ~ 1,
                                                                         solver.method="glm",
                                                                         solver.args=list("weights"=rep(1, n),
                                                                                               "family"="binomial"),
                                                                         predict.method="predict.glm",
                                                                         predict.args=list(type="response"))
    fit.propensity.score.glm.constant.only <- fit(object=one.ModelObj.glm.constant.only,
                                                                        data=data.df,
                                                                        response=data.df[, name.treatment])
    predicted.propensity.score.glm.constant.only <- as.numeric(predict(fit.propensity.score.glm.constant.only, newdata=data.df))
    original.outcome <- data.df[, name.outcome]
    treatment <- data.df[, name.treatment]
    outcome.redefined <- rep(NA, n)
    outcome.redefined[treatment == 0] <- original.outcome[treatment == 0] *
                                                         ((1 - predicted.propensity.score.glm.constant.only) / (1 - predicted.propensity.score.larger.model))[treatment == 0]
    outcome.redefined[treatment == 1] <- original.outcome[treatment == 1] *
                                                         (predicted.propensity.score.glm.constant.only / predicted.propensity.score.larger.model)[treatment == 1]
    return(list("outcome.redefined"=outcome.redefined,
             "one.ModelObj.glm.constant.only"=one.ModelObj.glm.constant.only))
}

TruncateProbability <- function(probability, threshold=0.05) {
    if (min(probability, na.rm=TRUE) < 0 | max(probability, na.rm=TRUE) > 1) {
        stop("the probability argument given to TruncateProbability is not a valid probability")
    }
    stopifnot(min(probability) > 0)
    stopifnot(is.numeric(threshold))
    if (threshold > 0.25 | threshold < 0) {
        stop("threshold should be between 0 and 0.25")
    }
    result <- probability
    result[probability < threshold] <- threshold
    result[probability > (1-threshold)] <- (1 - threshold)
    return(result)
}

DoPrediction <- function(data.matrix,
                                  data.df,
                                  name.response,
                                  type.response,
                                  names.features,
                                  observation.weights,
                                  method,
                                  k.cv.folds,
                                  lambda.choice,
                                  include.intercept,
                                  lambda.seq=NULL,
                                  exclude.A.from.penalty=FALSE) {
    stopifnot(is.logical(exclude.A.from.penalty))
    my.glm.formula <- NamesToGlmFormula(name.response=name.response,
                                                             names.features=names.features,
                                                             include.intercept=include.intercept)
    if (method %in% c("lasso", "ridge")) {
        features <- data.matrix[, names.features]
        response <- data.df[, name.response]
        if (type.response == "binary") {
            my.family <- "binomial"
            my.metric <- "auc"
        } else if (type.response == "continuous") {
            my.family <- "gaussian"
            my.metric <- "mse"
        }
        if (exclude.A.from.penalty == TRUE) {
            my.penalty.factor <- rep(NA, ncol(features))
            idx.treatment.neg.pos <- colnames(features) %in% "treatment.neg.pos"
            my.penalty.factor[idx.treatment.neg.pos] <- 0
            my.penalty.factor[!idx.treatment.neg.pos] <- 1
            stopifnot(all(!is.na(my.penalty.factor)))
            stopifnot(0 %in% my.penalty.factor)
        } else {
            my.penalty.factor <- rep(1, ncol(features))
        }
        one.cv.glmnet <- OneCvGlmnet(method=method, features=features, response=response, observation.weights=observation.weights,
                                                         my.family=my.family, my.metric=my.metric, lambda.seq=lambda.seq, k.cv.folds=k.cv.folds,
                                                         include.intercept=include.intercept, my.penalty.factor=my.penalty.factor)
        if (lambda.choice == "min") {
            optimal.lambda <- "lambda.min"
        } else if (lambda.choice == "1se") {
            optimal.lambda <- "lambda.1se"
        }
        if (type.response == "binary") {
            one.fit.predicted.probability <- as.numeric(predict.cv.glmnet(one.cv.glmnet, newx=features, s=optimal.lambda, type="response"))
            one.fit.predicted.class <- as.numeric(one.fit.predicted.probability >= 0.5)
            return(list("one.fit"=one.cv.glmnet,
                          "one.fit.predicted.probability"=one.fit.predicted.probability,
                          "one.fit.predicted.class"=one.fit.predicted.class))
        } else if (type.response == "continuous") {
            one.fit.predicted.response <- as.numeric(predict.cv.glmnet(one.cv.glmnet, newx=features, s=optimal.lambda, type="response"))
            return(list("one.fit"=one.cv.glmnet,
                          "one.fit.predicted.response"=one.fit.predicted.response))
        }
    } else if (method == "logistic.regression") {
        if (type.response != "binary") {
            stop("logistic regression is only appropriate for a binary resposne variable")
        }
        df.with.weights <- data.frame(data.df,
                                           "observation.weights"=observation.weights)
        #print(summary(df.with.weights))
        #print(my.glm.formula)
        one.fit <- stats::glm(my.glm.formula,
                            family="quasibinomial", 
                            data=df.with.weights,
                            weights=observation.weights)
        #print(coef(one.fit))
        one.fit.predicted.probability <- as.numeric(predict(one.fit, type="response"))
        one.fit.predicted.class <- as.numeric(one.fit.predicted.probability >= 0.5)
        return(list("one.fit"=one.fit,
                      "one.fit.predicted.probability"=one.fit.predicted.probability,
                      "one.fit.predicted.class"=one.fit.predicted.class))
#                       "df.with.weights"=df.with.weights))
    } else if (method == "linear.regression") {
        if (type.response == "binary") {
            warning("linear regression was used for a binary response variable")
        }
        df.with.weights <- data.frame(data.df,
                                           "observation.weights"=observation.weights)
        one.fit <- stats::glm(my.glm.formula,
                            family="gaussian",
                            data=df.with.weights,
                            weights=observation.weights)
        one.fit.predicted.response <- as.numeric(predict(one.fit, type="response"))
        return(list("one.fit"=one.fit,
                      "one.fit.predicted.response"=one.fit.predicted.response))
    } else {
        print(paste0("specified method was:", method))
        stop("invalid fitting method specified")
    }
}

FormatData <- function(data,
                                 name.outcome,
                                 name.treatment,
                                 type.outcome,
                                 names.influencing.treatment,
                                 names.influencing.rule) {
    n <- nrow(data)
    fac.treatment <- rep(NA, n)
    fac.treatment[data[, name.treatment] == 0] <- "no_treatment"
    fac.treatment[data[, name.treatment] == 1] <- "treatment"
    fac.treatment <- as.factor(fac.treatment)
    stopifnot(all(!is.na(fac.treatment)))
    treatment.neg.pos <- rep(NA, n)
    treatment.neg.pos[data[, name.treatment] == 0] <- -1
    treatment.neg.pos[data[, name.treatment] == 1] <- 1
    stopifnot(all(treatment.neg.pos %in% c(-1, 1)))
    fac.treatment.neg.pos <- rep(NA, n)
    fac.treatment.neg.pos[treatment.neg.pos == -1] <- "no_treatment"
    fac.treatment.neg.pos[treatment.neg.pos == 1] <- "treatment"
    fac.treatment.neg.pos <- as.factor(fac.treatment.neg.pos)
    if (type.outcome == "binary") {
        outcome.fac <- rep(NA, n)
        outcome.fac[data[, name.outcome] == 0] <- "no_outcome"
        outcome.fac[data[, name.outcome] == 1] <- "outcome"
        outcome.fac <- as.factor(outcome.fac)
    } else if (type.outcome == "continuous") {
        ## accomodating centered outcome
        outcome.centered <- data[, name.outcome] - mean(data[, name.outcome])
    }
    # Create model matrix objects
    ## L (variables influencing treatment)
    model.matrix.L <- stats::model.matrix(stats::as.formula(paste("~", names.influencing.treatment, collapse="+", sep=" ")), data=data)
    #print("orig names of model.matrix.L:")
    #print(names(model.matrix.L))
    if (ncol(model.matrix.L) == 2) {
        #names.model.matrix.L <- model.matrix.L
        model.matrix.L <- as.matrix(model.matrix.L[, -1, drop=FALSE]) 
        #colnames(model.matrix.L) <- names.model.matrix.L[2]
    } else {
        model.matrix.L <- model.matrix.L[, -1, drop=FALSE]
    }
    df.model.matrix.L <- as.data.frame(model.matrix.L)
    ## X (variables influencing rule)
    model.matrix.X <- stats::model.matrix(stats::as.formula(paste("~", names.influencing.rule, collapse="+", sep=" ")), data=data)
    #print("orig names of model.matrix.X:")
    #print(names(model.matrix.X))
    if (ncol(model.matrix.X) == 2) {
        model.matrix.X <- as.matrix(model.matrix.X[, -1, drop=FALSE]) 
        colnames(model.matrix.X) <- names.influencing.rule
    } else {
        model.matrix.X <- model.matrix.X[, -1, drop=FALSE] 
    }
    df.model.matrix.X <- as.data.frame(model.matrix.X)
    ## L and X (variables influencing treatment and rule)
    if (length(c(names.influencing.treatment, names.influencing.rule)) <= 1){
        stop("need at least two variables influencing treatment and rule combined")
    }
    model.matrix.L.and.X <- stats::model.matrix(stats::as.formula(paste("~", c(names.influencing.treatment, names.influencing.rule), collapse="+", sep=" ")), data=data)[, -1, drop=FALSE]
    df.model.matrix.L.and.X <- as.data.frame(model.matrix.L.and.X)
    model.matrix.all <- stats::model.matrix(stats::as.formula(paste("~", c(name.outcome, name.treatment, names.influencing.treatment, names.influencing.rule), collapse="+", sep=" ")), data=data)[, -1, drop=FALSE]
    if (type.outcome == "binary") {
        df.model.matrix.all <- data.frame("outcome"=data[, name.outcome], "outcome.fac"=outcome.fac,
                                                      "treatment"=data[, name.treatment], "fac.treatment"=fac.treatment,
                                                       df.model.matrix.L.and.X)
    } else {
        df.model.matrix.all <- data.frame("outcome"=data[, name.outcome], 
                                                       "outcome.centered"=outcome.centered, 
                                                       "treatment"=data[, name.treatment], "fac.treatment"=fac.treatment,
                                                       df.model.matrix.L.and.X)
    }
    model.matrix.X.times.A <- cbind(treatment.neg.pos, model.matrix.X * treatment.neg.pos)
    ## Need to handle this very carefulyl when a variable is a part of L and also a part of X; if i just combine with cbind(model.matrix.L, model.matrix.X.times.A), then the variable will be represented twice, and the earlier column will not be multiplied by A!
        #model.matrix.all.times.A <- cbind(model.matrix.L, model.matrix.X.times.A)
    #model.matrix.X.times.A.plus.AY.potential.dup <- cbind(model.matrix.X.times.A, model.matrix.L)
    #model.matrix.X.times.A.plus.AY <- model.matrix.X.times.A.plus.AY.potential.dup[, !duplicated(colnames(model.matrix.X.times.A.plus.AY.potential.dup))]
    df.model.matrix.L <- as.data.frame(model.matrix.L)
    df.model.matrix.X.times.A <- as.data.frame(model.matrix.X.times.A)
    if (type.outcome == "binary") {
        df.model.matrix.X.times.A.plus.AY <- data.frame("outcome"=data[, name.outcome], "outcome.fac"=outcome.fac,
                                                             "fac.treatment.neg.pos"=fac.treatment.neg.pos,
                                                              #df.model.matrix.L,
                                                              df.model.matrix.X.times.A)
        ## df.model.matrix.X.times.A.plus.AY <- data.frame("outcome"=data[, name.outcome], "outcome.fac"=outcome.fac,
        ##                                                      "fac.treatment.neg.pos"=fac.treatment.neg.pos,
        ##                                                       df.model.matrix.L,
        ##                                                       df.model.matrix.X.times.A)
    } else {
        df.model.matrix.X.times.A.plus.AY <- data.frame("outcome"=data[, name.outcome], "outcome.centered"=outcome.centered,
                                                                 "fac.treatment.neg.pos"=fac.treatment.neg.pos, 
                                                                #df.model.matrix.L,
                                                                 df.model.matrix.X.times.A)
        ## df.model.matrix.X.times.A.plus.AY <- data.frame("outcome"=data[, name.outcome], "outcome.centered"=outcome.centered,
        ##                                                          "fac.treatment.neg.pos"=fac.treatment.neg.pos, 
        ##                                                          df.model.matrix.L,
        ##                                                          df.model.matrix.X.times.A)
    }
    return(list("model.matrix.X.times.A"=model.matrix.X.times.A, "df.model.matrix.X.times.A"=df.model.matrix.X.times.A,
                  #"model.matrix.X.times.A.plus.AY"=model.matrix.X.times.A.plus.AY,
                  "df.model.matrix.X.times.A.plus.AY"=df.model.matrix.X.times.A.plus.AY,
                  "model.matrix.L"=model.matrix.L, "df.model.matrix.L"=df.model.matrix.L,
                  "model.matrix.X"=model.matrix.X, "df.model.matrix.X"=df.model.matrix.X,
                  "model.matrix.L.and.X"=model.matrix.L.and.X, "df.model.matrix.L.and.X"=df.model.matrix.L.and.X,
                  "model.matrix.all"=model.matrix.all, "df.model.matrix.all"=df.model.matrix.all))
}


EvaluateRuleOnce <- function(data,
                                  my.formatted.data,
                                  BuildRule.object=NULL,
                                  B=NULL, 
                                  study.design, #=c("RCT", "observational"),
                                  type.outcome,
                                  desirable.outcome,
                                  separate.propensity.estimation,
                                  clinical.threshold=0,
                                  names.influencing.treatment,
                                  names.influencing.rule,
                                  propensity.method="logistic.regression",
                                  observation.weights,
                                  additional.weights=rep(1, nrow(data)),
                                  lambda.choice=c("min", "1se"),
                                  bootstrap.CI=FALSE, 
                                  propensity.k.cv.folds=10,
                                  truncate.propensity.score,
                                  truncate.propensity.score.threshold) {
    # estimate propensity score in this dataset (same idea as when we estimated propensity score as a part of build_rule())
    data.df <- my.formatted.data$df.model.matrix.all
    n <- nrow(data.df)
    idx.control <- which(data.df[, "treatment"] == 0) # remember name of treatment variable is created by FormatData()
    idx.treatment <- which(data.df[, "treatment"] == 1)
    if (is.null(BuildRule.object) & is.null(B)) {
        stop("either BuildRule.object or B has to be specified")
    }
    stopifnot(propensity.method %in% c("logistic.regression"))
    if (is.null(BuildRule.object) == FALSE) {
        B <- PredictRule(BuildRule.object=BuildRule.object,
                          new.X=data[, names.influencing.rule, drop=FALSE],
                          desirable.outcome=desirable.outcome,
                          clinical.threshold=clinical.threshold,
                          return.predicted.response=FALSE)
    }
    idx.test.positives <- B==1
    n.positives <- sum(idx.test.positives)
    idx.test.negatives <- B==0
    n.negatives <- sum(idx.test.negatives)
    if (is.null(observation.weights) == TRUE) {
        if (study.design == "observational") {
            p.for.propensity <- ncol(my.formatted.data$df.model.matrix.L)
            if (p.for.propensity > n.positives & separate.propensity.estimation == TRUE) {
                if ((propensity.method %in% c("logistic.regression", "linear.regression")) & separate.propensity.estimation == TRUE) {
                    #warning("Within the test-positives subset, there are fewer observations than predictors in the propensity score model. Since the specified propensity.method does not perform variable selection, the separate.propensity.estimation argument has been changed to FALSE to estimate propensity score in the pooled sample of test-positives and test-negatives")
                    separate.propensity.estimation <- FALSE
                } else if ((propensity.method %in% c("lasso", "ridge")) & separate.propensity.estimation == TRUE) {
                    #warning("Within the test-positives subset, there are fewer observations than predictors in the propensity score model. The specified propensity.method does perform variable selection so estimation may still be possible, but if this is unexpected it may be better to set separate.propensity.estimation=FALSE to perform this estimation in the pooled sample of test-positives and test-negatives.")
                }
            }
            if (p.for.propensity > n.negatives & separate.propensity.estimation == TRUE) {
                if ((propensity.method %in% c("logistic.regression", "linear.regression")) & separate.propensity.estimation == TRUE) {
                    #warning("Within the test-negatives subset, there are fewer observations than predictors in the propensity score model. Since the specified propensity.method does not perform variable selection, the separate.propensity.estimation argument has been changed to FALSE to estimate propensity score in the pooled sample of test-negatives and test-negatives")
                    separate.propensity.estimation <- FALSE
                } else if ((propensity.method %in% c("lasso", "ridge")) & separate.propensity.estimation == TRUE) {
                    #warning("Within the test-negatives subset, there are fewer observations than predictors in the propensity score model. The specified propensity.method does perform variable selection so estimation is still possible, but if this is unexpected it may be better to set separate.propensity.estimation=FALSE to perform this estimation in the pooled sample of test-negatives and test-negatives.")
                }
            }
            n.positives.treatment <- sum(my.formatted.data$df.model.matrix.all[idx.test.positives, "treatment"] == 1)
            n.positives.control <- sum(my.formatted.data$df.model.matrix.all[idx.test.positives, "treatment"] == 0)
            n.negatives.treatment <- sum(my.formatted.data$df.model.matrix.all[idx.test.negatives, "treatment"] == 1)
            n.negatives.control <- sum(my.formatted.data$df.model.matrix.all[idx.test.negatives, "treatment"] == 0)
            if (n.positives.treatment < 2 | n.positives.control < 2 | n.negatives.treatment < 2 | n.negatives.control <= 2) {
                separate.propensity.estimation <- FALSE
                if (n.positives.treatment < 2) {
                    paste("there were fewer than two test-positive observations in the treatment group, so the separate.propensity.estimation argument has been changed to FALSE to estimate propensity score in the pooled sample of test-positives and test-negatives")
                }
                if (n.positives.control < 2) {
                    paste("there were fewer than two test-positive observations in the control group, so the separate.propensity.estimation argument has been changed to FALSE to estimate propensity score in the pooled sample of test-positives and test-negatives")
                }
                if (n.negatives.treatment < 2) {
                    paste("there were fewer than two test-negative observations in the treatment group, so the separate.propensity.estimation argument has been changed to FALSE to estimate propensity score in the pooled sample of test-positives and test-negatives")
                }
                if (n.negatives.control < 2) {
                    paste("there were fewer than two test-negative observations in the control group, so the separate.propensity.estimation argument has been changed to FALSE to estimate propensity score in the pooled sample of test-positives and test-negatives")
                }
            }
            if (separate.propensity.estimation == TRUE) {
                # Predict P(T=1 | L, B=1)
                #print("names for eval propensity score")
                #print(names(my.formatted.data$df.model.matrix.all[idx.test.positives, ]))
                propensity.score.L.object.test.positives <- DoPrediction(data.matrix=my.formatted.data$model.matrix.all[idx.test.positives, ],
                                                                          data.df=my.formatted.data$df.model.matrix.all[idx.test.positives, ],
                                                                          name.response="fac.treatment",
                                                                          type.response="binary",
                                                                          names.features=names(my.formatted.data$df.model.matrix.L),
                                                                          observation.weights=rep(1, n.positives),
                                                                          method=propensity.method,
                                                                          lambda.choice=lambda.choice,
                                                                          k.cv.folds=propensity.k.cv.folds,
                                                                          include.intercept=TRUE,
                                                                          exclude.A.from.penalty=FALSE)
                propensity.score.L.test.positives.probability <- TruncateProbability(probability=propensity.score.L.object.test.positives$one.fit.predicted.probability,
                                                                                                         threshold=truncate.propensity.score.threshold)
                idx.test.positives.treatment <- my.formatted.data$df.model.matrix.all[idx.test.positives, "treatment"] == 1
                idx.test.positives.control <- my.formatted.data$df.model.matrix.all[idx.test.positives, "treatment"] == 0
                
                # Predict P(T=1 | L, B=0)
                propensity.score.L.object.test.negatives <- DoPrediction(data.matrix=my.formatted.data$model.matrix.all[idx.test.negatives, ],
                                                                          data.df=my.formatted.data$df.model.matrix.all[idx.test.negatives, ],
                                                                          name.response="fac.treatment",
                                                                          type.response="binary",
                                                                          names.features=names(my.formatted.data$df.model.matrix.L),
                                                                          observation.weights=rep(1, n.negatives),
                                                                          method=propensity.method,
                                                                          lambda.choice=lambda.choice,
                                                                          k.cv.folds=propensity.k.cv.folds,
                                                                          include.intercept=TRUE,
                                                                          exclude.A.from.penalty=FALSE)
                propensity.score.L.test.negatives.probability <- TruncateProbability(probability=propensity.score.L.object.test.negatives$one.fit.predicted.probability,
                                                                                                          threshold=truncate.propensity.score.threshold)
                idx.test.negatives.treatment <- my.formatted.data$df.model.matrix.all[idx.test.negatives, "treatment"] == 1
                idx.test.negatives.control <- my.formatted.data$df.model.matrix.all[idx.test.negatives, "treatment"] == 0
                
                obs.weights <- rep(NA, n)
                obs.weights[idx.test.positives][idx.test.positives.treatment] <- (1 / propensity.score.L.test.positives.probability)[idx.test.positives.treatment]
                obs.weights[idx.test.positives][idx.test.positives.control] <- (1 / (1 - propensity.score.L.test.positives.probability))[idx.test.positives.control]
                obs.weights[idx.test.negatives][idx.test.negatives.treatment] <- (1 / propensity.score.L.test.negatives.probability)[idx.test.negatives.treatment]
                obs.weights[idx.test.negatives][idx.test.negatives.control] <- (1 / (1 - propensity.score.L.test.negatives.probability))[idx.test.negatives.control]
                obs.weights <- additional.weights * obs.weights
                fit.object <- list("one.fit.test.positives"=propensity.score.L.object.test.positives$one.fit, "one.fit.test.negatives"=propensity.score.L.object.test.negatives$one.fit)
            } else {
                # Predict P(T=1 | L)
                propensity.score.L.object <- DoPrediction(data.matrix=my.formatted.data$model.matrix.all,
                                                           data.df=my.formatted.data$df.model.matrix.all,
                                                           name.response="fac.treatment",
                                                           type.response="binary",
                                                           names.features=names(my.formatted.data$df.model.matrix.L),
                                                           observation.weights=rep(1, n),
                                                           method=propensity.method,
                                                           lambda.choice=lambda.choice,
                                                           k.cv.folds=propensity.k.cv.folds,
                                                           include.intercept=TRUE,
                                                           exclude.A.from.penalty=FALSE)
                propensity.score.L.probability <- TruncateProbability(probability=propensity.score.L.object$one.fit.predicted.probability,
                                                                                       threshold=truncate.propensity.score.threshold)
                obs.weights <- rep(NA, n)
                obs.weights[idx.treatment] <- (1 / propensity.score.L.probability)[idx.treatment]
                obs.weights[idx.control] <- (1 / (1 - propensity.score.L.probability))[idx.control]
                obs.weights <- additional.weights * obs.weights
                fit.object <- list("one.fit"=propensity.score.L.object$one.fit)
            }
        } else if (study.design == "RCT") {
            propensity.score.empirical <- mean(my.formatted.data$df.model.matrix.all[, "treatment"] == 1)
            obs.weights <- rep(NA, n)
            obs.weights[idx.treatment] <- 1 / propensity.score.empirical
            obs.weights[idx.control] <- 1 / (1 - propensity.score.empirical)
            obs.weights <- additional.weights * obs.weights
            print(paste0("constant being used as propensity score for all observations:", propensity.score.empirical))
            fit.object <- list("one.fit"=propensity.score.empirical)
        }
    } else {
        obs.weights <- additional.weights * observation.weights
        fit.object <- NULL
        print("user-specified observation weights being used to evaluate treatment rule")
    }
    ## ATE among test positives and test negatives (OWL always assumes higher values of the outcoem variable are better)
    weighted.outcome <- obs.weights * data.df[, "outcome"]
    if (n.positives > 0) {
        ATE.positives <- sum(weighted.outcome[B==1 & data.df[, "treatment"] == 1]) / n.positives -
                                      sum(weighted.outcome[B==1 & data.df[, "treatment"] == 0]) / n.positives
    } else {
        ATE.positives <- NA
    }
    if (n.negatives > 0) {
        ATE.negatives <- sum(weighted.outcome[B==0 & data.df[, "treatment"] == 1]) / n.negatives -
                                       sum(weighted.outcome[B==0 & data.df[, "treatment"] == 0]) / n.negatives
    } else {
        ATE.negatives <- NA
    }
    ABR <- ComputeABR(n.positives=n.positives, ATE.positives=ATE.positives,
                                   n.negatives=n.negatives, ATE.negatives=ATE.negatives)
    return(list("fit.object"=fit.object,
                 "recommended.treatment"=B,
                "ABR"=ABR,
                "n.positives"=n.positives,
                "ATE.positives"=ATE.positives,
                "n.negatives"=n.negatives,
                "ATE.negatives"=ATE.negatives))
}
