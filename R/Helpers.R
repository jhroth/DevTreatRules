Logit <- function(x) log(x / (1 - x))
Expit <- function(x) exp(x) / (1 + exp(x))
NamesToGlmFormula <- function(name.response, names.features, include.intercept) {
    if (include.intercept == TRUE) {
        result <- as.formula(paste0(name.response, " ~ ", paste0(names.features, collapse=" + ")))
    } else if (include.intercept == FALSE) {
        result <- as.formula(paste0(name.response, " ~ 0 +", paste0(names.features, collapse=" + ")))
    }
    return(result)
}

ComputeABR <- function(n.test.positives, ATE.test.positives,
                                     n.test.negatives, ATE.test.negatives) {
    if (n.test.positives == 0) {
        result <- -1 * ATE.test.negatives
    } else if (n.test.negatives == 0) {
        result <- ATE.test.positives
    } else {
        result <- (n.test.positives / (n.test.positives + n.test.negatives)) * ATE.test.positives +
                    (n.test.negatives / (n.test.positives + n.test.negatives)) * (-1 * ATE.test.negatives)
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
                                        n.lambda.seq=NULL,
                                        kernel=c("linear", "radial"),
                                        kparam.seq=NULL,
                                        cvFolds=10,
                                        OWL.verbose=TRUE) {
    kernel <- match.arg(kernel)
    if (is.null(n.lambda.seq)) {
        n.lambda.seq <- 10
    }
    if (is.null(lambda.seq)) {
        lambda.seq <- 2 ^ seq(from=-5, to=5, by=1) # from Yingqi; using from -10 to 10 gave me more errors in simulation scenario 
    }
    if (kernel == "linear") {
        kparam.seq <- NULL
    } else if (is.null(kparam.seq) & kernel == "radial") {
        kparam.seq <- 2 ^ seq(from=-10, to=10, by=1) # from Yingqi
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
                                  b.cv.repeats,
                                  lambda.choice,
                                  include.intercept,
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
        one.fit <- glm(my.glm.formula,
                            family="quasibinomial", 
                            data=df.with.weights,
                            weights=observation.weights)
        one.fit.predicted.probability <- as.numeric(predict(one.fit, type="response"))
        one.fit.predicted.class <- as.numeric(one.fit.predicted.probability >= 0.5)
        return(list("one.fit"=one.fit,
                      "one.fit.predicted.probability"=one.fit.predicted.probability,
                      "one.fit.predicted.class"=one.fit.predicted.class))
    } else if (method == "linear.regression") {
        if (type.response == "binary") {
            warning("linear regression was used for a binary response variable")
        }
        df.with.weights <- data.frame(data.df,
                                           "observation.weights"=observation.weights)
        one.fit <- glm(my.glm.formula,
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
    model.matrix.L <- model.matrix(as.formula(paste("~", names.influencing.treatment, collapse="+", sep=" ")), data=data)
    if (ncol(model.matrix.L) == 2) {
        model.matrix.L <- as.matrix(model.matrix.L[, -1, drop=FALSE]) 
        colnames(model.matrix.L) <- names.influencing.treatment
    } else {
        model.matrix.L <- model.matrix.L[, -1, drop=FALSE]
    }
    df.model.matrix.L <- as.data.frame(model.matrix.L)
    ## X (variables influencing rule)
    model.matrix.X <- model.matrix(as.formula(paste("~", names.influencing.rule, collapse="+", sep=" ")), data=data)
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
    model.matrix.L.and.X <- model.matrix(as.formula(paste("~", c(names.influencing.treatment, names.influencing.rule), collapse="+", sep=" ")), data=data)[, -1, drop=FALSE]
    df.model.matrix.L.and.X <- as.data.frame(model.matrix.L.and.X)
    model.matrix.all <- model.matrix(as.formula(paste("~", c(name.outcome, name.treatment, names.influencing.treatment, names.influencing.rule), collapse="+", sep=" ")), data=data)[, -1, drop=FALSE]
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
    model.matrix.all.times.A <- cbind(model.matrix.L, model.matrix.X.times.A)
    df.model.matrix.L <- as.data.frame(model.matrix.L)
    df.model.matrix.X.times.A <- as.data.frame(model.matrix.X.times.A)
    if (type.outcome == "binary") {
        df.model.matrix.all.times.A <- data.frame("outcome"=data[, name.outcome], "outcome.fac"=outcome.fac,
                                                             "fac.treatment.neg.pos"=fac.treatment.neg.pos,
                                                              df.model.matrix.L,
                                                              df.model.matrix.X.times.A)
    } else {
        df.model.matrix.all.times.A <- data.frame("outcome"=data[, name.outcome], "outcome.centered"=outcome.centered,
                                                                 "fac.treatment.neg.pos"=fac.treatment.neg.pos, 
                                                                 df.model.matrix.L,
                                                                 df.model.matrix.X.times.A)
    }
    return(list("model.matrix.X.times.A"=model.matrix.X.times.A, "df.model.matrix.X.times.A"=df.model.matrix.X.times.A,
                  "model.matrix.all.times.A"=model.matrix.all.times.A, "df.model.matrix.all.times.A"=df.model.matrix.all.times.A,
                  "model.matrix.L"=model.matrix.L, "df.model.matrix.L"=df.model.matrix.L,
                  "model.matrix.X"=model.matrix.X, "df.model.matrix.X"=df.model.matrix.X,
                  "model.matrix.L.and.X"=model.matrix.L.and.X, "df.model.matrix.L.and.X"=df.model.matrix.L.and.X,
                  "model.matrix.all"=model.matrix.all, "df.model.matrix.all"=df.model.matrix.all))
}


