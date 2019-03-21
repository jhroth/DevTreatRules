BuildRule <- function(data,
                             study.design, #=c("RCT", "observational"),
                             prediction.approach, #=c("OWL", "OWL.framework", "split.regression", "direct.interactions"), 
                             name.outcome,
                             type.outcome, #=c("binary", "continuous"),
                             name.treatment,
                             names.influencing.treatment,
                             names.influencing.rule,
                             additional.weights=rep(1, nrow(data)),
                             desirable.outcome,
                             propensity.method, 
                             truncate.propensity.score=TRUE,
                             truncate.propensity.score.threshold=0.05, 
                             rule.method=NULL, 
                             type.observation.weights=NULL,
                             propensity.k.cv.folds=10,
                             propensity.b.cv.repeats=1,
                             rule.k.cv.folds=10,
                             rule.b.cv.repeats=1,
                             lambda.choice=c("min", "1se"),
                             OWL.lambda.seq=NULL,
                             OWL.n.lambda.seq=NULL,
                             OWL.kernel="linear",
                             OWL.kparam.seq=NULL,
                             OWL.cvFolds=10,
                             OWL.verbose=TRUE,
                             OWL.framework.shift.by.min=TRUE,
                             direct.interactions.center.continuous.Y=TRUE,
                             direct.interactions.exclude.A.from.penalty=TRUE) {
    # LOTS of checks
    if (is.data.frame(data) == FALSE) {
        stop("dataset must be a data frame")
    }
    if (all(data[, name.treatment] %in% c(0, 1)) == FALSE) {
        stop("treatment needs to be coded as a binary indicator")
    }
    if (!(study.design %in% c("RCT", "observational", "naive"))) {
        stop("study design needs to be RCT, observational, or naive")
    }
    if (is.null(type.observation.weights) & study.design == "observational" & (prediction.approach %in% c("OWL", "OWL.framework"))) {
        type.observation.weights <- "IPW.L"
    }
    if (is.null(type.observation.weights) & study.design == "observational" & (prediction.approach %in% c("split.regression"))) {
        type.observation.weights <- "IPW.ratio"
    }
    if (is.null(type.observation.weights) & prediction.approach %in% c("direct.interactions")) {
        type.observation.weights <- "IPW.L"
    }
    lambda.choice <- match.arg(lambda.choice)
    if (min(data[, name.outcome], na.rm=TRUE) < 0 & prediction.approach %in% c("OWL", "OWL.framework")) {
        stop("negative values of the outcome are not allowed when fitting OWL or  OWL framework can accomodate this")
    }
    if (is.logical(desirable.outcome) == FALSE) {
        stop("desirable.outcome has to be TRUE or FALSE")
    }
    if (desirable.outcome == FALSE & prediction.approach %in% c("OWL", "OWL.framework")) {
        stop("the OWL approach assumes that larger values of the outcome variable are better")
    }
    stopifnot((type.outcome %in% c("binary", "continuous")) == TRUE)
    stopifnot(rule.method %in% c("glm.regression", "lasso", "ridge"))
    if (type.outcome == "binary" & rule.method == "glm.regression") {
        rule.method <- "logistic.regression"
    }
    if (type.outcome == "continuous" & rule.method == "glm.regression" & (prediction.approach %in% c("split.regression", "direct.interactions"))) {
        rule.method <- "linear.regression"
    }
    if (type.outcome == "continuous" & rule.method == "glm.regression" & (prediction.approach %in% c("OWL", "OWL.framework"))) {
        rule.method <- "logistic.regression"
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
    # format the data frame and return different model.matrix objects (as data frames and matrices)
    my.formatted.data <- FormatData(data=data,
                                                    name.outcome=name.outcome,
                                                    name.treatment=name.treatment,
                                                    type.outcome=type.outcome,
                                                    names.influencing.rule=names.influencing.rule,
                                                    names.influencing.treatment=names.influencing.treatment)
    n <- nrow(data)
    idx.control <- which(my.formatted.data$df.model.matrix.all[, "treatment"] == 0)
    idx.treatment <- which(my.formatted.data$df.model.matrix.all[, "treatment"] == 1)
    if (study.design == "observational") {
        if (prediction.approach %in% c("direct.interactions")) {
            response.for.propensity.score <- "fac.treatment.neg.pos"
            data.matrix.for.propensity.score <- my.formatted.data$model.matrix.all.times.A
            data.df.for.propensity.score <- my.formatted.data$df.model.matrix.all.times.A
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
                                                         b.cv.repeats=propensity.b.cv.repeats,
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
                                                         b.cv.repeats=propensity.b.cv.repeats,
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
                                                       b.cv.repeats=propensity.b.cv.repeats,
                                                       include.intercept=TRUE)
            propensity.score.X.probability <- TruncateProbability(probability=propensity.score.X.object$one.fit.predicted.probability,
                                                                                   threshold=truncate.propensity.score.threshold)
            obs.weights <- rep(NA, n)
            obs.weights[idx.treatment] <- (propensity.score.X.probability / propensity.score.L.and.X.probability)[idx.treatment]
            obs.weights[idx.control] <- ((1 - propensity.score.X.probability) / (1 - propensity.score.L.and.X.probability))[idx.control]
            obs.weights <- additional.weights * obs.weights
            propensity.score.object <- list(propensity.score.X.object, propensity.score.L.and.X.object)
        } else if (type.observation.weights == "IPW.L.and.X") {
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
                                                                     n.lambda.seq=OWL.n.lambda.seq,
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
                                                 k.cv.folds=propensity.k.cv.folds,
                                                 b.cv.repeats=propensity.b.cv.repeats,
                                                 include.intercept=TRUE)
        warning("the OWL framework approach assumes that larger values of the outcome variable are better")
        return(list("type.outcome"=type.outcome,
                       "propensity.score.object"=propensity.score.object,
                       "observation.weights"=OWL.weights,
                       "prediction.approach"=prediction.approach,
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
                                                         k.cv.folds=propensity.k.cv.folds,
                                                         b.cv.repeats=propensity.b.cv.repeats,
                                                         include.intercept=TRUE)
        predict.Y.with.X.object.treatment <- DoPrediction(data.matrix=my.formatted.data$model.matrix.all[idx.treatment, ],
                                                           data.df=my.formatted.data$df.model.matrix.all[idx.treatment, ],
                                                           name.response=name.outcome.for.split.regression,
                                                           type.response=type.outcome,
                                                           names.features=names(my.formatted.data$df.model.matrix.X),
                                                           observation.weights=obs.weights[idx.treatment],
                                                           method=rule.method,
                                                           lambda.choice=lambda.choice,
                                                           k.cv.folds=propensity.k.cv.folds,
                                                           b.cv.repeats=propensity.b.cv.repeats,
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
            predict.Y.with.X.times.A <- DoPrediction(data.matrix=my.formatted.data$model.matrix.all.times.A,
                                                      data.df=my.formatted.data$df.model.matrix.all.times.A,
                                                      name.response=name.outcome.for.interactions.approach,
                                                      type.response=type.outcome,
                                                      names.features=names(my.formatted.data$df.model.matrix.X.times.A), # already includes "treatment.neg.pos, which was handled in FormatData()
                                                      observation.weights=obs.weights,
                                                      method=rule.method,
                                                      lambda.choice=lambda.choice,
                                                      k.cv.folds=propensity.k.cv.folds,
                                                      b.cv.repeats=propensity.b.cv.repeats,
                                                      include.intercept=FALSE,
                                                      exclude.A.from.penalty=direct.interactions.exclude.A.from.penalty) 
            ## can NOT use fitted values from the above call, but i think that's OK. i'll handle them correctly in the PredictRule() function
            return(list("type.outcome"=type.outcome,
                        "propensity.score.object"=propensity.score.object,
                        "observation.weights"=obs.weights,
                        "prediction.approach"=prediction.approach,
                        "rule.method"=rule.method,
                        "lambda.choice"=lambda.choice,
                        "rule.object"=predict.Y.with.X.times.A$one.fit))
        }
}
