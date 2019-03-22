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
                                  propensity.b.cv.repeats=1,
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
                                        #new.X=data.df[, names.influencing.rule],
                          desirable.outcome=desirable.outcome,
                          clinical.threshold=clinical.threshold)
    }
    idx.test.positives <- B==1
    n.test.positives <- sum(idx.test.positives)
    idx.test.negatives <- B==0
    n.test.negatives <- sum(idx.test.negatives)
    if (is.null(observation.weights) == TRUE) {
        if (study.design == "observational") {
            p.for.propensity <- ncol(my.formatted.data$df.model.matrix.L)
            if (p.for.propensity > n.test.positives & separate.propensity.estimation == TRUE) {
                if ((propensity.method %in% c("logistic.regression", "linear.regression")) & separate.propensity.estimation == TRUE) {
                    warning("Within the test-positives subset, there are fewer observations than predictors in the propensity score model. Since the specified propensity.method does not perform variable selection, the separate.propensity.estimation argument has been changed to FALSE to estimate propensity score in the pooled sample of test-positives and test-negatives")
                    separate.propensity.estimation <- FALSE
                } else if ((propensity.method %in% c("lasso", "ridge")) & separate.propensity.estimation == TRUE) {
                    warning("Within the test-positives subset, there are fewer observations than predictors in the propensity score model. The specified propensity.method does perform variable selection so estimation may still be possible, but if this is unexpected it may be better to set separate.propensity.estimation=FALSE to perform this estimation in the pooled sample of test-positives and test-negatives.")
                }
            }
            if (p.for.propensity > n.test.negatives & separate.propensity.estimation == TRUE) {
                if ((propensity.method %in% c("logistic.regression", "linear.regression")) & separate.propensity.estimation == TRUE) {
                    warning("Within the test-negatives subset, there are fewer observations than predictors in the propensity score model. Since the specified propensity.method does not perform variable selection, the separate.propensity.estimation argument has been changed to FALSE to estimate propensity score in the pooled sample of test-negatives and test-negatives")
                    separate.propensity.estimation <- FALSE
                } else if ((propensity.method %in% c("lasso", "ridge")) & separate.propensity.estimation == TRUE) {
                    warning("Within the test-negatives subset, there are fewer observations than predictors in the propensity score model. The specified propensity.method does perform variable selection so estimation is still possible, but if this is unexpected it may be better to set separate.propensity.estimation=FALSE to perform this estimation in the pooled sample of test-negatives and test-negatives.")
                }
            }
            n.test.positives.treatment <- sum(my.formatted.data$df.model.matrix.all[idx.test.positives, "treatment"] == 1)
            n.test.positives.control <- sum(my.formatted.data$df.model.matrix.all[idx.test.positives, "treatment"] == 0)
            n.test.negatives.treatment <- sum(my.formatted.data$df.model.matrix.all[idx.test.negatives, "treatment"] == 1)
            n.test.negatives.control <- sum(my.formatted.data$df.model.matrix.all[idx.test.negatives, "treatment"] == 0)
            if (n.test.positives.treatment < 2 | n.test.positives.control < 2 | n.test.negatives.treatment < 2 | n.test.negatives.control <= 2) {
                separate.propensity.estimation <- FALSE
                if (n.test.positives.treatment < 2) {
                    paste("there were fewer than two test-positive observations in the treatment group, so the separate.propensity.estimation argument has been changed to FALSE to estimate propensity score in the pooled sample of test-positives and test-negatives")
                }
                if (n.test.positives.control < 2) {
                    paste("there were fewer than two test-positive observations in the control group, so the separate.propensity.estimation argument has been changed to FALSE to estimate propensity score in the pooled sample of test-positives and test-negatives")
                }
                if (n.test.negatives.treatment < 2) {
                    paste("there were fewer than two test-negative observations in the treatment group, so the separate.propensity.estimation argument has been changed to FALSE to estimate propensity score in the pooled sample of test-positives and test-negatives")
                }
                if (n.test.negatives.control < 2) {
                    paste("there were fewer than two test-negative observations in the control group, so the separate.propensity.estimation argument has been changed to FALSE to estimate propensity score in the pooled sample of test-positives and test-negatives")
                }
            }
            if (separate.propensity.estimation == TRUE) {
                # Predict P(T=1 | L, B=1)
                propensity.score.L.object.test.positives <- DoPrediction(data.matrix=my.formatted.data$model.matrix.all[idx.test.positives, ],
                                                                          data.df=my.formatted.data$df.model.matrix.all[idx.test.positives, ],
                                                                          name.response="fac.treatment",
                                                                          type.response="binary",
                                                                          names.features=names(my.formatted.data$df.model.matrix.L),
                                                                          observation.weights=rep(1, n.test.positives),
                                                                          method=propensity.method,
                                                                          lambda.choice=lambda.choice,
                                                                          k.cv.folds=propensity.k.cv.folds,
                                                                          b.cv.repeats=propensity.b.cv.repeats,
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
                                                                          observation.weights=rep(1, n.test.negatives),
                                                                          method=propensity.method,
                                                                          lambda.choice=lambda.choice,
                                                                          k.cv.folds=propensity.k.cv.folds,
                                                                          b.cv.repeats=propensity.b.cv.repeats,
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
                                                           b.cv.repeats=propensity.b.cv.repeats,
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
    if (n.test.positives > 0) {
        ATE.test.positives <- sum(weighted.outcome[B==1 & data.df[, "treatment"] == 1]) / n.test.positives -
                                      sum(weighted.outcome[B==1 & data.df[, "treatment"] == 0]) / n.test.positives
    } else {
        ATE.test.positives <- NA
    }
    if (n.test.negatives > 0) {
        ATE.test.negatives <- sum(weighted.outcome[B==0 & data.df[, "treatment"] == 1]) / n.test.negatives -
                                       sum(weighted.outcome[B==0 & data.df[, "treatment"] == 0]) / n.test.negatives
    } else {
        ATE.test.negatives <- NA
    }
    ABR <- ComputeABR(n.test.positives=n.test.positives, ATE.test.positives=ATE.test.positives,
                                   n.test.negatives=n.test.negatives, ATE.test.negatives=ATE.test.negatives)
    return(list("fit.object"=fit.object,
                 "recommended.treatment"=B,
                "ABR"=ABR,
                "n.test.positives"=n.test.positives,
                "ATE.test.positives"=ATE.test.positives,
                "n.test.negatives"=n.test.negatives,
                "ATE.test.negatives"=ATE.test.negatives))
}
