PredictRule <- function(BuildRule.object,
                                new.X,
                                desirable.outcome=NULL,
                                clinical.threshold=0,
                                return.predicted.response=FALSE) {
#    if (is.numeric(new.X)) {
    if (is.data.frame(new.X) == FALSE) {
        stop("new.X must be a data frame")
    }
    prediction.approach <- BuildRule.object$prediction.approach
    model.matrix.new.X <- model.matrix(as.formula(paste("~", colnames(new.X), collapse="+", sep=" ")), data=new.X)[, -1, drop=FALSE]
    df.model.matrix.new.X <- data.frame(model.matrix.new.X)
    rm(new.X)
    if (prediction.approach == "OWL") {
        one.result <- optTx(BuildRule.object$owl.object$one.fit,
                                   newdata=df.model.matrix.new.X)$optimalTx 
        recommended.treatment <- rep(NA, nrow(df.model.matrix.new.X))
        recommended.treatment[one.result == -1] <- 0
        recommended.treatment[one.result == 1] <- 1
        warning("the OWL approach assumes that larger values of the outcome variable are better")
    } else if (prediction.approach %in% c("split.regression", "direct.interactions")) {
        if (is.null(desirable.outcome)) {
            stop("desirable.outcome needs to be specified when using split.regression or direct.interactions approach")
        }
        stopifnot(clinical.threshold >=0)
        if (BuildRule.object$rule.method %in% c("lasso", "ridge")) {
            if (BuildRule.object$lambda.choice == "min") {
                optimal.lambda <- "lambda.min"
            } else if (BuildRule.object$lambda.choice == "1se") {
                optimal.lambda <- "lambda.1se"
            }
            if (prediction.approach %in% c("split.regression")) {
                predicted.outcome.under.treatment <- as.numeric(predict.cv.glmnet(BuildRule.object$rule.object.treatment, newx=model.matrix.new.X, s=optimal.lambda, type="response"))
                predicted.outcome.under.control <- as.numeric(predict.cv.glmnet(BuildRule.object$rule.object.control, newx=model.matrix.new.X, s=optimal.lambda, type="response"))
            } else if (prediction.approach %in% c("direct.interactions")) {
                my.mat.X <- cbind(1, model.matrix.new.X)
                colnames(my.mat.X)[1] <- "treatment.neg.pos"
                predicted.outcome <- as.numeric(predict.cv.glmnet(BuildRule.object$rule.object, newx=my.mat.X, s=optimal.lambda, type="response"))
            }
        } else if (BuildRule.object$rule.method %in% c("logistic.regression","linear.regression")) {
            if (prediction.approach %in% c("split.regression")) {
                predicted.outcome.under.control <- as.numeric(predict(BuildRule.object$rule.object.control, df.model.matrix.new.X, type="response"))
                predicted.outcome.under.treatment <- as.numeric(predict(BuildRule.object$rule.object.treatment, df.model.matrix.new.X, type="response"))
            } else if (prediction.approach %in% c("direct.interactions")) {
                predicted.outcome <- as.numeric(predict(BuildRule.object$rule.object,
                                                                       data.frame("treatment.neg.pos"=1, df.model.matrix.new.X),
                                                                       type="response"))
            }
        }  else {
            stop("if more methods are added I'll need to handle predicted values appropriately")
        }
        if (prediction.approach %in% c("split.regression")) {
            if (desirable.outcome == TRUE) {
                recommended.treatment <- as.numeric(predicted.outcome.under.treatment > predicted.outcome.under.control + clinical.threshold)
            } else {
                recommended.treatment <- as.numeric(predicted.outcome.under.treatment + clinical.threshold < predicted.outcome.under.control)
            }
        } else if (prediction.approach %in% c("direct.interactions")) {
            if (clinical.threshold != 0) {
                stop("non-zero clinical threshold is not supported for interactions.approach")
            }
            if (BuildRule.object$type.outcome == "binary") {
                cutoff <- 0.5
            } else if (BuildRule.object$type.outcome == "continuous") {
                cutoff <- 0
            }
            if (desirable.outcome == TRUE) {
                recommended.treatment <- as.numeric(predicted.outcome > cutoff + 0)
            } else {
                recommended.treatment <- as.numeric(predicted.outcome < cutoff + 0)
            }
        }
    } else if (prediction.approach == "OWL.framework") {
        if (BuildRule.object$rule.method %in% c("lasso", "ridge")) {
            OWL.framework.predicted.response.prob <- predict(BuildRule.object$rule.object, model.matrix.new.X, type="response")
            recommended.treatment <- as.numeric(OWL.framework.predicted.response.prob > 0.5)
        } else if (BuildRule.object$rule.method %in% c("logistic.regression","linear.regression")) {
            OWL.framework.predicted.response.prob <- predict(BuildRule.object$rule.object, df.model.matrix.new.X, type="response")
            recommended.treatment <- as.numeric(OWL.framework.predicted.response.prob > 0.5)
        } else {
            stop("if more methods are added I'll need to handle predicted values appropriately")
        }
        warning("the OWL approach assumes that larger values of the outcome variable are better")
    } else {
        stop("prediction approach needs to be one of split.regression, OWL, OWL.framework, and interactions.approach")
    }
    if (return.predicted.response == FALSE) {
        return(recommended.treatment)
    } else if (return.predicted.response == TRUE & prediction.approach == "direct.interactions") {
        return(list("recommended.treatment"=recommended.treatment,
                      "predicted.outcome"=predicted.outcome))
    } else if (return.predicted.response == TRUE & prediction.approach == "split.regression") {
        return(list("predicted.outcome.under.control"=predicted.outcome.under.control,
                    "predicted.outcome.under.treatment"=predicted.outcome.under.treatment,
                     "recommended.treatment"=recommended.treatment))
    } else if (return.predicted.response == TRUE & prediction.approach == "OWL.framework") {
        return(list("predicted.outcome"=OWL.framework.predicted.response.prob,
                      "recommended.treatment"=recommended.treatment))
    }
}
