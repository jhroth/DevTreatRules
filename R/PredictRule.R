#' Get the treatment rule implied by \code{BuildRule()}
#' 
#' Map the object returned by \code{BuildRule()} to the treatment rule corresponding to a particular dataset
#'
#' @param BuildRule.object The object returned by the \code{BuildRule()} function. 
#' @param new.X A data frame representing the dataset for which the treatment rule is desired.
#' @param desirable.outcome A logical equal to \code{TRUE} if higher values of the outcome are considered desirable (e.g. for a binary outcome, 1 is more desirable than 0). The \code{OWL.framework} and \code{OWL} approaches to treatment rule estimation require a desirable outcome.
#' @param clinical.threshold A numeric equal a positive number above which the predicted outcome under treatment must be superior to the predicted outcome under control for treatment to be recommended. Only used when \code{BuildRuleObject} was specified and derived from the split-regression or direct-interactions approach. Defaults to 0.
#' @param return.predicted.response logical indicating whether the predicted response variable (for \code{split.regression}, \code{OWL.framework}, and \code{OWL} approaches) or score function (for \code{direct.interactions}) should be returned in addition to its mapping to a binary treatment recommendation. Default is \code{FALSE}.
#' @return
#' \itemize{
#' \item If \code{return.predicted.response=FALSE} (the default), then the single object returned is a numeric vector of 0s and 1s, with length equal to the number of rows in \code{new.X}, where a 0 indicates treatment is not recommended and a 1 indicates treatment is recommended for the corresponding observation in \code{new.X}.
#' \item If \code{return.predicted.response=TRUE}, then the object returned is a list with some combination of the following components (depending on which prediction approach underlies the \code{BuildRule.object}).
#' \itemize{
#'   \item \code{recommended.treatment}: A numeric vector of 0s and 1s, with length equal to the number of rows in \code{new.X}, where a 0 indicates treatment is not recommended and a 1 indicates treatment is recommended for the corresponding observation in \code{new.X}.
#'   \item \code{predicted.outcome}: A numeric vector showing the predicted values of the score function mapped to \code{recommended.treatment}. Only returned if \code{return.predicted.response=TRUE} and the approach underlying \code{BuildRule.object}. was `direct.interactions'.
#'   \item \code{predicted.outcome.under.control}: A numeric vector showing the predicted values of the outcome under no treatment which, along with \code{predicted.outcome.under.treatment}, corresponds to \code{recommended.treatment}. Only returned if \code{return.predicted.response=TRUE} and the approach underlying \code{BuildRule.object}. was `split.regression'.
#'   \item \code{predicted.outcome.under.treatment}: A numeric vector showing the predicted values of the outcome under treatment which, along with \code{predicted.outcome.under.control}, corresponds to \code{recommended.treatment}. Only returned if \code{return.predicted.response=TRUE} and the approach underlying \code{BuildRule.object}. was `split.regression'.
#'   \item \code{predicted.treatment.prob}: A numeric vector showing the predicted treatment probability that corresponds to \code{recommended.treatment}. Only returned if \code{return.predicted.response=TRUE} and the approach underlying \code{BuildRule.object}. was `OWL.framework'.
#' }
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
#' one.prediction <- PredictRule(BuildRule.object=one.rule,
#'                                         new.X=validation.data[, c("age", paste0("gene_", 1:10))],
#'                                         desirable.outcome=TRUE,
#'                                         clinical.threshold=0)
#' table(one.prediction)
#' @import glmnet
#' @import DynTxRegime
#' @export

PredictRule <- function(BuildRule.object,
                                new.X,
                                desirable.outcome=NULL,
                                clinical.threshold=0,
                                return.predicted.response=FALSE) {
    if (is.data.frame(new.X) == FALSE) {
        stop("new.X must be a data frame")
    }
    prediction.approach <- BuildRule.object$prediction.approach
    model.matrix.new.X <- stats::model.matrix(stats::as.formula(paste("~", colnames(new.X), collapse="+", sep=" ")), data=new.X)[, -1, drop=FALSE]
    df.model.matrix.new.X <- data.frame(model.matrix.new.X)
    rm(new.X)
    if (prediction.approach == "OWL") {
        one.result <- optTx(BuildRule.object$owl.object$one.fit,
                                   newdata=df.model.matrix.new.X)$optimalTx 
        recommended.treatment <- rep(NA, nrow(df.model.matrix.new.X))
        recommended.treatment[one.result == -1] <- 0
        recommended.treatment[one.result == 1] <- 1
        #warning("the OWL approach assumes that larger values of the outcome variable are better")
    } else if (prediction.approach %in% c("split.regression", "direct.interactions")) {
        if (is.null(desirable.outcome)) {
            stop("desirable.outcome needs to be specified when using split.regression or direct.interactions approach")
        }
        stopifnot(clinical.threshold >= 0)
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
        #warning("the OWL.framework approach assumes that larger values of the outcome variable are better")
    } else {
        stop("prediction approach needs to be one of split.regression, OWL, OWL.framework, and interactions.approach")
    }
    if (return.predicted.response == FALSE) {
        return(as.numeric(recommended.treatment))
    } else if (return.predicted.response == TRUE & prediction.approach == "direct.interactions") {
        return(list("recommended.treatment"=recommended.treatment,
                      "predicted.outcome"=predicted.outcome))
    } else if (return.predicted.response == TRUE & prediction.approach == "split.regression") {
        return(list("predicted.outcome.under.control"=predicted.outcome.under.control,
                    "predicted.outcome.under.treatment"=predicted.outcome.under.treatment,
                     "recommended.treatment"=recommended.treatment))
    } else if (return.predicted.response == TRUE & prediction.approach == "OWL.framework") {
        return(list("predicted.treatment.prob"=OWL.framework.predicted.response.prob,
                      "recommended.treatment"=recommended.treatment))
    }
}
