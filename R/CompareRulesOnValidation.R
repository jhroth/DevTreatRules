#' Build treatment rules on a development dataset and evaluate performance on an independent validation dataset
#'
#' In many practical settings, \code{BuildRule()} has limited utility because it requires the specification of a single value in its  \code{prediction.approach} argument (even if there is no prior knowledge about which of the split-regression, OWL framework, and direct-interactions approaches will perform best) and a single value for the `propensity.score` and `rule.method` arguments (even if there is no prior knowledge about whether standard or penalized GLM will perform best). \code{CompareRulesOnValidation()} supports model selection in these settings by  essentially looping over calls to \code{BuildRule()} for different combinations of split-regression/OWL framework/direct-interactions and standard/lasso/ridge regression to simultaneously build the rules on a development dataset and evaluate them on an independent validation dataset.

#'
#' @param development.data A data frame representing the *development* dataset used to build treatment rules.
#' @param validation.data A data frame representing the independent *validation* dataset used to estimate the performance of treatment rules built on the development dataset.
#' @param vec.approaches A character vector (or element) indicating the values of the \code{prediction.approach} to be used for building the rule with \code{BuildRule()}. Default is \code{c(`split.regression', `OWL.framework', `direct.interactions')}.
#' @param vec.rule.methods A character vector (or element) indicating the values of the \code{rule.method} to be used for building the rule with \code{BuildRule()}. Default is \code{c(`glm.regression', `lasso', `ridge')}.
#' @param vec.propensity.methods A character vector (or element) indicating the values of \code{propensity.method} to be used for building the rule with \code{Build.Rule()}. Default is `logistic.regression' to allow for estimation of bootstrap-based CIs.
#' @param study.design.development Either `observational', `RCT', or `naive', representing the study design on the development dataset. For the \code{observational} design, the function will use inverse-probability-of-treatment observation weights (IPW) based on estimated propensity scores with predictors \code{names.influencing.treatment}; for the \code{RCT} design, the function will use IPW based on propensity scores equal to the observed sample proportions; for the \code{naive} design, all observation weights will be uniformly equal to 1.
#' @param study.design.validation Either `observational', `RCT', or `naive',representing the study design on the development dataset. Default is the value of \code{study.design.development}.
#' @param name.outcome.development A character indicating the name of the outcome variable in \code{development.data}.
#' @param type.outcome.development Either `binary' or `continuous', the form of \code{name.outcome.development}.
#' @param type.outcome.validation Either `binary' or `continuous', the form of \code{name.outcome.validation}. Default is the value of \code{type.outcome.development}.
#' @param name.outcome.validation A character indicating the name of the outcome variable in \code{validation.data}. Default is the value of \code{name.outcome.development}.
#' @param name.treatment.development A character indicating the name of the treatment variable in \code{development.data}.
#' @param name.treatment.validation A character indicating the name of the treatment variable in \code{validation.data}. Default is the value of \code{name.treatment.development}
#' @param names.influencing.treatment.development A character vector (or element) indicating the names of the variables in \code{development.data} that are expected to influence treatment assignment in the current dataset. Required for \code{study.design.development=}`observational'. 
#' @param names.influencing.treatment.validation A character vector (or element) indicating the names of the variables in \code{validation.data} that are expected to influence treatment assignment in \code{validation.data}. Required for Required for \code{study.design.validation=}`observational'. Default is the value of \code{names.influencing.treatment.development}.
#' @param names.influencing.rule.development A character vector (or element) indicating the names of the variables in \code{development.data} that may influence response to treatment and are expected to be observed in future clinical settings.
#' @param names.influencing.rule.validation A character vector (or element) indicating the names of the variables in \code{validation.data} that may influence response to treatment and are expected to be observed in future clinical settings. Default is the value of \code{names.influencing.rule.development}
#' @param desirable.outcome.development A logical equal to \code{TRUE} if higher values of the outcome on \code{development,data} are considered desirable (e.g. for a binary outcome, a 1 is more desirable than a 0). The \code{OWL.framework} and \code{OWL} prediction approaches require a desirable outcome.
#' @param desirable.outcome.validation A logical equal to \code{TRUE} if higher values of the outcome on \code{validation,data}  are considered desirable (e.g. for a binary outcome, a 1 is more desirable than a 0). The \code{OWL.framework} and \code{OWL} prediction approaches require a desirable outcome. Default is the value of \code{desirable.outcome.development}
#' @param clinical.threshold.validation A numeric equal to a positive number above which the predicted outcome under treatment must be superior to the predicted outcome under control for treatment to be recommended. Only used when \code{BuildRuleObject} was specified and derived from the split-regression or direct-interactions approach. Default is 0.
#' @param propensity.method.validation One of `logistic.regression', `lasso', or `ridge'. This is the underlying regression model used to estimate propensity scores (for \code{study.design=}`observational' on \code{validation.data}. If \code{bootstrap.CI=TRUE}, then \code{propensity.method} must be `logistic.regression'. Default is `logistic.regression' to allow for estimation of bootstrap-based CIs.
#' @param additional.weights.development A numeric vector of observation weights that will be multiplied by IPW weights in the rule development stage, with length equal to the number of rows in \code{development.data}. This can be used, for example, to account for a non-representative sampling design or an IPW adjustment for missingness. The default is a vector of 1s.
#' @param additional.weights.validation A numeric vector of observation weights that will be multiplied by IPW weights in the rule evaluation stage, with length equal to the number of rows in \code{validation.data}. This can be used, for example, to account for a non-representative sampling design or an IPW adjustment for missingness. The default is a vector of 1s.
#' @param truncate.propensity.score A logical variable dictating whether estimated propensity scores less than \code{truncate.propensity.score.threshold} away from 0 or 1 should be truncated to be \code{truncate.propensity.score.threshold} away from 0 or 1.
#' @param truncate.propensity.score.threshold A numeric value between 0 and 0.25.
#' @param type.observation.weights Default is NULL, but other choices are `IPW.L', `IPW.L.and.X', and `IPW.ratio', where L indicates the \code{names.influencing.treatment} variables, X indicates the \code{names.influencing.rule} variables. The default behavior is to use the `IPW.ratio' observation weights (propensity score based on X divided by propensity score based on L and X) for \code{prediction.approach=}`split.regression' and to use `IPW.L' observation weights (inverse of propensity score based on L) for the `direct.interactions', `OWL', and `OWL.framework' prediction approaches.
#' @param propensity.k.cv.folds An integer specifying how many folds to use for K-fold cross-validation that chooses the tuning parameter when \code{propensity.method} is `lasso' or `ridge'. Default is 10.
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
#' @param bootstrap.CI Logical indicating whether the ATE/ABR estimates on the validation set should be accompanied by 95\% confidence intervals based on the bootstrap. Default is \code{FALSE}.
#' @param bootstrap.CI.replications An integer specifying how many bootstrap replications should underlie the computed CIs. Default is 1000.
#' @return A list with components:
#' \itemize{
#'   \item \code{list.summaries}: A list with number of elements equal to the length of \code{vec.approaches}. Each element is a matrix that, for a given prediction approach, shows estimated rule performance with 5 columns if \code{bootstrap.CI=FALSE}  (number of test-positives, number of test-negatives, ATE in test-positives, ATE in test-negatives, ABR) for the different combinations of \code{vec.rule.methods} or 9 columns if \code{bootstrap.CI=TRUE} (those same 5 summaries plus the bounds for 95\% CIs for ATE in test-positives and ATE in test-negatives) and, in the rows, the \code{vec.propensity.methods} in addition to the two naive rules (treating all observations and treating no observations).
#'   \item \code{list.rules}: A list with number of elements equal to the length of \code{vec.approaches}. Each element is another list that, for a given prediction approach, stores the object returned by \code{BuildRule()} for the different combinations of \code{vec.rule.methods} and \code{vec.propensity.methods} in the rows.
#' }
#' @examples
#' set.seed(123)
#' example.split <- SplitData(data=obsStudyGeneExpressions,
#'                                     n.sets=3, split.proportions=c(0.5, 0.25, 0.25))
#' development.data <- example.split[example.split$partition == "development", ]
#' validation.data <- example.split[example.split$partition == "validation", ]
#' model.selection <- CompareRulesOnValidation(development.data=development.data,
#'                validation.data=validation.data,
#'                study.design.development="observational",
#'                vec.approaches=c("split.regression", "OWL.framework", "direct.interactions"),
#'                vec.rule.methods=c("glm.regression", "lasso"),
#'                vec.propensity.methods="logistic.regression",
#'                name.outcome.development="no_relapse",
#'                type.outcome.development="binary",
#'                name.treatment.development="intervention",
#'                names.influencing.treatment.development=c("prognosis", "clinic", "age"),
#'                names.influencing.rule.development=c("age", paste0("gene_", 1:10)),
#'                desirable.outcome.development=TRUE)
#' model.selection$list.summaries$split.regression
#' @export
CompareRulesOnValidation <- function(development.data,
                                                 validation.data,
                                                 vec.approaches=c("split.regression", "OWL.framework", "direct.interactions"),
                                                 vec.rule.methods=c("glm.regression", "lasso", "ridge"),
                                                 vec.propensity.methods="logistic.regression",
                                                 study.design.development, 
                                                 name.outcome.development,
                                                 type.outcome.development, 
                                                 name.treatment.development,
                                                 names.influencing.treatment.development,
                                                 names.influencing.rule.development,
                                                 desirable.outcome.development,
                                                 additional.weights.development=rep(1, nrow(development.data)),
                                                 study.design.validation=study.design.development, 
                                                 name.outcome.validation=name.outcome.development,
                                                 type.outcome.validation=type.outcome.development, 
                                                 name.treatment.validation=name.treatment.development,
                                                 names.influencing.treatment.validation=names.influencing.treatment.development,
                                                 names.influencing.rule.validation=names.influencing.rule.development,
                                                 desirable.outcome.validation=desirable.outcome.development,
                                                 clinical.threshold.validation=0,
                                                 propensity.method.validation="logistic.regression",
                                                 additional.weights.validation=rep(1, nrow(validation.data)),
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
                                                 direct.interactions.exclude.A.from.penalty=TRUE,
                                                 bootstrap.CI=FALSE,
                                                 bootstrap.CI.replications=100) {
    stopifnot(type.outcome.development == type.outcome.validation)
    stopifnot(is.logical(bootstrap.CI))
    list.summaries <- vector("list", length(vec.approaches))
    #names(list.summaries) <- vec.approaches
    list.rules <- vector("list", length(vec.approaches))
    for (a in 1:length(vec.approaches)) {
        one.approach <- vec.approaches[a]
        # for storing summaries in matrices
        if (bootstrap.CI == FALSE) {
            mat.summaries <- matrix(NA, nrow=length(vec.rule.methods)*length(vec.propensity.methods) + 2, ncol=5)
            colnames(mat.summaries) <- c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")
            rownames(mat.summaries) <- c(rep(NA, length(vec.rule.methods)*length(vec.propensity.methods)), "treat.all", "treat.none")
        } else {
            mat.summaries <- matrix(NA, nrow=length(vec.rule.methods)*length(vec.propensity.methods) + 2, ncol=9)
            colnames(mat.summaries) <- c("Positives", "Negatives",
                                                       "ATE in Positives", "ATE in Positives - CI lower", "ATE in Positives - CI upper",
                                                       "ATE in Negatives", "ATE in Negatives - CI lower", "ATE in Negatives - CI upper", 
                                                       "ABR")
            rownames(mat.summaries) <- c(rep(NA, length(vec.rule.methods)*length(vec.propensity.methods)), "treat.all", "treat.none")
        }
        list.rules.one.approach <- vector("list", length(vec.rule.methods)*length(vec.propensity.methods))
        row.number <- 1
        if (one.approach != "OWL") {
            for (r in 1:length(vec.rule.methods)) {
                one.rule.method <- vec.rule.methods[r]
                for (p in 1:length(vec.propensity.methods)) {
                    one.propensity.method <- vec.propensity.methods[p]
                    # build rule on development dataset for each combination of approach/rule.method/propensity/method 
                    build.one <- BuildRule(development.data=development.data,
                                           study.design=study.design.development,
                                           prediction.approach=one.approach,
                                           name.outcome=name.outcome.development,
                                           type.outcome=type.outcome.development,
                                           name.treatment=name.treatment.development,
                                           names.influencing.treatment=names.influencing.treatment.development,
                                           names.influencing.rule=names.influencing.rule.development,
                                           desirable.outcome=desirable.outcome.development,
                                           propensity.method=one.propensity.method,
                                           additional.weights=additional.weights.development,
                                           rule.method=one.rule.method, 
                                           truncate.propensity.score=truncate.propensity.score,
                                           truncate.propensity.score.threshold=truncate.propensity.score.threshold,
                                           type.observation.weights=type.observation.weights,
                                           propensity.k.cv.folds=propensity.k.cv.folds,
                                           rule.k.cv.folds=rule.k.cv.folds,
                                           lambda.choice=lambda.choice,
                                           OWL.lambda.seq=OWL.lambda.seq,
                                           OWL.kernel=OWL.kernel,
                                           OWL.kparam.seq=OWL.kparam.seq,
                                           OWL.cvFolds=OWL.cvFolds,
                                           OWL.verbose=OWL.verbose,
                                           OWL.framework.shift.by.min=OWL.framework.shift.by.min,
                                           direct.interactions.center.continuous.Y=direct.interactions.center.continuous.Y,
                                           direct.interactions.exclude.A.from.penalty=direct.interactions.exclude.A.from.penalty)
                    
                    # evaluate the rule returned by BuildRule() on the validation data
                    evaluate.one <- EvaluateRule(evaluation.data=validation.data,
                                                 BuildRule.object=build.one,
                                                 study.design=study.design.validation,
                                                 name.outcome=name.outcome.validation,
                                                 type.outcome=type.outcome.validation,
                                                 desirable.outcome=desirable.outcome.validation,
                                                 clinical.threshold=clinical.threshold.validation,
                                                 name.treatment=name.treatment.validation,
                                                 names.influencing.treatment=names.influencing.treatment.validation,
                                                 names.influencing.rule=names.influencing.rule.validation,
                                                 propensity.method=propensity.method.validation,
                                                 show.treat.all=FALSE,
                                                 show.treat.none=FALSE,
                                                 additional.weights=additional.weights.validation,
                                                 bootstrap.CI=bootstrap.CI,
                                                 bootstrap.CI.replications=bootstrap.CI.replications)
                    # store results
                    one.rowname <- paste0("propensity_", one.propensity.method, "_rule_", one.rule.method)
                    rownames(mat.summaries)[row.number] <- one.rowname
                    if (bootstrap.CI == FALSE) {
                        mat.summaries[one.rowname, c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")] <-
                            as.numeric(evaluate.one$summaries["estimated.rule", c("n.positives", "n.negatives", "ATE.positives", "ATE.negatives", "ABR")])
                    } else {
                        mat.summaries[one.rowname, c("Positives", "Negatives",
                                                                    "ATE in Positives", "ATE in Positives - CI lower", "ATE in Positives - CI upper",
                                                                    "ATE in Negatives", "ATE in Negatives - CI lower", "ATE in Negatives - CI upper", 
                                                                    "ABR")] <-
                            as.numeric(evaluate.one$summaries["estimated.rule", c("n.positives", "n.negatives",
                                                                                                       "ATE.positives", "ATE.positives.CI.LB", "ATE.positives.CI.UB",
                                                                                                       "ATE.negatives", "ATE.negatives.CI.LB", "ATE.negatives.CI.UB",
                                                                                                       "ABR")])
                    }
                    list.rules.one.approach[[row.number]] <- build.one
                    names(list.rules.one.approach)[[row.number]] <- one.rowname
                    row.number <- row.number + 1
                }
            }
            evaluate.one.with.naive <- EvaluateRule(evaluation.data=validation.data,
                                                    BuildRule.object=build.one,
                                                    study.design=study.design.validation,
                                                    name.outcome=name.outcome.validation,
                                                    type.outcome=type.outcome.validation,
                                                    desirable.outcome=desirable.outcome.validation,
                                                    clinical.threshold=clinical.threshold.validation,
                                                    name.treatment=name.treatment.validation,
                                                    names.influencing.treatment=names.influencing.treatment.validation,
                                                    names.influencing.rule=names.influencing.rule.validation,
                                                    propensity.method=propensity.method.validation,
                                                    show.treat.all=TRUE,
                                                    show.treat.none=TRUE,
                                                    additional.weights=additional.weights.validation,
                                                    bootstrap.CI=bootstrap.CI,
                                                    bootstrap.CI.replications=bootstrap.CI.replications)
            mat.summaries["treat.all", c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")] <-
                as.numeric(evaluate.one.with.naive$summaries["treat.all", c("n.positives", "n.negatives", "ATE.positives", "ATE.negatives", "ABR")])
            mat.summaries["treat.none", c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")] <-
                as.numeric(evaluate.one.with.naive$summaries["treat.none", c("n.positives", "n.negatives", "ATE.positives", "ATE.negatives", "ABR")])
        } else {
            for (p in 1:length(vec.propensity.methods)) {
                one.propensity.method <- vec.propensity.methods[p]
                # build rule on development dataset for each combination of approach/rule.method/propensity/method 
                build.one <- BuildRule(development.data=development.data,
                                       study.design=study.design.development,
                                       prediction.approach=one.approach,
                                       name.outcome=name.outcome.development,
                                       type.outcome=type.outcome.development,
                                       name.treatment=name.treatment.development,
                                       names.influencing.treatment=names.influencing.treatment.development,
                                       names.influencing.rule=names.influencing.rule.development,
                                       desirable.outcome=desirable.outcome.development,
                                       propensity.method=one.propensity.method,
                                       additional.weights=additional.weights.development,
                                       rule.method=one.rule.method, 
                                       truncate.propensity.score=truncate.propensity.score,
                                       truncate.propensity.score.threshold=truncate.propensity.score.threshold,
                                       type.observation.weights=type.observation.weights,
                                       propensity.k.cv.folds=propensity.k.cv.folds,
                                       rule.k.cv.folds=rule.k.cv.folds,
                                       lambda.choice=lambda.choice,
                                       OWL.lambda.seq=OWL.lambda.seq,
                                       OWL.kernel=OWL.kernel,
                                       OWL.kparam.seq=OWL.kparam.seq,
                                       OWL.cvFolds=OWL.cvFolds,
                                       OWL.verbose=OWL.verbose,
                                       OWL.framework.shift.by.min=OWL.framework.shift.by.min,
                                       direct.interactions.center.continuous.Y=direct.interactions.center.continuous.Y,
                                       direct.interactions.exclude.A.from.penalty=direct.interactions.exclude.A.from.penalty)
                
                # evaluate the rule returned by BuildRule() on the validation data
                evaluate.one <- EvaluateRule(evaluation.data=validation.data,
                                             BuildRule.object=build.one,
                                             study.design=study.design.validation,
                                             name.outcome=name.outcome.validation,
                                             type.outcome=type.outcome.validation,
                                             desirable.outcome=desirable.outcome.validation,
                                             clinical.threshold=clinical.threshold.validation,
                                             name.treatment=name.treatment.validation,
                                             names.influencing.treatment=names.influencing.treatment.validation,
                                             names.influencing.rule=names.influencing.rule.validation,
                                             propensity.method=propensity.method.validation,
                                             show.treat.all=FALSE,
                                             show.treat.none=FALSE,
                                             additional.weights=additional.weights.validation,
                                             bootstrap.CI=bootstrap.CI,
                                             bootstrap.CI.replications=bootstrap.CI.replications)
                # store results
                one.rowname <- paste0("propensity_", one.propensity.method)
                rownames(mat.summaries)[row.number] <- one.rowname
                mat.summaries[one.rowname, c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")] <-
                    as.numeric(evaluate.one$summaries["estimated.rule", c("n.positives", "n.negatives", "ATE.positives", "ATE.negatives", "ABR")])
                list.rules.one.approach[[row.number]] <- build.one
                names(list.rules.one.approach)[[row.number]] <- one.rowname
                row.number <- row.number + 1
            }
            mat.summaries["treat.all", c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")] <-
                as.numeric(evaluate.one.with.naive$summaries["treat.all", c("n.positives", "n.negatives", "ATE.positives", "ATE.negatives", "ABR")])
            mat.summaries["treat.none", c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")] <-
                as.numeric(evaluate.one.with.naive$summaries["treat.none", c("n.positives", "n.negatives", "ATE.positives", "ATE.negatives", "ABR")])
        }
        list.summaries[[a]] <- mat.summaries
        names(list.summaries)[[a]] <- one.approach
        list.rules[[a]] <- list.rules.one.approach
        names(list.rules)[[a]] <- one.approach
    }
    return(list("list.summaries"=list.summaries,
                  "list.rules"=list.rules))
}




