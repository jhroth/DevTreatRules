CompareRulesOnValidation <- function(development.data,
                                                 validation.data,
                                                 study.design.development, 
                                                 vec.approaches,
                                                 vec.rule.methods,
                                                 vec.propensity.methods,
                                                 name.outcome.development,
                                                 type.outcome.development, 
                                                 name.treatment.development,
                                                 names.influencing.treatment.development,
                                                 names.influencing.rule.development,
                                                 desirable.outcome.development,
                                                 additional.weights.development=rep(1, nrow(development.data)),
                                                 # idk if setting defaults like this will work or if i should do =NULL and then stuff later
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
                                                 # we'll see
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
                                                 direct.interactions.exclude.A.from.penalty=TRUE,
                                                 bootstrap.CI=FALSE,
                                                 bootstrap.CI.replications=100) {
    stopifnot(type.outcome.development == type.outcome.validation)
    list.summaries <- vector("list", length(vec.approaches))
    #names(list.summaries) <- vec.approaches
    list.rules <- vector("list", length(vec.approaches))
    for (a in 1:length(vec.approaches)) {
        one.approach <- vec.approaches[a]
        # for storing summaries in matrices
        mat.summaries <- matrix(NA, nrow=length(vec.rule.methods)*length(vec.propensity.methods), ncol=5)
        colnames(mat.summaries) <- c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")
        rownames(mat.summaries) <- rep(NA, length(vec.rule.methods)*length(vec.propensity.methods))
        list.rules.one.approach <- vector("list", length(vec.rule.methods)*length(vec.propensity.methods))
        row.number <- 1
        for (r in 1:length(vec.rule.methods)) {
            one.rule.method <- vec.rule.methods[r]
            for (p in 1:length(vec.propensity.methods)) {
                one.propensity.method <- vec.propensity.methods[p]
                # build rule on development dataset for each combination of approach/rule.method/propensity/method 
                build.one <- BuildRule(data=development.data,
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
                                       propensity.b.cv.repeats=propensity.b.cv.repeats,
                                       rule.k.cv.folds=rule.k.cv.folds,
                                       rule.b.cv.repeats=rule.b.cv.repeats,
                                       lambda.choice=lambda.choice,
                                       OWL.lambda.seq=OWL.lambda.seq,
                                       OWL.n.lambda.seq=OWL.n.lambda.seq,
                                       OWL.kernel=OWL.kernel,
                                       OWL.kparam.seq=OWL.kparam.seq,
                                       OWL.cvFolds=OWL.cvFolds,
                                       OWL.verbose=OWL.verbose,
                                       OWL.framework.shift.by.min=OWL.framework.shift.by.min,
                                       direct.interactions.center.continuous.Y=direct.interactions.center.continuous.Y,
                                       direct.interactions.exclude.A.from.penalty=direct.interactions.exclude.A.from.penalty)
                
                # evaluate the rule returned by BuildRule() on the validation data
                evaluate.one <- EvaluateRule(data=validation.data,
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
                                             additional.weights=additional.weights.validation,
                                             bootstrap.CI=bootstrap.CI,
                                             bootstrap.CI.replications=bootstrap.CI.replications)
                
                # store results
                one.rowname <- paste0("propensity_", one.propensity.method, "_rule_", one.rule.method)
                rownames(mat.summaries)[row.number] <- one.rowname
                mat.summaries[one.rowname, c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")] <- 
                as.numeric(evaluate.one[c("n.test.positives", "n.test.negatives", "ATE.test.positives", "ATE.test.negatives", "ABR")])
                list.rules.one.approach[[row.number]] <- build.one
                names(list.rules.one.approach)[[row.number]] <- one.rowname
                row.number <- row.number + 1
            }
        }
        list.summaries[[a]] <- mat.summaries
        names(list.summaries)[[a]] <- one.approach
        list.rules[[a]] <- list.rules.one.approach
        names(list.rules)[[a]] <- one.approach
    }
    return(list("list.summaries"=list.summaries,
                  "list.rules"=list.rules))
}




