GetComparisonRules <- function(scenario,
                                     B.optimal,
                                     test.df,
                                     study.design,
                                     name.outcome, type.outcome, desirable.outcome,
                                     name.treatment, names.influencing.treatment, names.influencing.rule,
                                     propensity.method,
                                     clinical.threshold=0, 
                                     bootstrap.CI=FALSE) {
    # optimal rule
    evaluate.optimal.rule <- EvaluateRule(evaluation.data=test.df,
                                      BuildRule.object=NULL,
                                      B=B.optimal,
                                      study.design=study.design,
                                      name.outcome=name.outcome,
                                      type.outcome=type.outcome,
                                      desirable.outcome=desirable.outcome,
                                      clinical.threshold=0, 
                                      name.treatment=name.treatment,
                                      names.influencing.treatment=names.influencing.treatment,
                                      names.influencing.rule=names.influencing.rule,
                                      propensity.method=propensity.method,
                                      bootstrap.CI=bootstrap.CI)
    # treating all
    B.treat.all <- rep(1, nrow(test.df))
    evaluate.treat.all <- EvaluateRule(evaluation.data=test.df,
                                        BuildRule.object=NULL,
                                        B=B.treat.all,
                                        study.design=study.design,
                                        name.outcome=name.outcome,
                                        type.outcome=type.outcome,
                                        desirable.outcome=desirable.outcome,
                                        clinical.threshold=0, 
                                        name.treatment=name.treatment,
                                        names.influencing.treatment=names.influencing.treatment,
                                        names.influencing.rule=names.influencing.rule,
                                        propensity.method=propensity.method,
                                        bootstrap.CI=bootstrap.CI)
    # treating none
    B.treat.none <- rep(0, nrow(test.df))
    evaluate.treat.none <- EvaluateRule(evaluation.data=test.df,
                                        BuildRule.object=NULL,
                                        B=B.treat.none,
                                        study.design=study.design,
                                        name.outcome=name.outcome,
                                        type.outcome=type.outcome,
                                        desirable.outcome=desirable.outcome,
                                        clinical.threshold=0, 
                                        name.treatment=name.treatment,
                                        names.influencing.treatment=names.influencing.treatment,
                                        names.influencing.rule=names.influencing.rule,
                                        propensity.method=propensity.method,
                                        bootstrap.CI=bootstrap.CI)
    return(list("evaluate.optimal.rule"=evaluate.optimal.rule,
                  "evaluate.treat.all"=evaluate.treat.all,
                   "evaluate.treat.none"=evaluate.treat.none))
}
