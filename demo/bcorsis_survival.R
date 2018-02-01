library(Ball)

data("genlung")
predictor <- genlung[["covariate"]]
surv_status <- genlung[["survival"]]
result <- bcorsis(x = predictor,
                  y = surv_status,
                  d = "small", method = "survival")
top_gene <- colnames(predictor)[result[["ix"]]]
head(top_gene, n = 1)
