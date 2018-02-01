library(Ball)

load("genlung")
result <- bcorsis(x = genlung[["covariate"]],
                  y = genlung[["survival"]],
                  d = "small", method = "survival")
top_gene <- colnames(genlung[["covariate"]])[result[["ix"]]]
head(top_gene, n = 1)
