# Lung cancer genomic data

Publicly available lung cancer genomic data from the Chemores Cohort
Study, containing the expression levels of mRNA, miRNA, artificial noise
variables as well as clinical variables and response.

## Format

- `genlung$survival`: A data.frame containing \\n=123\\ complete
  observations. The first column is disease-free survival time and the
  second column is censoring status.

  `genlung$covariate`: A data.frame containing \\p=2000\\ covariates.

## Details

Tissue samples were analysed from a cohort of 123 patients, who
underwent complete surgical resection at the Institut Mutualiste
Montsouris (Paris, France) between 30 January 2002 and 26 June 2006. The
studied outcome was the "Disease-Free Survival Time". Patients were
followed until the first relapse occurred or administrative censoring.
In this genomic dataset, the expression levels of Agilent miRNA probes
(\\p=939\\) were included from the \\n=123\\ cohort samples. The miRNA
data contains normalized expression levels. See below the paper by Lazar
et al. (2013) and Array Express data repository for the complete
description of the samples, tissue preparation, Agilent array
technology, and data normalization. In addition to the genomic data,
five clinical variables, also evaluated on the cohort samples, are
included as continuous variable ('Age') and nominal variables
('Type','KRAS.status','EGFR.status','P53.status'). See Lazar et al.
(2013) for more details. Moreover, we add 1056 standard gaussian
variables which are independent with the censored response as noise
covariates. This dataset represents a situation where the number of
covariates dominates the number of complete observations or \\p \>\> n\\
case.

## References

Lazar V. et al. (2013). Integrated molecular portrait of non-small cell
lung cancers. BMC Medical Genomics 6:53-65.
