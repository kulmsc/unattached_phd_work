
author=$1

Rscript simple_report.R $author

Rscript extreme_compare.R $author
Rscript extreme_timeline.R $author

Rscript enet_residual_pheno.R $author
Rscript enet_extreme.R $author
