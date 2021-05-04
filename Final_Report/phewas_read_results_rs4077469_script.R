

citation("survival")
citation("dplyr")
citation()
citation("tableone")


require(data.table)
require(PheWAS)
require(dplyr)
require(tableone)

# Checking PheWAS Results for rs4077469

results_rs4077469 = fread("phewasresults_rs4077469_covariates_final_data.tsv", stringsAsFactors=F, header=T, na.strings=c(""," ","NA"))

results_d = results_rs4077469 

#List the significant results
results_d[results_d$bonferroni&!is.na(results_d$p),]
#List the top 10 results
results_d[order(results_d$p)[1:10],]
#List the top 3 results
results_d[order(results_d$p)[1:3],]