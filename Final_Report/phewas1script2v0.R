#V2 R script for running phewas. 

require(data.table)
require(PheWAS)
require(dplyr)


gdata1 = fread("ukb_chr1_rs4077469.raw", stringsAsFactors=F, header=T, na.strings=c(""," ","NA"))
gdata1[1:10,]
gdata1 = rename(gdata1, id = IID) %>% select( id, rs4077469_T)
gdata1[1:10,]
genotypes = gdata1

phedata = fread("icd10_data_with_phecodes2.txt", stringsAsFactors=F, header=T, na.strings=c(""," ","NA"))
head(phedata)
str(phedata)
summary(phedata)

phedata0 = phedata %>% rename(phecode = code) %>% select(id, phecode, count) 
head(phedata0)
str(phedata0)
summary(phedata0)
dim(phedata0)

phenotypes=createPhewasTable(phedata0, min.code.count = 1, add.exclusions = F, translate = F)
# Recommended to turn off exlucusions: add.exclusions = F
phenotypes[1:10, 1:10]

#Run the PheWAS
results=phewas(phenotypes,genotypes,cores=1,
               significance.threshold=c("bonferroni"))

#Plot the results
pdf("phewasplotcount1.pdf")
phewasManhattan(results, annotate.angle=0,
                title="My Example PheWAS Manhattan Plot for rs4077469")
dev.off()

#Add PheWAS descriptions
results_d=addPhecodeInfo(results)
#List the significant results
results_d[results_d$bonferroni&!is.na(results_d$p),]
#List the top 10 results
results_d[order(results_d$p)[1:10],]

write.table(results_d, "phewas1count1run1.tsv"
            ,quote = F, row.names = F
            , col.names = T, sep = "\t")

