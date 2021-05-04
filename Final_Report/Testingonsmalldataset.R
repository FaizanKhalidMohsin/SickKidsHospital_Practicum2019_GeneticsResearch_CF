#V2 R script for running phewas. 

require(data.table)
require(PheWAS)
require(dplyr)
require(tableone)

gdata0 = fread("ukb_chr1_rs4077469.raw", stringsAsFactors=F, header=T, na.strings=c(""," ","NA"))
gdata0[1:10,]
gdata1 = rename(gdata0, id = IID) %>% select( id, rs4077469_T)
gdata1[1:10,]
covariates = rename(gdata0, id = IID) %>% select( id, SEX)
covariates[1:10,]

gdata2 = rename(gdata0, id = IID) %>% select( id, SEX, rs4077469_T)

phedata = fread("icd10_data_with_phecodes2.txt", stringsAsFactors=F, header=T, na.strings=c(""," ","NA"))
head(phedata)
str(phedata)
summary(phedata)
dim(phedata)

phedata0 = phedata %>% rename(phecode = code) %>% select(id, phecode, count) 
head(phedata0)
str(phedata0)
dim(phedata0)

######### Testing small sample first id's = 100 #############################
first_100 = unique(phedata0$id)[1:100]
test_genotypes = filter(gdata1, id %in% first_100)
genotypes = test_genotypes
dim(genotypes)
head(genotypes)
############################################################

length(gdata1$rs4077469_T)
sum(is.na(gdata1$rs4077469_T))
table(gdata1$rs4077469_T)

# phenotypes
phenotypes = createPhewasTable(phedata0, min.code.count = 1, add.exclusions = F, translate = F)
# Recommended to turn off exlucusions: add.exclusions = F
phenotypes[1:10, 1:10]
dim(phenotypes)
str(phenotypes[1:10, 1:10])
rs4077469_T_data = left_join(phenotypes, gdata2, by="id")
dim(rs4077469_T_data)
str(rs4077469_T_data[,1:10])
#rs4077469_T_data$"218.1"


rs4077469_T_data_table = select(rs4077469_T_data, id, '218.1', '218', SEX, rs4077469_T)
rs4077469_T_data_table[1:10,]

# Distribution of cases and controls by allele type
table(rs4077469_T_data_table$`218.1`, 
      rs4077469_T_data_table$rs4077469_T)

table(rs4077469_T_data_table$`218`, 
      rs4077469_T_data_table$rs4077469_T)

table(rs4077469_T_data_table$SEX, 
      rs4077469_T_data_table$rs4077469_T)


# Distribution of cases and controls by gender and allele type
table(filter(rs4077469_T_data_table, SEX == 1)$`218.1`, 
      filter(rs4077469_T_data_table, SEX == 1)$rs4077469_T)

table(filter(rs4077469_T_data_table, SEX == 2)$`218.1`, 
      filter(rs4077469_T_data_table, SEX == 2)$rs4077469_T)

variables_names = names(rs4077469_T_data_table)

# table1 = CreateTableOne(vars = variables_names
#                , factorVars = variables_names
#                , data = rs4077469_T_data_table
#                , strata = "SEX")
# table_1 = print(table1, showAllLevels = T)
# write.csv(table_1, file = "Table1_rs4077469_T_gender.csv")

# rs4077469_T_data1 = left_join(gdata2, phenotypes, by="id")
# rs4077469_T_data_table1 = select(rs4077469_T_data1, id, '218.1', '218', SEX, rs4077469_T)

#Run the PheWAS
results=phewas(phenotypes, genotypes
               , covariates = covariates 
               , cores=1
               , significance.threshold=c("bonferroni"))

#Plot the results
pdf("Testphewasplot_rs4077469_gender.pdf")
phewasManhattan(results, annotate.angle=0,
                title="Test PheWAS Manhattan Plot for rs4077469 with Gender")
dev.off()

#Add PheWAS descriptions
results_d=addPhecodeInfo(results)
#List the significant results
results_d[results_d$bonferroni&!is.na(results_d$p),]
#List the top 10 results
results_d[order(results_d$p)[1:10],]

write.table(results_d, "phewasresults_rs4077469_gender.tsv"
            ,quote = F, row.names = F
            , col.names = T, sep = "\t")

