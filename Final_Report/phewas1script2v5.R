#V2 R script for running phewas. 

require(data.table)
require(PheWAS)
require(dplyr)
require(tableone)


gdata0 = fread("ukb_chr1_rs4077469.raw", stringsAsFactors=F, header=T, na.strings=c(""," ","NA"))
gdata0[1:10,]
gdata1 = rename(gdata0, id = IID) %>% select( id, rs4077469_T)
gdata1[1:10,]
genotypes = gdata1
genotypes[1:10,]
covariate_sex = rename(gdata0, id = IID) %>% select( id, SEX)


gdata2 = rename(gdata0, id = IID) %>% select( id, SEX, rs4077469_T)

covariates_age_gender = fread("pheno_age_sex.tab",  stringsAsFactors=F, header=T, na.strings=c(""," ","NA"))
colnames(covariates_age_gender) = c("id", "age", "SEX")
covariates = mutate(covariates_age_gender, age2 = age**2)
print(covariates[1:10,])


phedata = fread("icd10_data_with_phecodes2.txt", stringsAsFactors=F, header=T, na.strings=c(""," ","NA"))
head(phedata)
str(phedata)
summary(phedata)
dim(phedata)

phedata0 = phedata %>% rename(phecode = code) %>% select(id, phecode, count) 
head(phedata0)
str(phedata0)
print(dim(phedata0))

phenotypes=createPhewasTable(phedata0, min.code.count = 1, add.exclusions = F, translate = F)
# Recommended to turn off exlucusions: add.exclusions = F 
summary(phedata0)
phenotypes[1:10, 1:10]
print(dim(phenotypes))
str(phenotypes[1:10, 1:10])
# rs4077469_T_data0 = left_join(phenotypes, gdata1, by="id")
# rs4077469_T_data = left_join(rs4077469_T_data0, covariates, by="id")

rs4077469_T_data0 = left_join(gdata1, phenotypes,  by="id")
rs4077469_T_data = left_join(rs4077469_T_data0, covariates, by="id")

dim(rs4077469_T_data)
str(rs4077469_T_data[,1:10])
head(rs4077469_T_data$"218.1")
# write.table(rs4077469_T_data, "data_rs4077469_T_phecodes_covariates.tsv"
#             ,quote = F, row.names = F
#             , col.names = T, sep = "\t")

rs4077469_T_data_table1 = select(rs4077469_T_data, '218.1', '218', age, age2, SEX, rs4077469_T)
rs4077469_T_data_table1[1:10,]

# Distribution of alleles
print(length(gdata2$rs4077469_T))
print(sum(is.na(gdata2$rs4077469_T)))
print(table(gdata2$rs4077469_T))


print(length(rs4077469_T_data_table1$rs4077469_T))
print(sum(is.na(rs4077469_T_data_table1$rs4077469_T)))
print(table(rs4077469_T_data_table1$rs4077469_T))

# # Distribution of cases and controls by allele type
# table(rs4077469_T_data_table$`218.1`, 
#       rs4077469_T_data_table$rs4077469_T)
# 
# table(rs4077469_T_data_table$`218`, 
#       rs4077469_T_data_table$rs4077469_T)
# 
# table(rs4077469_T_data_table$SEX, 
#       rs4077469_T_data_table$rs4077469_T)
# 
# 
# # Distribution of cases and controls by gender and allele type
# table(filter(rs4077469_T_data_table, SEX == 1)$`218.1`, 
#       filter(rs4077469_T_data_table, SEX == 1)$rs4077469_T)
# 
# table(filter(rs4077469_T_data_table, SEX == 2)$`218.1`, 
#       filter(rs4077469_T_data_table, SEX == 2)$rs4077469_T)

# Create and save table 1

variables_names = names(rs4077469_T_data_table1)
print(variables_names)

factor_names = c("218.1",
                 "218",
                 "SEX",
                 "rs4077469_T")
print(factor_names)


table1 = CreateTableOne(vars = variables_names
                        , factorVars = factor_names
                        , data = rs4077469_T_data_table
                        #, strata = "SEX"
                        )
table_1 = print(table1, showAllLevels = T)
#write.csv(table_1, file = "Table1_rs4077469_T.csv")

table1 = CreateTableOne(vars = variables_names
               , factorVars = factor_names
               , data = rs4077469_T_data_table
               , strata = "SEX")
table_1 = print(table1, showAllLevels = T)
#write.csv(table_1, file = "Table1_rs4077469_T_gender.csv")


table1 = CreateTableOne(vars = variables_names
                        , factorVars = factor_names
                        , data = rs4077469_T_data_table
                        , strata = "rs4077469_T")
table_1 = print(table1, showAllLevels = T)
write.csv(table_1, file = "Table1_rs4077469_T_bySNP.csv")

table1 = CreateTableOne(vars = variables_names
                        , factorVars = factor_names
                        , data = rs4077469_T_data_table
                        , strata = c("218.1", "SEX"))
table_1 = print(table1, showAllLevels = T)
#write.csv(table_1, file = "Table1_rs4077469_T_gender_Uterine_leiomyoma.csv")



table1 = CreateTableOne(vars = variables_names
                        , factorVars = factor_names
                        , data = rs4077469_T_data_table
                        , strata = c("218", "SEX"))
table_1 = print(table1, showAllLevels = T)
#write.csv(table_1, file = "Table1_rs4077469_T_gender_BenignNeoplasm_of_Uterus.csv")



#Run the PheWAS
results=phewas(phenotypes, genotypes
               , covariates = covariates 
               , cores=1
               , significance.threshold=c("bonferroni")
               )

#Plot the results
pdf("phewasplot_rs4077469_covariates.pdf")
phewasManhattan(results, annotate.angle=0,
                title="Manhattan Plot for rs4077469 with Gender, Age & Age-sqrd")
dev.off()

#Add PheWAS descriptions
results_d=addPhecodeInfo(results)
#List the significant results
results_d[results_d$bonferroni&!is.na(results_d$p),]
#List the top 10 results
results_d[order(results_d$p)[1:10],]

write.table(results_d, "phewasresults_rs4077469_covariates.tsv"
            ,quote = F, row.names = F
            , col.names = T, sep = "\t")

