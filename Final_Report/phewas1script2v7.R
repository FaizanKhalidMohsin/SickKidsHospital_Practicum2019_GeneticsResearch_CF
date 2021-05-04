#V2 R script for running phewas. 

require(data.table)
require(PheWAS)
require(dplyr)
require(tableone)


unrelated = fread("09-kinship-degree2-unrelatedunrelated.txt", col.names = c("fid", "id"), stringsAsFactors=F, header=F, na.strings=c(""," ","NA"))
#unrelated_tb_removed = fread("09-kinship-degree2-unrelatedunrelated_toberemoved.txt", col.names = c("fid", "id") , stringsAsFactors=F, header=F, na.strings=c(""," ","NA"))
#454029+ 36100  
#488377- 36100

print(str(unrelated))
unrelated$id = as.numeric(unrelated$id)
unrelated1 = unrelated  %>% filter(!is.na(id)) %>% filter(id > 0 )
print(dim(unrelated1))
print(str(unrelated1))

gdata0 = fread("ukb_chr1_rs4077469.raw", stringsAsFactors=F, header=T, na.strings=c(""," ","NA"))
gdata0[1:10,]
gdata1 = rename(gdata0, gender = SEX, id = IID) 
gdata1[1:10,]
print(str(gdata1))

# Creating the independent data set for the genotypic data. 
gdata2 = inner_join(unrelated1, gdata1, by ="id")
dim(gdata2)
gdata2[1:10,]
print(str(gdata2))
summary(gdata2)


# Reading in genetic ethnic group data ( caucasian = T/F)
ethnic_data0 = fread("ukb24727_22006_genetic_ethnic_groups.tab", 
                    col.names = c("id", "caucasian"), stringsAsFactors=F, 
                    header=T, na.strings=c(""," ","NA"))
ethnic_data = filter(ethnic_data0, !is.na(caucasian))
dim(ethnic_data)
print(str(ethnic_data))
summary(ethnic_data)


gdata3 = inner_join(gdata2, ethnic_data, by = "id")
print(dim(gdata3))
summary(gdata3)
str(gdata3)
genotypes = gdata3 %>% select( id, rs4077469_T)
genotypes[1:10,]


covariates_age_gender = fread("pheno_age_sex.tab",  stringsAsFactors=F, 
                              col.names = c("id", "age", "SEX"), 
                              header=T, na.strings=c(""," ","NA"))
covariates1 = mutate(covariates_age_gender, age2 = age**2)
covariates = semi_join(covariates1, genotypes, by = "id")
summary(covariates)
print(str(covariates))
print(dim(covariates))
print(covariates[1:10,])


phedata = fread("icd10_data_with_phecodes2.txt", stringsAsFactors=F, header=T, na.strings=c(""," ","NA"))
head(phedata)
print(str(phedata))
summary(phedata)
dim(phedata)

phedata0 = phedata %>% rename(phecode = code) %>% select(id, phecode, count) 
head(phedata0)
str(phedata0)
print(dim(phedata0))
summary(phedata0)

phenotypes=createPhewasTable(phedata0, min.code.count = 1, add.exclusions = F, translate = F)
# Recommended to turn off exlucusions: add.exclusions = F 
phenotypes[1:10, 1:10]
print(dim(phenotypes))
str(phenotypes[1:10, 1:10])

# Do not use the below one. Use "phenotypes" data.  
phenotypes1 = semi_join(phenotypes, genotypes, by = "id")
print(dim(phenotypes1))
print(phenotypes1[1:10, 1:10])

genotypes_1 = semi_join(genotypes, phenotypes,  by = "id")
print(dim(genotypes_1))
print(genotypes_1[1:10,])


# Missing id's
phenotypes_missing_id = anti_join(phenotypes, genotypes, by = "id")
print(dim(phenotypes_missing_id))
print(phenotypes_missing_id[1:10, 1:10])

genotypes_missing_id = anti_join(genotypes, phenotypes, by = "id")
print(dim(genotypes_missing_id))
print(phenotypes_missing_id[1:10, 1:10])

missing_ids_from_gdata1 = anti_join(phenotypes, gdata1, by = "id")
print(dim(missing_ids_from_gdata1))
print(missing_ids_from_gdata1[1:10, 1:10])




# Also, do not use this genotypes1. Use "genotypes".
genotypes1 = semi_join(genotypes, phenotypes1, by = "id")
summary(genotypes1)
str(genotypes1)
dim(genotypes1)


print(dim(genotypes))
print(dim(phenotypes))
print(dim(covariates))

#Run the PheWAS
results=phewas(phenotypes, genotypes
               , covariates = covariates 
               , cores=1
               , significance.threshold=c("bonferroni")
               )

#Plot the results
pdf("phewasplot_rs4077469_covariates_caucasians.pdf")
phewasManhattan(results, annotate.angle=0,
                title="Manhattan Plot for rs4077469 and Caucasians with Covariates")
dev.off()

#Add PheWAS descriptions
results_d=addPhecodeInfo(results)
#List the significant results
results_d[results_d$bonferroni&!is.na(results_d$p),]
#List the top 10 results
results_d[order(results_d$p)[1:10],]

r = results_d[order(results_d$p)[1:10],]
write.csv(r, "Results_Top10_rs4077469_covariates_final_data.csv")

# Save the entire PheWAS Study results
write.table(results_d, "phewasresults_rs4077469_covariates_final_data.tsv"
            ,quote = F, row.names = F
            , col.names = T, sep = "\t")

phecode1 = as.character(results_d$phecode[1])
phecode2 = as.character(results_d$phecode[2])

#gdata3 = rename(gdata3, gender = SEX)
rs4077469_T_data0 = inner_join(gdata3, phenotypes,  by="id") %>% filter(!is.na(rs4077469_T))
rs4077469_T_data = left_join(rs4077469_T_data0, covariates, by="id")
print(dim(rs4077469_T_data))

#str(rs4077469_T_data[,1:10])
#head(rs4077469_T_data$"218.1")
# write.table(rs4077469_T_data, "data_rs4077469_T_phecodes_covariates_final_data.tsv"
#             ,quote = F, row.names = F
#             , col.names = T, sep = "\t")

#rs4077469_T_data_table1 = select(rs4077469_T_data, '218.1', '218', age, age2, gender, SEX, rs4077469_T)
# Select the data variables for Creating Table1, including the top two phecode results.
rs4077469_T_data_table1 = select(rs4077469_T_data, phecode1, phecode2, age, age2, gender, SEX, rs4077469_T)
rs4077469_T_data_table1[1:10,]
print(summary(rs4077469_T_data_table1))
print(str(rs4077469_T_data_table1))


# Distribution of alleles
print(length(gdata3$rs4077469_T))
print(sum(is.na(gdata3$rs4077469_T)))
print(table(gdata3$rs4077469_T))

print(length(rs4077469_T_data_table1$rs4077469_T))
print(sum(is.na(rs4077469_T_data_table1$rs4077469_T)))
print(table(rs4077469_T_data_table1$rs4077469_T))


# Create and save table 1

variables_names = names(rs4077469_T_data_table1)
print(variables_names)

factor_names = c(phecode1,
                 phecode2,
                 "gender",
                 "SEX",
                 "rs4077469_T")
print(factor_names)


table1 = CreateTableOne(vars = variables_names
                        , factorVars = factor_names
                        , data = rs4077469_T_data_table1
                        #, strata = "SEX"
                        )

table_1 = print(table1, showAllLevels = T)
write.csv(table_1, file = "Table1_final_v2_rs4077469_T.csv")

table1 = CreateTableOne(vars = variables_names
                        , factorVars = factor_names
                        , data = rs4077469_T_data_table1
                        , strata = "SEX")
table_1 = print(table1, showAllLevels = T)
write.csv(table_1, file = "Table1_final_v2_rs4077469_T_gender.csv")


table1 = CreateTableOne(vars = variables_names
                        , factorVars = factor_names
                        , data = rs4077469_T_data_table1
                        , strata = "rs4077469_T")
table_1 = print(table1, showAllLevels = T)
write.csv(table_1, file = "Table1_final_v2_rs4077469_T_bySNP.csv")

table1 = CreateTableOne(vars = variables_names
                        , factorVars = factor_names
                        , data = rs4077469_T_data_table1
                        , strata = c(phecode1, "SEX"))
table_1 = print(table1, showAllLevels = T)
write.csv(table_1, file = "Table1_final_v2_rs4077469_T_gender_Uterine_leiomyoma.csv")



table1 = CreateTableOne(vars = variables_names
                        , factorVars = factor_names
                        , data = rs4077469_T_data_table1
                        , strata = c(phecode2, "SEX"))
table_1 = print(table1, showAllLevels = T)
write.csv(table_1, file = "Table1_final_v2_rs4077469_T_gender_BenignNeoplasm_of_Uterus.csv")


