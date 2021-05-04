# R script for running phewas for genes in the X-Chromosome. 
# This is for Males
#
# Things to change for new Gene and SNP (not an exhaustive list)
# ukb_chrX_rs5905176.raw
# SLC6A14 
# rs5905176
# rs5905176_G
# ChX

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

gdata0 = fread("ukb_chrX_rs5905176.raw", stringsAsFactors=F, header=T, na.strings=c(""," ","NA"))
gdata0[1:10,]
print(dim(gdata0))
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
dim(ethnic_data0)
ethnic_data = filter(ethnic_data0, !is.na(caucasian))
dim(ethnic_data)
print(str(ethnic_data))
summary(ethnic_data)


gdata3 = inner_join(gdata2, ethnic_data, by = "id")
print(dim(gdata3))
summary(gdata3)
str(gdata3)
genotypes = gdata3 %>% select( id, rs5905176_G)
genotypes[1:10,]


covariates_age_gender = fread("pheno_age_sex.tab",  stringsAsFactors=F, 
                              col.names = c("id", "age", "SEX"), 
                              header=T, na.strings=c(""," ","NA"))
print(dim(covariates_age_gender))
covariates1 = mutate(covariates_age_gender, age2 = age**2)
covariates = semi_join(covariates1, genotypes, by = "id")
print(dim(covariates))
summary(covariates)
print(str(covariates))
print(covariates[1:10,])

females = filter(covariates, SEX == 0)
print(dim(females))
print(str(females))

covariates = females %>% select(-SEX)
print(dim(covariates))
print(str(covariates))

genotypes = semi_join(genotypes, females, by = "id")

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

phenotypes = semi_join(phenotypes, females, by = "id")

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
print(sum(is.na(genotypes1$rs5905176_G)))

print(dim(genotypes))
print(dim(phenotypes))
print(dim(covariates))


#Run the unadjusted PheWAS (unvariable analysis)

results_uni=phewas(phenotypes, genotypes
               , cores=1
               , significance.threshold=c("bonferroni")
               )


#Plot the results
pdf("phewasplot_uni_SLC6A14_rs5905176_G_females.pdf")
phewasManhattan(results_uni, annotate.angle=0,
                title="Manhattan Plot for SLC6A14 & rs5905176 for Females")
dev.off()

#Add PheWAS descriptions
results_uni_d=addPhecodeInfo(results_uni)
#List the significant results
results_uni_d[results_uni_d$bonferroni&!is.na(results_uni_d$p),]
#List the top 10 results
results_uni_d[order(results_uni_d$p)[1:10],]

# Save the top 10 results
r = results_uni_d[order(results_uni_d$p)[1:10],]
write.csv(r, "Results_uni_Top10_SLC6A14_rs5905176_G_females.csv")

# Save the entire PheWAS Study results
write.table(results_uni_d, "phewasresults_uni_SLC6A14_rs5905176_G_females.tsv"
            ,quote = F, row.names = F
            , col.names = T, sep = "\t")


# Save "results" to plot the PheWAS Study later again if need be.
write.table(results_uni, "phewasresults_uni_SLC6A14_rs5905176_G_data_forplotting_females.tsv"
            ,quote = F, row.names = F
            , col.names = T, sep = "\t")




#Run the adjusted PheWAS with covariates
results=phewas(phenotypes, genotypes
               , covariates = covariates 
               , cores=1
               , significance.threshold=c("bonferroni")
               )


#Plot the results
pdf("phewasplot_SLC6A14_rs5905176_G_covariates_females.pdf")
phewasManhattan(results, annotate.angle=0,
                title="Manhattan Plot for SLC6A14 & rs5905176 with Covariates for Females")
dev.off()

#Add PheWAS descriptions
results_d=addPhecodeInfo(results)
#List the significant results
results_d[results_d$bonferroni&!is.na(results_d$p),]
#List the top 10 results
results_d[order(results_d$p)[1:10],]

# Save the top 10 results
r = results_d[order(results_d$p)[1:10],]
write.csv(r, "Results_Top10_SLC6A14_rs5905176_G_covariates_females.csv")

# Save the entire PheWAS Study results
write.table(results_d, "phewasresults_SLC6A14_rs5905176_G_covariates_females.tsv"
            ,quote = F, row.names = F
            , col.names = T, sep = "\t")


# Save "results" to plot the PheWAS Study later again if need be.
write.table(results, "phewasresults_SLC6A14_rs5905176_G_data_forplotting_females.tsv"
            ,quote = F, row.names = F
            , col.names = T, sep = "\t")


# Store the first two phecodes. 
phecode1 = as.character(r$phecode[1])
phecode2 = as.character(r$phecode[2])

#gdata3 = rename(gdata3, gender = SEX)
rs5905176_G_data0 = inner_join(gdata3, phenotypes,  by="id") %>% filter(!is.na(rs5905176_G))
rs5905176_G_data = inner_join(rs5905176_G_data0, covariates, by="id")
print(dim(rs5905176_G_data))

#str(rs5905176_G_data[,1:10])
#head(rs5905176_G_data$"218.1")
# write.table(rs5905176_G_data, "data_rs5905176_G_phecodes_covariates_final_data.tsv"
#             ,quote = F, row.names = F
#             , col.names = T, sep = "\t")

#rs5905176_G_data_table1 = select(rs5905176_G_data, '218.1', '218', age, age2, gender, SEX, rs5905176_G)
# Select the data variables for Creating Table1, including the top two phecode results.
rs5905176_G_data_table1 = select(rs5905176_G_data, phecode1, phecode2, age, age2, rs5905176_G)
rs5905176_G_data_table1[1:10,]
print(summary(rs5905176_G_data_table1))
print(str(rs5905176_G_data_table1))


# Distribution of alleles
print(length(gdata3$rs5905176_G))
print(sum(is.na(gdata3$rs5905176_G)))
print(table(gdata3$rs5905176_G))

print(length(rs5905176_G_data_table1$rs5905176_G))
print(sum(is.na(rs5905176_G_data_table1$rs5905176_G)))
print(table(rs5905176_G_data_table1$rs5905176_G))


# Create and save table 1's

variables_names = names(rs5905176_G_data_table1)
print(variables_names)

factor_names = c(phecode1,
                 phecode2,
                 "rs5905176_G")
print(factor_names)


table1 = CreateTableOne(vars = variables_names
                        , factorVars = factor_names
                        , data = rs5905176_G_data_table1
                        #, strata = "SEX"
                        )

table_1 = print(table1, showAllLevels = T)
write.csv(table_1, file = "Table1_final_SLC6A14_rs5905176_G.csv")

# table1 = CreateTableOne(vars = variables_names
#                         , factorVars = factor_names
#                         , data = rs5905176_G_data_table1
#                         , strata = "gender")
# table_1 = print(table1, showAllLevels = T)
# write.csv(table_1, file = "Table1_final_SLC6A14_rs5905176_G_gender.csv")


table1 = CreateTableOne(vars = variables_names
                        , factorVars = factor_names
                        , data = rs5905176_G_data_table1
                        , strata = "rs5905176_G")
table_1 = print(table1, showAllLevels = T)
write.csv(table_1, file = "Table1_final_SLC6A14_rs5905176_G_bySNP.csv")

table1 = CreateTableOne(vars = variables_names
                        , factorVars = factor_names
                        , data = rs5905176_G_data_table1
                        , strata = c(phecode1))
table_1 = print(table1, showAllLevels = T)
write.csv(table_1, file = "Table1_final_SLC6A14_rs5905176_G_phecode1_females.csv")


table1 = CreateTableOne(vars = variables_names
                        , factorVars = factor_names
                        , data = rs5905176_G_data_table1
                        , strata = c(phecode2))
table_1 = print(table1, showAllLevels = T)
write.csv(table_1, file = "Table1_final_SLC6A14_rs5905176_G_phecode2_females.csv")


# Create input for Flow Chart.

sink('SLC6A14_rs5905176_G_ChX_females_flow_chart.txt')

print("dim(gdata0)")
print(dim(gdata0))

print("dim(unrelated)")
print(dim(unrelated))

print("dim(unrelated1)")
print(dim(unrelated1))

print("gdata2 = inner_join(unrelated1, gdata1)")
print("dim(gdata2)")
print(dim(gdata2))

print("dim(ethnic_data0)")
print(dim(ethnic_data0))

print("ethnic_data = filter(ethnic_data0, !is.na(caucasian))")
print("dim(ethnic_data)")
print(dim(ethnic_data))

print("gdata3 = inner_join(gdata2, ethnic_data)")
print(dim(gdata3))
print(dim(gdata3))

print("dim(covariates_age_gender)")
print(dim(covariates_age_gender))

print("dim(females)")
print(dim(females))

print("covariates = semi_join(covariates1, genotypes)")
print("dim(covariates)")
print(dim(covariates))

print("dim(phenotypes)")
print(dim(phenotypes))

pritn("genotypes_1 = semi_join(genotypes, phenotypes)")
print("dim(genotypes_1)")
print(dim(genotypes_1))

print("sum(is.na(genotypes1$rs5905176_G))")
print(sum(is.na(genotypes1$rs5905176_G)))

sink()







