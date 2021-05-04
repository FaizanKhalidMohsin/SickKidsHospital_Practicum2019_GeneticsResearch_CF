#install.packages("devtools")

#It may be necessary to install required as not all package dependencies are installed by devtools:
#install.packages(c("dplyr","tidyr","ggplot2","MASS","meta","ggrepel","DT"))

#library(devtools)
#install_github("PheWAS/PheWAS")

#devtools::install_github("PheWAS/PheWAS")

#install.packages(c("dplyr","tidyr","ggplot2","MASS","meta","ggrepel","DT"))
library(PheWAS)

#Set the random seed so it is replicable
set.seed(1)
#Generate some example data
ex=generateExample()

?generateExample

#Extract the two parts from the returned list
id.icd9.count=ex$id.icd9.count
genotypes=ex$genotypes
id.icd9.count$count = 1
dim(genotypes)
SEX = sample(c(1, 2), 5000, replace= T)
covariates = data.frame(id = genotypes$id, SEX)
dim(genotypes)
?sample
#Create the phecode table- translates the icd9s, adds 
#exclusions, and reshapes to a wide format
?createPhewasTable
phenotypes=createPhewasTable(id.icd9.count, min.code.count = 1)
phenotypes1 = phenotypes[1:500, ]
phemap
PheWAS::phecode_map()
?phe_map
?createPhewasTable
t0core1 = Sys.time()
?phewas
#Run the PheWAS
results=phewas(phenotypes,genotypes
               #, covariates = covariates
               , cores=2,
               significance.threshold=c("bonferroni"))
t1core1 = sys.time() - t0core1
#Plot the results
pdf("Phewasplotmincount1.pdf")
phewasManhattan(results, annotate.angle=0,
                title="My Example PheWAS Manhattan Plot")
dev.off()
#Add PheWAS descriptions
results_d=addPhecodeInfo(results)
names(results_d)
unique(results_d$group)
unique(results_d$description)
#List the significant results
results_d[results_d$bonferroni&!is.na(results_d$p),]
#List the top 10 results
results_d[order(results_d$p)[1:10],]
names(results_d)
#*********************************** FUNCTIONS *****************************************#
phewas=function (phenotypes, genotypes, data, covariates = c(NA), adjustments = list(NA), 
          outcomes, predictors, cores = 1, additive.genotypes = T, 
          significance.threshold=c("bonferroni"), alpha = 0.05, unadjusted = F, return.models = F, 
          min.records = 20, MASS.confint.level = NA, quick.confint.level) 
{
  if (missing(phenotypes)) {
    if (!missing(outcomes)) 
      phenotypes = outcomes
    else stop("Either phenotypes or outcomes must be passed in.")
  }
  if (missing(genotypes)) {
    if (!missing(predictors)) 
      genotypes = predictors
    else stop("Either genotypes or predictors must be passed in.")
  }
  association_method = phe_as
  if (unadjusted) {
    association_method = phe_as_unadjusted
    if (!is.na(covariates) | !is.na(adjustments)) 
      warning("Covariates and adjustments are ignored in unadjusted mode.")
  }
  if (missing(data)) {
    phe = phenotypes
    gen = genotypes
    cov = covariates
    adjustment = adjustments
    id = intersect(names(phenotypes), names(genotypes))
    if (length(id) == 0) {
      stop("There is no shared column to merge phenotypes and genotypes!")
    }
    message(paste("Merging data using these shared columns: ", 
                  id))
    phenotypes = names(phenotypes)
    phenotypes = phenotypes[!(phenotypes %in% id)]
    genotypes = names(genotypes)
    genotypes = genotypes[!(genotypes %in% id)]
    if (length(phenotypes) < 1 || length(genotypes) < 1) {
      stop("Either phenotypes or genotypes contained no non-shared columns, yielding no variables for analysis after the merge.")
    }
    data = merge(phe, gen, by = id)
    if (!is.null(names(covariates)[-1])) {
      covariates = names(covariates)
      if (sum(id %in% covariates) != length(id)) {
        stop(paste("The shared ID column(s) do not all exist in covariates: ", 
                   id))
      }
      covariates = covariates[!(covariates %in% id)]
      data = merge(data, cov, by = id)
    }
    if (!is.null(names(adjustments)[-1])) {
      adjustments = names(adjustments)
      if (sum(id %in% adjustments) != length(id)) {
        stop(paste("The shared ID column(s) do not all exist in adjustments: ", 
                   id))
      }
      adjustments = as.list(c(NA, adjustments[!(adjustments %in% 
                                                  id)]))
      data = merge(data, adjustment, by = id)
    }
  }
  para = (cores > 1)
  if (length(phenotypes) < 1 || length(genotypes) < 1) {
    stop("You must provide at least one genotype/predictor and one phenotype/outcome for analysis.")
  }
  full_list = data.frame(t(expand.grid(phenotypes, genotypes, 
                                       adjustments, stringsAsFactors = F)), stringsAsFactors = F)
  if (para) {
    if (exists("phewas.cluster.handle")) {
      message("Old cluster detected (phewas.cluster.handle), removing...")
      try(stopCluster(phewas.cluster.handle), silent = T)
      rm(phewas.cluster.handle, envir = .GlobalEnv)
    }
    message("Starting cluster...")
    assign("phewas.cluster.handle", makeCluster(cores), envir = .GlobalEnv)
    message("Cluster created, finding associations...")
    clusterExport(phewas.cluster.handle, c("data", "covariates"), 
                  envir = environment())
    result <- parLapplyLB(phewas.cluster.handle, full_list, 
                          association_method, additive.genotypes, confint.level = MASS.confint.level, 
                          min.records, return.models)
    stopCluster(phewas.cluster.handle)
    rm(phewas.cluster.handle, envir = .GlobalEnv)
  }
  else {
    message("Finding associations...")
    result = lapply(full_list, FUN = association_method, 
                    additive.genotypes, min.records, return.models, confint.level = MASS.confint.level, 
                    data, covariates)
  }
  if (return.models) {
    message("Collecting models...")
    models = lapply(result, function(x) {
      attributes(x)$model
    })
    names(models) = sapply(models, function(x) {
      paste0(as.character(terms(x))[c(2, 1, 3)], collapse = " ")
    })
  }
  message("Compiling results...")
  successful.phenotypes = na.omit(sapply(result, function(x) {
    attributes(x)$successful.phenotype
  }))
  n.tests = length(successful.phenotypes)
  successful.phenotypes = unique(successful.phenotypes)
  successful.genotypes = unique(na.omit(sapply(result, function(x) {
    attributes(x)$successful.genotype
  })))
  sig = bind_rows(result)
  if (max(grepl(pattern = "[Error: The model did not converge]", 
                sig$note, fixed = TRUE))) {
    warning("Not all models converged, check the notes column for details.")
  }
  message("Cleaning up...")
  attributes(sig)$alpha = alpha
  attributes(sig)$n.tests = n.tests
  if (!missing(significance.threshold)) {
    message("Finding significance thresholds...")
    thresh = match(c("p-value", "bonferroni", "fdr", "simplem-genotype", 
                     "simplem-phenotype", "simplem-product"), significance.threshold)
    sm.g = 1
    sm.p = 1
    if (!is.na(thresh[1])) {
      sig$p.value = sig$p <= alpha
    }
    if (!is.na(thresh[2])) {
      sig$bonferroni = sig$p <= alpha/n.tests
      attributes(sig)$bonferroni = alpha/n.tests
    }
    if (!is.na(thresh[3])) {
      sig$fdr = p.adjust(sig$p, method = "fdr") <= alpha
    }
    if (!is.na(thresh[4]) | !is.na(thresh[6])) {
      if (length(successful.genotypes) > 1) {
        eigs = eigen(cor(data[, genotypes], use = "pairwise.complete.obs", 
                         method = "spearman"))[[1]]
        max.eig = sum(eigs)
        sm.g = which.max(cumsum(eigs) > 0.995 * max.eig)
      }
      else {
        sm.g = 1
      }
      sig$simplem.genotype = sig$p <= alpha/sm.g
      attributes(sig)$simplem.genotype = alpha/sm.g
      attributes(sig)$simplem.genotype.meff = sm.g
    }
    if (!is.na(thresh[5]) | !is.na(thresh[6])) {
      if (length(successful.phenotypes > 1)) {
        eigs = try(cor(data[, successful.phenotypes], 
                       use = "pairwise.complete.obs", method = "spearman"), 
                   silent = T)
        if (class(eigs) != "try-error") {
          eigs[is.na(eigs)] = 0
          eigs = eigen(eigs)[[1]]
          max.eig = sum(eigs)
          sm.p = which.max(cumsum(eigs) > 0.995 * max.eig)
        }
        else {
          warning("Phentoype correlation generation failed; this is typically due to sparse phenotypes.")
          sm.p = NA
        }
      }
      else {
        sm.p = 1
      }
      sig$simplem.phenotype = sig$p <= alpha/sm.p
      attributes(sig)$simplem.phenotype = alpha/sm.p
      attributes(sig)$simplem.phenotype.meff = sm.p
    }
    if (!is.na(thresh[6])) {
      sm = sm.g * sm.p
      sig$simplem.product = sig$p <= alpha/sm
      attributes(sig)$simplem.product = alpha/sm
      attributes(sig)$simplem.product.meff = sm
    }
  }
  if (!missing(outcomes)) 
    names(sig)[names(sig) == "phenotype"] = "outcome"
  if (!missing(predictors)) 
    names(sig)[names(sig) == "snp"] = "predictor"
  if (return.models) {
    sig = list(results = sig, models = models)
  }
  if (!missing(quick.confint.level)) {
    if (quick.confint.level >= 1 | quick.confint.level <= 
        0) {
      warning("Quick confidence interval requested, but a value in the range (0,1) was not supplied")
    }
    else {
      sig.names = names(sig)
      two.sided = (1 - quick.confint.level)/2
      sig = sig %>% mutate(lower.q = beta + se * qnorm(two.sided), 
                           upper.q = beta + se * qnorm(two.sided, lower.tail = F))
      sig = sig %>% mutate(lower.q = ifelse(sig$type == 
                                              "logistic", exp(lower.q), lower.q), upper.q = ifelse(sig$type == 
                                                                                                     "logistic", exp(upper.q), upper.q))
      sig = sig[, c(sig.names[1:5], "lower.q", "upper.q", 
                    sig.names[6:length(sig.names)])]
    }
  }
  return(sig)
}
phewasManhattan <-
  function(d, annotate.phenotype.description=T, ...) {
    if(sum(c("phenotype","p") %in% names(d))<2 ) stop("Data input must contain columns phenotype and p.")
    if(class(d$phenotype)!="character") {
      if(class(d$phenotype)=="factor") {
        warning("Factor phenotype input mapped to characters")
        d$phenotype=as.character(d$phenotype)
      } else {
        stop("Non-character or non-factor phenotypes passed in, so an accurate phecode mapping is not possible.")
      }
    }
    #Check to see if it looks 0-padded
    if(min(nchar(d$phenotype))<3) warning("Phenotypes with length <3 observed, ensure they are are 0-padded (e.g., \"008\")")
    
    #Add the groups and phecode descriptions as requested
    d=addPhecodeInfo(d,groupnums =T, groupcolors = T) %>% rename(phenotype=phecode)
    
    phenotypeManhattan(d, annotate.phenotype.description=annotate.phenotype.description, ...)
  }


#Multiple Sclerosis
which(names(phe)==335)
549

ww=549
Y=phe[,ww]
X=gen[,2]
  

dat=na.omit(data.frame(cbind(Y,X)))

X=dat$X

model1=glm(dat[,1]~X,family=binomial)

summary(model1)
result$X548

jpeg('X3.0.0boxplot.jpg')
dev.off()

