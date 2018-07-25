
### Script to validate transcripts in Cavalli dataset
### Yura Grabovska and Marion Mateos July 17 2018

# source("https://bioconductor.org/biocLite.R") ### use http if https does not work
# biocLite()


library(car)
library(stats)
library(survival)
library(BiocInstaller)

biocLite(pkgs = "affy")
biocLite("survival")

library(affy)
library(biomaRt)
# install.packages("affy")

# sessionInfo()

require(affy)

### unsure if need to run parallel but may help 

# library(foreach)
# library(tictoc)
library(parallel)
library(doParallel)
registerDoParallel(16)

##expression array data
eset <- readRDS("/home/yuri/eSet_for_Dan_GSE85218.rds")

##obtain a dataframe of the phenotype factors from the paper # (I omitted the copy number data because it wasn’t straightforward to tabulate from the paper supplements) 

pd <- pData(eset)

eset.exprs <- exprs(eset) ## obtain the expression values as a matrix

rownames(eset.exprs) <- gsub("_at", "", rownames(eset.exprs)) ## gets rid of the "_at" after the ENGS IDs ### length = 763

head(eset.exprs)

# “/home/dan/mygit/consensus_scripts/taylor_data”

# taylor <- View("/home/dan/mygit/consensus_scripts/taylor_data")

## the newer versions of the R biomaRt package run a series of smaller queries which are less likely to timeout internally without error install.packages("/home/yuri/biomaRt_2.37.0.tar.gz", repos = NULL, type="source") ## devel level version

require(biomaRt)
useMart()
## load the 'mart' and tell the program we want to use Ensembl and specifically the human annotations 

# listMarts ### sometimes temperamental and have to run the useMart or listMarts more than once, then rerun ensembl.mart line


ensembl.mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

## set up the query: attributes are all the things we want to pull and filters are the IDs we want to filter by, takes 5-10 minutes without doParallel
hugene.anno <- getBM(attributes = c("ensembl_gene_id",
"hgnc_symbol",
"external_gene_name",
"description",
"chromosome_name",
"start_position",
"end_position",
"strand",
"band",
"gene_biotype",
"go_id",
"name_1006",
"definition_1006",
"go_linkage_type"),
filters = "ensembl_gene_id",
values = rownames(eset.exprs),
mart = ensembl.mart,
verbose = FALSE) ## setting verbose to TRUE will throw out all the raw API interface calls

## collapse data to unique gene-level, with all information concatenated by a separator; it's a messy call but it works and is relatively fast
require(plyr)
hugene.anno2 <- ddply(hugene.anno, .(ensembl_gene_id), summarize,
                     external_gene_name_u = paste(unique(external_gene_name), collapse=","),
                     hgnc_symbol_u = paste(unique(hgnc_symbol), collapse=","),
                     description_u = paste(unique(description), collapse=","),
                     chromosome_name_u = paste(unique(chromosome_name), collapse=","),
                     start_position_u = paste(unique(start_position), collapse=","),
                     end_position_u = paste(unique(end_position), collapse=","),
                     strand_u = paste(unique(strand), collapse=","),
                     band_u = paste(unique(band), collapse=","),
                     gene_biotype_u = paste(unique(gene_biotype), collapse=","),
                     go_id_u = paste(unique(go_id), collapse=","),
                     name_1006_u = paste(unique(name_1006), collapse=","),
                     definition_1006_u = paste(unique(definition_1006), collapse=","),
                     go_linkage_type_u = paste(unique(go_linkage_type), collapse=","))

# write.csv(hugene.anno, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/hugene.anno.20180717.csv")

### now to validate the transcripts in survival analysis
View(pData(eset))

eset$OS <- eset$OS_.years.

### goi of interest
### create categorical variable
### then run survival analysis

# goi <- "ENSG00000165304"
# goi.df <- eset.exprs[rowname = goi,] 
# goi.cat <- ifelse(goi >median(goi, na.rm = T), "high","low")


### match in data
### this worked below 17/7/18 for MELK

MELK <- eset.exprs["ENSG00000165304", ]
# plot(MELK)
# qqnorm(MELK) ### demonstrates that it is normally distributed

MELK.cat <- ifelse(MELK>median(MELK, na.rm = T), "high", "low")

index <- match(names(MELK.cat), rownames(eset))
matched.eset <- eset[index[!is.na(index)],]

# View(matched.eset)

### add in MELK categorical variable into the dataset directly to compare with OS outcomes, and look at G3G4 separately

matched.eset$MELK <- MELK.cat

# km.OS<- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
km.OS.MELK <- survfit(Surv(eset$OS, eset$Dead)~MELK.cat, type = "kaplan-meier", conf.type= "log")
km.OS.MELK <- survfit(Surv(matched.eset$OS, matched.eset$Dead)~matched.eset$MELK, type = "kaplan-meier", conf.type= "log") ### generates same curve as above code
summary(km.OS.MELK)
plot(km.OS.MELK)

### create relevant variables to then use to adjust in multivariable analysis

matched.eset$LCA  <- ifelse((matched.eset$histology=="LCA"), "LCA", "non-LCA")
matched.eset$mets <- ifelse(matched.eset$Met.status_.1.Met._0.M0.=="1","metastatic","non-metastatic")


### include in multivariable cox regression analysis

cox.OS.MELK <- coxph(Surv(matched.eset$OS, matched.eset$Dead)~matched.eset$MELK + matched.eset$LCA + matched.eset$Gender + matched.eset$Subgroup + matched.eset$mets, data = matched.eset)

#str(summary(cox.OS.MELK))

cox.n <- summary(cox.OS.MELK)[[4]] ### n             ###cox.OS.MELK$n          ### cox.OS.MELK[[11]] 
cox.nevent <- summary(cox.OS.MELK)[[6]] ### nevent   ###cox.OS.MELK$nevent     ### cox.OS.MELK[[12]] 

summary(cox.OS.MELK)[[7]] ### this is the table of relevance p value
cox.pval.MELK <- summary(cox.OS.MELK)[[7]][1,5] ### this accesses the p value for MELK (row 1, position 5)
cox.HR.MELK <- summary(cox.OS.MELK)[[7]][1,2]
cox.lower.95CI.MELK <- summary(cox.OS.MELK)[[8]][1,3]
cox.upper.95CI.MELK <- summary (cox.OS.MELK)[[8]][1,4]
summary.cox.MELK <- list(pval = cox.pval.MELK, HR = cox.HR.MELK, L95CI = cox.lower.95CI.MELK, U95CI =cox.upper.95CI.MELK, n = cox.n, nevent = cox.nevent, table = summary(cox.OS.MELK)[[7]])  
  

### example from cox multivariate in clinical_data_functions

# cox.multivar.surv_8 <- function (time, event, marker, FacA, FacB, FacC, FacD, FacE, FacF, FacG, FacH, strata = NULL, data) {
  # if(is.null(strata)){
  #   cox.temp <- coxph(Surv(time, event)~marker + FacA + FacB + FacC + FacD + FacE + FacF + FacG + FacH, data=data)
  # }else {
 #    cox.temp <- coxph(Surv(time, event)~marker + FacA + FacB +FacC +FacD + FacE + FacF + FacG + FacH, data=data)
 #  }  
 #  cox.p.val <- summary(cox.temp)$coefficients[1,5] ### updated 14/11
 #  cox.HR <- summary(cox.temp)$coefficients[1,2] ### updated 14/11
 # cox.lower.95CI <- summary(cox.temp)$conf.int[1,3] ### as now multivariate, therefore need to access 1st row results
 # cox.upper.95CI <- summary(cox.temp)$conf.int[1,4]
 # cox.Zscore <- summary(cox.temp)$coefficients[1,4] ### added this in to access Z score
 # cox.n <-summary(cox.temp)$n
 # cox.nevent <-summary(cox.temp)$nevent
 # summary.cox <- list(cox.pval = cox.p.val,cox.HR = cox.HR, cox.lower.95CI = cox.lower.95CI, cox.upper.95CI =cox.upper.95CI, cox.Zscore = cox.Zscore, n = cox.n, n.event = cox.nevent)
 # return (summary.cox)
}

### restrict to G3G4 cohort
View(matched.eset)
----------------------------------------------------------------------------------------------

### example of code that was trialled
# rownames(eset) ### ensemblID
# names(MELK.cat)### sample names
# colnames(eset.exprs) ### sample names ### also colnames(eset)
#colnames(eset) ### sample names

# eset.expr.high<- which(eset.exprs["ENSG00000165304", ] < median(eset.exprs["ENSG00000165304", ]))
#index <- match(colnames(eset), names(MELK.cat)) ### did not work

# matched.eset <- eset[index[!is.na(index)], ] ### did not work 
