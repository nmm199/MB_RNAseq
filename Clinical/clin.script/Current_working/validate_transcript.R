
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
biocLite(pkgs = "dplyr")

library(affy)
library(biomaRt)
library (ggplot2)
library(dplyr)

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


### match in data with MELK
### this worked below 17/7/18 - 8/8/18 for MELK

MELK <- eset.exprs["ENSG00000165304", ]
# plot(MELK)
# qqnorm(MELK) ### demonstrates that it is normally distributed

MELK.cat <- ifelse(MELK>median(MELK, na.rm = T), "high", "low")

index <- match(names(MELK.cat), rownames(eset))
matched.eset <- eset[index[!is.na(index)],]

### can add in MYC and MYCN data here

### add in MELK categorical variable into the matched dataset directly to compare with OS outcomes

matched.eset$MELK <- MELK.cat
matched.eset$MELKexp <- MELK

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
  
### Plot transcript expression in Cavalli dataset, unfiltered, with median expression as cutoff

plot(MELK, xlab = "individual samples", ylab = "MELK expression", main = "Expression of MELK in validation cohort")
abline(h=median(MELK), lty = 1, col = "red") ### v for vertical ie x axis ### h, for horizontal
  
MELK.high <- which (MELK > median(MELK))
length(MELK.high)
MELK.low <- which (MELK < median (MELK))
length (MELK.low)
 
### look at expression in Group 3/Group 4. Generate dataframe for all expression data (rather than just matched.eset for MELK)
### rename the matched.eset here
matched.eset.test <- eset 

# sub <- matched.eset.all$Subgroup

# matched.eset$LCA  <- ifelse((matched.eset$histology=="LCA"), "LCA", "non-LCA")
# matched.eset$mets <- ifelse((matched.eset$Met.status_.1.Met._0.M0.=="1"),"metastatic","non-metastatic")
sub <- matched.eset.test@phenoData@data$Subgroup ### an alternate way to access the expression data subcolumns
Group3.exp <- matched.eset.test[,which(sub=="Group3")]
fdata <- featureData(matched.eset.test)
Group3.exp@featureData <- fdata ### upload feature data so that it contains data
pData(Group3.exp) ### 144 samples, all expression features (21641)


pData(Group3.Expression)
### Group 3 and 4
G3G4.exp <- matched.eset.test[, which(sub=="Group3"|sub=="Group4")]
pData(G3G4.exp)
### creating G3G4 subgroup

#sub <- matched.eset@phenoData@data$Subgroup
#sel <- which(sub=="Group3")
#subset <- sub[sel]
#sub_expression <- matched.eset[,subset]
#g3 <- exprs(matched.eset)[,matched.eset@phenoData@data$Subgroup=="Group3"]
#matched.eset [matched.eset$Subgroup == "Group3", ]

### matched.eset.all is a dataframe that can be used to match any expression set in (goi)

## goi 

### create matched G3G4 dataframe with MELK expression data

matched.eset.G3G4 <- matched.eset.test[index[!is.na(index)],]
matched.eset.G3G4$MELK   <- eset.exprs["ENSG00000165304", ]
plot(matched.eset.G3G4$MELK)

plot(matched.eset.G3G4$MELK, xlab = "individual samples", ylab = "MELK expression", main = "Expression of MELK in G3G4 validation cohort")
abline(h=median(matched.eset.G3G4$MELK), lty = 1, col = "red") ### v for vertical ie x axis ### h, for horizontal

MELK.cat.G3G4 <- ifelse(matched.eset.G3G4$MELK>median(matched.eset.G3G4$MELK, na.rm = T), "high", "low")
index <- match(names(MELK.cat.G3G4), rownames(G3G4.exp))
matched.eset.G3G4$MELK.cat <- MELK.cat.G3G4
# matched.eset.G3G4$MELKexp <- matched.eset.G3G4$MELK


# km.OS<- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
km.OS.G3G4.MELK <- survfit(Surv(matched.eset.G3G4$OS, matched.eset.G3G4$Dead)~matched.eset.G3G4$MELK.cat, type = "kaplan-meier", conf.type= "log") ### generates same curve as above code
summary(km.OS.G3G4.MELK)
plot(km.OS.G3G4.MELK)

### create relevant variables to then use to adjust in multivariable analysis

matched.eset.G3G4$LCA  <- ifelse((matched.eset.G3G4$histology=="LCA"), "LCA", "non-LCA")
matched.eset.G3G4$mets <- ifelse(matched.eset.G3G4$Met.status_.1.Met._0.M0.=="1","metastatic","non-metastatic")


### include in multivariable cox regression analysis

cox.OS.G3G4.MELK <- coxph(Surv(matched.eset.G3G4$OS, matched.eset.G3G4$Dead)~matched.eset.G3G4$MELK + matched.eset.G3G4$LCA + matched.eset.G3G4$Gender + matched.eset.G3G4$S + matched.eset.G3G4$mets, data = matched.eset.G3G4)

#str(summary(cox.OS.MELK))

cox.n.G3G4 <- summary(cox.OS.G3G4.MELK)[[4]] ### n             ###cox.OS.MELK$n          ### cox.OS.MELK[[11]] 
cox.nevent.G3G4 <- summary(cox.OS.G3G4.MELK)[[6]] ### nevent   ###cox.OS.MELK$nevent     ### cox.OS.MELK[[12]] 

summary(cox.OS.G3G4.MELK)[[7]] ### this is the table of relevance p value
cox.pval.G3G4.MELK <- summary(cox.OS.G3G4.MELK)[[7]][1,5] ### this accesses the p value for MELK (row 1, position 5)
cox.HR.G3G4.MELK <- summary(cox.OS.G3G4.MELK)[[7]][1,2]
cox.lower.95CI.G3G4.MELK <- summary(cox.OS.G3G4.MELK)[[8]][1,3]
cox.upper.95CI.G3G4.MELK <- summary (cox.OS.G3G4.MELK)[[8]][1,4]
summary.cox.G3G4.MELK <- list(pval = cox.pval.MELK, HR = cox.HR.MELK, L95CI = cox.lower.95CI.MELK, U95CI =cox.upper.95CI.MELK, n = cox.n, nevent = cox.nevent, table = summary(cox.OS.MELK)[[7]])  


### need to get high risk low risk G3G4 status and also to pull in MYC data
----------------------------------------------------------------------------------------------

### example of code that was trialled
# rownames(eset) ### ensemblID
# names(MELK.cat)### sample names
# colnames(eset.exprs) ### sample names ### also colnames(eset)
#colnames(eset) ### sample names

# eset.expr.high<- which(eset.exprs["ENSG00000165304", ] < median(eset.exprs["ENSG00000165304", ]))
#index <- match(colnames(eset), names(MELK.cat)) ### did not work

# matched.eset <- eset[index[!is.na(index)], ] ### did not work 


