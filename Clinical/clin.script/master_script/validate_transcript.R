
### Script to validate transcripts in Cavalli dataset
### Yura Grabovska and Marion Mateos July 17 2018

# source("https://bioconductor.org/biocLite.R") ### use http if https does not work
# biocLite() ### to install program


library(car)
library(stats)
library(survival)
library(BiocInstaller)

#biocLite(pkgs = "affy")
#biocLite("survival")
#biocLite(pkgs = "dplyr")

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
# eset <- readRDS("/home/yuri/eSet_for_Dan_GSE85218.rds")
eset <- readRDS ("/home/yuri/eSet_GSE85218_26Sep2018.rds") 

##obtain a dataframe of the phenotype factors from the paper # (Yura omitted the copy number data because it wasn’t straightforward to tabulate from the paper supplements) 

pd <- pData(eset)

eset.exprs <- exprs(eset) ## obtain the expression values as a matrix

rownames(eset.exprs) <- gsub("_at", "", rownames(eset.exprs)) ## gets rid of the "_at" after the ENGS IDs ### length = 763

head(eset.exprs)

# “/home/dan/mygit/consensus_scripts/taylor_data”

# taylor <- View("/home/dan/mygit/consensus_scripts/taylor_data")

## the newer versions of the R biomaRt package run a series of smaller queries which are less likely to timeout internally without error install.packages("/home/yuri/biomaRt_2.37.0.tar.gz", repos = NULL, type="source") ## devel level version

require(biomaRt)
useMart("ENSEMBL_MART_ENSEMBL") ### used to be useMart()
## load the 'mart' and tell the program we want to use Ensembl and specifically the human annotations 

listMarts ### sometimes temperamental and have to run the useMart or listMarts more than once, then rerun ensembl.mart line. It worked after I ran the listMarts command then listMarts() then reran ensembl.mart (18/9/18)
listMarts() ### added 18/9/18

# listMarts(mart = NULL, host="www.ensembl.org", path="/biomart/martservice",
         # port=80, includeHosts = FALSE, archive = FALSE, ssl.verifypeer = TRUE, 
         # ensemblRedirect = TRUE, verbose = FALSE)

# ?listMarts
ensembl.mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# ensembl.mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")


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

### can add in MYC and MYCN data here

eset_match <- read.csv(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/eset_master_20181015.csv", header = TRUE, sep = ",", quote = "\"", row.names = 1)
head(eset_match)

#rownames(eset_match)<- eset_match[,1]

names(eset_match)
eset_match$Sample_Name
head(eset_match)
View(pData(eset))

# colnames(eset_match)
### eset_match$Sample_Name this is the MB_SubtypeStudy_number

### next step is to move MYC and MYCN columns from eset_match into the eset that uses pData
### then can create this as main dataframe from which to subset expression datasets/ match in expression data. 

eset$MYC <- eset_match$MYC  ### this is OK as the sample numbers are in same order
eset$MYCN <- eset_match$MYCN
eset$meth7 <- eset_match$meth7
eset$q13loss <- eset_match$q13loss_YN
identical(rownames(pData(eset)),rownames(eset_match)) ### TRUE indicates that rownames are identical therefore valid to use above column additions

View(pData(eset)) 


##############################################################################################################
##############################################################################################################
### look to validate the goi in the Cavalli dataset

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/Validate_function.R")

### when ran the function above ### 

#### 
### need to add in the validation script for overall cohort to the validate_function.R

##########################################################################################################   

### running for individual goi
  
list.goi.MELK <- goi.validate(goi = "ENSG00000165304", data = eset )
list.goi.NQO2 <- goi.validate(goi =  "ENSG00000124588", data = eset )
# list.goi.LOC256880 <- goi.validate (goi = "ENSG00000245322", data = eset) ### error subscript out of bounds
list.goi.MRPS12 <- goi.validate (goi = "ENSG00000128626", data = eset)
list.goi.MLYCD <- goi.validate (goi = "ENSG00000103150", data = eset)
list.goi.AKTIP <- goi.validate(goi = "ENSG00000166971", data = eset)
list.goi.SLC25A34 <- goi.validate (goi ="ENSG00000162461", data = eset )
list.goi.POU6F1 <- goi.validate (goi = "ENSG00000184271", data = eset)
list.goi.CXXC4 <- goi.validate (goi ="ENSG00000168772", data = eset )
list.goi.CHEK1 <- goi.validate (goi = "ENSG00000149554", data = eset)
list.goi.CDK4 <- goi.validate (goi = "ENSG00000135446", data = eset)
list.goi.ROGDI <- goi.validate (goi ="ENSG00000067836", data = eset )

write.csv(list.goi.NQO2$cox_summary$summary_nogender, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/NQO2_cox_nogender.csv")
write.csv(list.goi.NQO2$cox_summary$summary_gender, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/NQO2_cox_gender.csv")
write.csv (list.goi.MRPS12[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/Validation/MRPS12_cox_nogender.csv")
write.csv (list.goi.MLYCD[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/Validation/MLYCD_cox_nogender.csv")
write.csv (list.goi.AKTIP[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/Validation/AKTIP_cox_nogender.csv")
write.csv (list.goi.SLC25A34[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/Validation/SLC25A34_cox_nogender.csv")
write.csv (list.goi.POU6F1[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/Validation/POU6F1_cox_nogender.csv")
write.csv (list.goi.CXXC4[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/Validation/CXXC4_cox_nogender.csv") 
write.csv (list.goi.CHEK1[[4]], file =  "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/Validation/CHEK1_cox_nogender.csv")
write.csv (list.goi.CDK4[[4]], file =  "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/Validation/CDK4_cox_nogender.csv")
write.csv (list.goi.ROGDI[[4]], file =  "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/Validation/ROGD1_cox_nogender.csv")

### saveRDS 

saveRDS(list.goi.NQO2, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/NQ02.rds")
saveRDS(list.goi.MELK, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/MELK.rds")

### write out specific .csv based on the files in the rds object

NQO2 <- readRDS(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/NQ02.rds")
MELK <- readRDS(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/MELK.rds")
# NQO2$pval



#########################################################
#########################################################

### see if data below can be used to improve the validate_function.R

# km.OS<- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
### keep this until determine km curves within the goi.validate function

km.OS.MELK <- survfit(Surv(eset$OS, eset$Dead)~MELK.cat, type = "kaplan-meier", conf.type= "log")
km.OS.MELK <- survfit(Surv(matched.eset$OS, matched.eset$Dead)~matched.eset$MELK, type = "kaplan-meier", conf.type= "log") ### generates same curve as above code
summary(km.OS.MELK)
plot(km.OS.MELK)

plot(matched.eset.G3G4$MELKexp, xlab = "individual samples", ylab = "MELK expression", main = "Expression of MELK in G3G4 validation cohort")
abline(h=median(matched.eset.G3G4$MELKexp), lty = 1, col = "red") ### v for vertical ie x axis ### h, for horizontal

km.OS.G3G4.MELK <- survfit(Surv(matched.eset.G3G4.incl$OS, matched.eset.G3G4.incl$Dead)~matched.eset.G3G4.incl$MELK, type = "kaplan-meier", conf.type= "log") ### generates same curve as above code
summary(km.OS.G3G4.MELK)
plot(km.OS.G3G4.MELK)


### multivariable cox regression analysis for 
### i) Lancet
### ii) PNET 5

cox.OS.G3G4.MELK.Lancet <- coxph(Surv(matched.eset.G3G4.incl$OS, matched.eset.G3G4.incl$Dead)~matched.eset.G3G4.incl$MELK + matched.eset.G3G4.incl$q13loss + matched.eset.G3G4.incl$meth7_HR +  matched.eset.G3G4.incl$Gender + matched.eset.G3G4.incl$MYC, data = matched.eset.G3G4.incl)

#str(summary(cox.OS.MELK))

cox.n.G3G4.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[4]] ### n             ###cox.OS.MELK$n          ### cox.OS.MELK[[11]] 
cox.nevent.G3G4.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[6]] ### nevent   ###cox.OS.MELK$nevent     ### cox.OS.MELK[[12]] 

summary(cox.OS.G3G4.MELK.Lancet)[[7]] ### this is the table of relevance p value
cox.pval.G3G4.MELK.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[7]][1,5] ### this accesses the p value for MELK (row 1, position 5)
cox.HR.G3G4.MELK.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[7]][1,2]
cox.lower.95CI.G3G4.MELK.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[8]][1,3]
cox.upper.95CI.G3G4.MELK.Lancet <- summary (cox.OS.G3G4.MELK.Lancet)[[8]][1,4]
summary.cox.G3G4.MELK.Lancet <- list(pval = cox.pval.G3G4.MELK.Lancet, HR = cox.HR.G3G4.MELK.Lancet, L95CI = cox.lower.95CI.G3G4.MELK.Lancet, U95CI =cox.upper.95CI.G3G4.MELK.Lancet, n = cox.n.G3G4.Lancet, nevent = cox.nevent.G3G4.Lancet, table = summary(cox.OS.G3G4.MELK.Lancet)[[7]], HR_table = summary(cox.OS.G3G4.MELK.Lancet)$conf.int)  

write.csv (summary.cox.G3G4.MELK.Lancet, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/MELK.multivar.G3G4.csv")

#####################################################################################################################

### ORIGINAL DATA FOR MELK, CAN RUN FOR VISUALISATION HOWEVER CAN GENERATE COX MODEL FROM VALIDATE FUNCTION ABOVE
### CAVEAT: the G3G4 subsetting appears to give different results as only 369 results were included in the manual function below

### match in data with MELK
### this worked below 17/7/18 - 8/8/18 for MELK

MELK <- eset.exprs["ENSG00000165304", ]
plot(MELK)
qqnorm(MELK) ### demonstrates that it is normally distributed
summary(MELK)
MELK.cat <- ifelse(MELK>median(MELK, na.rm = T), "high", "low")

index <- match(names(MELK.cat), rownames(eset))
matched.eset <- eset[index[!is.na(index)],]

summary_cavalli <- summary(pData(matched.eset))
# write.csv (summary_cavalli, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/summary.cavalli.csv")

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

### include in multivariable cox regression analysis adding in MYC or MYCN

cox.OS.MELK <- coxph(Surv(matched.eset$OS, matched.eset$Dead)~matched.eset$MELK + matched.eset$LCA + matched.eset$Gender + matched.eset$Subgroup + matched.eset$mets + matched.eset$MYC + matched.eset$MYCN, data = matched.eset)

#str(summary(cox.OS.MELK))

cox.n <- summary(cox.OS.MELK)[[4]] ### n             ###cox.OS.MELK$n          ### cox.OS.MELK[[11]] 
cox.nevent <- summary(cox.OS.MELK)[[6]] ### nevent   ###cox.OS.MELK$nevent     ### cox.OS.MELK[[12]] 

summary(cox.OS.MELK)[[7]] ### this is the table of relevance p value
cox.pval.MELK <- summary(cox.OS.MELK)[[7]][1,5] ### this accesses the p value for MELK (row 1, position 5)
cox.HR.MELK <- summary(cox.OS.MELK)[[7]][1,2]
cox.lower.95CI.MELK <- summary(cox.OS.MELK)[[8]][1,3]
cox.upper.95CI.MELK <- summary (cox.OS.MELK)[[8]][1,4]
summary.cox.MELK <- list(pval = cox.pval.MELK, HR = cox.HR.MELK, L95CI = cox.lower.95CI.MELK, U95CI =cox.upper.95CI.MELK, n = cox.n, nevent = cox.nevent, table = summary(cox.OS.MELK)[[7]])  

### putting in the Lancet Oncology factors including q13 loss (added 15/10/18)

cox.OS.MELK.Lancet <- coxph(Surv(matched.eset$OS, matched.eset$Dead)~matched.eset$MELK + matched.eset$meth7 + matched.eset$MYC + matched.eset$q13loss + matched.eset$Gender)

summary(cox.OS.MELK.Lancet)

### create function to output these p value and HR characteristics

### Plot transcript expression in Cavalli dataset, unfiltered, with median expression as cutoff

plot(MELK, xlab = "individual samples", ylab = "MELK expression", main = "Expression of MELK in validation cohort")
abline(h=median(MELK), lty = 1, col = "red") ### v for vertical ie x axis ### h, for horizontal

MELK.high <- which (MELK > median(MELK))
length(MELK.high)
MELK.low <- which (MELK < median (MELK))
length (MELK.low)

#########################################################################################################
### look at expression in Group 3/Group 4. Generate dataframe for all expression data (rather than just matched.eset for MELK)
### rename the matched.eset here
matched.eset.all <- matched.eset 

# sub <- matched.eset.all$Subgroup
sub <- matched.eset.all@phenoData@data$Subgroup ### an alternate way to access the expression data subcolumns

Group3.exp <- matched.eset.all[,which(sub=="Group3")]
fdata <- featureData(matched.eset.all)
Group3.exp@featureData <- fdata ### upload feature data so that it contains data
pData(Group3.exp) ### 144 samples, all expression features (21641)


### Group 3 and 4
G3G4.exp <- matched.eset.all[, which(sub=="Group3"|sub=="Group4")]
pData(G3G4.exp)    ### 470 samples, expression set retained

### matched.eset.all is an expression set that can be used to match any expression set in (goi), pData(eset) converts to dataframe
View(pData(matched.eset.all))

### create matched G3G4 dataframe with MELK expression data
### note that MELK is called in from earlier, therefore already in G3G4.exp dataset

View(pData(G3G4.exp))

matched.eset.G3G4 <- pData(G3G4.exp)
# class(matched.eset.G3G4)


###############################

### Exploring characteristics of the dataset, note that expression dataset needs the prefix(pData(eset)) to convert to dataframe
# dim(G3G4.exp)
# names(pData(G3G4.exp))
# summary(pData(G3G4.exp))
# str(pData(G3G4.exp))
# class(pData(G3G4.exp))
# rownames(pData(G3G4.exp))
# plot(G3G4.exp$MELKexp)
# G3G4.exp$MELK

### MELK-specific script #################################################################

plot(matched.eset.G3G4$MELKexp, xlab = "individual samples", ylab = "MELK expression", main = "Expression of MELK in G3G4 validation cohort")
abline(h=median(matched.eset.G3G4$MELKexp), lty = 1, col = "red") ### v for vertical ie x axis ### h, for horizontal

# MELK.cat.G3G4 <- ifelse(matched.eset.G3G4$MELK>median(matched.eset.G3G4$MELK, na.rm = T), "high", "low")
# index <- match(names(MELK.cat.G3G4), rownames(G3G4.exp))
# matched.eset.G3G4$MELK.cat <- MELK.cat.G3G4
# matched.eset.G3G4$MELKexp <- matched.eset.G3G4$MELK

### create new dataframe to remove the two samples which are coded as SHH_Inf and WNT

matched.eset.G3G4.incl <- matched.eset.G3G4[(matched.eset.G3G4$meth7=="Grp3_LowRisk"|matched.eset.G3G4$meth7=="Grp3_HighRisk"|matched.eset.G3G4$meth7=="Grp4_LowRisk"|matched.eset.G3G4$meth7=="Grp4_HighRisk"),  ]
# matched.eset.G3G4.excl <- matched.eset.G3G4[(matched.eset.G3G4$meth7=="SHH_Inf"|matched.eset.G3G4$meth7=="WNT"),  ] ### note can also use != nomenclature (does not equal)


matched.eset.G3G4.incl$meth7_HR <- matched.eset.G3G4.incl$meth7=="Grp4_HighRisk"|matched.eset.G3G4.incl$meth7 =="Grp3_HighRisk"

# km.OS<- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
km.OS.G3G4.MELK <- survfit(Surv(matched.eset.G3G4.incl$OS, matched.eset.G3G4.incl$Dead)~matched.eset.G3G4.incl$MELK, type = "kaplan-meier", conf.type= "log") ### generates same curve as above code
summary(km.OS.G3G4.MELK)
plot(km.OS.G3G4.MELK)


### multivariable cox regression analysis for 
### i) Lancet
### ii) PNET 5

cox.OS.G3G4.MELK.Lancet <- coxph(Surv(matched.eset.G3G4.incl$OS, matched.eset.G3G4.incl$Dead)~matched.eset.G3G4.incl$MELK + matched.eset.G3G4.incl$q13loss + matched.eset.G3G4.incl$meth7_HR  + matched.eset.G3G4.incl$MYC + matched.eset.G3G4.incl$Gender, data = matched.eset.G3G4.incl)
cox.OS.G3G4.Lancet <- coxph(Surv(matched.eset.G3G4.incl$OS, matched.eset.G3G4.incl$Dead)~ matched.eset.G3G4.incl$q13loss + matched.eset.G3G4.incl$meth7_HR  + matched.eset.G3G4.incl$MYC, data = matched.eset.G3G4.incl )

cox.OS.G3G4.Lancet.output <- summary(cox.OS.G3G4.Lancet)$conf.int

# write.csv(cox.OS.G3G4.Lancet.output, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/Lancet.G3G4.csv")

#str(summary(cox.OS.MELK))

#####################################################################################################################################

cox.n.G3G4.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[4]] ### n             ###cox.OS.MELK$n          ### cox.OS.MELK[[11]] 
cox.nevent.G3G4.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[6]] ### nevent   ###cox.OS.MELK$nevent     ### cox.OS.MELK[[12]] 

summary(cox.OS.G3G4.MELK.Lancet)[[7]] ### this is the table of relevance p value
cox.pval.G3G4.MELK.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[7]][1,5] ### this accesses the p value for MELK (row 1, position 5)
cox.HR.G3G4.MELK.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[7]][1,2]
cox.lower.95CI.G3G4.MELK.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[8]][1,3]
cox.upper.95CI.G3G4.MELK.Lancet <- summary (cox.OS.G3G4.MELK.Lancet)[[8]][1,4]
summary.cox.G3G4.MELK.Lancet <- list(pval = cox.pval.G3G4.MELK.Lancet, HR = cox.HR.G3G4.MELK.Lancet, L95CI = cox.lower.95CI.G3G4.MELK.Lancet, U95CI =cox.upper.95CI.G3G4.MELK.Lancet, n = cox.n.G3G4.Lancet, nevent = cox.nevent.G3G4.Lancet, table = summary(cox.OS.G3G4.MELK.Lancet)[[7]], HR_table = summary(cox.OS.G3G4.MELK.Lancet)$conf.int)  

write.csv (summary.cox.G3G4.MELK.Lancet, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/MELK.multivar.G3G4.csv")

# cox.list.95CI <- summary(cox.OS.G3G4.MELK.Lancet)$conf.int  ### this brings up all the hazard ratios and 95CI range



