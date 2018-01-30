
### Script for analysing the output from univariate analysis of RNA expression data compared to survival outcomes
### Date: September 25 2017
### Author: Dr Marion Mateos

### file input

results.master <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.master.allgenes.20180104.rds") 
# results.master <- readRDS(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/superceded/results.master.allgenes.rds")
# results.master <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.10.051017.rds") ### has cox Z score
# results.master <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.novel.12100.12150.rds")
# results.master <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.novel.rds")  ### currently error as file has error

### read in functions file

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_master.R")
source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_extract_master.R")

 ### file output 

### dataframes with p value, adjusted p value
### all.survival.p.bothgroups: includes OS, EFS and PFS for overall group and G3G4 
### significant adjusted p values for OS, EFS, PFS separately in dataframes for overall group and for G3G4
### graphical output for adjusted p values compared to p values
### univariate cox regression p values, Z scores, HR, 95 CI
### multivariate cox regression p values Z scores, HR, 95 CI
### working on logistic regression p values, OR, 95 CI
### working on chi squared p values, adjusted p values
### graphical output to include survival curves, relationship to subgroup


##########################################################################

### Extracted dataframes

### KM Survival overall
### see extraction functions for KM

######################
### OS p values
km.OS.all.results <- extract.km.OS.pval(results.master, name = "surv.km.OS.all")
km.OS.G3G4.results <- extract.km.OS.pval (results.master, name = "surv.km.OS.G3G4")
km.OS.SHH.results <- extract.km.OS.pval(results.master, name = "surv.km.OS.SHH") ### this is an example to put into the clinical_data_extract_DW.R file
km.OS.SHH.old.results <- extract.km.OS.pval(results.master, name = "surv.km.OS.SHH.old") 

### EFS p values
km.EFS.all.results <- extract.km.EFS.pval(results.master, name = "surv.km.EFS.all")
km.EFS.G3G4.results <- extract.km.EFS.pval(results.master, name = "surv.km.EFS.G3G4")

### PFS p values
km.PFS.all.results <- extract.km.PFS.pval(results.master, name = "surv.km.PFS.all")
km.PFS.G3G4.results <- extract.km.PFS.pval(results.master, name = "surv.km.PFS.G3G4")
km.PFS.SHH.results <- extract.km.PFS.pval(results.master, name = "surv.km.PFS.SHH")
km.PFS.SHH.old.results <- extract.km.PFS.pval(results.master, name = "surv.km.PFS.SHH.old")


### significant dataframes for adjusted p values

significant.km.EFS.all <- km.EFS.all.results[which(km.EFS.all.results[, 2]<0.05),]
significant.km.EFS.G3G4 <- km.EFS.G3G4.results [which(km.EFS.G3G4.results[, 2]<0.05),]

significant.km.OS.all <- km.OS.all.results [which (km.OS.all.results[,2] <0.05),]
significant.km.OS.G3G4 <- km.OS.G3G4.results[which(km.OS.G3G4.results [, 2]<0.05), ]
significant.km.OS.SHH <- km.OS.SHH.results [which(km.OS.SHH.results [, 2]<0.05), ]
significant.km.OS.SHH.old <- km.OS.SHH.old.results[which(km.OS.SHH.old.results[, 2]<0.05),]

significant.km.PFS.all <- km.PFS.all.results[which(km.PFS.all.results[, 2]<0.05),]              ### 4/12/17 nrow = 452 (for mb.vsd)
significant.km.PFS.G3G4 <- km.PFS.G3G4.results[which(km.PFS.G3G4.results[, 2]<0.05),]           ### 4/12/17 nrow = 379 (for mb.vsd)
significant.km.PFS.SHH <- km.PFS.SHH.results[which(km.PFS.SHH.results[, 2]<0.05),]              ### 4/12/17 nrow = 223 (for mb.vsd)
significant.km.PFS.SHH.old <- km.PFS.SHH.old.results[which(km.PFS.SHH.old.results[, 2]<0.05),]  ### 4/12/17 nrow = 570 (for mb.vsd)


##########################################################################################################################
##########################################################################################################################

### cox PFS for continuous variable, overall category
### this section has been updated below , specify results.master then subset.index

cox.PFS.cat.all.df <- extract.cox (results.master, 3)

cox.PFS.cont.all.df <- extract.cox(results.master, 4)

cox.PFS.cat.G3G4.df <- extract.cox(results.master, 5)

cox.PFS.cont.G3G4.df <- extract.cox(results.master, 6)

cox.PFS.cat.SHH.df <- extract.cox (results.master, 7)

cox.PFS.cont.SHH.df <- extract.cox (results.master, 8)

cox.PFS.cat.SHH.old.df <- extract.cox.SHH.old (results.master, 9) 
                                       
cox.PFS.cont.SHH.old.df <- extract.cox.SHH.old (results.master, 10) 


### then to filter the lists by a defined threshold e.g adjusted p<0.05

### create significant dataframes  

sig.cox.PFS.cat.all <- cox.PFS.cat.all.df [which(cox.PFS.cat.all.df[, 2]<0.05),]  ### nrow = 81

sig.cox.PFS.cont.all <- cox.PFS.cont.all.df[which(cox.PFS.cont.all.df[, 2]<0.05),] ### nrow = 698

sig.cox.PFS.cat.G3G4 <- cox.PFS.cat.G3G4.df[which(cox.PFS.cat.G3G4.df[, 2]<0.05),]

sig.cox.PFS.cont.G3G4 <- cox.PFS.cont.G3G4.df[which(cox.PFS.cont.G3G4.df[, 2]<0.05),]

sig.cox.PFS.cat.SHH <- cox.PFS.cat.SHH.df[which(cox.PFS.cat.SHH.df[, 2]<0.05),]

sig.cox.PFS.cont.SHH <- cox.PFS.cont.SHH.df[which(cox.PFS.cont.SHH.df[, 2]<0.05)]

sig.cox.PFS.cat.SHH.old <- cox.PFS.cat.SHH.old.df[which(cox.PFS.cat.SHH.old.df[, 2]<0.05),] ### no significant results 11/10/17

sig.cox.PFS.cont.SHH.old <- cox.PFS.cont.SHH.old.df[which(cox.PFS.cont.SHH.old.df[, 2]<0.05),] ### makes more sense as continuous variable 


###########################################################################################

### annotate those with ensembl gene IDs, removing those with NA

try(annot.cox.PFS.cont.all <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cont.all)), silent = T)
try(annot.sig.cox.PFS.cont.all <- cbind (annot.cox.PFS.cont.all, sig.cox.PFS.cont.all), silent = T)


try(annot.cox.PFS.cat.all <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cat.all)), silent = T)
try(annot.sig.cox.PFS.cat.all <- cbind(annot.cox.PFS.cat.all, sig.cox.PFS.cat.all), silent = T)


try(annot.cox.PFS.cont.G3G4 <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cont.G3G4)), silent = T)
annot.sig.cox.PFS.cont.G3G4 <- cbind(annot.cox.PFS.cont.G3G4, sig.cox.PFS.cont.G3G4)
# annot.sig.cox.PFS.cont.G3G4.clean <- annot.cox.PFS.cont.G3G4[complete.cases(annot.sig.cox.PFS.cont.G3G4),]

try(annot.cox.PFS.cat.G3G4 <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cat.G3G4)), silent = T)
annot.sig.cox.PFS.cat.G3G4 <- cbind(annot.cox.PFS.cat.G3G4, sig.cox.PFS.cat.G3G4)


# try(annot.cox.PFS.cont.SHH <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cont.SHH)), silent = T) ### error as does not contain data
# annot.sig.cox.PFS.cont.SHH <- cbind(annot.cox.PFS.cont.SHH, sig.cox.PFS.cont.SHH)
# try(annot.cox.PFS.cat.SHH <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cat.SHH)), silent = T)  ### error as does not contain data

try(annot.cox.PFS.cont.SHH.old <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cont.SHH.old)), silent = T)
annot.sig.cox.PFS.cont.SHH.old <- cbind (annot.cox.PFS.cont.SHH.old, sig.cox.PFS.cont.SHH.old)
try(annot.cox.PFS.cat.SHH.old <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cat.SHH.old)), silent = T)
annot.sig.cox.PFS.cat.SHH.old <- cbind(annot.cox.PFS.cat.SHH.old, sig.cox.PFS.cat.SHH.old)

### save annotated files
### unhash SHH cat and SHH contin if useful
# clean.annot.sig.cox.PFS.cont.all <- annot.sig.cox.PFS.cont.all[complete.cases(annot.sig.cox.PFS.cont.all),] ### complete.cases removes NAs

write.csv(annot.sig.cox.PFS.cont.all, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/annot.significant.cox.PFS.cont.allgroups.csv")
write.csv(annot.sig.cox.PFS.cat.all, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/annot.significant.cox.PFS.cat.allgroups.csv")

write.csv(annot.sig.cox.PFS.cont.G3G4, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/annot.sig.cox.PFS.cont.G3G4.csv")
write.csv(annot.sig.cox.PFS.cat.G3G4, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/annot.sig.cox.PFS.cat.G3G4.csv")


# write.csv(annot.sig.cox.PFS.cont.SHH, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant.cox.PFS.cont.SHH.csv")
# write.csv (annot.sig.cox.PFS.cat.SHH, file = ""/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant.cox.PFS.cat.SHH.csv")

try(write.csv (annot.sig.cox.PFS.cont.SHH.old, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/annot.sig.cox.PFS.cont.SHH.old.csv"), silent = T)
try(write.csv(annot.sig.cox.PFS.cat.SHH.old, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/annot.sig.cox.PFS.cat.SHH.old.csv"), silent = T)

########################################################################################
########################################################################################

### Cox OS overall for categorical variable

cox.OS.cat.all.df <- extract.cox.OS (results.master, 11) 

cox.OS.cont.all.df <- extract.cox.OS (results.master, 12)

cox.OS.cat.G3G4.df <- extract.cox.OS (results.master, 13)

cox.OS.cont.G3G4.df <- extract.cox.OS (results.master, 14)

cox.OS.cat.SHH.df <- extract.cox.SHH.old (results.master, 15)

cox.OS.cont.SHH.df <- extract.cox.SHH.old (results.master, 16)

cox.OS.cat.SHH.old.df <- extract.cox.SHH.old (results.master, 17)

cox.OS.cont.SHH.old.df <- extract.cox.SHH.old (results.master, 18)

##################

### significant dataframes & annotation
### can work on this more later for SHH, SHH.old if needed

sig.cox.OS.cat.all <- cox.OS.cat.all.df[which(cox.OS.cat.all.df[, 2]<0.05),]

try(annot.sig.cox.OS.cat.all <- annotate.HTseq.IDs(rownames(sig.cox.OS.cat.all)), silent = T)
# write.csv(annot.sig.cox.OS.cat.all,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/annot.sig.cox.OS.cat.allgroups.csv")

###

sig.cox.OS.cont.all <- cox.OS.cont.all.df[which(cox.OS.cont.all.df[, 2]<0.05),]
try(annot.sig.cox.OS.cont.all <- annotate.HTseq.IDs(rownames(sig.cox.OS.cont.all)), silent = T)
# write.csv(annot.sig.cox.OS.cont.all,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/annot.sig.cox.OS.cont.allgroups.csv")

###

sig.cox.OS.cat.G3G4 <- cox.OS.cat.G3G4.df[which(cox.OS.cat.G3G4.df[, 2]<0.05),]
try(annot.sig.cox.OS.cat.G3G4 <- annotate.HTseq.IDs(rownames(sig.cox.OS.cat.G3G4)), silent = T)
write.csv(annot.sig.cox.OS.cat.G3G4,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/annot.sig.cox.OS.cat.G3G4.csv")

###

sig.cox.OS.cont.G3G4 <- cox.OS.cont.G3G4.df[which(cox.OS.cont.G3G4.df[, 2]<0.05),]
try(annot.sig.cox.OS.cat.G3G4 <- annotate.HTseq.IDs(rownames(sig.cox.OS.cont.G3G4)), silent = T)
write.csv(sig.cox.OS.cont.G3G4,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/annot.sig.cox.OS.cont.G3G4.csv")

###

sig.cox.OS.cat.SHH <- cox.OS.cat.SHH.df[which(cox.OS.cat.SHH.df[, 2]<0.05),] ### no results

sig.cox.OS.cont.SHH.df <- cox.OS.cont.SHH.df[which(cox.OS.cont.SHH.df[, 2]<0.05), ]

sig.cox.OS.cat.SHH.old <- cox.OS.cat.SHH.old.df[which(cox.OS.cat.SHH.old.df[, 2]<0.05),]

sig.cox.OS.cont.SHH.old.df <- cox.OS.cont.SHH.old.df [which(cox.OS.cont.SHH.old.df[,2]<0.05),]

########################################################################

### Cox EFS for all - these are all categorical expression data


cox.EFS.cat.all.df <- extract.cox.OS (results.master, 1)

cox.EFS.cat.G3G4.df <- extract.cox.OS (results.master, 2)

# cox.EFS.cat.G3G4.v2.df <- extract.cox.SHH.old (results.master, 2) ### this is to compare the output when subsetting is increased to x[[4]]<17


sig.cox.EFS.cat.all <- cox.EFS.cat.all.df[which(cox.EFS.cat.all.df[, 2]<0.05),]

sig.cox.EFS.cat.G3G4 <- cox.EFS.cat.G3G4.df[which(cox.EFS.cat.G3G4.df[,2]<0.05), ]

# sig.cox.EFS.cat.G3G4.v2 <- cox.EFS.cat.G3G4.v2.df[which(cox.EFS.cat.G3G4.v2.df[, 2]<0.05), ] 

### explore difference from sig.cox.EFS.cat.G3G4, I wonder if it is related to the significant cox dataframe creation
### seems to be that the more accurate subsetting, the more relevant targets are found within the dataset i.e x[[4]]< 6 for a more complete dataset is better than generic x[[4]]<17
### there are fewer candidates found when have more liberal NA rule e.g for x[[4]]<17 (liberal, fewer candidates, less NAs) compared to x[[4]]<6 (strict)
### other difference is the adjusted p value (perhaps with fewer candidates, there is less correction)


try(write.csv(sig.cox.EFS.cat.all,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/sig.cox.EFS.all.csv"), silent =T)

try(write.csv(sig.cox.EFS.cat.G3G4,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/sig.cox.EFS.G3G4.csv"), silent = T)

#################################################################################################
#################################################################################################

### extract logistic regression p value 

### decision made to hard code as difficulty generating function, in hindsight could replace with "reg.log.list" like in chi square function (extract.chi.all, name = "....")

logreg.LCA.pval <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.LCA$LR.pval)})  ### x$reg.log.list$log.reg.LCA[[1]]
logreg.LCA.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.LCA$LR.OR.val)})
logreg.LCA.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.LCA$lower.95CI)})
logreg.LCA.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.LCA$upper.95CI)})
extract.logreg.LCA.df <- log.reg.dataframe(pval = logreg.LCA.pval, OR = logreg.LCA.OR, L95CI = logreg.LCA.L95CI, U95CI= logreg.LCA.U95CI)

###

logreg.relapse.pval <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.relapse$LR.pval)})
logreg.relapse.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.relapse$LR.OR.val)})
logreg.relapse.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.relapse$lower.95CI)})
logreg.relapse.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.relapse$upper.95CI)})
extract.logreg.relapse.df <- log.reg.dataframe(pval = logreg.relapse.pval, OR = logreg.relapse.OR, L95CI = logreg.relapse.L95CI, U95CI = logreg.relapse.U95CI)  

###
logreg.mstatus.pval <-lapply(results.master, function(x){return(x$reg.log.list$log.reg.mstatus$LR.pval)}) 
logreg.mstatus.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.mstatus$LR.OR.val)}) 
logreg.mstatus.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.mstatus$lower.95CI)}) 
logreg.mstatus.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.mstatus$upper.95CI)})
extract.logreg.mstatus.df <- log.reg.dataframe(pval = logreg.mstatus.pval, OR = logreg.mstatus.OR, L95CI = logreg.mstatus.L95CI, U95CI = logreg.mstatus.U95CI)

###
logreg.age.cat.pval <-lapply(results.master, function(x){return(x$reg.log.list$log.reg.age.cat$LR.pval)}) 
logreg.age.cat.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.age.cat$LR.OR.val)}) 
logreg.age.cat.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.age.cat$lower.95CI)}) 
logreg.age.cat.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.age.cat$upper.95CI)})
extract.logreg.age.cat.df <- log.reg.dataframe(pval = logreg.age.cat.pval , OR = logreg.age.cat.OR , L95CI = logreg.age.cat.L95CI, U95CI = logreg.age.cat.U95CI)

###
logreg.meth.pval <-lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth$LR.pval)}) 
logreg.meth.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth$LR.OR.val)}) 
logreg.meth.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth$lower.95CI)}) 
logreg.meth.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth$upper.95CI)})
extract.logreg.meth.df <- log.reg.dataframe(pval = logreg.meth.pval , OR = logreg.meth.OR , L95CI = logreg.meth.L95CI, U95CI = logreg.meth.U95CI)

###

logreg.meth7.pval <-lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth7$LR.pval)}) 
logreg.meth7.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth7$LR.OR.val)}) 
logreg.meth7.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth7$lower.95CI)}) 
logreg.meth7.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth7$upper.95CI)})
extract.logreg.meth7.df <- log.reg.dataframe(pval = logreg.meth7.pval , OR = logreg.meth7.OR , L95CI = logreg.meth7.L95CI, U95CI = logreg.meth7.U95CI)


###

logreg.MYC.pval <-lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYC$LR.pval)}) 
logreg.MYC.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYC$LR.OR.val)}) 
logreg.MYC.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYC$lower.95CI)}) 
logreg.MYC.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYC$upper.95CI)})
extract.logreg.MYC.df <- log.reg.dataframe(pval = logreg.MYC.pval , OR = logreg.MYC.OR , L95CI = logreg.MYC.L95CI, U95CI = logreg.MYC.U95CI)


###

logreg.MYCN.pval <-lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCN$LR.pval)}) 
logreg.MYCN.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCN$LR.OR.val)}) 
logreg.MYCN.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCN$lower.95CI)}) 
logreg.MYCN.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCN$upper.95CI)})
extract.logreg.MYCN.df <- log.reg.dataframe(pval = logreg.MYCN.pval , OR = logreg.MYCN.OR , L95CI = logreg.MYCN.L95CI, U95CI = logreg.MYCN.U95CI)

###
logreg.MYCMYCN.pval <-lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCMYCN$LR.pval)}) 
logreg.MYCMYCN.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCMYCN$LR.OR.val)}) 
logreg.MYCMYCN.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCMYCN$lower.95CI)}) 
logreg.MYCMYCN.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCMYCN$upper.95CI)})
extract.logreg.MYCMYCN.df <- log.reg.dataframe(pval = logreg.MYCMYCN.pval , OR = logreg.MYCMYCN.OR , L95CI = logreg.MYCMYCN.L95CI, U95CI = logreg.MYCMYCN.U95CI)

###

logreg.resection.pval <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.resection$LR.pval)})
logreg.resection.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.resection$LR.OR.val)})
logreg.resection.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.resection$lower.95CI)})
logreg.resection.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.resection$upper.95CI)})
extract.logreg.resection.df <- log.reg.dataframe(pval = logreg.resection.pval, OR = logreg.resection.OR, L95CI = logreg.resection.L95CI, U95CI = logreg.resection.U95CI) 

###
logreg.sex.pval <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.sex$LR.pval)})
logreg.sex.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.sex$LR.OR.val)})
logreg.sex.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.sex$lower.95CI)})
logreg.sex.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.sex$upper.95CI)})
extract.logreg.sex.df <- log.reg.dataframe(pval = logreg.sex.pval, OR = logreg.sex.OR, L95CI = logreg.sex.L95CI, U95CI = logreg.sex.U95CI) 

###

logreg.TERT.pval <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TERT$LR.pval)})
logreg.TERT.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TERT$LR.OR.val)})
logreg.TERT.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TERT$lower.95CI)})
logreg.TERT.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TERT$upper.95CI)})
extract.logreg.TERT.df <- log.reg.dataframe(pval = logreg.TERT.pval, OR = logreg.TERT.OR, L95CI = logreg.TERT.L95CI, U95CI = logreg.TERT.U95CI) 


###

logreg.TP53.pval <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TP53$LR.pval)})
logreg.TP53.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TP53$LR.OR.val)})
logreg.TP53.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TP53$lower.95CI)})
logreg.TP53.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TP53$upper.95CI)})
extract.logreg.TP53.df <- log.reg.dataframe(pval = logreg.TP53.pval, OR = logreg.TP53.OR, L95CI = logreg.TP53.L95CI, U95CI = logreg.TP53.U95CI) 


######################################################################
### if wish to make a list of all the logistic regression results 
logistic.reg.results <- as.list(mget(ls(pattern="extract.logreg"))) 

######################################################################

### creating log reg list for all significant variables
# significant.logreg.df.all <- logistic.reg.results[which(logistic.reg.results[, 2]<0.05),] ### this does not work because logistic.reg.results object is a list


#######################################################################

### extract chi square p value

chi.age.cat.infant.result <- extract.chi.all(results.master, name = "list.age.cat.infant")
chi.CSI.result <- extract.chi.all(results.master, name = "list.CSI")
chi.LCA.result <- extract.chi.all(results.master, name = "list.LCA")
chi.meth4.result <-extract.chi.all(results.master, name = "list.meth.4")
chi.meth7.result <- extract.chi.all(results.master, name = "list.meth.7")
chi.mstatus.result <- extract.chi.all(results.master, name = "list.mstatus")
chi.MYC.result <- extract.chi.all(results.master, name = "list.MYC")
chi.MYCMYCN.result <- extract.chi.all(results.master, name = "list.MYCMYCN")
chi.MYCN.result <- extract.chi.all(results.master, name = "list.MYCN")
chi.q13loss.result <- extract.chi.all(results.master, name = "list.q13loss")
chi.relapse.result <- extract.chi.all(results.master, name = "list.relapse")
chi.resection.result <- extract.chi.all(results.master, name = "list.resection")
chi.RTX.result <- extract.chi.all(results.master, name = "list.RTX")
chi.sex.result <- extract.chi.all(results.master, name = "list.sex")
chi.TERT.result <- extract.chi.all(results.master, name = "list.TERT")
chi.TP53.result <- extract.chi.all(results.master, name = "list.TP53")


### later if wish to extract the chi squared statistic, will need to rerun results.master with the updated naming for the chi squared function 4/12/17

### extract adj p <0.05 for relapse, mstatus, MYC, MYCN, MYCMYCN

significant.chi.relapse <- chi.relapse.result[which(chi.relapse.result[,2]<0.05), ]  ### n=4388 4/12/17 for mb.vsd
significant.chi.mstatus <- chi.mstatus.result[which(chi.mstatus.result[,2]<0.05), ]  ### n=3875, 4/12/17 for mb.vsd
significant.chi.MYC <- chi.MYC.result [which(chi.MYC.result[,2]<0.05), ]            ### n=4640, 4/12/17 for mb.vsd
significant.chi.MYCN <- chi.MYCN.result[which(chi.MYCN.result[,2]<0.05),]
significant.chi.MYCMYCN <- chi.MYCMYCN.result[which(chi.MYCMYCN.result[,2]<0.05), ]  ### n=214 4/12/17 for mb.vsd


########################################################################

###  multivariate cox, looking for transcripts that are significant beyond either the current PNET5, the Lancet oncology paper (Schwalbe et al 2017) or a combined model taking both models together


multivar.cox.OS.combined.cat.df <- extract.multivar.cox(results.master, 1)  ### updated so that p value is for biomarker not overall model p val 21/11/17

multivar.cox.OS.combined.cont.df <- extract.multivar.cox(results.master,2) 

multivar.cox.OS.lancetG3G4.cat.df <- extract.multivar.cox(results.master, 3) 

multivar.cox.OS.lancetG3G4.cont.df <- extract.multivar.cox(results.master, 4)

multivar.cox.OS.PNET5.cat.df <- extract.multivar.cox(results.master, 5)

multivar.cox.OS.PNET5.cont.df <- extract.multivar.cox(results.master, 6)

multivar.cox.OS.SHHold.cat.df <- extract.multivar.cox(results.master, 7)

multivar.cox.OS.SHHold.cont.df <- extract.multivar.cox(results.master, 8)

multivar.cox.PFS.combined.cat.df <- extract.multivar.cox.PFS(results.master, 9) 

multivar.cox.PFS.combined.cont.df <- extract.multivar.cox.PFS (results.master, 10)  ### subsetting worked with x[[5]]<10

multivar.cox.PFS.lancetG3G4.cat.df <- extract.multivar.cox.PFS (results.master, 11) ### had to increase subset to x[[5]]< 12

multivar.cox.PFS.lancetG3G4.cont.df <- extract.multivar.cox.PFS(results.master, 12) ### worked with x[[5]]<12

multivar.cox.PFS.PNET5.cat.df <- extract.multivar.cox.PFS (results.master, 13) ### worked with x[[5]]<13

multivar.cox.PFS.PNET5.cont.df <- extract.multivar.cox.PFS(results.master, 14) ### worked with x[[5]]<14

multivar.cox.PFS.SHHold.cat.df <- extract.multivar.cox.PFS.SHH (results.master, 15) ### worked with x[[5]]<15 

multivar.cox.PFS.SHHold.cont.df <- extract.multivar.cox.PFS.SHH (results.master, 16) ### worked with x[[5]]<16




### generating significant dataframes for the multivariate cox modelling, ie transcripts that perform about and beyond current clinical risk models

significant.multivar.cox.OS.combined.cat <- multivar.cox.OS.combined.cat.df [which(multivar.cox.OS.combined.cat.df[,2]<0.05), ]

significant.multivar.cox.OS.combined.cont <- multivar.cox.OS.combined.cont.df [which(multivar.cox.OS.combined.cont.df[,2]<0.05), ] ###n=13 21/12/17 for all transcripts

significant.multivar.cox.OS.lancetG3G4.cat <- multivar.cox.OS.lancetG3G4.cat.df[which(multivar.cox.OS.lancetG3G4.cat.df[,2]<0.05),]

significant.multivar.cox.OS.lancetG3G4.cont <- multivar.cox.OS.lancetG3G4.cont.df[which(multivar.cox.OS.lancetG3G4.cont.df[,2]<0.05),]

significant.multivar.cox.OS.PNET5.cat <- multivar.cox.OS.PNET5.cat.df[which(multivar.cox.OS.PNET5.cat.df[,2]<0.05),]  ### n=43

significant.multivar.cox.OS.PNET5.cont <- multivar.cox.OS.PNET5.cont.df[which(multivar.cox.OS.PNET5.cont.df[,2]<0.05),] ### n=30

significant.multivar.cox.OS.SHHold.cat <- multivar.cox.OS.SHHold.cat.df[which(multivar.cox.OS.SHHold.cat.df[,2]<0.05), ]

significant.multivar.cox.OS.SHHold.cont <- multivar.cox.OS.SHHold.cont.df[which(multivar.cox.OS.SHHold.cont.df[,2]<0.05),]

###need to check the PFS data Jan 2018, as the adjusted p value for the PFS dataframes are mostly =1, therefore the transcripts being extracted for the significant dataframes are those with infinite hazard ratio 95CI

significant.multivar.cox.PFS.combined.cat <- multivar.cox.PFS.combined.cat.df [which (multivar.cox.PFS.combined.cat.df[,2]<0.05), ] ###need to check this Jan 2018

significant.multivar.cox.PFS.combined.cont <- multivar.cox.PFS.combined.cont.df [which (multivar.cox.PFS.combined.cont.df[,2]<0.05), ]  ###need to check this Jan 2018

significant.multivar.cox.PFS.lancetG3G4.cat <- multivar.cox.PFS.lancetG3G4.cat.df [which (multivar.cox.PFS.lancetG3G4.cat.df[,2]<0.05), ]  ###need to check this Jan 2018

significant.multivar.cox.PFS.lancetG3G4.cont <- multivar.cox.PFS.lancetG3G4.cont.df[which (multivar.cox.PFS.lancetG3G4.cont.df[,2]<0.05),]

significant.multivar.cox.PFS.PNET5.cat <- multivar.cox.PFS.PNET5.cat.df[which(multivar.cox.PFS.PNET5.cat.df[,2]<0.05),]

significant.multivar.cox.PFS.PNET5.cont <- multivar.cox.PFS.PNET5.cont.df[which(multivar.cox.PFS.PNET5.cont.df[,2]<0.05),]

significant.multivar.cox.PFS.SHHold.cat <- multivar.cox.PFS.SHHold.cat.df[which(multivar.cox.PFS.SHHold.cat.df[,2]<0.05),]

significant.multivar.cox.PFS.SHHold.cont <- multivar.cox.PFS.SHHold.cont.df[which(multivar.cox.PFS.SHHold.cont.df[,2]<0.05),]


### annotated dataframes for significant multivar cox  ### note that the rownames (ENSG id) have been truncated so may need to be optimised 4/1/18


try(annot.sig.multi.cox.OS.combined.cat <- annotate.HTseq.IDs(rownames(significant.multivar.cox.OS.combined.cat)),silent = T)
try(annot.sig.multi.cox.OS.combined.cont <- annotate.HTseq.IDs(rownames(significant.multivar.cox.OS.combined.cont)), silent =T) ### error that given dataset hsapiens_gene_ensembl is not valid, use listDatasets() function to check
try(annot.sig.multi.cox.OS.lancetG3G4.cont <- annotate.HTseq.IDs(rownames(significant.multivar.cox.OS.lancetG3G4.cont)), silent = T) ###note that sometimes need to run the annotate separately outside of the try statement for it to work
try(annot.sig.multi.cox.OS.PNET5.cont <- annotate.HTseq.IDs(rownames(significant.multivar.cox.OS.PNET5.cont)),silent = T)
try(annot.sig.multivar.cox.OS.SHHold.cont <- annotate.HTseq.IDs(rownames(significant.multivar.cox.OS.SHHold.cont)), silent = T) ###  may not exist


# write.csv(annot.sig.multi.cox.OS.lancetG3G4.cont, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.annot.sig.multi.cox.OS.lancetG3G4.cont.csv")
          
### can also add in explicitly where to make changes if new results.master file is being utilised as the basis for this clinical_data_extract_DW.R file
############################################################################
############################################################################

### schema when integrating new subsetting into a function
### 1. determine relevant formulae e. extract.cox.OS, cox.dataframe
### 2. can use hardcoding first to determine the cut off subset required e.g x[[5]]<15, see example for extract.multivar.cox.PFS.comb.cat.pval below
### 3. then create function using x, subset.index; and keep this function within the current file to tweak until all subset indices are working, or generate new functions as required
### 4. then move new function into main "clinical_data_functions_extract_master.R" file
### 5. OR can name the index, such as $p.val rather than [[subset.index]] described as a number
  
########################################

#########################################################
#########################################################

### graphical depiction of p values against adjusted p values, may wish to add in abline
### may wish to alter graphics later

histo.p.adj.km.EFS.all <- hist(adjusted.p.km.EFS.all)

# library(density)
plot(ecdf(adjusted.p.km.EFS.all))
plot(density(adjusted.p.km.EFS.all))
hist(km.EFS.p.extract.assembled.all)
lines(density(km.EFS.p.extract.assembled.all), col = "red")

cox.PFS.cat.G3G4.df


head(cox.PFS.cat.G3G4.df)

hist(cox.PFS.cat.G3G4.df[,1])
hist(cox.PFS.cat.G3G4.df[,2])

cox.PFS.cat.G3G4.df[,1] -> x
"Cox PFS categorical G3/G4" -> test.name
breaks = 100

plotHist <- function(x, test.name, breaks = 100, xlab = "p-value", cutoff = 0.05, text.pos = 0.9){
  hist.res <- hist(x, breaks = breaks, plot = F)
  max(hist.res$counts) -> temp.height
  length(which(x<0.05)) -> temp.no.sig
  
  if(length(cutoff)==1){
  ifelse(hist.res$breaks<cutoff,"red", "grey") -> hist.cols
  }else{
    ifelse(hist.res$breaks<cutoff[1]|hist.res$breaks>cutoff[2],"red", "grey") -> hist.cols
  }
  
    hist(x, breaks = breaks,  xlab = xlab, main = paste("Histogram of", test.name), col = hist.cols)
    
    if(length(cutoff)==1){
      abline(v= cutoff, lty = 2 , col = "red")
      text(text.pos, temp.height-(temp.height*0.1), paste("Number Genes p <", cutoff, "=", temp.no.sig), pos = 2)
    }else{
      abline(v= cutoff, lty = 2 , col = "red")
      text(text.pos, temp.height-(temp.height*0.1), paste("Number Genes p <", cutoff[1],"or p >", cutoff[2],"=", temp.no.sig), pos = 2)
    }
}

cox.PFS.cat.G3G4.df[,3] -> x

plotEcdf <- function(x, y = NULL, test.name, xlab = "z-score", cutoff=c(-2,2)){
  cdf.x <- ecdf(x)
  plot(ecdf(x), xlab = xlab, main = paste("cumulative density plot of", test.name), col = "red")
  abline(h = 0.5, v = 0)
  abline(v = cutoff, lty = 2)
  min(x, na.rm = T) -> min.x
  temp.no.dn.x <- length(which(x<cutoff[1]))
  max(x, na.rm = T) -> max.x
  temp.no.up.x <- length(which(x>cutoff[2]))
  text(min.x-(0.1*min.x), 0.9, paste("Number Genes z <", cutoff[1],temp.no.dn.x), pos = 4)
  text(max.x-(0.1*max.x), 0.2, paste("Number Genes z >", cutoff[2],temp.no.up.x), pos = 2)
  if(!is.null(y)){
    cdf.y <- ecdf(y)
  }  
  
}



plotHist(cox.PFS.cat.G3G4.df[,1], "Cox PFS categorical G3/G4", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(cox.PFS.cat.G3G4.df[,2], "Cox PFS categorical G3/G4", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)
plotHist(cox.PFS.cat.G3G4.df[,3], "Cox PFS categorical G3/G4", breaks = 100, xlab = "Z-score", cutoff = c(-2, 2))



plot(ecdf(cox.PFS.cat.G3G4.df[,3]))
plot(density(x, na.rm = "T"))
hist(km.EFS.p.extract.assembled.all)
lines(density(km.EFS.p.extract.assembled.all), col = "red")



