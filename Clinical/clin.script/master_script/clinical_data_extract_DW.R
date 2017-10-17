
### Script for analysing the output from univariate analysis of RNA expression data compared to survival outcomes
### Date: September 25 2017
### Author: Dr Marion Mateos

### file input

results.master <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.master.allgenes.rds") 
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
### this section is to be updated as of Oct 11 2017

###
function(results.master){

### OS p value for G3G4
  cbind(
extracted.km.OS.pval = extractKM(results.master, 3),
extracted.km.OS.pval.G3G4 = extractKM(results.master, 4)
)

######

### put all the functions you want to run here


###
extracted.results.list <-list(
  extracted.km.OS.pval=extracted.km.OS.pval,
  extracted.km.OS.pval.G3G4=extracted.km.OS.pval.G3G4
)
return(extracted.results.list)
}





### OS p value for SHH

extracted.km.OS.pval.SHH <- lapply(results.master, function(x){return(x[[1]][[5]][[1]])}) 
km.OS.p.extract.assembled.SHH <- do.call(rbind, extracted.km.OS.pval.SHH)
adjusted.p.km.OS.SHH <-p.adjust(km.OS.p.extract.assembled.SHH, method = "BH") 
OS.pvalue.SHH <- cbind(km.OS.p.extract.assembled.SHH, adjusted.p.km.OS.SHH)
colnames(OS.pvalue.SHH) <- c("OS.p.value.SHH", "OS.adjusted.pval.SHH")

### OS p value for SHH.old
extracted.km.OS.pval.SHH.old <- lapply(results.master, function(x){return(x[[1]][[6]][[1]])}) 
km.OS.p.extract.assembled.SHH.old <- do.call(rbind, extracted.km.OS.pval.SHH.old)
adjusted.p.km.OS.SHH.old <-p.adjust(km.OS.p.extract.assembled.SHH.old, method = "BH") 
OS.pvalue.SHH.old <- cbind(km.OS.p.extract.assembled.SHH.old, adjusted.p.km.OS.SHH.old)
colnames(OS.pvalue.SHH.old) <- c("OS.p.value.SHH.old", "OS.adjusted.pval.SHH.old")

### EFS p value for overall group

extracted.km.EFS.pval.all <- lapply(results.master, function(x){return(x[[1]][[1]][[1]])})  ### note that pulls out same as x[[1]][[1]][[1]] but able to combine dataframes and do p adjust on this element
km.EFS.p.extract.assembled.all <- do.call(rbind, extracted.km.EFS.pval.all)
adjusted.p.km.EFS.all <-p.adjust(km.EFS.p.extract.assembled.all, method = "BH") 
EFS.pvalue.all.combined <- cbind(km.EFS.p.extract.assembled.all, adjusted.p.km.EFS.all)
colnames(EFS.pvalue.all.combined) <- c("EFS.p.value.all", "EFS.adjusted.pval.all")

### EFS p value for G3G4

extracted.km.EFS.pval.G3G4 <- lapply(results.master, function(x){return(x[[1]][[2]][[1]])}) 
km.EFS.p.extract.assembled.G3G4 <- do.call(rbind, extracted.km.EFS.pval.G3G4)
adjusted.p.km.EFS.G3G4 <-p.adjust(km.EFS.p.extract.assembled.G3G4, method = "BH") 
EFS.pvalue.G3G4.combined <- cbind(km.EFS.p.extract.assembled.G3G4, adjusted.p.km.EFS.G3G4)
colnames(EFS.pvalue.G3G4.combined) <- c("EFS.p.value.G3G4", "EFS.adjusted.pval.G3G4")

### PFS p value for overall 
extracted.km.PFS.pval.all <- lapply(results.master, function(x){return(x[[1]][[7]][[1]])}) ### changed now that SHH and SHH.old analysis have been added
km.PFS.p.extract.assembled.all <- do.call(rbind, extracted.km.PFS.pval.all)
adjusted.p.km.PFS.all <-p.adjust(km.PFS.p.extract.assembled.all, method = "BH") 
PFS.pvalue.all.combined <- cbind(km.PFS.p.extract.assembled.all, adjusted.p.km.PFS.all)
colnames(PFS.pvalue.all.combined) <- c("PFS.p.value.all", "PFS.adjusted.pval.all")

### PFS p value for G3G4
extracted.km.PFS.pval.G3G4 <- lapply(results.master, function(x){return(x[[1]][[8]][[1]])})  ### changed as above
km.PFS.p.extract.assembled.G3G4 <- do.call(rbind, extracted.km.PFS.pval.G3G4)
adjusted.p.km.PFS.G3G4 <-p.adjust(km.PFS.p.extract.assembled.G3G4, method = "BH") 
PFS.pvalue.G3G4.combined <- cbind(km.PFS.p.extract.assembled.G3G4, adjusted.p.km.PFS.G3G4)
colnames(PFS.pvalue.G3G4.combined) <- c("PFS.p.value.G3G4", "PFS.adjusted.pval.G3G4")

### PFS p value for SHH
extracted.km.PFS.p.val.SHH <- lapply(results.master, function(x){return(x[[1]][[9]][[1]])})
km.PFS.p.extract.assembled.SHH <- do.call(rbind, extracted.km.PFS.p.val.SHH)
adjusted.p.km.PFS.SHH <- p.adjust(km.PFS.p.extract.assembled.SHH, method = "BH")
PFS.pvalue.SHH.combined <- cbind(km.PFS.p.extract.assembled.SHH, adjusted.p.km.PFS.SHH)
colnames(PFS.pvalue.SHH.combined)<- c("PFS.p.value.SHH", "PFS.adjusted.pval.SHH")

### PFS p value for SHH_child
extracted.km.PFS.p.val.SHH.old <- lapply(results.master, function(x){return(x[[1]][[10]][[1]])})
km.PFS.p.extract.assembled.SHH.old <- do.call(rbind, extracted.km.PFS.p.val.SHH.old)
adjusted.p.km.PFS.SHH.old <- p.adjust(km.PFS.p.extract.assembled.SHH.old, method = "BH")
PFS.pvalue.SHH.combined.old <- cbind(km.PFS.p.extract.assembled.SHH.old, adjusted.p.km.PFS.SHH.old)
colnames(PFS.pvalue.SHH.combined.old)<- c("PFS.p.value.SHH.old", "PFS.adjusted.pval.SHH.old")


### combined dataframe with OS, PFS, EFS results for overall and G3G4, unadjusted and adjusted p values

OS.pvalues.bothgroups <- cbind(OS.pvalue.all, OS.pvalue.G3G4,  OS.pvalue.SHH, OS.pvalue.SHH.old)
EFS.pvalues.bothgroups <- cbind(EFS.pvalue.all.combined, EFS.pvalue.G3G4.combined)
PFS.pvalues.bothgroups <- cbind(PFS.pvalue.all.combined, PFS.pvalue.G3G4.combined, PFS.pvalue.SHH.combined, PFS.pvalue.SHH.combined.old )

all.survival.p.bothgroups <- cbind(OS.pvalues.bothgroups, EFS.pvalues.bothgroups, PFS.pvalues.bothgroups)

### extract those goi with p<0.05 in adjusted p values for survival

significant.p.EFS.all <- EFS.pvalue.all.combined[which(EFS.pvalue.all.combined[, 2]<0.05),]
significant.p.EFS.G3G4 <- EFS.pvalue.G3G4.combined[which(EFS.pvalue.all.combined[, 2]<0.05),]
significant.p.OS.all <- OS.pvalue.all[which(OS.pvalue.all[, 2]<0.05), ]
significant.p.OS.G3G4 <- OS.pvalue.G3G4[which(OS.pvalue.G3G4 [, 2]<0.05), ]
significant.p.OS.SHH <- OS.pvalue.SHH[which(OS.pvalue.SHH[, 2]<0.05), ]
significant.p.OS.SHH.old <- OS.pvalue.SHH.old[which(OS.pvalue.SHH.old[, 2]<0.05),]
significant.p.PFS.all <- PFS.pvalue.all.combined[which(PFS.pvalue.all.combined[, 2]<0.05),]
significant.p.PFS.G3G4 <- PFS.pvalue.G3G4.combined[which(PFS.pvalue.G3G4.combined[, 2]<0.05),]
significant.p.PFS.SHH <- PFS.pvalue.SHH.combined[which(PFS.pvalue.SHH.combined[, 2]<0.05),]
significant.p.PFS.SHH.old <- PFS.pvalue.SHH.combined.old[which(PFS.pvalue.SHH.combined.old[, 2]<0.05),]



### need to update up until here 11/10/17

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


### DW suggests to create list of outputs, rather than significant dataframes
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

write.csv (annot.sig.cox.PFS.cont.SHH.old, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/annot.sig.cox.PFS.cont.SHH.old.csv")
write.csv(annot.sig.cox.PFS.cat.SHH.old, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/annot.sig.cox.PFS.cat.SHH.old.csv")

########################################################################################
########################################################################################

### Cox OS overall for categorical variable
### Then work on EFS then work on km survival


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
### can work on this more later

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



########################################################################
### checking dimensions in each dataframe

# nrow(sig.cox.OS.cat.all)
# nrow(sig.cox.OS.cont.all)
# nrow (sig.cox.OS.cat.G3G4)
# nrow(sig.cox.OS.cont.G3G4)


#########################################################################

### Cox EFS for all - these are all categorical expression data

try(cox.EFS.pval.all <- lapply(results.master, function(x){return(x[[4]][[1]][[1]])}), silent = T)
try (cox.EFS.Zscore.all <- lapply(results.master, function(x){return(x[[4]][[1]][[5]])}), silent = T)
try(cox.EFS.HR.all <- lapply(results.master, function(x){return(x[[4]][[1]][[2]])}), silent = T)
try(cox.U95CI.EFS.HR.all <- lapply(results.master, function(x){return(x[[4]][[1]][[4]])}), silent = T)  
try(cox.L95CI.EFS.HR.all <- lapply(results.master, function(x){return(x[[4]][[1]][[3]])}), silent = T)

try(cox.EFS.all.df <- cox.dataframe(pval = cox.EFS.pval.all, Zscore = cox.EFS.Zscore.all, HR = cox.EFS.HR.all, L95CI = cox.L95CI.EFS.HR.all, U95CI = cox.U95CI.EFS.HR.all), silent = T)
try(colnames(cox.EFS.all.df) <- c("cox.EFS.pval.all", "cox.EFS.adj.pval.all", "cox.EFS.Zscore.all", "cox.EFS.HR.all", "cox.U95CI.EFS.HR.all", "cox.L95CI.EFS.HR.all"), silent = T)

try(significant.cox.EFS.all <- cox.EFS.all.df[which(cox.EFS.all.df[, 2]<0.05),], silent = T)
try(write.csv(significant.cox.EFS.all,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant.cox.EFS.all.csv"), silent =T)

### cox EFS for G3G4

try(cox.EFS.pval.G3G4  <- lapply(results.master, function(x){return(x[[4]][[2]][[1]])}), silent = T)
try(cox.EFS.Zscore.G3G4 <- lapply(results.master, function(x){return(x[[4]][[2]][[5]])}), silent = T)
try(cox.EFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[2]][[2]])}), silent = T)
try(cox.U95CI.EFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[2]][[4]])}), silent = T)  
try(cox.L95CI.EFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[2]][[3]])}), silent = T)

try(cox.EFS.G3G4.df <- cox.dataframe(pval = cox.EFS.pval.G3G4, Zscore = cox.EFS.Zscore.G3G4, HR = cox.EFS.HR.G3G4, L95CI = cox.L95CI.EFS.HR.G3G4, U95CI = cox.U95CI.EFS.HR.G3G4), silent = T)
try(colnames(cox.EFS.G3G4.df) <- c("cox.EFS.pval.G3G4", "cox.EFS.adj.pval.G3G4", "cox.EFS.Zscore.G3G4", "cox.EFS.HR.G3G4", "cox.U95CI.EFS.HR.G3G4", "cox.L95CI.EFS.HR.G3G4"), silent = T)

try(significant.cox.EFS.G3G4 <- cox.EFS.G3G4.df[which(cox.EFS.G3G4.df[, 2]<0.05),], silent = T)
try(write.csv(significant.cox.EFS.G3G4,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant.cox.EFS.G3G4.csv"), silent = T)

########################################################################

#################################################################################################

### extract logistic regression p value 
### create function like for cox regression dataframe (clinical_data_functions_master.R)

### example input for logistic regression relapse

log.reg.dataframe <- function(pval, OR, L95CI, U95CI){
logreg.pval.assembled <- do.call(rbind, pval) ### check if recurrent error here when more than one input goi
logreg.adj.pval <- p.adjust(logreg.pval.assembled, method = "BH")
logreg.HR.assembled <- do.call(rbind, HR)
logreg.L95CI.assembled <- do.call(rbind, L95CI)
logreg.U95CI.assembled <- do.call(rbind, U95CI)
logreg.allresults.df <- cbind(logreg.pval.assembled, logreg.adj.pval, logreg.HR.assembled, logreg.HR.assembled, logreg.L95CI.assembled, logreg.U95CI.assembled)
colnames(logreg.allresults.df)<- c("logreg.pval", "logreg.adj.pval","logreg.HR", "logreg.HR.L95CI", "logreg.HR.U95CI" )
return (logreg.allresults.df)
}


### potentially just return p values and adjusted p values for each logistic regression category currently

logreg.LCA.list <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.LCA)})
logreg.LCA.pval <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.LCA[[1]])})
logreg.LCA.df <- do.call(rbind, logreg.LCA.pval) ### possibly do not need this
logreg.LCA.adj.pval <- p.adjust(logreg.LCA.df, method = "BH")
logreg.LCA.combined.pval<- cbind (logreg.LCA.pval, logreg.LCA.df, logreg.LCA.adj.pval)



logreg.relapse.list <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.relapse)})
logreg.relapse.pval <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.relapse[[1]])})
logreg.relapse.adj.pval <- p.adjust(logreg.relapse.pval, method = "BH")
logreg.relapse.combined.pval <- cbind(logreg.relapse.pval, logreg.relapse.adj.pval)


### could create a foreach loop

#logreg.allresults <- foreach (i = 1:length(results.master[[i]]$reg.log.list)$dopar% {
 # logreg.pval <- lapply(results.master, function(x){return(x$reg.log.list[[i]][[1]])})
 # logreg.adj.pval <- p.adjust(logreg.pval, method = "BH")
 # logreg.combined.pval <- cbind (logreg.pval, logreg.adj.pval)
 # return(logreg.combined.pval)
## }
# )

######################################################################

### creating log reg list for all variables

#significant.logreg.df.all <- logreg.allresults.df[which(cox.allresults.df[, 2]<0.05)]


#######################################################################

### extract chi square p value


########################################################################

### examples for how to then extract lists you are interested in

#lapply(results.master, function(x){return(x[[3]][[2]])}) -> extracted.results

#do.call(rbind, extracted.results) -> compiled.results

#p.adjust(compiled.results[,1], method = "BH") -> adjusted.p.values

#hist(adjusted.p.values)



#############################################################################

############################################################################
### function to return files:
extracted.dataframes <- list(cox.PFS.cat.all.df, 
                             cox.PFS.cont.all.df,
                             cox.PFS.cat.G3G4.df,
                             cox.PFS.cont.G3G4.df,
                             cox.PFS.cat.SHH.df,
                             cox.PFS.cont.SHH.df,
                             cox.PFS.cat.SHH.old.df,
                             cox.PFS.cont.SHH.old.df,
                             cox.OS.cat.all.df,
                             cox.OS.cont.all.df,
                             cox.OS.cat.G3G4.df,
                             cox.OS.cont.G3G4.df,
                             cox.OS.cat.SHH.df,
                             cox.OS.cont.SHH.df,
                             cox.OS.cat.SHH.old.df,
                             cox.OS.cont.SHH.old.df
                             )
# return (extracted.dataframes)



### script that worked previously for the novel transcript file that has now been superceded 28/9/17


#########
### examples based on previous dataframes

# cox.PFS.pval.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[1]])})
# cox.PFS.Zscore.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[5]])})
# cox.PFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[2]])})
# cox.U95CI.PFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[4]])})
# cox.L95CI.PFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[3]])})

#########
### older script for COX OS all

# cox.OS.pval.all <- lapply(results.master, extract.coxpval.OS.all)
# cox.OS.pval.all <- lapply (results.master, function(x){return(x[[4]][[5]][[1]])})
# cox.OS.Zscore.all <- lapply(results.master, function(x){return(x[[4]][[5]][[5]])})
# cox.OS.HR.all <- lapply(results.master, function(x){return(x[[4]][[5]][[2]])})
# cox.U95CI.OS.HR.all <- lapply(results.master, function(x){return(x[[4]][[5]][[4]])})
# cox.L95CI.OS.HR.all <- lapply(results.master, function(x){return(x[[4]][[5]][[3]])})

# cox.OS.all.df <- cox.dataframe(pval = cox.OS.pval.all, Zscore = cox.OS.Zscore.all, HR = cox.OS.HR.all, L95CI = cox.L95CI.OS.HR.all, U95CI = cox.U95CI.OS.HR.all)
# colnames(cox.OS.all.df) <- c("cox.OS.pval.all", "cox.OS.adj.pval.all", "cox.OS.Zscore.all", "cox.OS.HR.all","cox.L95CI.OS.HR.all", "cox.U95CI.OS.HR.all")

# significant.cox.OS.all <- cox.OS.all.df[which(cox.OS.all.df[, 2]<0.05),]

# write.csv(significant.cox.OS.all,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant.cox.OS.all.csv")






###########################

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
