
### Script for analysing the output from univariate analysis of RNA expression data compared to survival outcomes
### Date: September 25 2017
### Author: Dr Marion Mateos

### example file input

results.file <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.5.rds") ### has cox Z score
results.master <- results.file

### file output 

### dataframes with p value, adjusted p value

### all.survival.p.bothgroups: includes OS, EFS and PFS for overall group and G3G4 

### significant adjusted p values for OS, EFS, PFS separately in dataframes for overall group and for G3G4

### graphical output for adjusted p values compared to p values

### cox regression p values, Z scores, HR, 95 CI

### working on logistic regression p values, OR, 95 CI

### working on chi squared p values, adjusted p values

### graphical output to include survival curves, relationship to subgroup


##########################################################################

### Extracted dataframes
### OS p value for overall group

extracted.km.OS.pval <- lapply(results.master, function(x){return(x[[1]][[3]][[1]])}) 
km.OS.p.extract.assembled <- do.call(rbind, extracted.km.OS.pval)
adjusted.p.km.OS <-p.adjust(km.OS.p.extract.assembled, method = "BH") 
OS.pvalue.all <- cbind(km.OS.p.extract.assembled, adjusted.p.km.OS)
colnames(OS.pvalue.all) <- c("OS.p.value.all", "OS.adjusted.pval.all")

### OS p value for G3G4

extracted.km.OS.pval.G3G4 <- lapply(results.master, function(x){return(x[[1]][[4]][[1]])}) 
km.OS.p.extract.assembled.G3G4 <- do.call(rbind, extracted.km.OS.pval.G3G4)
adjusted.p.km.OS.G3G4 <-p.adjust(km.OS.p.extract.assembled.G3G4, method = "BH") 
OS.pvalue.G3G4 <- cbind(km.OS.p.extract.assembled.G3G4, adjusted.p.km.OS.G3G4)
colnames(OS.pvalue.G3G4) <- c("OS.p.value.G3G4", "OS.adjusted.pval.G3G4")


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
extracted.km.PFS.pval.all <- lapply(results.master, function(x){return(x[[1]][[5]][[1]])})
km.PFS.p.extract.assembled.all <- do.call(rbind, extracted.km.PFS.pval.all)
adjusted.p.km.PFS.all <-p.adjust(km.PFS.p.extract.assembled.all, method = "BH") 
PFS.pvalue.all.combined <- cbind(km.PFS.p.extract.assembled.all, adjusted.p.km.PFS.all)
colnames(PFS.pvalue.all.combined) <- c("PFS.p.value.all", "PFS.adjusted.pval.all")

### PFS p value for G3G4
extracted.km.PFS.pval.G3G4 <- lapply(results.master, function(x){return(x[[1]][[6]][[1]])})
km.PFS.p.extract.assembled.G3G4 <- do.call(rbind, extracted.km.PFS.pval.G3G4)
adjusted.p.km.PFS.G3G4 <-p.adjust(km.PFS.p.extract.assembled.G3G4, method = "BH") 
PFS.pvalue.G3G4.combined <- cbind(km.PFS.p.extract.assembled.G3G4, adjusted.p.km.PFS.G3G4)
colnames(PFS.pvalue.G3G4.combined) <- c("PFS.p.value.G3G4", "PFS.adjusted.pval.G3G4")

### combined dataframe with OS, PFS, EFS results for overall and G3G4, unadjusted and adjusted p values

OS.pvalues.bothgroups <- cbind(OS.pvalue.all, OS.pvalue.G3G4)
EFS.pvalues.bothgroups <- cbind(EFS.pvalue.all.combined, EFS.pvalue.G3G4.combined)
PFS.pvalues.bothgroups <- cbind(PFS.pvalue.all.combined, PFS.pvalue.G3G4.combined)

all.survival.p.bothgroups <- cbind(OS.pvalues.bothgroups, EFS.pvalues.bothgroups, PFS.pvalues.bothgroups)

### extract those goi with p<0.05 in adjusted p values for survival

significant.p.EFS.all <- EFS.pvalue.all.combined[which(EFS.pvalue.all.combined[, 2]<0.05),]
significant.p.EFS.G3G4 <- EFS.pvalue.G3G4.combined[which(EFS.pvalue.all.combined[, 2]<0.05),]
significant.p.OS.all <- OS.pvalue.all[which(OS.pvalue.all[, 2]<0.05), ]
significant.p.OS.G3G4 <- OS.pvalue.G3G4[which(OS.pvalue.G3G4 [, 2]<0.05), ]
significant.p.PFS.all <- PFS.pvalue.all.combined[which(PFS.pvalue.all.combined[, 2]<0.05),]
significant.p.PFS.G3G4 <- PFS.pvalue.G3G4.combined[which(PFS.pvalue.G3G4.combined[, 2]<0.05),]

### graphical depiction of p values against adjusted p values, may wish to add in abline

histo.p.adj.km.EFS.all <- hist(adjusted.p.km.EFS.all)

library(density)
plot(ecdf(adjusted.p.km.EFS.all))
plot(density(adjusted.p.km.EFS.all))
hist(km.EFS.p.extract.assembled.all)
lines(density(km.EFS.p.extract.assembled.all), col = "red")

histo.p.adj.km.EFS.G3G4 <- hist(adjusted.p.km.EFS.G3G4)
histo.p.adj.km.OS.all <- hist(adjusted.p.km.OS)
histo.p.adj.km.OS.G3G4 <- hist(adjusted.p.km.OS.G3G4)
histo.p.adj.km.PFS.all <- hist(adjusted.p.km.PFS.all)
histo.p.adj.km.PFS.G3G4 <- hist(adjusted.p.km.PFS.G3G4) ### all same value, therefore no histogram generated

#################################################

### extract cox regression p value and Z score

### generating EFS dataframe

extracted.cox.EFS.pval <- lapply(results.master, function(x){return(x[[4]][[1]][[1]])})
cox.EFS.pval.assembled.all <- do.call(rbind, extracted.cox.EFS.pval)
adjusted.p.cox.EFS.all <-p.adjust(cox.EFS.pval.assembled.all, method = "BH") 
EFS.cox.p.all.combined <- cbind(cox.EFS.pval.assembled.all, adjusted.p.cox.EFS.all)
colnames(EFS.cox.p.all.combined) <- c("EFS.cox.pval.all", "EFS.cox.adj.pval.all")

### extract cox Z score
extracted.cox.EFS.Zscore.all <- lapply(results.master, function(x){return(x[[4]][[1]][[5]])})
cox.EFS.extract.Zscore.assembled.all <- do.call(rbind, extracted.cox.EFS.Zscore.all)
colnames(cox.EFS.extract.Zscore.assembled.all)<- c("EFS.cox.Zscore.all")

### extract Hazard ratio cox
extracted.EFS.cox.HR.all <- lapply(results.master, function(x){return(x[[4]][[1]][[2]])})
cox.EFS.HR.assembled.all <- do.call(rbind,extracted.EFS.cox.HR.all )
colnames(cox.EFS.HR.assembled.all)<- c("EFS.cox.HR.all")

### extract 95 Confidence intervals
### upper 95 CI

U95CI.EFS.cox.HR.all <- lapply(results.master, function(x){return(x[[4]][[1]][[4]])})
cox.EFS.HR.U95CI.assembled.all <- do.call(rbind, U95CI.EFS.cox.HR.all)
colnames(cox.EFS.HR.U95CI.assembled.all)<-c("EFS.cox.HR.U95CI")

### lower 95 CI
L95CI.EFS.cox.HR.all <- lapply(results.master, function(x){return(x[[4]][[1]][[3]])})
cox.EFS.HR.L95CI.assembled.all <- do.call(rbind, L95CI.EFS.cox.HR.all)
colnames(cox.EFS.HR.L95CI.assembled.all)<-c("EFS.cox.HR.L95CI")

### construct cox regression dataframe for all parameters of interest

cox.EFS.HR.95CI.all.df <- cbind(cox.EFS.HR.assembled.all,cox.EFS.HR.L95CI.assembled.all, cox.EFS.HR.U95CI.assembled.all)
cox.EFS.all.df <- cbind(EFS.cox.p.all.combined, cox.EFS.extract.Zscore.assembled.all, cox.EFS.HR.95CI.all.df)


############################
### replicate cox dataframe for OS, PFS then EFS G3G4

### PFS for overall

extracted.cox.PFS.pval <- lapply(results.master, function(x){return(x[[4]][[3]][[1]])})
cox.PFS.pval.assembled.all <- do.call(rbind, extracted.cox.PFS.pval)
adjusted.p.cox.PFS.all <-p.adjust(cox.PFS.pval.assembled.all, method = "BH") 
PFS.cox.p.all.combined <- cbind(cox.PFS.pval.assembled.all, adjusted.p.cox.PFS.all)
colnames(PFS.cox.p.all.combined) <- c("PFS.cox.pval.all", "PFS.cox.adj.pval.all")

### extract cox Z score
extracted.cox.PFS.Zscore.all <- lapply(results.master, function(x){return(x[[4]][[3]][[5]])})
cox.PFS.extract.Zscore.assembled.all <- do.call(rbind, extracted.cox.PFS.Zscore.all)
colnames(cox.PFS.extract.Zscore.assembled.all)<- c("PFS.cox.Zscore.all")

### extract Hazard ratio cox
extracted.PFS.cox.HR.all <- lapply(results.master, function(x){return(x[[4]][[3]][[2]])})
cox.PFS.HR.assembled.all <- do.call(rbind,extracted.PFS.cox.HR.all )
colnames(cox.PFS.HR.assembled.all)<- c("PFS.cox.HR.all")

### extract 95 Confidence intervals
### upper 95 CI

U95CI.PFS.cox.HR.all <- lapply(results.master, function(x){return(x[[4]][[3]][[4]])})
cox.PFS.HR.U95CI.assembled.all <- do.call(rbind, U95CI.PFS.cox.HR.all)
colnames(cox.PFS.HR.U95CI.assembled.all)<-c("PFS.cox.HR.U95CI")

### lower 95 CI
L95CI.PFS.cox.HR.all <- lapply(results.master, function(x){return(x[[4]][[3]][[3]])})
cox.PFS.HR.L95CI.assembled.all <- do.call(rbind, L95CI.PFS.cox.HR.all)
colnames(cox.PFS.HR.L95CI.assembled.all)<-c("PFS.cox.HR.L95CI")

### construct cox regression dataframe for all parameters of interest

cox.PFS.HR.95CI.all.df <- cbind(cox.PFS.HR.assembled.all,cox.PFS.HR.L95CI.assembled.all, cox.PFS.HR.U95CI.assembled.all)
cox.PFS.all.df <- cbind(PFS.cox.p.all.combined, cox.PFS.extract.Zscore.assembled.all, cox.PFS.HR.95CI.all.df)

### extract those goi with p<0.05 in adjusted p values for survival

significant.cox.PFS.all <- cox.PFS.all.df[which(cox.PFS.all.df[, 2]<0.05),]

################## 

### PFS for G3G4

extracted.cox.PFS.pval.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[1]])})
cox.PFS.pval.assembled.G3G4 <- do.call(rbind, extracted.cox.PFS.pval.G3G4)
adjusted.p.cox.PFS.G3G4 <-p.adjust(cox.PFS.pval.assembled.G3G4, method = "BH") 
PFS.cox.p.G3G4.combined <- cbind(cox.PFS.pval.assembled.G3G4, adjusted.p.cox.PFS.G3G4)
colnames(PFS.cox.p.G3G4.combined) <- c("PFS.cox.pval.G3G4", "PFS.cox.adj.pval.G3G4")


### extract cox Z score
extracted.cox.PFS.Zscore.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[5]])})
cox.PFS.extract.Zscore.assembled.G3G4 <- do.call(rbind, extracted.cox.PFS.Zscore.G3G4)
colnames(cox.PFS.extract.Zscore.assembled.G3G4)<- c("PFS.cox.Zscore.G3G4")

### extract Hazard ratio cox
extracted.PFS.cox.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[2]])})
cox.PFS.HR.assembled.G3G4 <- do.call(rbind,extracted.PFS.cox.HR.G3G4 )
colnames(cox.PFS.HR.assembled.G3G4)<- c("PFS.cox.HR.G3G4")

### extract 95 Confidence intervals
### upper 95 CI

U95CI.PFS.cox.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[4]])})
cox.PFS.HR.U95CI.assembled.G3G4 <- do.call(rbind, U95CI.PFS.cox.HR.G3G4)
colnames(cox.PFS.HR.U95CI.assembled.G3G4)<-c("PFS.cox.HR.U95CI.G3G4")

### lower 95 CI
L95CI.PFS.cox.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[3]])})
cox.PFS.HR.L95CI.assembled.G3G4 <- do.call(rbind, L95CI.PFS.cox.HR.G3G4)
colnames(cox.PFS.HR.L95CI.assembled.G3G4)<-c("PFS.cox.HR.L95CI.G3G4")

### construct cox regression dataframe for all parameters of interest

cox.PFS.HR.95CI.G3G4.df <- cbind(cox.PFS.HR.assembled.G3G4,cox.PFS.HR.L95CI.assembled.G3G4, cox.PFS.HR.U95CI.assembled.G3G4)
cox.PFS.G3G4.df <- cbind(PFS.cox.p.G3G4.combined, cox.PFS.extract.Zscore.assembled.G3G4, cox.PFS.HR.95CI.G3G4.df)

### extract those goi with p<0.05 in adjusted p values for survival

significant.cox.PFS.G3G4 <- cox.PFS.all.df[which(cox.PFS.G3G4.df[, 2]<0.05),]

#####################################################################
### to replicate script above for EFS G3G4, OS overall and OS G3G4



########################################################################

### extract logistic regression p value 

#######################################################################

### extract chi square p value


########################################################################

### examples for how to then extract lists you are interested in

#lapply(results.master, function(x){return(x[[3]][[2]])}) -> extracted.results

#do.call(rbind, extracted.results) -> compiled.results

#p.adjust(compiled.results[,1], method = "BH") -> adjusted.p.values

#hist(adjusted.p.values)

