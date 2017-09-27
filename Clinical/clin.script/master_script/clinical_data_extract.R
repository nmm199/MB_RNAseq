
### Script for analysing the output from univariate analysis of RNA expression data compared to survival outcomes
### Date: September 25 2017
### Author: Dr Marion Mateos

### example file input

results.file <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.1000.rds") 
results.file <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.5.rds") ### has cox Z score
results.master <- results.file

### need to read in functions file

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_master.R")

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

# library(density)
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

### extract cox regression p value and Z scores into individual dataframes

#################################################

### Cox PFS overall

cox.PFS.pval.all <- lapply(results.master, function(x){return(x[[4]][[3]][[1]])})
cox.PFS.Zscore.all <- lapply(results.master, function(x){return(x[[4]][[3]][[5]])})
cox.PFS.HR.all <- lapply(results.master, function(x){return(x[[4]][[3]][[2]])})
cox.L95CI.PFS.HR.all <- lapply(results.master, function(x){return(x[[4]][[3]][[3]])})
cox.U95CI.PFS.HR.all <- lapply(results.master, function(x){return(x[[4]][[3]][[4]])})

cox.PFS.all.df <- cox.dataframe(pval = cox.PFS.pval.all, Zscore = cox.PFS.Zscore.all, HR = cox.PFS.HR.all, L95CI = cox.L95CI.PFS.HR.all, U95CI = cox.U95CI.PFS.HR.all )

colnames(cox.PFS.all.df) <- c("cox.PFS.pval.all", "cox.PFS.adj.pval.all", "cox.PFS.Zscore.all", "cox.PFS.HR.all", "cox.PFS.L95CI.all", "cox.PFS.U95CI.all")

significant.cox.PFS.all <- cox.PFS.all.df[which(cox.PFS.all.df[, 2]<0.05),]

write.csv(significant.cox.PFS.all, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant_cox_PFS_all.csv")


### Cox PFS for G3G4

cox.PFS.pval.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[1]])})
cox.PFS.Zscore.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[5]])})
cox.PFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[2]])})
cox.U95CI.PFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[4]])})
cox.L95CI.PFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[4]][[3]])})

cox.PFS.G3G4.df <- cox.dataframe (pval = cox.PFS.pval.G3G4, Zscore = cox.PFS.Zscore.G3G4, HR = cox.PFS.HR.G3G4, L95CI = cox.L95CI.PFS.HR.G3G4 , U95CI = cox.U95CI.PFS.HR.G3G4)
colnames(cox.PFS.G3G4.df) <- c("cox.PFS.pval.G3G4", "cox.PFS.adj.pval.G3G4", "cox.PFS.Zscore.G3G4", "cox.PFS.HR.G3G4", "cox.PFS.L95CI.G3G4", "cox.PFS.U95CI.G3G4")

significant.cox.PFS.G3G4 <- cox.PFS.G3G4.df[which(cox.PFS.G3G4.df[, 2]<0.05),]

write.csv(significant.cox.PFS.G3G4, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant_cox_PFS_G3G4.csv")


### Cox OS overall

cox.OS.pval.all <- lapply (results.master, function(x){return(x[[4]][[5]][[1]])})
cox.OS.Zscore.all <- lapply(results.master, function(x){return(x[[4]][[5]][[5]])})
cox.OS.HR.all <- lapply(results.master, function(x){return(x[[4]][[5]][[2]])})
cox.U95CI.OS.HR.all <- lapply(results.master, function(x){return(x[[4]][[5]][[4]])})
cox.L95CI.OS.HR.all <- lapply(results.master, function(x){return(x[[4]][[5]][[3]])})

cox.OS.all.df <- cox.dataframe(pval = cox.OS.pval.all, Zscore = cox.OS.Zscore.all, HR = cox.OS.HR.all, L95CI = cox.L95CI.OS.HR.all, U95CI = cox.U95CI.OS.HR.all)
colnames(cox.OS.all.df) <- c("cox.OS.pval.all", "cox.OS.adj.pval.all", "cox.OS.Zscore.all", "cox.OS.HR.all","cox.L95CI.OS.HR.all", "cox.U95CI.OS.HR.all")
  
significant.cox.OS.all <- cox.OS.all.df[which(cox.OS.all.df[, 2]<0.05),]

write.csv(significant.cox.OS.all,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant.cox.OS.all.csv")


### Cox OS for G3G4

cox.OS.pval.G3G4 <- lapply (results.master, function(x){return(x[[4]][[6]][[1]])})
cox.OS.Zscore.G3G4 <- lapply(results.master, function(x){return(x[[4]][[6]][[5]])})
cox.OS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[6]][[2]])})
cox.U95CI.OS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[6]][[4]])})
cox.L95CI.OS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[6]][[3]])})

cox.OS.G3G4.df <- cox.dataframe(pval = cox.OS.pval.G3G4, Zscore = cox.OS.Zscore.G3G4, HR = cox.OS.HR.G3G4, L95CI = cox.L95CI.OS.HR.G3G4, U95CI = cox.U95CI.OS.HR.G3G4)
colnames(cox.OS.G3G4.df) <- c("cox.OS.pval.G3G4", "cox.OS.adj.pval.G3G4", "cox.OS.Zscore.G3G4", "cox.OS.HR.G3G4","cox.L95CI.OS.HR.G3G4", "cox.U95CI.OS.HR.G3G4")

significant.cox.OS.G3G4 <- cox.OS.G3G4.df[which(cox.OS.G3G4.df[, 2]<0.05),]
write.csv(significant.cox.OS.G3G4,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant.cox.OS.G3G4.csv")


### Cox EFS for all

cox.EFS.pval.all <- lapply(results.master, function(x){return(x[[4]][[1]][[1]])})
cox.EFS.Zscore.all <- lapply(results.master, function(x){return(x[[4]][[1]][[5]])})
cox.EFS.HR.all <- lapply(results.master, function(x){return(x[[4]][[1]][[2]])})
cox.U95CI.EFS.HR.all <- lapply(results.master, function(x){return(x[[4]][[1]][[4]])})  
cox.L95CI.EFS.HR.all <- lapply(results.master, function(x){return(x[[4]][[1]][[3]])})

cox.EFS.all.df <- cox.dataframe(pval = cox.EFS.pval.all, Zscore = cox.EFS.Zscore.all, HR = cox.EFS.HR.all, L95CI = cox.L95CI.EFS.HR.all, U95CI = cox.U95CI.EFS.HR.all)
colnames(cox.EFS.all.df) <- c("cox.EFS.pval.all", "cox.EFS.adj.pval.all", "cox.EFS.Zscore.all", "cox.EFS.HR.all", "cox.U95CI.EFS.HR.all", "cox.L95CI.EFS.HR.all")

significant.cox.EFS.all <- cox.EFS.all.df[which(cox.EFS.all.df[, 2]<0.05),]
write.csv(significant.cox.EFS.all,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant.cox.EFS.all.csv")

### cox EFS for G3G4

cox.EFS.pval.G3G4  <- lapply(results.master, function(x){return(x[[4]][[2]][[1]])})
cox.EFS.Zscore.G3G4 <- lapply(results.master, function(x){return(x[[4]][[2]][[5]])})
cox.EFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[2]][[2]])})
cox.U95CI.EFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[2]][[4]])})  
cox.L95CI.EFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[2]][[3]])})

cox.EFS.G3G4.df <- cox.dataframe(pval = cox.EFS.pval.G3G4, Zscore = cox.EFS.Zscore.G3G4, HR = cox.EFS.HR.G3G4, L95CI = cox.L95CI.EFS.HR.G3G4, U95CI = cox.U95CI.EFS.HR.G3G4)
colnames(cox.EFS.G3G4.df) <- c("cox.EFS.pval.G3G4", "cox.EFS.adj.pval.G3G4", "cox.EFS.Zscore.G3G4", "cox.EFS.HR.G3G4", "cox.U95CI.EFS.HR.G3G4", "cox.L95CI.EFS.HR.G3G4")

significant.cox.EFS.G3G4 <- cox.EFS.G3G4.df[which(cox.EFS.G3G4.df[, 2]<0.05),]
write.csv(significant.cox.EFS.G3G4,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant.cox.EFS.G3G4.csv")

########################################################################

### extract logistic regression p value 
### create function like for cox regression dataframe (clinical_data_functions_master.R)

### example input for logistic regression relapse

log.reg.dataframe <- function(pval, OR, L95CI, U95CI){
logreg.pval.assembled <- do.call(rbind, pval)
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







