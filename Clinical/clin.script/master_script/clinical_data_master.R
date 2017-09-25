
### File introduction
### File name: clinical_data_master.R  based on clinical_data_110917.R, clinical_data_4.R with excerpts from clinical_data_7.R and clinPathAssess_backup.R

### Aim of file is to 
# 1. Run basic descriptive statistics on a cohort of children treated for medulloblastoma, whose details are contained within the local clinical database
# 2. Analyse genes of interest in relation to univariate and multivariate risk prediction models for survival (overall, event-free and progression-free)
# 3. This script covers analysis up to and including univariate risk factor analysis
# 4. Multivariate analysis / AUC will be covered by a separate script


### Author: Dr Marion Mateos
### Date: Sep 11 2017


### R version 3.4.0 (2017-04-21)
### Platform: x86_64-pc-linux-gnu (64-bit)
### Running under: Ubuntu 16.04.2 LTS


# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

### Packages and version numbers
## [1] survival_2.41-3     
## RColorBrewer_1.1-2  
## car_2.1-4           
## gplots_3.0.1        
## NMF_0.20.6          
## Biobase_2.36.2      
## BiocGenerics_0.22.0
## cluster_2.0.6       
## rngtools_1.2.4      
## pkgmaker_0.22       
## registry_0.3       


### Libraries to be used
# install.packages('gplots')
# install.packages('survival')


### Functions used

### names of functions for info on function see source file
### "chi.sq"
### "cor.result"
### "lin.reg"
### "km.log.test"
### "km.log.test.OS"
### "cox.result.OS"
### "km.log.test.EFS"
### updatepData

### External files required

### clinical database "x.data"
### 7 molecular group data "meth.data"
### cytogenetic arm data "cytogen.data"
### RNA expression data "RNA.data"


###############################################################################
### deals with making a GOI.vsd vector this is the part that is going to change
### at the moment this is just any old expression data
### you will plug in your own goi - isoforms, novel genes, etc
### just for demonstration purposes at the moment
###############################################################################

# cat ("reading in expression data", sep ="\n")
### unhash to use the following file for all the expression data 
# RNA.data <- "/home/dan/mygit/rna_seq_mb/paper/MB.vsd.txt" ### run it first on this 070917, then update to the novel vsd and the foreach loop if it is working on the single goi

### unhash to use the following file for novel RNA transcripts
RNA.data <- "/home/dan/mygit/rna_seq_mb/paper/vsd.novel.txt"  ### updated060917

# mb.vsd <- read.delim(file="/home/dan/mygit/rna_seq_mb/paper/vsd.novel.txt") ### or use this as the whole command in one line

mb.vsd <- read.delim(RNA.data)


##############################################################################
source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_master.R")


### loading in clinical data object = test.pData
load("/home/nmm199/R/MB_RNAseq/Clinical/test.pData")

### set file for pdf output
#pdf.file <- "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/marker.results.pdf"
pdf.file <- "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/marker.results.novel.pdf"

### set file for log output
#log.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/pDatalog.txt"
log.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/pDatalog.novel.txt"

### need to define goi.vsd in order to use the ClinPathAsess function
goi.vsd <- as.numeric(mb.vsd[1,])
names(mb.vsd) -> names(goi.vsd)

## run of one
results <- clinPathAssess(test.pData,
                          goi.vsd,
                          pdf.file = pdf.file,
                          log.file = log.file)


############################
### to run multiple genes
library(foreach)
library(tictoc)
library(parallel)
library(doParallel)
registerDoParallel(10)

tic()

results.master <- foreach(i = 1:nrow(mb.vsd))%dopar%{
#i=1515
as.numeric(mb.vsd[i,]) -> x
names(x) <- colnames(mb.vsd)
names(goi.vsd) <- gsub("T","",names(mb.vsd))
#clinPathAssess(test.pData,x)
return(clinPathAssess(test.pData,x))
}

#results.master <- foreach(i = 1:nrow(mb.vsd))%dopar%{
  ###i=25
  #as.numeric(mb.vsd[i,]) -> x
  #names(x) <- colnames(mb.vsd)
  #names(goi.vsd) <- gsub("T","",names(mb.vsd))
  ###clinPathAssess(test.pData,x)
  #return(clinPathAssess(test.pData,x))
#}
toc()
names(results.master) <- row.names(mb.vsd)[1:nrow(mb.vsd)]
#names(results.master) <- row.names(mb.vsd)[1:5]

##################################################################
### save outputted results
### update name according to input file

#saveRDS (results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.5.rds")
#saveRDS (results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.rds")
saveRDS (results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.novel.rds")

### then reload this when examining the results

#results.file <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.1000.rds") ### generated before cox Z score extracted 
results.file <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.5.rds") ### has cox Z score
results.master <- results.file


############################################

### examples for how to then extract lists you are interested in

#lapply(results.master, function(x){return(x[[3]][[2]])}) -> extracted.results

#do.call(rbind, extracted.results) -> compiled.results

#p.adjust(compiled.results[,1], method = "BH") -> adjusted.p.values

#hist(adjusted.p.values)

############################################

### delete the section below once the "clinical_data_extract.R" file is finalised. 

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

### 

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



### extract logistic regression p value 

### extract chi square p value

