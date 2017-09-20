
### File introduction
### File name: clinical_data_110917.R  based on clinical_data_4.R with excerpts from clinical_data_7.R and clinPathAssess_backup.R

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
RNA.data <- "/home/dan/mygit/rna_seq_mb/paper/MB.vsd.txt" ### run it first on this 070917, then update to the novel vsd and the foreach loop if it is working on the single goi
# RNA.data <- "/home/dan/mygit/rna_seq_mb/paper/vsd.novel.txt"  ### updated060917
#mb.vsd <- read.delim(file="/home/dan/mygit/rna_seq_mb/paper/vsd.novel.txt")
mb.vsd <- read.delim(RNA.data)


##############################################################################
source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_master.R")


### loading in clinical data object = test.pData
load("/home/nmm199/R/MB_RNAseq/Clinical/test.pData")

### set file for pdf output
pdf.file <- "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/marker.results.pdf"

### set file for log output
log.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/pDatalog.txt"

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
results.master <- foreach(i = 1:5)%dopar%{
  #i=25
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
# names(results.master) <- row.names(mb.vsd)[1:nrow(mb.vsd)]
names(results.master) <- row.names(mb.vsd)[1:5]

### save outputted results
### update name according to input file
### curren naming is for first 1000 goi for the mb.vsd file

# saveRDS (results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.1000.rds")
 
### then reload this when examining the results
# results.file <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.1000.rds")


############################################

### examples for how to then extract lists you are interested in

lapply(results.master, function(x){return(x[[3]][[2]])}) -> extracted.results

do.call(rbind, extracted.results) -> compiled.results

p.adjust(compiled.results[,1], method = "BH") -> adjusted.p.values

hist(adjusted.p.values)

############################################

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


### do the same for PFS for overall and G3G4




### combined dataframe with OS, PFS, EFS results for overall and G3G4, unadjusted and adjusted p values

EFS.pvalues.bothgroups <- cbind(EFS.pvalue.all.combined, EFS.pvalue.G3G4.combined)
OS.pvalues.bothgroups <- cbind(OS.pvalue.all, OS.pvalue.G3G4)


### graphical depiction of p values against adjusted p values, may wish to add in abline

histo.p.adj.km.EFS.all <- hist(adjusted.p.km.EFS)
histo.p.adj.km.EFS.G3G4 <- hist(adjusted.p.km.EFS.G3G4)
histo.p.adj.km.OS.all <- hist(adjusted.p.km.OS)
histo.p.adj.km.OS.G3G4 <- hist(adjusted.p.km.OS.G3G4)

### work out how to depict both on same graph
### extracted.km.EFS.pval

#################################################


