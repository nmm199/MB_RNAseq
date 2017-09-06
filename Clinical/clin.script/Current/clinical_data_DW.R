
### File introduction
### File name: clinical_data_4.R

### Aim of file is to 
# 1. Run basic descriptive statistics on a cohort of children treated for medulloblastoma, whose details are contained within the local clinical database
# 2. Analyse genes of interest in relation to univariate and multivariate risk prediction models for survival (overall, event-free and progression-free)
# 3. This script covers analysis up to and including univariate risk factor analysis
# 4. Multivariate analysis / AUC will be covered by a separate script


### Author: Dr Marion Mateos
### Date: July 3 2017


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

library(NMF)
library(gplots)
library(car)
library(stats)
library(survival)

### Functions used

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/Current/clinical_data_functions_DW.R")

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
cat ("reading in expression data", sep ="\n")
# RNA.data <- "/home/dan/mygit/rna_seq_mb/paper/MB.vsd.txt"
RNA.data <- "/home/dan/mygit/rna_seq_mb/paper/vsd.novel.txt"  ### updated060917
#mb.vsd <- read.delim(file="/home/dan/mygit/rna_seq_mb/paper/vsd.novel.txt")
mb.vsd <- read.delim(RNA.data)

####  providing a putative biomarker
goi <- "ENSG00000136997"
goi.vsd <- as.numeric(mb.vsd[goi,]) 
### the output would be a vector with a continuous variable, names equal NMB numbers
names(goi.vsd) <- gsub("T","",names(mb.vsd))
#####################################################################################
### update your pData object

x.data <- "/home/nmm199/R/MB_RNAseq/Input_data/database270617_060917.csv" ### changed as previous file database270617_310817 had PFS_R with non-numerical values
cat ("reading in clinical database", sep ="\n")
### add in row names to original database file
pData <- read.csv(x.data, row.names = 1)

meth.data <- "/home/nmm199/R/MB_RNAseq/Input_data/all7subgroupCalls.csv"
meth7 <- read.csv(meth.data, header=TRUE, sep=",", quote="\"", dec=".", row.names=1)

cytogen.data <- "/home/nmm199/R/MB_RNAseq/Input_data/arm_calls_clean280617.txt"
cytogen <- read.table (cytogen.data, header=T, sep="\t")

test.pData <- updatepData(pData, meth7, cytogen, pdf.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/temp.pdf", log.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/temp.log.txt")
save(test.pData, file = "/home/nmm199/R/MB_RNAseq/Clinical/test.pData") 

log.file = "pDatalog.txt"

################################################################################

pdf.file <- "marker.results.pdf"

results <- clinPathAssess(test.pData,goi.vsd,pdf.file = "marker.results.pdf",log.file = "marker.results.txt")


library(foreach)
library(tictoc)
library(parallel)
library(doParallel)
registerDoParallel(10)

tic()
out <- foreach(i = 1:10)%dopar%{
  as.numeric(mb.vsd[i,]) -> x
  names(x) <- colnames(mb.vsd)
  return(clinPathAssess(test.pData,x))
}
toc()

### some script from earlier that may help, to look into 060917

#head(mb.vsd)
#apply(mb.vsd,1,function(x){clinPathAssess(test.pData,x)}) -> results.multiple

#library(foreach)
#library(doParallel)
#registerDoParallel(10)
#out <- foreach(i = 1:10)%dopar%{
#  x <- as.numeric(mb.vsd[i,])
#  names(x) <- colnames(mb.vsd)
#  clinPathAssess(test.pData, x )
#}

#library(foreach)
#library(doParallel)
#registerDoParallel(10)
#out <- foreach(i = 1:nrow(mb.vsd))%dopar%{
#  x <- as.numeric(mb.vsd[i,])
#  names(x) <- colnames(mb.vsd)
#  clinPathAssess(test.pData, x )
#}
