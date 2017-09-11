
### File introduction
### File name: clinical_data_070917.R  based on clinical_data_4.R with excerpts from clinical_data_7.R and clinPathAssess_backup.R

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

cat ("reading in expression data", sep ="\n")
RNA.data <- "/home/dan/mygit/rna_seq_mb/paper/MB.vsd.txt" ### run it first on this 070917, then update to the novel vsd and the foreach loop if it is working on the single goi
# RNA.data <- "/home/dan/mygit/rna_seq_mb/paper/vsd.novel.txt"  ### updated060917
#mb.vsd <- read.delim(file="/home/dan/mygit/rna_seq_mb/paper/vsd.novel.txt")
mb.vsd <- read.delim(RNA.data)


##############################################################################
source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/Current/clinical_data_functions_110917.R")


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
out <- foreach(i = 1:1000)%dopar%{
 as.numeric(mb.vsd[i,]) -> x
 names(x) <- colnames(mb.vsd)
names(goi.vsd) <- gsub("T","",names(mb.vsd)) 
return(clinPathAssess(test.pData,x))
}
toc()


lapply(out, function(x){return(x[[3]][[2]])}) -> extracted.results

do.call(rbind, extracted.results) -> compiled.results

p.adjust(compiled.results[,1], method = "BH") -> adjusted.p.values

hist(adjusted.p.values)


