
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
### clinPathAssess

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

### unhash to use the following file for all the expression data 
RNA.data <- "/home/dan/mygit/rna_seq_mb/paper/MB.vsd.txt" ### run it first on this 070917, then update to the novel vsd and the foreach loop if it is working on the single goi

mb.vsd.novel <- read.delim(file="/home/dan/mygit/rna_seq_mb/paper/vsd.novel.txt") ### novel transcripts, updated 060917

mb.vsd <- read.delim(RNA.data)

##############################################################################
source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_master.R")

### loading in clinical data object = test.pData
load("/home/nmm199/R/MB_RNAseq/Clinical/test.pData")

### set file for pdf output
pdf.file <- "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/marker.results.pdf"
#pdf.file <- "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/marker.results.novel.pdf"

### set file for log output
log.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/pDatalog.txt"
#log.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/pDatalog.novel.txt"

################################################################################

### Run of one: need to define goi.vsd in order to use the ClinPathAsess function
### unhash this section (lines 97 - 115) when running one goi

# goi.vsd <- as.numeric(mb.vsd[1,]) ### can choose a specific row, or can specify a goi within inverted commas

# goi <- "ENSG00000008196"                
# goi.vsd <- as.numeric(mb.vsd[goi, ])    
# names(mb.vsd) -> names(goi.vsd)        

## results for run of one
# results <- clinPathAssess(test.pData,
                         # goi.vsd,
                         # pdf.file = pdf.file,
                         # log.file = log.file)


# names(results)<- row.names(goi.vsd)

# saveRDS (results, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.ENSG00000008196.rds")
# readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.ENSG00000008196.rds")

###############################################################################

### To run multiple genes in the clinPathAssess function 

library(foreach)
library(tictoc)
library(parallel)
library(doParallel)
registerDoParallel(10)

tic()

### unhash here when trouble shooting
# i=17
# as.numeric(mb.vsd.novel[i,]) -> x
# names(x) <- colnames(mb.vsd.novel)
# names(x) <- gsub("T","",names(mb.vsd.novel))
# clinPathAssess(test.pData,x)      ### unhash here when trouble shooting


### unhash when running the complete transcript set

# results.master <- foreach(i = 1:nrow(mb.vsd))%dopar%{
# as.numeric(mb.vsd[i,]) -> x
# names(x) <- colnames(mb.vsd)
# names(x) <- gsub("T","",names(mb.vsd))
# return(clinPathAssess(test.pData,x)) 
# }

### unhash when running the novel transcript set

# results.master <- foreach(i = 1:nrow(mb.vsd.novel))%dopar%{
# as.numeric(mb.vsd.novel[i,]) -> x
# names(x) <- colnames(mb.vsd.novel)
# names(x) <- gsub("T","",names(mb.vsd.novel)) ### check that this is correct
# return(clinPathAssess(test.pData,x)) 
# }


### script for [1:10] ie isolated set of transcripts. have changed names(goi.vsd) to names(x), goi.vsd is specified as "x" in script below:
results.master <- foreach(i = 1:10)%dopar%{
  as.numeric(mb.vsd[i,]) -> x
  names(x) <- colnames(mb.vsd)
  names(x) <- gsub("T","",names(mb.vsd)) 
  return(clinPathAssess(test.pData,x))
}


### then unhash the relevant outputs

# names(results.master) <- row.names(mb.vsd)
# names(results.master) <- row.names(mb.vsd.novel)
# names(results.master) <- row.names(mb.vsd)[1:nrow(mb.vsd)]
  names(results.master) <- row.names(mb.vsd)[1:10]

###########################################################
  
  ### if ongoing errors try this:
  
  # tryCatch{ #### This catches the error and outputs it to screen but allows the program to continue running
  # results.master <- foreach(i = 1:nrow(mb.vsd))%dopar%{
  #  as.numeric(mb.vsd [i,]) -> x
  # names(x) <- colnames(mb.vsd)
  #  names(x) <- gsub("T","",names(mb.vsd)) ### check that this is correct
  #  error=function(e){cat("ERROR :",conditionMessage(e), "\n")} 
  #  return(clinPathAssess(test.pData,x)) 
  # }
  # }  #### prints the error message to screen
  
  
############################################################
  
### Annotate with known gene sets

 annot.results <- annotate.HTseq.IDs(row.names(mb.vsd))

# annot.novel <- annotate.HTseq.IDs(row.names(mb.vsd.novel)) 

toc()

##################################################################
### save outputted results

### update name according to input file

saveRDS (results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.10.rds")
#saveRDS (results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.rds")
# saveRDS (results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.rds")
# saveRDS (results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.novel.rds")

### then reload this when examining the results

#results.file <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.1000.rds") ### generated before cox Z score extracted 
# results.file <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.5.rds") ### has cox Z score
# results.master <- results.file

############################################

### examples for how to then extract lists you are interested in: see separate script "clinical_data_extract.R" 

#lapply(results.master, function(x){return(x[[3]][[2]])}) -> extracted.results
#do.call(rbind, extracted.results) -> compiled.results
#p.adjust(compiled.results[,1], method = "BH") -> adjusted.p.values
#hist(adjusted.p.values)

############################################



