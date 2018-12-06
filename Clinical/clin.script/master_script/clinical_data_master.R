
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
library(NMF)

### Functions used

### names of functions for info on function see source file
### "chi.sq"
### "cor.result"
### "lin.reg"
### "km.log.test"
### "km.log.test.OS"
### "cox.result.OS"
### "km.log.test.EFS"
### "cox.result.surv"
### "updatepData"
### "clinPathAssess"
### "cox.dataframe"
### "annotate.HTseq.IDs"

### input files required

### clinical database "x.data"
### 7 molecular group data "meth.data"
### cytogenetic arm data "cytogen.data"
### RNA expression data "RNA.data"
### test.pData file (this has been generated via pData_input_master.R.   If the clinical database file is updated then need to rerun pData_input_master.R)

### output files
### a defined RDS object from which other dataframes can be extracted using clinical_data_extract.R

########################################################################################################################################################
########################################################################################################################################################

### deals with making a GOI.vsd vector this is the part that is going to change
### at the moment this is just any old expression data
### you will plug in your own goi - isoforms, novel genes, etc

###############################################################################

### unhash to use the following file for all the expression data 

RNA.data <- "/home/dan/mygit/rna_seq_mb/paper/MB.vsd.txt" ### run it first on this 070917, then update to the novel vsd and the foreach loop if it is working on the single goi

mb.vsd.novel <- read.delim(file="/home/dan/mygit/rna_seq_mb/paper/vsd.novel.txt") ### novel transcripts, updated 060917

mb.vsd <- read.delim(RNA.data)

mb.vsd.random <- randomize(mb.vsd) ### generate this first then run the clinPathAssess function on this.
rownames(mb.vsd.random) <- rownames(mb.vsd)


##############################################################################

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_master.R")

### loading in clinical data object = test.pData
load("/home/nmm199/R/MB_RNAseq/Clinical/test.pData")

### set file for pdf output
pdf.file <- "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/marker.results.pdf" ### this is the name of the file generated when run clinPathAssess
# pdf.file <- "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/marker.results.novel.pdf"

### set file for log output
log.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/pDatalog.txt"
# log.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/pDatalog.novel.txt"

### run up to here if wish to interrogate clinPathAssess function, by directly nominating input variables within clinical_data_functions_master.R
################################################################################

### If running on gene filtered file, use gp.filt.mb.vsd section 

###############################################################################

### To run multiple genes in the clinPathAssess function 

library(foreach)
library(tictoc)
library(parallel)
library(doParallel)
registerDoParallel(16) ### best to do in multiples of 8 or 16 as divides easier over the NICR compsvr (32/64 cores)

tic()

###############################################################################
################################################################################
### creating filtered files 30/1/18
### filtering out patients
filt.mb.vsd <- mb.vsd[-grep ("T", names(mb.vsd ))] 
filt.mb.vsd.random <- mb.vsd.random [, -grep ("T", colnames(mb.vsd.random))] ### or names (mb.vsd.random) 27/2/18
# rownames(filt.mb.vsd.random)<- rownames(mb.vsd)
filt.mb.vsd.novel <- mb.vsd.novel[-grep("T", names(mb.vsd.novel))]

### then can insert filt.mb.vsd directly into clinPathAssess function below

### further filtering of samples based on preset filters of variance stabilised reads
### filters include: change in expression levels, proportion (remove top and bottom x percent), base expression i.e lowest expression level, proportion of samples that need to express the minimum base level)

gp.index <- apply(2^filt.mb.vsd,1,gp.style.filter, fold.change = 3, delta = 300, prop = 0.05, base = 30, prop.base = 0.05)
gp.index.novel <- apply(2^filt.mb.vsd.novel,1,gp.style.filter, fold.change = 3, delta = 300, prop = 0.05, base = 30, prop.base = 0.05)
# gp.index.random <- apply(2^filt.mb.vsd.random,1,gp.style.filter, fold.change = 3, delta = 300, prop = 0.05, base = 30, prop.base = 0.05)
### the mb vsd expression set is transformed from log(2) to exponential, to allow characterisation of delta change (in absolute terms)
### PROP = 0.05, removes top 5% of extreme data before filter 
### base = minimum expression (based on reads)
### PROP.BASE has to be expressed in at least 5% of patients, have room to increase these



### create a dataframe based on these transcripts that meet filter criteria, which will be passed through to the clinPathAssess function below

gp.filt.mb.vsd <- filt.mb.vsd[gp.index,]
gp.filt.mb.vsd.novel <- filt.mb.vsd.novel[gp.index.novel, ]

# gp.filt.mb.vsd.random <- filt.mb.vsd.random[gp.index.random, ]
gp.filt.mb.vsd.random <- randomize(gp.filt.mb.vsd) ### previous gp.filt.mb.vsd.random did not filter from the larger file, therefore randomise on different file
rownames(gp.filt.mb.vsd.random) <- rownames(gp.filt.mb.vsd)

############################################################################################################################################
### there are some additional options for filtering, added 6/12/18
### rowvsd function
n <- 3000
apply(vsd.matrix, 1, sd) -> sd.genes
names(head(sd.genes[order(sd.genes, decreasing = T)], n)) -> most.variable
vsd.matrix[most.variable, ] -> filt.vsd.matrix

library(genefilter)

cvfun <- cv(a=2, b=Inf, na.rm=TRUE)
cvfun(vsd.matrix) -> index.cv
vsd.matrix[index.cv,] -> filt.vsd.matrix

pOverAfun <- pOverA(p=0.05, A=100, na.rm=TRUE)
pOverAfun(vsd.matrix) -> index.pOverA
vsd.matrix[index.pOverA,] -> filt.vsd.matrix

### or 

vsd.matrix[index.pOverA & index.cv,] -> filt.vsd.matrix

#### remove genes not independently prognostic???
coxfilter()


######################################################################################################################################################################
######################################################################################################################################################################
### can run specific input transcript files, then choose the relevant input and output below, with results.master name and destination, annotated file for known genes
######################################################################################################################################################################
### for run of n=1
### define goi

# goi <- "ENSG00000166971" ### AKTIP
# goi <- "ENSG00000162461"  ### SLC25A34
# goi <- "ENSG00000184271" ### POU6F1
# goi <- "ENSG00000161664" ### ASB16
# goi <- "ENSG00000260440" ### LINC01544
# goi <- "ENSG00000119042" ### SATB2
# goi <- "ENSG00000128626" ### MRPS12
# goi <- "ENSG00000067836" ### ROGDI
# goi <- "ENSG00000103150" ### MLYCD  
# goi <- "ENSG00000124588"   ### NQO2
# goi <- "ENSG00000168772" ### CXXC4
# goi <- "ENSG00000173818" ### ENDOV
# goi <- "ENSG00000165304" ### MELK
# goi <-  "ENSG00000178980" ### SELENOW
# goi <- "ENSG00000245322"  ### LOC256880
# goi <- "ENSG00000266872"

### create dataframe
# goi.vsd <- as.numeric (gp.filt.mb.vsd[goi,])
# names(goi.vsd) <- names(gp.filt.mb.vsd)

### to determine distribution of the biomarker using qqnorm (qqplot) to assess value as a biomarker
# qqnorm(goi.vsd)

### unhash when running results for run of one
# results.master <- clinPathAssess(test.pData,
#  goi.vsd,
#  pdf.file = pdf.file,
#  log.file = log.file)

# names(results.master)<- row.names(goi.vsd)

# saveRDS (results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/results.master.transcript.filt.mb.vsd.rds")
### change name of transcript output
### saveRDS (results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/results.master.ENSG00000178980.gp.filt.mb.vsd.rds")
 
### read results back in
 
# readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/results.master.transcript.filt.mb.vsd.rds")

### Graphical output: marker.results.pdf graphs generated when do run of 1 in clinPathAssess
### OR see clinical_data_functions_master.R line 703 onwards for instructions on generating single KM curves within ClinPathAssess function
######################################################################################################################################################################

### unhash when running the filtered transcript set (filt.mb.vsd) (removed duplicates from mb.vsd with "NMBXXXT")

# results.master <- foreach(i = 1:nrow(filt.mb.vsd))%dopar%{
  ### results.master <- foreach(i = 1:100)%dopar%{
  # as.numeric(filt.mb.vsd[i,]) -> x
  # names(x) <- colnames(filt.mb.vsd)
  # return(clinPathAssess(test.pData,x)) 
# }

################################################################################
### unhash when running the filtered transcript set that removes sample duplicates AND pre-filters based on expression features of the transcripts
### USE THIS SCRIPT FOR MAIN ANALYSIS 25/10/18

# results.master <- foreach(i = 1:nrow(gp.filt.mb.vsd))%dopar%{
 ### results.master <- foreach(i = 1:100)%dopar%{
# as.numeric(gp.filt.mb.vsd[i,]) -> x
# names(x) <- colnames(gp.filt.mb.vsd)
# return(clinPathAssess(test.pData,x)) 
# }


# results.master <- foreach(i = 1:10)%dopar%{
  ### results.master <- foreach(i = 1:100)%dopar%{
#  as.numeric(gp.filt.mb.vsd[i,]) -> x
#  names(x) <- colnames(gp.filt.mb.vsd)
#  return(clinPathAssess(test.pData,x)) 
# }

################################################################################
### unhash when running the novel transcript analysis
### interchange gp.filt.mb.vsd.novel with filt.mb.vsd.novel

  results.master <- foreach(i = 1:nrow(gp.filt.mb.vsd.novel))%dopar%{
  ### results.master <- foreach(i = 1:5)%dopar%{
  as.numeric(gp.filt.mb.vsd.novel[i,]) -> x
  names(x) <- colnames(gp.filt.mb.vsd.novel)
   return(clinPathAssess(test.pData,x)) 
  }

### unhash when running the novel transcript set

# results.master <- foreach(i = 1:nrow(gp.filt.mb.vsd.novel))%dopar%{
  ### results.master <- foreach(i = 1:10)%dopar%{
 # as.numeric(mb.vsd.novel[i,]) -> x
 # names(x) <- colnames(mb.vsd.novel)
 # names(x) <- gsub("T","",names(mb.vsd.novel)) ### check that this is correct
 # return(clinPathAssess(test.pData,x)) 
# }

### compare to output using above 1st command for novel transcripts

################################################################################
### unhash when running the randomised dataset  ### can interchange with gp.filt.mb.vsd.random 

# results.master <- foreach(i = 1:nrow(gp.filt.mb.vsd.random))%dopar%{
#  as.numeric(gp.filt.mb.vsd.random [i,]) -> x
#  names(x) <- colnames(gp.filt.mb.vsd.random)
#  return(clinPathAssess(test.pData,x)) 
# }

################################################################################
################################################################################
### Unhash this section below when running one goi within clinPathAssess function on complete transcripts, however has been replaced by using filt.mb.vsd files above. 

# goi <- "ENSG00000124588"   ### NQO2
# goi <- "ENSG00000173818.16"

### goi.vsd <- as.numeric(mb.vsd[1,]) ### can choose a specific row, OR can specify a goi within inverted commas 
# goi.vsd <- as.numeric(mb.vsd[goi, ])   

# names(mb.vsd) -> names(goi.vsd)        

## results for run of one
# results.master <- clinPathAssess(test.pData,
#  goi.vsd,
#  pdf.file = pdf.file,
#  log.file = log.file)

# names(results.master)<- row.names(goi.vsd)

# saveRDS (results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/results.master.ENSG00000124588.rds")

# readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/results.master.ENSG00000124588.rds")


######################################################################################################################################################################
### unhash the relevant name for the output 
######################################################################################################################################################################
### note files: unfiltered complete transcripts (mb.vsd), randomised complete (mb.vsd.random), complete novel (mb.vsd.novel)
### sample filtered: will need to use these now with sample filtration (filt.mb.vsd)  
### gene and sample filtered (gp.filt.mb.vsd)

# names (results.master) <- row.names (filt.mb.vsd)
# names (results.master)<- row.names(gp.filt.mb.vsd) 
 names (results.master)<- row.names(gp.filt.mb.vsd.novel) ### interchange (gp.)filt.mb.vsd.novel
# names (results.master) <- row.names(gp.filt.mb.vsd.random) ### interchange (gp.)filt.mb.vsd.random

### superceded
# names(results.master) <- row.names(mb.vsd)
# names (results.master) <- row.names(mb.vsd.random)
# names(results.master) <- row.names(mb.vsd.novel)
# names(results.master) <- row.names(mb.vsd)[1:10]

toc()

######################################################################################################################################################################
### run up to here for the clinPathAssess function 
######################################################################################################################################################################################################################################################
### save RDS

### update name according to input file

# saveRDS (results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.allgenes.10.051017.rds")
# saveRDS (results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.master.allgenes.20180104.rds")
# saveRDS(results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.master.allgenes.random.20180104.rds")
# saveRDS (results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.master.allgenes.novel.20180104.rds")

# saveRDS(results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.filt.complete.20180220.rds") ### this is the filtered file for samples, contains > 60000 transcripts (filt.mb.vsd)
# saveRDS(results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.filt.genefilter.20181031.rds") ### this is the filtered file for both samples and genes, ~9000 transcripts (gp.filt.mb.vsd)
 saveRDS(results.master, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.filt.genefilter.novel.20181031.rds") ### interchange with complete.novel and genefilter.novel
# saveRDS (results.master, file =  "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.filt.genefilter.random.20180327.rds") ### randomised file, interchange with genefilter.random

### to examine results, reload relevant results master file and see clinical_data_extract_DW.R script file

######################################################################################################################################################################
### run up to here for the clinPathAssess function 
######################################################################################################################################################################

### Annotate with known gene sets for known transcripts (not novel transcripts)

# annot.results <- annotate.HTseq.IDs(row.names(mb.vsd))
# annot.novel <- annotate.HTseq.IDs(row.names(mb.vsd.novel)) ### not useful unless known transcript

 annot.filt.results <- annotate.HTseq.IDs(row.names(gp.filt.mb.vsd)) ### interchange filt.mb.vsd and  gp.filt.mb.vsd
# annot.filt.novel <- annotate.HTseq.IDs(row.names(gp.filt.mb.vsd.novel)) ### interchange filt.mb.vsd.novel and gp.filt.mb.vsd.novel
# annot.filt.random <- annotate.HTseq.IDs(row.names(gp.filt.mb.vsd.random))

# write.csv(annot.filt.results, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.annot.filt.complete.20180220.csv") ### filt.mb.vsd this is the filtered file for samples only, > 60000 transcripts
 write.csv(annot.filt.results, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.filt.genefilter.20181031.csv" ) ### this is the filtered file for both samples and genes, ~9000 transcripts 

# write.csv(annot.filt.random, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.filt.genefilter.random.20180227.csv") ### interchange genefilter.random and complete.random
# write.csv(annot.filt.novel, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.filt.genefilter.novel.20181026.csv") ### interchange genefilter.novel and complete.novel

###############################################################################
###############################################################################

### script for mb.vsd dataset for Janet
 
guiltByAssociation <-function(data, associated.gene, cores = 10){
  library(foreach)
  library(tictoc)
  library(parallel)
  library(doParallel)
  registerDoParallel(cores)
  
  index <- match(names(data), names(associated.gene))
  matched.data <- data[,!is.na(index)]
  matched.associated.gene <- associated.gene[index[!is.na(index)]] 
  
  guiltAssociation <- function(x,y){
    return(c(cor = cor.test(x,y)$estimate,
             p.val = cor.test(x,y)$p.value))
  }
  
  res <- foreach(i = 1:nrow(data), .combine = rbind)%dopar%{
    x <- as.numeric(matched.data[i,])
    return(guiltAssociation(x,matched.associated.gene))
  }
  
  rownames(res) <- rownames(data)
  return(res)
}

 
MYC <- as.numeric(mb.vsd["ENSG00000136997.15_1",])
names(MYC) <- colnames(mb.vsd)
guilt.res.MYC <- guiltByAssociation(mb.vsd, MYC)

annot <- annotate.HTseq.IDs(rownames(guilt.res.MYC)) ### note this worked when I loaded the annotate.HTseq.IDs function again (clinical_data_functions_master.R)
cbind(guilt.res.MYC, adj.p.val=p.adjust(guilt.res.MYC[,2], method = "BH"), annot) -> guilt.res.MYC
guilt.res.MYC[!is.na(guilt.res.MYC[,1]),] -> guilt.res.MYC
head(guilt.res.MYC[order(guilt.res.MYC[,1]),],20) ### get the first 20 associated transcripts 
tail(guilt.res.MYC[order(guilt.res.MYC[,1]),],20) ### get the last 20 associated transcripts



#####################################################################
#####################################################################

### earlier commands for clinical data
# library(readxl)
# database <- read_excel("~/R/Data/database.xlsx")

### installing packages for survival analysis
### choose mirror e.g Cambridge
# install.packages('survival')
# install.packages ('pwr')
# install.packages ('powerSurvEpi')
library(survival)
library(pwr)
library(powerSurvEpi)

### example for calculating power relative to exposure/control, survival and clinical trials
# powerCT.default(nE=, nC=, pE=, pC=, RR=, alpha=0.05) 
### where E=experimental group, C=control group, pE= probability of event in the experimental group, RR is relative risk
### if calculating ratios between E:C group, where k=ratio E:C, m= total number of expected events over both groups
# powerCT.default0(k=, m=, RR=, alpha=)

################################################################################
################################################################################
### troubleshooting
################################################################################

### script for  isolated set of transcripts to see that function is working. Changed names(goi.vsd) to names(x), goi.vsd is specified as "x" in script below:

### this is for the main expression dataset
# i = 1
# results.master <- foreach(i = 1:25)%dopar%{
# as.numeric(filt.mb.vsd[i,]) -> x
# names(x) <- colnames(filt.mb.vsd)
# names(x) <- gsub("T","",names(filt.mb.vsd)) 
# return(clinPathAssess(test.pData,x))
# }

# results.master <- foreach(i = 12200:12250)%dopar%{      ### investigating range to determine location of error
# as.numeric(mb.vsd.novel[i,]) -> x
#  names(x) <- colnames(mb.vsd.novel)
# names(x) <- gsub("T","",names(mb.vsd.novel)) 
#  return(clinPathAssess(test.pData,x))
# }

### other troubleshooting

# i=25 
# as.numeric(mb.vsd.novel[i,]) -> x
# names(x) <- colnames(mb.vsd.novel)
#  names(x) <- gsub("T","",names(mb.vsd.novel))
#  clinPathAssess(test.pData,x)      ### unhash here when trouble shooting
### additional graphics for distribution plots in primary dataset

# fig1 = ggplot(matched.test.incl.pData, aes(x=matched.test.incl.pData$Followup, y=matched.goi.vsd.incl)) ### need to update with additional graphics as per below
### define goi then run through clinPathAssess function  manually 
fig2 = ggplot(matched.test.incl.pData, aes(x=matched.test.incl.pData$OS.cat, y=matched.goi.vsd.incl)) +
                 geom_boxplot() +
                 geom_point(colour = "red") +
                 xlab ("Survival") +
                 ylab ("Gene expression") +
                 ggtitle ("Gene expression and survivorship")

                 
               




