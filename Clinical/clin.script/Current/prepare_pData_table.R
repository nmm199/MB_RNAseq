### File introduction
### File name: prepare_PData_table.R

### Aim of file is to 
# 1. Run basic descriptive statistics on a cohort of children treated for medulloblastoma, whose details are contained within the local clinical database
# 2. Analyse genes of interest in relation to univariate and multivariate risk prediction models for survival (overall, event-free and progression-free)
# 3. This script covers analysis up to and including univariate risk factor analysis
# 4. Multivariate analysis / AUC will be covered by a separate script


### Author: Dr Marion Mateos
### Date: July 19 2017
### this file was created by modifying clinical_data_4.R and subsequent 'for loops' added by Dr Louise Pease
### this file needs to be run before the biomarker_discovery_functions.R file ("nmm199/home/R/MB_RNAseq/Clinical/clin.script/Current/biomarker_discovery_functions.R")
### both this pData file and the biomarker discovery function file need to be run prior to the Biomarker_assess.R file ("nmm199/home/R/MB_RNAseq/Clinical/clin.script/Current/Biomarker_assess.R")


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
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("gplots", "car"))
library('gplots')
library(car)
library(stats)
library('survival')
# library(scales)



### Functions used

#source (file = "/home/nlp71/Marion/biomarker_discovery_functions.R")

### names of functions
### "chi.sq"
### "cor.result"
### "lin.reg"
### "logisticRegression"
### "km.log.test"
### "km.log.test.sub"
### "km.log.test.OS"
### "km.log.test.OS.sub"
### "cox.result.OS"
### "km.log.test.EFS"
### "km.log.test.EFS.sub"

### External files required

### clinical database "x.data"
### 7 molecular group data "meth.data"
### cytogenetic arm data "cytogen.data"
### RNA expression data "RNA.data"


#############################################
### setting up output files for log, messages and PDF files

#sink("pDatalog.txt")

# cat ()

#pdf("marker.results.pdf", width = 10, height = 10) 


### Read in pData 
### add in row names to original database file
cat ("reading in clinical database", sep ="\n")
#pData <- read.csv(file="/home/nlp71/Marion/database270617.csv", row.names = 1)
pData <- read.csv(file = "/home/nmm199/R/MB_RNAseq/Input data/database270617.csv", row.names = 1)
#pData <- read.csv(file = "/home/nmm199/R/MB_RNAseq/Input data/database270617_v2.csv", row.names = 1)
### updated the file to remove relapse to death as example to see if improves chi_square output

#meth.data <- "/home/nmm199/R/MB_RNAseq/Input data/all7subgroupCalls.csv"
#cytogen.data <- "/home/nmm199/R/MB_RNAseq/Input data/arm_calls_clean280617.txt"
#RNA.data <- "/home/dan/mygit/rna_seq_mb/paper/MB.vsd.txt"

### Read in methylation data
#meth.data <- "/home/nlp71/Marion/all7subgroupCalls.csv"
meth.data<- "/home/nmm199/R/MB_RNAseq/Input data/all7subgroupCalls.csv"
#cytogen.data <- "/home/nlp71/Marion/arm_calls_clean280617.txt"
cytogen.data <- "/home/nmm199/R/MB_RNAseq/Input data/arm_calls_clean280617.txt"

### reading in files
cat ("reading in expression data", sep ="\n")
mb.vsd <- read.delim(file="/home/dan/mygit/rna_seq_mb/paper/MB.vsd.txt")

#############################################
### setting up output files for log, messages and PDF files

#sink("pDatalog.txt")
# cat ()
#pdf("marker.results.pdf", width = 10, height = 10) 

####  providing a putative biomarker
### hashed out as now in main script to set one goi.vsd per gene (unhashed 7/8)
#goi <- "ENSG00000136997"

### unhash the line below if wish to run all goi results in RNAseq data (7/8/17)
#goi.vsd <- as.numeric(mb.vsd[goi,]) 

### the output would be a vector with a continuous variable, names equal NMB numbers
#names(goi.vsd) <- gsub("T","",names(mb.vsd))

###################################################
#### convert pData object columns into categorical variables

### identifier
NMB <- pData$NMB.NO
### sex
sex <- pData$Sex_R
sex <- ifelse(sex==1, "male","female")

### These are still used in the calculations and data frame but the age and group factors set above are included also as these could provide more insight
age.cont <- pData$Age_R
age.cat.infant <- pData$Age_R<3
age.cat.adult.16 <- pData$Age_R>16 
age.cat.adult.21 <- pData$Age_R>21 
age.filter <- !age.cat.infant&!age.cat.adult.16

### I have updated most of these so that the objects are set in one line only
### metastatic
mstatus <- pData$Mstatus_R

### resection status
resection <- ifelse(pData$Rstatus_R==1, "Gross total resection", "subtotal resection")

### RTX
RTX <- ifelse(pData$RTX_R == "Yes", "RTX", "No RTX")

### CSI
CSI <- ifelse(pData$RTXCSI_R =="Yes", "CSI", "No CSI")

### molecular subgroups
meth <-pData$X450K_R
meth7 <- read.csv(meth.data, header=TRUE, sep=",", quote="\"", dec=".", row.names=1)
meth7.cat <- meth7[c(1)]

### cytogenic data 
cytogen <- read.table (cytogen.data, header=T, sep="\t")
cytogen.q13.cat <- ifelse(cytogen$q13 == "Loss", "q13 loss", "no q13 loss") 

### pathology
histopath <-pData$path_R
LCA <- ifelse (pData$path_R =="LCA", "LCA", "non LCA")

### MYC
MYC.cat <-ifelse(pData$MYC_R==0, "MYC non ampl", "MYC ampl")

### MYCN
MYCN.cat <- ifelse(pData$MYCN_R, "MYCN non ampl", "MYCN ampl")

### MYC/MYCN ampl
MYCN.cat.ampl <- MYCN.cat =="MYCN ampl"
MYC.cat.ampl <- MYC.cat =="MYC ampl"
MYCMYCN.cat <- ifelse(MYC.cat.ampl =="TRUE"|MYCN.cat.ampl =="TRUE", "MYC MYCN ampl", "MYC MYCN non ampl")

MYC.df <- data.frame(MYC.cat, 
                     MYCN.cat, 
                     MYCMYCN.cat
)

###TP53
TP53.cat <- ifelse(pData$TP53_R==0, "TP53 WT", "TP53 mut")

###TERT
TERT.cat <- ifelse(pData$TERT_R==0, "TERT neg", "TERT pos")

### curative group
curative <- ifelse(pData$Curative_R=="Yes", "curative", "non curative")


### categorical survival PFS (defined as relapse, which includes relapse/progression), OS. 
relapse <- ifelse(pData$PFS_R==0, "non-relapse", "relapse")

OS.cat <-pData$OS_R

### continuous survival times PFS (years) and Followup (OS in years)
Followup <- pData$Followup_R
PFS <-pData$PFS_yr_R
#Relapsetodeath <- pData$relapsetodeath_R
Event <- pData$Event_R
EFS <- pData$EFS_R

####################################################

### visualisation of the data

cat ("reading out the summary statistics for the entire cohort n=802", sep ="\n")

### assessing the data, parametric etc
qqnorm(age.cont)
plot(age.cont, col='red') 
#text(age.cont)

### Getting summary characteristics
summary(age.cont)
mean(age.cont, na.rm = T)
median(age.cont, na.rm = T)

cat ("age summary statistics n=802", sep ="\n")
t.test(age.cont, conf.level=0.95)

boxplot(age.cont~meth, col=c("yellow","green","red","blue"))


### install.packages('gplots')

age.means <- plotmeans(age.cont~meth, data=pData)


###########################################
### Combined dataframes 

### need to add the 7 group methylation status using dataframes
test.meth7 <- data.frame(meth7.cat)

test.meth7$Sample_Name <- rownames(test.meth7)

cytogen.q13.cat <- cytogen [c("SampleID", "q13")]

### need to convert q13 loss into loss, and rest into "no loss"

cytogen.q13 <- ifelse(cytogen.q13.cat$q13 =="Loss", "q13 Loss", "No q13 loss")

### make q13 loss dataframe

cytogen.q13.df <- data.frame(cytogen.q13.cat[,-1], 
                             row.names=cytogen.q13.cat [,1],
                             cytogen.q13)

#### creation of groups for the analysis of the data 
#### the group factor is designed to take account of age, gender and sex differences in development
data=pData
agegrpfac=rep(NA,length(data))
#sexfac=rep(NA,length(data$Sex_R))
agegrpfac[data$Age_R<=3&data$Sex_R==1]="Infant.M"
agegrpfac[data$Age_R<=3&data$Sex_R==2]="Infant.F"
agegrpfac[data$Age_R>3&data$Age_R<=10&data$Sex_R==2]="Junior.F"
agegrpfac[data$Age_R>3&data$Age_R<=12&data$Sex_R==1]="Junior.M"
agegrpfac[data$Age_R>3&data$Age_R<=10&data$Sex_R==2]="Junior.F"
agegrpfac[data$Age_R>3&data$Age_R<=12&data$Sex_R==1]="Junior.M"
agegrpfac[data$Age_R>10&data$Age_R<=16&data$Sex_R==2]="Teenager.F"
agegrpfac[data$Age_R>12&data$Age_R<=16&data$Sex_R==1]="Teenager.M"
agegrpfac[data$Age_R>16&data$Age_R<=22&data$Sex_R==2]="Adolescent.F"
agegrpfac[data$Age_R>16&data$Age_R<=22&data$Sex_R==1]="Adolescent.M"
agegrpfac[data$Age_R>22&data$Sex_R==2]="Adult.F"
agegrpfac[data$Age_R>22&data$Sex_R==1]="Adult.M"
data=pData
childfac=rep(NA,length(data))
childfac[data$Age_R>3&data$Age_R<=16&data$Sex_R==2]="Child.F"
childfac[data$Age_R>3&data$Age_R<=16&data$Sex_R==1]="Child.M"

#### the age factor is designed to assign samples to groups purely based on age not taking account of sex differences in development
data=pData$Age_R
# age factor
agefac=rep(NA,length(data))
agefac[data<=3]="Infant"
agefac[data>3&data<=10]="Junior"
agefac[data>10&data<=16]="Teenager"
agefac[data>16&data<=21]="Adult_16_21"
agefac[data>22]="Adult_over_22"

#### The sex factor determines the gender of the samples
data=pData$Sex_R
sexfac=rep(NA,length(data))
sexfac[data==1]="M"
sexfac[data==2]="F"
#### Set the subgroup based on methylation data 4 subgroups
data=pData$X450K_R
subgroup4fac=rep(NA,length(data))
subgroup4fac[data=="G4"]="G4"
subgroup4fac[data=="G3"]="G3"
subgroup4fac[data=="SHH"]="SHH"
subgroup4fac[data=="WNT"]="WNT"
data=pData$X450K_R
group3or4fac=rep(NA,length(data))
group3or4fac[data=="G3" | data=="G4"]= "G3G4"

#### This command generates a new variable facs, this variable combines the age and gender factors such that they do not take account of sex differences in developmental stage 
#facs = paste(agefac, sexfac, subgroup, sep=".")
#f = factor(facs)
#### This command generates a new variable facs, this variable combines the age and gender factors such that they DO take account of sex differences in developmental stage 
#facs = paste(agefac, sexfac, subgroup, sep=".")
#facs = paste(agegrpfac, subgroupfac, sep=".")
#f = factor(facs)
#### We can get the summary statistics from the data using the table function
#### When puberty ages are set as different for males and females it becomes apparent Juniors are mostly affected and far more males than females
#table(agegrpfac)
#### Again Juniors mostly affected
#table(agefac)
#### Nearly twice as many males as females
#table(sexfac)
#### Only the child samples have subgroups, grp4 high risk much more common in junior males 
#table(facs)
#### Group 4 most common
#table(subgroup4fac)
#table(subgroup7fac)

test.pData <- data.frame(NMB, 
                         age.cont,
                         age.cat.infant,
                         age.cat.adult.16,
                         age.cat.adult.21,
                         age.filter,
                         agefac,
                         sex,
                         sexfac,
                         agegrpfac,
                         subgroup4fac,
                         group3or4fac,
                         childfac,
                         mstatus,
                         resection,
                         RTX,
                         CSI,
                         curative,
                         meth,
                         histopath,
                         LCA,
                         MYC.cat,
                         MYCN.cat,
                         MYCMYCN.cat,
                         TP53.cat,
                         TERT.cat,
                         relapse,
                         OS.cat,
                         PFS,
                         Followup,
                         #Relapsetodeath,
                         Event,
                         EFS
)
rownames(test.pData) <- rownames(pData)

test.pData$meth7 <- test.meth7[match(rownames(test.pData), rownames(test.meth7)), ]$all.calls

test.pData$q13loss <- cytogen.q13.df[match(rownames(test.pData), rownames(cytogen.q13.df)),]$cytogen.q13
data=test.pData$meth7
subgroup7fac=rep(NA,length(data))
subgroup7fac[data=="Grp4_HighRisk"]="4HR"
subgroup7fac[data=="Grp4_LowRisk"]="4LR"
subgroup7fac[data=="Grp3_LowRisk"]="3LR"
subgroup7fac[data=="Grp3_HighRisk"]="3HR"
subgroup7fac[data=="SHH_Inf"]="SHH_Inf"
subgroup7fac[data=="SHH_Old"]="SHH_Old"
subgroup7fac[data=="WNT"]="WNT"
test.pData$subgroup7fac <- subgroup7fac


index <- match(names(goi.vsd), rownames(test.pData)) 
matched.test.pData <- test.pData[index[!is.na(index)],] 
is.vector(matched.test.pData)
matched.goi.vsd <- goi.vsd[!is.na(index)] 
matched.goi.vsd.cat <- ifelse(matched.goi.vsd>median(goi.vsd, na.rm = T), "high","low") 


### add cytogenetic data to existing data frame, NMB650 duplicate (650, 650p for paraffin, however not included in the survival cohort)

cytogen.q13.cat <- cytogen [c("SampleID", "q13")]

### need to convert q13 loss into loss, and rest into "no loss"

cytogen.q13 <- ifelse(cytogen.q13.cat$q13 =="Loss", "q13 Loss", "No q13 loss")

### make q13 loss dataframe

cytogen.q13.df <- data.frame(cytogen.q13.cat[,-1], 
                             row.names=cytogen.q13.cat [,1],
                             cytogen.q13)

# View(cytogen.q13.df)

test.pData$q13loss <- cytogen.q13.df[match(rownames(test.pData), rownames(cytogen.q13.df)),]$cytogen.q13

age.mean.goi <- plotmeans(test.pData$age.cont ~ test.pData$meth, data=test.pData)

#matched.test.pData$q13loss <- cytogen.q13.df[match(rownames(matched.test.pData), rownames(cytogen.q13.df)),]$cytogen.q13

#age.mean.goi <- plotmeans(matched.test.pData$age.cont ~ matched.test.pData$meth, data=matched.test.pData)

### summary data for all variables within the age 3-16 yo group, to construct survival cohort
age.df <- data.frame(NMB, age.filter, age.cont)
summary(age.df)

save.image("pData_saved.RData")

