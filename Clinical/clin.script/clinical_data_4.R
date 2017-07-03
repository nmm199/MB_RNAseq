
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
library('gplots')
library(car)
library('survival')



### Functions used

source (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/clinical_data_functions.R")
### names of functions
### "chi.sq"
### "cor.result"
### "km.log.test"
### "km.log.test.OS"



### External files required

### clinical database "x.data"
### 7 molecular group data "meth.data"
### cytogenetic arm data "cytogen.data"
### RNA expression data "RNA.data"


x.data <- "/home/nmm199/R/MB_RNAseq/Input data/database160617.csv"
meth.data <- "/home/nmm199/R/MB_RNAseq/Input data/all7subgroupCalls.csv"
cytogen.data <- "/home/nmm199/R/MB_RNAseq/Input data/arm_calls_clean280617.txt"
RNA.data <- "/home/dan/mygit/rna_seq_mb/paper/MB.vsd.txt"



#############################################
### setting up output files for log, messages and PDF files

sink("pDatalog.txt")

# cat ()

pdf("marker.results.pdf", width = 10, height = 10) 



##################################################
### reading in files

cat ("reading in expression data")

mb.vsd <- read.delim(RNA.data)
  
cat ("reading in clinical database")

### add in row names to original database file

pData <- read.csv(x.data, row.names = 1)


####  providing a putative biomarker

goi <- "ENSG00000136997"

goi.vsd <- as.numeric(mb.vsd[goi,]) 


### the output would be a vector with a continuous variable, names equal NMB numbers

names(goi.vsd) <- gsub("T","",names(mb.vsd))



###################################################


#### convert pData object columns into categorical variables

### identifier
NMB <- pData$NMB.NO


### sex
sex <- pData$Sex_R
sex <- ifelse(sex==1, "male","female")


### age 
age.cont <- pData$Age_R
age.cat.infant <- pData$Age_R<3
age.cat.adult.16 <- pData$Age_R>16 
age.cat.adult.21 <- pData$Age_R>21 
age.filter <- !age.cat.infant&!age.cat.adult.16

### metastatic
mstatus <- pData$Mstatus_R

### resection status
resection <- pData$Rstatus_R
resection <- ifelse(resection==1, "Gross total resection", "subtotal resection")

### RTX
RTX <- pData$RTX_value

### CSI
CSI <-pData$RTXCSI_R

### molecular subgroups
meth <-pData$X450K_R
meth7 <- read.csv(meth.data, header=TRUE, sep=",", quote="\"", dec=".", row.names=1)
meth7.cat <- meth7[c(1)]

cytogen <- read.table (cytogen.data, header=T, sep="\t")

cytogen.q13.cat <- ifelse(cytogen$q13 == "Loss", "Loss", "no loss") 

### pathology
histopath <-pData$path_R
LCA <- ifelse (pData$path_R =="LCA", "LCA", "non LCA")

### MYC
MYC.cat <-pData$MYC_R
MYC.cat <-ifelse(MYC.cat==0, "MYC non ampl", "MYC ampl")

### MYCN
MYCN.cat <-pData$MYCN_R
MYCN.cat <- ifelse(MYCN.cat==0, "MYCN non ampl", "MYCN ampl")

### MYC/MYCN ampl

MYCN.cat.ampl <- MYCN.cat =="MYCN ampl"
MYC.cat.ampl <- MYC.cat =="MYC ampl"
MYCMYCN.cat <- ifelse(MYC.cat.ampl =="TRUE"|MYCN.cat.ampl =="TRUE", "MYC MYCN ampl", "MYC MYCN non ampl")

MYC.df <- data.frame(MYC.cat, 
                     MYCN.cat, 
                     MYCMYCN.cat
                     )


###TP53
TP53.cat<- pData$TP53_R
TP53.cat <- ifelse(TP53.cat==0, "TP53 WT", "TP53 mut")

###TERT
TERT.cat <- pData$TERT_R
TERT.cat <- ifelse(TERT.cat==0, "TERT neg", "TERT pos")

### categorical survival PFS (defined as relapse, which includes relapse/progression), OS. 
relapse <- pData$PFS_R
relapse <- ifelse(relapse==0, "non-relapse", "relapse")

OS.cat <-pData$OS_R

### continuous survival times PFS (years) and Followup (OS in years)
Followup <- pData$Followup_R
PFS <-pData$PFS_yr_R
Relapsetodeath <- pData$relapsetodeath_R
Event <- pData$Event_R
EFS <- pData$EFS_R


####################################################

### visualisation of the data

cat ("reading out the summary statistics for the entire cohort n=802")

### assessing the data, parametric etc
qqnorm(age.cont)
plot(age.cont, col='red')


### Getting summary characteristics
summary(age.cont)
mean(age.cont, na.rm = T)
median(age.cont, na.rm = T)


cat ("age summary statistics n=802")
t.test(age.cont, conf.level=0.95)
boxplot(age.cont~meth, col=c("yellow","green","red","blue"))

### install.packages('gplots')

age.means <- plotmeans(age.cont~meth, data=pData)


###########################################
### Combined dataframes 

test.pData <- data.frame(NMB, 
                    age.cont,
                    age.cat.infant,
                    age.cat.adult.16,
                    age.cat.adult.21,
                    age.filter,
                    sex,
                    mstatus,
                    resection,
                    RTX,
                    CSI,
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
                    Relapsetodeath,
                    Event,
                    EFS
                    )
rownames(test.pData) <- rownames(pData)

View(test.pData)

### need to add the 7 group methylation status using dataframes
test.meth7 <- data.frame(meth7.cat)

test.meth7$Sample_Name <- rownames(test.meth7)
test.pData$meth7 <- test.meth7[match(rownames(test.pData), rownames(test.meth7)), ]$all.calls

View(test.pData)
save(test.pData, file = "/home/nmm199/R/MB_RNAseq/Clinical/test.pData") 



#############################################


### goi.vsd #### name of a putative biomarker
### matching the expression data (goi data) with the initially compiled clinical data with important variables

index <- match(names(goi.vsd), rownames(test.pData)) 
matched.test.pData <- test.pData[index[!is.na(index)],] 
is.vector(matched.test.pData)
matched.goi.vsd <- goi.vsd[!is.na(index)] 
matched.goi.vsd.cat <- ifelse(matched.goi.vsd>median(goi.vsd, na.rm = T), "high","low") 



#############################################

### summary data 

### visualising means of the goi in the matched dataset

age.mean.goi <- plotmeans(matched.test.pData$age.cont ~ matched.test.pData$meth, data=matched.test.pData)


### summary data for all variables within the age 3-16 yo group, to construct survival cohort
age.df <- data.frame(NMB, age.filter, age.cont)
summary(age.df)



##############################################

### Chi squared analysis


list.age.cat.infant <- chi.sq(x = matched.test.pData$age.cat.infant, y = matched.goi.vsd.cat)

list.sex <- chi.sq (x = matched.test.pData$sex, y= matched.goi.vsd.cat)

list.mstatus <- chi.sq(x = matched.test.pData$mstatus, y=matched.goi.vsd.cat) 

list.relapse <- chi.sq(x = matched.test.pData$relapse, y = matched.goi.vsd.cat) 

list.resection <- chi.sq (x = matched.test.pData$resection, y = matched.goi.vsd.cat)

list.meth.4 <- chi.sq (x= matched.test.pData$meth, y = matched.goi.vsd.cat)

list.meth.7 <- chi.sq (x = matched.test.pData$meth7, y = matched.goi.vsd.cat)

list.path <- chi.sq (x = matched.test.pData$histopath, y = matched.goi.vsd.cat)

list.MYC <- chi.sq (x = matched.test.pData$MYC.cat, y = matched.goi.vsd.cat)

list.MYCN <- chi.sq (x = matched.test.pData$MYCN.cat, y = matched.goi.vsd.cat)

list.MYCMYCN <- chi.sq(x= matched.test.pData$MYCMYCN.cat, y= matched.goi.vsd.cat)

list.TP53 <- chi.sq (x = matched.test.pData$TP53.cat, y = matched.goi.vsd.cat)

list.TERT <- chi.sq (x = matched.test.pData$TERT.cat, y = matched.goi.vsd.cat)


### is the biomarker overrepresented in poor prognostic groups or those who received different therapy

list.RTX <- chi.sq (x = matched.test.pData$RTX, y = matched.goi.vsd.cat)

list.CSI <- chi.sq (x = matched.test.pData$CSI, y = matched.goi.vsd.cat)


### run Fisher's exact test on those where count < 5 e.g pathology "other"

histopath.table <- table(as.factor(matched.test.pData$histopath), as.factor(matched.goi.vsd.cat))

histopath.result <- fisher.test(histopath.table)
histopath.result



###################################

### Correlation coefficients

cat ("correlation coefficients for association between variables")


list.cor.age <- cor.result(x = matched.test.pData$age.cont, y = matched.goi.vsd)
list.cor.age


##################################

###  linear regression

lin.reg.age <- lin.reg(x= matched.test.pData$age.cont, y= matched.goi.vsd)
lin.reg.age

#################################

### Mann-whitney U (aka wilcoxon rank sum test) for non-parametric data

### age continuous 
age.cont.wilcox <- wilcox.test(matched.test.pData$age.cont ~ matched.goi.vsd.cat, exact = F, correct = F)
age.cont.wilcox
# age.cont.wilcox$p.value

### comparison to t-test for understanding data
# t.test(matched.test.pData$age.cont ~ matched.goi.vsd.cat)



##################################
### logistic regression

message("processing logistic regression individually")

### age categorical

log.reg.age.cat <- glm(age.cat.infant ~ matched.goi.vsd, family = binomial(link= 'logit'), data=matched.test.pData)
summary(log.reg.age.cat)
age.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$age.cat.infant, col = c("red", "blue"), xlab = "Infant", ylab = "Biomarker expression", main = "Correlation between biomarker and age (infant vs non infant)")


### sex 

log.reg.sex <- glm (matched.test.pData$sex ~ matched.goi.vsd, family = binomial (link = 'logit'), data= matched.test.pData)
summary(log.reg.sex)
sex.boxplot <- boxplot (matched.goi.vsd ~ matched.test.pData$sex, col = c("red", "blue"), xlab = "Gender", ylab = "Expression of biomarker", main = "Biomarker expression and gender")

### metastatic status

log.reg.mstatus <- glm(matched.test.pData$mstatus ~ matched.goi.vsd,  family = binomial(link= 'logit'), data=matched.test.pData)
summary(log.reg.mstatus)

mstatus.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$mstatus, col = c("red", "blue"), xlab = "M status", ylab = "Biomarker expression", main = "Correlation between biomarker and metastatic status")

### this also works for metastatic status
# log.reg.mstatus2 <- glm(mstatus ~ matched.goi.vsd,  family = binomial(link= 'logit'), data=matched.test.pData)
# summary(log.reg.mstatus2)


### relapse 

log.reg.relapse <- glm(matched.test.pData$relapse ~ matched.goi.vsd, family = binomial(link= 'logit'), data=matched.test.pData)
summary(log.reg.relapse)
relapse.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$relapse, col = c("red", "blue"), xlab = "Relapse status", ylab = "Biomarker expression",  main = "Correlation between biomarker and relapse")


### resection

log.reg.resection <- glm (matched.test.pData$resection ~ matched.goi.vsd, family = binomial(link= 'logit'), data=matched.test.pData)
summary(log.reg.resection)
resection.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$resection, col = c("red", "blue"), xlab = "Resection status", ylab = "Biomarker expression", main = "Correlation between biomarker and resection status")


##########################################


### logistic regression and molecular subgroups
### includes some more detailed assumption testing and post-hoc analysis
message("logistic regression and association with molecular subgroups")

log.reg.meth <- glm (matched.test.pData$meth ~ matched.goi.vsd, family = binomial(link='logit'), data = matched.test.pData)
summary(log.reg.meth)                  
str(log.reg.meth)

log.reg.meth7 <- glm (matched.test.pData$meth7 ~ matched.goi.vsd, family = binomial(link='logit'), data = matched.test.pData)
summary(log.reg.meth7) 

  
### visualise distribution of biomarker in cohort

message ("visualisation of biomarker and relationship with methylation")

qqnorm(matched.goi.vsd)

### visualise relationship between biomarker and methylation groups

meth.boxplot <- boxplot(matched.goi.vsd~matched.test.pData$meth, col=c("yellow","green","red","blue"), xlab = "Methylation subgroup", ylab = "Biomarker expression", main = "Correlation between biomarker and 4 molecular subgroups")

meth7.boxplot <- boxplot(matched.goi.vsd~matched.test.pData$meth7, col=c("yellow","green","red","blue"),  xlab = "Methylation subgroup", ylab = "Biomarker expression", main = "Correlation between biomarker and 7 molecular subgroups")

### if biomarker is normally distributed, can use ANOVA (one-way ANOVA)

meth.aov <- aov (matched.goi.vsd ~ matched.test.pData$meth, data=matched.test.pData)
summary(meth.aov)
plot(meth.aov)

meth7.aov <- aov (matched.goi.vsd ~ matched.test.pData$meth7, data = matched.test.pData)
summary(meth7.aov)
plot(meth7.aov)

### AOV depends on Levene's test for homogeneity of variances, H0= variances are equal, HA= variances are not equal. Therefore wish for H0 to be true ie cannot reject H0 , i.e variances are equal (homogenous)
### check variances are equal (i.e p>0.05, to confirm assumption upon which one-way ANOVA is created

cat ("confirm that variances are equal, ie. leveneTest is not significant")
leveneTest(matched.goi.vsd ~ matched.test.pData$meth7)

# plot(matched.goi.vsd, matched.test.pData$meth)

# qqplot(matched.goi.vsd, matched.test.pData$meth7)


### post-hoc tests for methylation

### pairwise t-test to determine where the difference lies between the groups
meth.pw <- pairwise.t.test(matched.goi.vsd, matched.test.pData$meth)

meth.pw

### alternative coding, creates same boxplots as above

# plot(matched.test.pData$meth, matched.goi.vsd, col = c("Yellow", "Green", "red","blue"),  xlab = "methylation groups", ylab= "Expression of biomarker", main = "Correlation with methylation and biomarker")


# can specify bonferroni correction, which may be over-conservative
# meth.pw.bon <- pairwise.t.test(matched.goi.vsd, matched.test.pData$meth, p.adj='bonf')
# meth.pw.bon

###################################

### logistic regression continued

message ("logistic regression continued")

### histopath

log.reg.histopath <- glm (matched.test.pData$histopath ~ matched.goi.vsd, family = binomial(link='logit'), data=matched.test.pData)
summary(log.reg.histopath)

histopath.pw <- pairwise.t.test(matched.goi.vsd, matched.test.pData$histopath)
histopath.pw
histopath.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$histopath, col = c("red", "blue"), xlab = "Histopathology subtype", ylab = "Biomarker expression", main = "Correlation between biomarker and histopathology")


### MYC.cat

log.reg.MYC <- glm (matched.test.pData$MYC.cat ~ matched.goi.vsd, family = binomial(link='logit'), data=matched.test.pData)
summary(log.reg.MYC)
MYC.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$MYC.cat, col = c("Red", "Blue"), xlab = "MYCN expression", ylab = "Biomarker expression", main = "Correlation between biomarker and MYC expression")

### MYCN.cat

log.reg.MYCN <- glm (matched.test.pData$MYCN.cat ~ matched.goi.vsd, family = binomial (link='logit'), data=matched.test.pData)
summary(log.reg.MYCN)
MYCN.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$MYCN.cat, col = c("Red", "Blue"), xlab = "MYCN expression", ylab = "Biomarker expression", main = "Correlation between biomarker and MYCN expression")

### TP53

log.reg.TP53 <- glm (matched.test.pData$TP53.cat ~ matched.goi.vsd,family = binomial (link='logit'), data=matched.test.pData)
summary(log.reg.TP53)
TP53.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$TP53.cat, col = c("Red", "Blue"), xlab = "TP53 mutational status", ylab = "Biomarker expression", main = "Correlation between biomarker and TP53 mutational status")


### additional subgroup specific tests

### TERT (may need to specific subgroup)

log.reg.TERT <- glm (matched.test.pData$TERT.cat ~ matched.goi.vsd, family = binomial (link = 'logit'), data = matched.test.pData)
summary(log.reg.TP53)
        
TERT.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$TERT.cat, col = c("Red", "Blue"), xlab = "TERT status", ylab = "Biomarker expression", main = "Correlation between biomarker and TERT status")

####################
### summary data for each methylation subgroup (n=4), no age restriction for first analyses

message ("processing summary stats for 4 molecular subgroups, no age restriction")

### G3

G3 <- matched.test.pData$meth =="G3" 
View(matched.test.pData$meth)
G3.group <- matched.test.pData [G3, ]
summary (G3.group)

## G4

G4 <- matched.test.pData$meth =="G4"
G4.group <- matched.test.pData [G4, ]
summary(G4.group)

### confirmed that the above is valid through subsetting the object and defining removal of na
# str(G4.group)
# median(G4.group$age.cont, na.rm=T)
# summary(G4.group$age.cont)
# View(G4.group)

### SHH

SHH <- matched.test.pData$meth =="SHH"
SHH.group <- matched.test.pData [SHH, ]
summary(SHH.group)

### WNT

WNT <- matched.test.pData$meth == "WNT"
WNT.group <- matched.test.pData [WNT, ]
summary(WNT.group)

message ("processing summary stats for 7 molecular subgroups, no age restriction")

### G3 subgroups

G3.high.group <- matched.test.pData$meth7 =="Grp3_HighRisk"
G3.high.data <- matched.test.pData [G3.high.group, ]
summary(G3.high.data)

G3.low.group <- matched.test.pData$meth7 =="Grp3_LowRisk"
G3.low.data <- matched.test.pData [G3.low.group, ]
summary(G3.low.data)


### G4 subgroups

G4.high.group <- matched.test.pData$meth7 =="Grp4_HighRisk"
G4.high.data <- matched.test.pData [G4.high.group, ]
summary(G4.high.data)

G4.low.group <- matched.test.pData$meth7 =="Grp4_LowRisk"
G4.low.data <- matched.test.pData [G4.low.group, ]
summary(G4.low.data)

### SHH
SHH.inf.group <- matched.test.pData$meth7 == "SHH_Inf"
SHH.inf.data <- matched.test.pData [SHH.inf.group, ]
summary (SHH.inf.data)

SHH.old.group <- matched.test.pData$meth7 =="SHH_Old"
SHH.old.data <- matched.test.pData [SHH.old.group, ]
summary (SHH.old.data)


################

### restrict analysis to age 3-16 yo 

message("restrict analysis to age 3-16 years")
message ("create dataframe for age 3-16 years")

Age.incl <-matched.test.pData$age.filter =="TRUE"
Age.incl.df <- matched.test.pData [Age.incl, ]
summary(Age.incl.df)
View(Age.incl.df)

### compare to prior dataframes to check accuracy of new dataframe

message ("comparing with previous data frames for accuracy")
summary(test.pData)
summary (matched.test.pData)

message ("processing summary stats for 3-16 yo for 4 molecular subgroups")

### G3, aged 3-16 years
G3.incl <- Age.incl.df$meth =="G3"
G3.group.incl <- Age.incl.df [G3.incl, ]
summary(G3.group.incl)
nrow(G3.group.incl)

### G4, aged 3-16 years
G4.incl <- Age.incl.df$meth =="G4"
G4.group.incl <- Age.incl.df [G4.incl, ]
summary(G4.group.incl)
nrow(G4.group.incl)

### SHH, aged 3-16 years
SHH.incl <- Age.incl.df$meth == "SHH"
SHH.group.incl <- Age.incl.df [SHH.incl, ]
summary (SHH.group.incl)
nrow(SHH.group.incl)

### WNT, aged 3-16 years
WNT.incl <- Age.incl.df$meth == "WNT"
WNT.group.incl <- Age.incl.df [WNT.incl, ]
summary (WNT.group.incl)
nrow (WNT.group.incl)


### defining features of 7 molecular groups, group 3 and 4 subgroups

message ("processing summary stats for 7 molecular subgroup data, age 3-16 years")

G3.high <- Age.incl.df$meth7 == "Grp3_HighRisk"
G3.high.incl <- Age.incl.df [G3.high, ]
summary(G3.high.incl)

G3.low <- Age.incl.df$meth7 == "Grp3_LowRisk"
G3.low.incl <- Age.incl.df [G3.low, ]
summary(G3.low.incl)


G4.high <- Age.incl.df$meth7 =="Grp4_HighRisk"
G4.high.incl <- Age.incl.df [G4.high, ]
summary (G4.high.incl)

G4.low <- Age.incl.df$meth7 =="Grp4_LowRisk"
G4.low.incl <- Age.incl.df [G4.low, ]
summary (G4.low.incl)

SHH.inf <- Age.incl.df$meth7 == "SHH_Inf"
SHH.inf.incl <- Age.incl.df [SHH.inf, ]
summary (SHH.inf.incl)

SHH.old <- Age.incl.df$meth7 =="SHH_Old"
SHH.old.incl <- Age.incl.df [SHH.old, ]
summary (SHH.old.incl)




########################

### add cytogenetic data to existing data frame, NMB650 to be resolved but will not be included in survival analysis therefore add cytogenetic data at this stage
### duplicate NMB650 data emailed Ed 

cytogen.q13.cat <- cytogen [c("SampleID", "q13")]

### need to convert q13 loss into loss, and rest into "no loss"

cytogen.q13 <- ifelse(cytogen.q13.cat$q13 =="Loss", "Loss", "No loss")

### make q13 loss dataframe

cytogen.q13.df <- data.frame(cytogen.q13.cat[,-1], 
                             row.names=cytogen.q13.cat [,1],
                             cytogen.q13)

View(cytogen.q13.df)

matched.test.pData$q13loss <- cytogen.q13.df[match(rownames(matched.test.pData), rownames(cytogen.q13.df)),]$cytogen.q13

View(matched.test.pData)









######## other code to consider using
### creating survival dataframe for age 3-16
### incorporating cytogenetic data
### then can run the chi.sq analysis
### transform c/s 13 loss yes/no as a variable 


#cytogen.univar <- foreach(i = 1:length(cytogen))%do%{
# temp.starts <- start(cytogen[[i]]) 
# temp.chisq <- chi.sq(as.factor(cytogen[[i])], matched.test.incl.pData$relapse)
# temp.vector <- rep(NA, max(temp.starts)) 
# temp.vector[temp.starts] <-cytogen[[i]]$score
#  return(temp.vector)
#}



# index.cytogen.incl <- match(names(goi.vsd), rownames(Age.incl.df)) 
# matched.test.incl.pData <- Age.incl.df[index.incl[!is.na(index.incl)],] 
# is.vector(matched.test.incl.pData)
# matched.goi.vsd.incl <- goi.vsd[!is.na(index.incl)] 
# matched.goi.vsd.cat.incl <- ifelse(matched.goi.vsd.incl>median(goi.vsd, na.rm = T), "high","low")





#######################


### visualisation options to explore

# plot(ks.test(x[y],x[!y]))

# plot(ecdf(x[y]), col = "red")

# lines(ecdf(x[!y]), col = "blue")

# plot(density(x[y], na.rm = T), col = "red")
# lines(density(x[!y], na.rm = T), col = "blue")
# rug(x[y], col = "red")
# rug(x[!y], col = "blue")
# temp.pval <- 0.01
# text(0.06, 0.06, paste("P.val  = ", temp.pval))
# abline(v = median(x[y], na.rm = T), lty = 2, col = "red")
# abline(v = median(x[!y], na.rm = T), lty = 2, col = "blue")


########################
### survival analysis



### evaluate effect of biomarker on PFS
message ("evaluate effect of biomarker on PFS")

### change relapse to binary
relapse.bin <- ifelse(matched.test.pData$relapse == "relapse", 1, 0)


message ("restrict survival analysis for age 3-16 years")
message ("creating matched data frame")

index.incl <- match(names(goi.vsd), rownames(Age.incl.df)) 
matched.test.incl.pData <- Age.incl.df[index.incl[!is.na(index.incl)],] 
is.vector(matched.test.incl.pData)
matched.goi.vsd.incl <- goi.vsd[!is.na(index.incl)] 
matched.goi.vsd.cat.incl <- ifelse(matched.goi.vsd.incl>median(goi.vsd, na.rm = T), "high","low")


relapse.bin.incl <- ifelse(matched.test.incl.pData$relapse == "relapse", 1,0)


# matched.test.incl.pData$PFS -> time
# relapse.bin.incl -> event
# matched.goi.vsd.cat.incl -> marker

km.log.test <- function(time, event, marker, out.file = "none"){
if(out.file!="none"){
  pdf(out.file)
}
  km.PFS.incl <- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
plot(km.PFS.incl, col = c("red", "blue"),xlab = "time to progression/relapse (years)", ylab = "PFS", xlim = c(0,10), main = "Biomarker expression and progression-free survival (PFS)",  lty = 1:2)
PFS.names <- c("biomarker - high", "biomarker - low")
legend (x="topright", PFS.names,  lty= 1:2, col = c("red","blue"))
PFS.incl.logrank <- survdiff(Surv(time, event) ~ marker)
1 - pchisq(PFS.incl.logrank$chisq, length(PFS.incl.logrank$obs)-1) -> surv.p.val
text(4,0.1,paste("p =",round(surv.p.val, 3)), pos = 4, cex = 1)
if(out.file!="none"){
  dev.off()
}
}

km.log.test(matched.test.incl.pData$PFS,relapse.bin.incl, matched.goi.vsd.cat.incl )







### cox regression analysis
### exp(coef) is the hazard ratio, with 95% CI. p value is listed above against the variable
# help ("coxph")

message("cox regression analysis, PFS and biomarker")

cox.relapse <- coxph (Surv(matched.test.pData$PFS, relapse.bin) ~ matched.goi.vsd.cat)
summary(cox.relapse)$logtest
summary(cox.relapse)


### 

message("cox regression analysis, PFS, biomarker, age 3-16 years")

cox.relapse.incl <- coxph (Surv(matched.test.incl.pData$PFS, relapse.bin.incl) ~ matched.goi.vsd.cat.incl)
summary(cox.relapse.incl)$logtest
summary(cox.relapse.incl)



  

##########################################

### OS

### check alive status is binary
matched.test.pData$OS.cat

# OS.cat.bin <- ifelse(matched.test.pData$OS.cat == "Dead", 1, 0)
# OS.cat.bin

OS.cat.bin.incl <- ifelse(matched.test.incl.pData$OS.cat == "Dead", 1,0)
km.OS.incl <- survfit(Surv(matched.test.incl.pData$Followup, OS.cat.bin.incl)~matched.goi.vsd.cat.incl)
km.OS.incl


km.log.test.OS <- function(time, event, marker, out.file = "none"){
  if(out.file!="none"){
    pdf(out.file)
  }
  km.OS.incl <- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
  plot(km.OS.incl, col = c("red", "blue"),xlab = "overall survival (years)", ylab = "OS", xlim = c(0,10), main = "Biomarker expression and overall survival (OS)",  lty = 1:2)
  OS.names <- c("biomarker - high", "biomarker - low")
  legend (x="topright", OS.names,  lty= 1:2, col = c("red","blue"))
  OS.incl.logrank <- survdiff(Surv(time, event) ~ marker)
  1 - pchisq(OS.incl.logrank$chisq, length(OS.incl.logrank$obs)-1) -> surv.p.val
  text(4,0.1,paste("p =",round(surv.p.val, 3)), pos = 4, cex = 1)
  if(out.file!="none"){
    dev.off()
  }
}

km.log.test.OS(matched.test.incl.pData$Followup, OS.cat.bin.incl, matched.goi.vsd.cat.incl )






### cox regression analysis
#message("cox regression for OS on entire cohort")

#cox.OS <- coxph (Surv(matched.test.pData$Followup, OS.cat.bin) ~ matched.goi.vsd.cat)
#summary(cox.OS)$logtest
#summary(cox.OS)

### 

message("cox regression on age 3-16 years")

cox.OS.incl <- coxph (Surv(matched.test.incl.pData$Followup, OS.cat.bin.incl) ~ matched.goi.vsd.cat.incl)
summary(cox.OS.incl)$logtest
summary(cox.OS.incl)


cox.result.OS <- function (time, event, marker, strata = NULL)  
{
 if(is.null(strata)){
   cox.temp <- coxph (Surv(time, event)~marker, data= matched.test.incl.pData)
 }else{
   cox.temp <- coxph (Surv(time, event)~marker, data= matched.test.incl.pData)
}
     summary.cox <- c(rownames(summary(cox.temp)$coefficients),summary(cox.temp)$coefficients,
     summary(cox.temp)$n,
     summary(cox.temp)$nevent)
     names(summary.cox) <- c("marker_name",colnames(summary(cox.temp)$coefficients), "n","nevents")
  return (summary.cox)
}


cox.result.OS.1 <- cox.result.OS (time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.cat.incl)
cox.result.OS.1 
        
        t# cox.OS.incl <- survfit(cox.result.OS(matched.test.incl.pData, OS.cat.bin.incl, matched.goi.vsd.cat.incl))





### example below to delete after worked out cox graph function
km.log.test.OS <- function(time, event, marker, out.file = "none"){
  if(out.file!="none"){
    pdf(out.file)
  }
  km.OS.incl <- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
  plot(km.OS.incl, col = c("red", "blue"),xlab = "overall survival (years)", ylab = "OS", xlim = c(0,10), main = "Biomarker expression and overall survival (OS)",  lty = 1:2)
  OS.names <- c("biomarker - high", "biomarker - low")
  legend (x="topright", OS.names,  lty= 1:2, col = c("red","blue"))
  OS.incl.logrank <- survdiff(Surv(time, event) ~ marker)
  1 - pchisq(OS.incl.logrank$chisq, length(OS.incl.logrank$obs)-1) -> surv.p.val
  text(4,0.1,paste("p =",round(surv.p.val, 3)), pos = 4, cex = 1)
  if(out.file!="none"){
    dev.off()
  }
}



#####################################

### EFS: 

### change event-free status to binary
matched.test.pData$Event

EFS.cat.bin <- ifelse(matched.test.pData$Event == "Event", 1, 0)
EFS.cat.bin

chi.sq(as.factor(EFS.cat.bin), as.factor(matched.goi.vsd.cat))
### creating a table as km curve suggests that all in the biomarker high group had event. Ask Dan

# table.EFS<- table(as.factor(EFS.cat.bin), as.factor(matched.goi.vsd.cat))
# table.EFS.perc <- prop.table(table.EFS)*100


km.EFS <- survfit(Surv(matched.test.pData$EFS, EFS.cat.bin)~matched.goi.vsd.cat)

# km.EFS <- survfit(Surv(matched.test.pData$EFS, EFS.cat.bin==1)~matched.goi.vsd.cat) is the same as when not specifying status==1
# km.EFS

plot(km.EFS, col = c("red", "blue"), xlab = "time to event (years)", ylab = "EFS", main = "Biomarker expression and event-free survival (EFS)",  lty = 1:2)
EFS.names <- c("biomarker - high", "biomarker - low")
legend (15,1, EFS.names,  lty= 1:2)
# temp.logrank.pval <- 0.262
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))

EFS.logrank <- survdiff(Surv(matched.test.pData$EFS, EFS.cat.bin) ~ matched.goi.vsd.cat)
EFS.logrank
summary(EFS.logrank)

message ("will need to manually adjust the graph based on the log-rank p value")

plot(km.EFS, col = c("red", "blue"), xlab = "time to event (years)", ylab = "PFS", main = "Biomarker expression and event-free survival (EFS)",  lty = 1:2)
EFS.names <- c("biomarker - high", "biomarker - low")
legend (15,1, EFS.names,  lty= 1:2)
temp.logrank.pval <- 0.262
text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))


message ("restrict survival analysis for age 3-16 years")

EFS.cat.bin.incl <- ifelse(matched.test.incl.pData$Event == "Event", 1,0)
summary(EFS.cat.bin.incl)
km.EFS.incl <- survfit(Surv(matched.test.incl.pData$EFS, EFS.cat.bin.incl)~matched.goi.vsd.cat.incl)
km.EFS.incl

plot(km.EFS.incl, col = c("red", "blue"), xlab = "event-free survival (years)", ylab = "EFS", main = "Biomarker expression and event-free survival (EFS)",  lty = 1:2)
EFS.names.incl <- c("biomarker - high", "biomarker - low")
legend (15,1, EFS.names.incl,  lty= 1:2)
# temp.logrank.pval <- 0.67
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))

EFS.incl.logrank <- survdiff(Surv(matched.test.incl.pData$EFS, OS.cat.bin.incl) ~ matched.goi.vsd.cat.incl)
EFS.incl.logrank

### now go back in and add p value from EFS.incl.logrank into graph
message ("need to add p value from EFS.incl.logrank into graph")

plot(km.EFS.incl, col = c("red", "blue"), xlab = "event-free survival (years)", ylab = "EFS", main = "Biomarker expression and event-free survival (EFS)",  lty = 1:2)
EFS.names.incl <- c("biomarker - high", "biomarker - low")
legend (15,1, EFS.names.incl,  lty= 1:2)
temp.logrank.pval <- 0.543
text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))


### cox regression analysis
message("cox regression for EFS on entire cohort")

cox.EFS <- coxph (Surv(matched.test.pData$EFS, EFS.cat.bin) ~ matched.goi.vsd.cat)
summary(cox.EFS)$logtest
summary(cox.EFS)

### 

message("cox regression for EFS on age 3-16 years")

cox.EFS.incl <- coxph (Surv(matched.test.incl.pData$EFS, OS.cat.bin.incl) ~ matched.goi.vsd.cat.incl)
summary(cox.EFS.incl)$logtest
summary(cox.EFS.incl)





###################################

### multivariate analysis
### to use R+, M+, LCA, sex, MYCMYCN.cat, TP53 across entire group, restricted to age 3-16 yo
### cross reference with dataset of curative intent (CSI and max resection)


### for relapse

log.reg.multi <- glm(matched.test.incl.pData$relapse ~ matched.test.incl.pData$resection 
                    + matched.test.incl.pData$mstatus 
                    + matched.test.incl.pData$sex 
                    + matched.test.incl.pData$resection
                    + matched.test.incl.pData$histopath
                    + matched.test.incl.pData$MYC.cat
                    + matched.test.incl.pData$MYCN.cat
                    + matched.test.incl.pData$TP53.cat, family=binomial (link = 'logit'),
                    data = matched.test.incl.pData)

summary(log.reg.multi)

### compare to model with inclusion of biomarker

log.reg.multi.goi <- glm(matched.test.incl.pData$relapse ~ matched.test.incl.pData$resection 
                     + matched.test.incl.pData$mstatus 
                     + matched.test.incl.pData$sex 
                     + matched.test.incl.pData$resection
                     + matched.test.incl.pData$histopath
                     + matched.test.incl.pData$MYC.cat
                     + matched.test.incl.pData$MYCN.cat
                     + matched.test.incl.pData$TP53.cat
                     + matched.goi.vsd.cat.incl, 
                     family=binomial (link = 'logit'),
                     data = matched.test.incl.pData)

summary(log.reg.multi.goi)

### consider breaking down into subgroups, as the above may overwhelm the dataset with missing values

message ("processing subgroup specific biomarkers")

### Group 3 and Group 4
message ("create combined Group 3 Group 4 dataframe for age 3-16 yo")

G3.match <- matched.test.incl.pData$meth=="G3" 
G4.match <- matched.test.incl.pData$meth=="G4"

G3.match.df <- matched.test.incl.pData [G3.match, ]
G4.match.df <- matched.test.incl.pData [G4.match, ]

G3G4.match.df <- rbind(G3.match.df, G4.match.df)
View(G3G4.match.df)
nrow(G3G4.match.df)

index.incl <- match(names(goi.vsd), rownames(G3G4.match.df)) 
matched.G3G4.incl.pData <- G3G4.match.df[index.incl[!is.na(index.incl)],] 
is.vector(matched.G3G4.incl.pData)
matched.goi.vsd.G3G4.incl <- goi.vsd[!is.na(index.incl)] 
matched.goi.vsd.cat.G3G4.incl <- ifelse(matched.goi.vsd.G3G4.incl>median(goi.vsd, na.rm = T), "high","low")

# View(matched.G3G4.incl.pData)
       
### how would you grep a value in a dataframe
# G3G4 <- grep("G[34]", "G3G4")  




##########################################

### backward elimination to remove least significant variable at each stage, model with G3G4
### assign levels to relevant factors, below function does not work

MYC.G3G4.cat <- relevel (matched.G3G4.incl.pData$MYC.cat, ref="MYC non ampl")
mstatus.G3G4 <- relevel (matched.G3G4.incl.pData$mstatus, ref = "M0/M1")
resection.G3G4 <- relevel (matched.G3G4.incl.pData$resection, ref = "Gross total resection")
histopath.G3G4 <- relevel (matched.G3G4.incl.pData$histopath, ref = "FAV")
sex.G3G4 <- relevel (matched.G3G4.incl.pData$sex, ref = "female")
MYCN.G3G4.cat <- relevel (matched.G3G4.incl.pData$MYCN.cat, ref = "MYCN non ampl")

# relevel.G3G4 <- function (x,y){
# relevel.temp <- relevel (x, ref = "y")
# return (relevel.temp)
# }

# MYC.G3G4.cat <- relevel.G3G4(x = matched.G3G4.incl.pData$MYC.cat, y = "MYC non ampl")


log.reg.multi.G3G4.1 <- glm(matched.G3G4.incl.pData$relapse ~ resection.G3G4 
                              + mstatus.G3G4
                              + sex.G3G4
                              + resection.G3G4
                              + histopath.G3G4
                              + MYC.G3G4.cat
                              + MYCN.G3G4.cat,
                              family=binomial (link = 'logit'),
                              data = matched.G3G4.incl.pData)

summary(log.reg.multi.G3G4.1)
confint (log.reg.multi.G3G4.1)
exp(coef(log.reg.multi.G3G4.1))
exp(cbind(OR = coef(log.reg.multi.G3G4.1), confint(log.reg.multi.G3G4.1)))
    

### question: proceed in this manner, and do manually at the time? 


### then compare to the addition of the biomarker, perform backward elimination


log.reg.multi.G3G4.1A <- glm(matched.G3G4.incl.pData$relapse ~ resection.G3G4 
                            + mstatus.G3G4
                            + sex.G3G4
                            + resection.G3G4
                            + histopath.G3G4
                            + MYC.G3G4.cat
                            family=binomial (link = 'logit'),
                            data = matched.G3G4.incl.pData)

summary(log.reg.multi.G3G4.1A)
confint (log.reg.multi.G3G4.1A)
exp(coef(log.reg.multi.G3G4.1A))
exp(cbind(OR = coef(log.reg.multi.G3G4.1A), confint(log.reg.multi.G3G4.1A)))



log.reg.multi.G3G4.1B <- glm(matched.G3G4.incl.pData$relapse ~ matched.G3G4.incl.pData$resection 
                             + matched.G3G4.incl.pData$mstatus 
                             + matched.G3G4.incl.pData$sex 
                             + matched.G3G4.incl.pData$resection
                             + matched.G3G4.incl.pData$histopath
                             + matched.G3G4.incl.pData$MYC.cat,
                             family=binomial (link = 'logit'),
                             data = matched.G3G4.incl.pData)


summary(log.reg.multi.G3G4.1B)


### compare with addition of biomarker

log.reg.multi.goi.G3G4 <- glm(matched.G3G4.incl.pData$relapse ~ matched.G3G4.incl.pData$resection 
                         + matched.G3G4.incl.pData$mstatus 
                         + matched.G3G4.incl.pData$sex 
                         + matched.G3G4.incl.pData$resection
                         + matched.G3G4.incl.pData$histopath
                         + matched.G3G4.incl.pData$MYC.cat
                         + matched.goi.vsd.cat.G3G4.incl, 
                         family=binomial (link = 'logit'),
                         data = matched.G3G4.incl.pData)


### other subgroups e.g SHH TP53








############################
### multivariate cox regression analysis


time <- matched.test.incl.pData$PFS
# event <- ifelse(matched.test.incl.pData$relapse=="relapse", 1, 0)
event <- relapse.bin.incl
y <-matched.test.incl.pData$histopath
a <-matched.test.incl.pData$mstatus
b <-matched.test.incl.pData$resection
c <-matched.test.incl.pData$MYC.cat 
d <-matched.test.incl.pData$MYCN.cat

View(test.pData)

cox.result <- function (time, event, y, a, b, c, d)  
{
  cox.temp <- coxph (Surv(time, event)~y + a + b + c + d, data= matched.test.incl.pData)
  summary.cox <- summary(cox.temp)$logtest
  return (summary(summary.cox))
}



# summary(summary.cox)

### function not currently running
  

cox.result.1 <- cox.result (time = matched.test.incl.pData$PFS, event = relapse.bin.incl, 
                           y = matched.test.incl.pData$histopath,
                           a = matched.test.incl.pData$mstatus,
                           b = matched.test.incl.pData$resection,
                           c = matched.test.incl.pData$MYC.cat,
                           d = matched.test.incl.pData$MYCN.cat
                            #  )

#summary(cox.result.1)

### only input values below pre-defined p threshold for multivariate model
### let's use p<0.2 from the univariate logistic regression

cox.multi <- coxph (Surv(matched.test.incl.pData$PFS, relapse.bin.incl) ~ )
summary(cox.multi)$logtest
summary(cox.multi)


cox.multi.marker <- coxph (Surv(matched.test.incl.pData$PFS, relapse.bin.incl) ~ matched.goi.vsd.cat.incl)
summary(cox.multi)$logtest
summary(cox.multi)


### also compare the model to that with the addition of a single biomarker






###################################
### AUC and ROC analysis







sink()
####################################

### code from Dan below

matched.test.pData$PFS
matched.test.pData$relapse

View(matched.test.pData)
y <- matched.goi.vsd
y.cat <- matched.goi.vsd.cat
time <- matched.test.pData$PFS
event <-  matched.test.pData$relapse=="relapse"


km <- survfit(Surv(time, event)~y.cat)
temp.survdiff <- survdiff(Surv(time, event)~y.cat)
str(temp.survdiff)
plot(km, col = c("red","blue"), xlab = "time (years)", ylab = "PFS")
legend
summary(km)
summary(temp.survdiff)
cox.result <- coxph (Surv(time, event)~y.cat)
summary(cox.result)$logtest
summary(cox.result)




### to look at relationship between variables

### multivariate logistic regression and comparison to existing models

sink()
dev.off()




##############################
### troubleshooting tips

### example codes below

# pData <- read.table(file= "/home/nmm199/test/data/Test_data.txt",  header=T)
# pData <- read.table(file= "/home/nmm199/test/data/Test2.txt",  header=T, sep = "\t")
# read.table(file="/home/nmm199/test/data/Test2.txt", header = TRUE, sep="", quote="'", dec=".")



### if unsure of how function works, can try specifying the variables, and explore the temporary output.

# x = matched.test.pData$age.cat.infant
# y = matched.goi.vsd.cat
# chi.sq <- function(x,y){
#  table.temp <- table(x, y) 
#  table.temp.perc <- prop.table(table.temp)*100
#  summary.table(table.temp)
#  chi.test.temp <- chisq.test(table.temp) 
#  chi.test.temp.stat <- c(stat=chi.test.temp$statistic, p.value=chi.test.temp$p.value) 
#  chi.test.temp.res <- chi.test.temp$residuals
#  aheatmap(chi.test.temp.res, Rowv=NA, Colv = NA)
# list.temp <- list  (table.temp, 
#                     table.temp.perc,
#                     chi.test.temp.stat,
#                     chi.test.temp.res
#                      )

#  return(list.temp)
#}

# list.temp



### exploring the commands within the chi square function

# table.age.cat <- table(matched.test.pData$age.cat.infant, matched.goi.vsd.cat) 
# table.age.cat.perc <- prop.table(table.age.cat)*100
# summary.table(table.age.cat)
# chi.test.age.cat <- chisq.test(table.age.cat) 
# chi.test.age.cat.res <- c(stat=chi.test.age.cat$statistic, p.value=chi.test.age.cat$p.value) 



### relabelling heatmap if needed for output

# aheatmap(chi.test.age.cat$residuals, Rowv=NA, Colv = NA, labRow = "Infant", labCol="Expression", main="Correlations")
# help(aheatmap)

# aheatmap(chi.test.age.cat$residuals, Rowv=NA, Colv = NA, main = "Expression correlated to Infant status") 


### generate list, without the function above
# list.age <- list(age= table.age.cat, 
#             perc=table.age.cat.perc,
#            residual=chi.test.age.cat.res,
#             chi=chi.test.age.cat
#             )
# list.age

### subsetting list results for equivalent slice of object

# list.age[[1]]
# list.age$age



### other examples for t.test
# t.test(x~y, data=data, var.equal=T)
# chisq.test(tab, correct=F)
# fisher.test(tab)
# glm(formula = meth ~ age, family = binomial(link = "logit"), 
#  data = data)


########################################

### record of code that does not work

### error in trying to use age.cont with categorical biomarker or is use biomarker is dependent and age (categorical) as independent factor
# log.reg.age.cat2<- glm(matched.goi.vsd ~ matched.test.pData$age.cat.infant, family = binomial(link='logit'), data=matched.test.pData)

# log.reg.age.cont <- glm(age.cont ~ matched.goi.vsd.cat, family = binomial(link= 'logit'), data=matched.test.pData)
# matched.goi.vsd.cat.bin <- ifelse(matched.goi.vsd.cat=="high", 1, 0)
# matched.goi.vsd.cat.bin.2 <- matched.goi.vsd.cat.bin, na.rm=T
# also na.action = omit as another option to incorporate

# log.reg.age.cont <- glm (formula = matched.test.pData$age.cont ~ matched.goi.vsd.cat.bin , family = binomial(link='logit'), data=matched.test.pData)
# matched.goi.vsd.cat.bin <- as.factor(matched.goi.vsd.cat.bin)

# is.factor(matched.goi.vsd.cat.bin)


# log.reg <- function(x,y){
# glm(formula = x~y, family = binomial(link = 'logit'), data = matched.test.pData)
# summary.log.temp <- summary(log.reg.temp)
# aov.temp <- aov(x~y)
# sum.aov.temp <- summary(aov.temp) 
# list.log.reg <- list (summary.log.temp, 
#                       anova.temp,
#                       aov.temp)
# return(list.log.reg)
#}

# x <- age.cat.infant.bin
# y <- matched.goi.vsd

# age.cat.infant.bin<- ifelse(matched.test.pData$age.cat.infant, TRUE==1, 0)
# log.reg.age <- log.reg(x = matched.test.pData$age.cat.infant, y = matched.goi.vsd)
# log.reg.age.bin <- log.reg (x = age.cat.infant.bin, y = matched.goi.vsd)

### attempt with Dan below for records
# x <- matched.test.pData$age.cont
# y <- matched.goi.vsd.cat=="high"


# cat.variable.cont <- function(x, y){
# data <- data.frame(x,y)
#log.reg <- function(x,y){
#glm(formula = x~y, family = binomial(link = "logit"), data = data)
#anova(model, test="Chisq")

#aov(x~y) -> aov.out
#str(summary(aov.out))
#summary(aov.out) -> sum.aov.out
#sum.aov.out


# subgroup.meth<- function (x,y){
# temp.sub <- x =="temp.group"
# temp.group <- matched.test.pData [temp.sub, ]
# summary(temp.group)
# return(summary(temp.group))
# }

# G3.subgroup <- subgroup.meth(x = matched.test.pData$meth )
# G3.subgroup

### cox ph

### i think this code below is redundant 

# cox.relapse.df <- data.frame(matched.goi.vsd.cat, matched.test.pData$relapse, matched.test.pData$PFS)
# cox.relapse.df <- data.frame(matched.goi.vsd.cat, relapse.bin, matched.test.pData$PFS)
# surv.cox.rel <- survfit (cox.relapse, newdata=cox.relapse.df)
# summary(surv.cox.rel)                

######################

# message ("adding in cytogenetic 13q data") - superceded because 13q loss data was added earlier on to test.pData

# index.incl.cytogen <- match(rownames (Age.incl.df),row.names(cytogen.q13.cat [,1]))
# matched.test.incl.pData <- Age.incl.df[index.incl.cytogen[!is.na(index.incl.cytogen)],]



#######################
### redundant PFS code 



# km.PFS <- survfit(Surv(matched.test.pData$PFS, relapse.bin)~matched.goi.vsd.cat)
# km.PFS

# plot(km.PFS, col = c("red", "blue"), xlab = "time to progression/relapse (years)", ylab = "PFS", main = "Biomarker expression and progression-free survival (PFS)",  lty = 1:2)
# PFS.names <- c("biomarker - high", "biomarker - low")
# legend (12,1, PFS.names,  lty= 1:2)
# temp.logrank.pval <- 0.494
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))

 #PFS.logrank <- survdiff(Surv(matched.test.pData$PFS, relapse.bin) ~ matched.goi.vsd.cat)
# PFS.logrank

# message ("will need to manually adjust the graph based on the log-rank p value")

# plot(km.PFS, col = c("red", "blue"), xlab = "time to progression/relapse (years)", ylab = "PFS", main = "Biomarker expression and progression-free survival (PFS)",  lty = 1:2)
# PFS.names <- c("biomarker - high", "biomarker - low")
# legend (15,1, PFS.names,  lty= 1:2)
# temp.logrank.pval <- 0.494
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))
### now go back in and add p value from PFS.incl.logrank into graph
# message ("need to add p value from PFS.incl.logrank into graph")

#plot(km.PFS.incl, col = c("red", "blue"), xlab = "time to progression/relapse (years)", ylab = "PFS", main = "Biomarker expression and progression-free survival (PFS)",  lty = 1:2)
#PFS.names <- c("biomarker - high", "biomarker - low")
#legend (15,1, PFS.names,  lty= 1:2)
# temp.logrank.pval <- 0.016
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))



###############################
### redundant OS code

# km.OS <- survfit(Surv(matched.test.pData$Followup, OS.cat.bin)~matched.goi.vsd.cat)
# km.OS

# plot(km.OS, col = c("red", "blue"), xlab = "time to progression/relapse (years)", ylab = "OS", main = "Biomarker expression and overall survival (OS)",  lty = 1:2)
# OS.names <- c("biomarker - high", "biomarker - low")
# legend (15,1, OS.names,  lty= 1:2)
# temp.logrank.pval <- 0.494
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))


# km.log.test(matched.test.incl.pData$Followup,OS.cat.bin, matched.goi.vsd.cat.incl )

# OS.logrank <- survdiff(Surv(matched.test.pData$Followup, OS.cat.bin) ~ matched.goi.vsd.cat)
# OS.logrank
# summary(OS.logrank)

# message ("will need to manually adjust the graph based on the log-rank p value")

# plot(km.OS, col = c("red", "blue"), xlab = "time to progression/relapse (years)", ylab = "PFS", main = "Biomarker expression and progression-free survival (PFS)",  lty = 1:2)
# OS.names <- c("biomarker - high", "biomarker - low")
 #legend (15,1, PFS.names,  lty= 1:2)
#temp.logrank.pval <- 0.331
#text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))


 #message ("restrict survival analysis for age 3-16 years")



 #plot(km.OS.incl, col = c("red", "blue"), xlab = "overall survival (years)", ylab = "OS", main = "Biomarker expression and overall survival (PFS)",  lty = 1:2)
#OS.names.incl <- c("biomarker - high", "biomarker - low")
#legend (15,1, OS.names.incl,  lty= 1:2)
# temp.logrank.pval <- 0.67
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))

 #OS.incl.logrank <- survdiff(Surv(matched.test.incl.pData$Followup, OS.cat.bin.incl) ~ matched.goi.vsd.cat.incl)
 #OS.incl.logrank

### now go back in and add p value from OS.incl.logrank into graph
# message ("need to add p value from OS.incl.logrank into graph")

# plot(km.OS.incl, col = c("red", "blue"), xlab = "overall survival (years)", ylab = "OS", main = "Biomarker expression and overall survival (PFS)",  lty = 1:2)
# OS.names.incl <- c("biomarker - high", "biomarker - low")
# legend (15,1, OS.names.incl,  lty= 1:2)
# temp.logrank.pval <- 0.625
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))
