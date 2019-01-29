### Run univariable logistic and cox regression in test.pData to determine which factors are significant in univariate analysis
### input: test.pData
### output: list of factors and p values and adjusted p values

### Author: Dr Marion Mateos
### Date: October 26 2018

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_master.R")

### loading in clinical data object = test.pData
load("/home/nmm199/R/MB_RNAseq/Clinical/test.pData")

### inputs
### goi (MYC) 
### gp.filt.mb.vsd object ### need to generate this from the clinical_data_master.R file
### run the clinPath Assess function within clinical_data_functions_master.R to generate matched.test.pData or use the script below

goi <- "ENSG00000136997" ### MYC

goi.vsd <- as.numeric(gp.filt.mb.vsd[goi,]) 

names(goi.vsd) <- names(gp.filt.mb.vsd) 


#################################################################################################

index <- match(names(goi.vsd), rownames(test.pData)) 
matched.test.pData <- test.pData[index[!is.na(index)],] 
is.vector(matched.test.pData)
#as.data.frame(matched.test.pData) ### added 070917 in attempt to avoid downstream" $ atomic in a vector" error
matched.goi.vsd <- goi.vsd[!is.na(index)] 
matched.goi.vsd.cat <- ifelse(matched.goi.vsd>median(goi.vsd, na.rm = T), "high","low") 

#################################################################################################

names(matched.test.pData)
summary.217 <- as.data.frame(summary(matched.test.pData)) ### clinical characteristics of the 217 patients with matched clinical data and RNA

# write.csv (summary.217, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/summary.217.csv")

#################################################################################################
### can read in the matched.test.pData object using MYC as the goi if have not generated it above

# load (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/summary.217.csv")

### generate WNT vs non WNT

matched.test.pData$WNT_PNET5<- ifelse(matched.test.pData$meth == "WNT", "WNT", "non-WNT")

### in order to create levels in logistic regression, to set 'reference level' for risk (odds ratio), variables must be factors ### Jan 18 2019
as.factor(names(matched.test.pData))
matched.test.pData$WNT <- as.factor(matched.test.pData$WNT_PNET5)

# head(matched.test.pData$CSI) ### check relevant levels

### define reference level for each
matched.test.pData$mstatus <- relevel (matched.test.pData$mstatus, ref = "M0/M1" )
matched.test.pData$LCA <- relevel (matched.test.pData$LCA, ref = "non LCA")
matched.test.pData$MYC.cat <- relevel (matched.test.pData$MYC.cat, ref = "MYC non ampl")
matched.test.pData$MYCN.cat<- relevel (matched.test.pData$MYCN.cat, ref = "MYCN non ampl" )
matched.test.pData$MYCMYCN.cat <- relevel (matched.test.pData$MYCMYCN.cat, ref = "MYC MYCN non ampl" )
matched.test.pData$WNT <- relevel (matched.test.pData$WNT, ref = "WNT")
matched.test.pData$TP53.cat <- relevel (matched.test.pData$TP53.cat, ref = "TP53 WT")

#### run logistic regression for OS and PFS

logreg.pData.sex <- logisticRegression(x=matched.test.pData$OS.cat, y = matched.test.pData$sex, data = matched.test.pData)

logreg.pData.sex.PFS <- logisticRegression(x =matched.test.pData$relapse, y= matched.test.pData$sex, data = matched.test.pData)


try(logreg.pData.mstatus <- logisticRegression(x= matched.test.pData$OS.cat, y= matched.test.pData$mstatus, data=matched.test.pData), silent= T)

try(logreg.pData.mstatus.PFS <- logisticRegression(x= matched.test.pData$relapse, y= matched.test.pData$mstatus, data=matched.test.pData), silent= T)


try(logreg.pData.resection <- logisticRegression(x =matched.test.pData$OS.cat, y = matched.test.pData$resection, data=matched.test.pData), silent= T)

try(logreg.pData.resection.PFS <- logisticRegression(x =matched.test.pData$relapse, y = matched.test.pData$resection, data=matched.test.pData), silent= T)


try(logreg.pData.LCA <- logisticRegression(x= matched.test.pData$OS.cat, y = matched.test.pData$LCA, data=matched.test.pData), silent= T)

try(logreg.pData.LCA.PFS <- logisticRegression(x= matched.test.pData$relapse, y = matched.test.pData$LCA, data=matched.test.pData), silent= T)


try(logreg.pData.MYC <- logisticRegression (x= matched.test.pData$OS.cat, y= matched.test.pData$MYC.cat,  data=matched.test.pData), silent= T)

try(logreg.pData.MYC.PFS <- logisticRegression (x= matched.test.pData$relapse, y= matched.test.pData$MYC.cat,  data=matched.test.pData), silent= T)



try(logreg.pData.MYCN <- logisticRegression (x=matched.test.pData$OS.cat , y=matched.test.pData$MYCN.cat, data=matched.test.pData), silent= T)

try(logreg.pData.MYCN.PFS <- logisticRegression (x=matched.test.pData$relapse , y=matched.test.pData$MYCN.cat, data=matched.test.pData), silent= T)


try(logreg.pData.MYCMYCN <- logisticRegression (x=matched.test.pData$OS.cat, y= matched.test.pData$MYCMYCN.cat, data=matched.test.pData), silent= T)

try(logreg.pData.MYCMYCN.PFS <- logisticRegression (x=matched.test.pData$relapse, y= matched.test.pData$MYCMYCN.cat, data=matched.test.pData), silent= T)


try(logreg.pData.meth <- logisticRegression (x= matched.test.pData$OS.cat, y= matched.test.pData$meth, data = matched.test.pData), silent= T)

try(logreg.pData.meth.PFS <- logisticRegression (x= matched.test.pData$relapse, y= matched.test.pData$meth, data = matched.test.pData), silent= T)


try(logreg.pData.meth7 <- logisticRegression (x = matched.test.pData$OS.cat, y= matched.test.pData$meth7,  data = matched.test.pData), silent= T)

try(logreg.pData.meth7.PFS <- logisticRegression (x = matched.test.pData$relapse, y= matched.test.pData$meth7,  data = matched.test.pData), silent= T)


logreg.pData.q13loss <- logisticRegression(x=matched.test.pData$OS.cat, y= matched.test.pData$q13loss, data = matched.test.pData)

logreg.pData.q13loss.PFS <- logisticRegression(x= matched.test.pData$relapse, y= matched.test.pData$q13loss, data = matched.test.pData)

  
try(logreg.pData.TP53 <- logisticRegression (x = matched.test.pData$OS.cat , y= matched.test.pData$TP53.cat, data=matched.test.pData), silent= T)

try(logreg.pData.TP53.PFS <- logisticRegression (x = matched.test.pData$relapse , y= matched.test.pData$TP53.cat, data=matched.test.pData), silent= T)


try(logreg.pData.TERT <- logisticRegression (x= matched.test.pData$OS.cat, y= matched.test.pData$TERT.cat, data = matched.test.pData), silent= T)

try(logreg.pData.TERT.PFS <- logisticRegression (x= matched.test.pData$relapse, y= matched.test.pData$TERT.cat, data = matched.test.pData), silent= T)


logreg.pData.WNT <- logisticRegression(x = matched.test.pData$OS.cat, y = matched.test.pData$WNT, data = matched.test.pData)

logreg.pData.WNT.PFS <- logisticRegression(x = matched.test.pData$relapse, y = matched.test.pData$WNT, data = matched.test.pData)

### generate list of results

logistic.reg.results <- as.list(mget(ls(pattern="logreg.pData"))) 


#################################################################################################

### Determine which results are significant in G3G4
### make dataframe

summary(matched.test.pData$meth)

G3G4.matched.pData <- matched.test.pData[(matched.test.pData$meth=="G3"| matched.test.pData$meth =="G4"),  ] 

G3G4.matched.pData <- G3G4.matched.pData[!is.na(G3G4.matched.pData$meth), ] ### this is how to remove NA, n=136 G3G4 patients

### generate G3G4_HR variable

G3G4.matched.pData$G3G4_HR <- ifelse (G3G4.matched.pData$meth7.cat=="Grp4_HighRisk"|G3G4.matched.pData$meth7.cat=="Grp3_HighRisk", "G3G4_HR", "G3G4_LR")

# G3G4.matched.pData$G3G4_HR
# summary(G3G4.matched.pData$meth)
# summary(G3G4.matched.pData$meth7.cat)

### transform variables into factors to allow levels to be stipulated
## need to run above relevel commands on matched.test.pData to allow follow through of reference variables 

as.factor(names(G3G4.matched.pData))
G3G4.matched.pData$G3G4_HR <- relevel(as.factor(G3G4.matched.pData$G3G4_HR), ref = "G3G4_LR")
# G3G4.matched.pData$q13loss ### default is "no q13 loss" as reference


### univariate logistic regression associations

logreg.G3G4.sex <- logisticRegression(x=G3G4.matched.pData$OS.cat, y = G3G4.matched.pData$sex, data = G3G4.matched.pData)

logreg.G3G4.sex.PFS <- logisticRegression(x = G3G4.matched.pData$relapse, y= G3G4.matched.pData$sex, data = G3G4.matched.pData)

try(logreg.G3G4.mstatus <- logisticRegression(x= G3G4.matched.pData$OS.cat, y= G3G4.matched.pData$mstatus, data=G3G4.matched.pData), silent= T)

try(logreg.G3G4.mstatus.PFS <- logisticRegression(x= G3G4.matched.pData$relapse, y= G3G4.matched.pData$mstatus, data=G3G4.matched.pData), silent= T)


try(logreg.G3G4.resection <- logisticRegression(x =G3G4.matched.pData$OS.cat, y = G3G4.matched.pData$resection, data=G3G4.matched.pData), silent= T)

try(logreg.G3G4.resection.PFS <- logisticRegression(x =G3G4.matched.pData$relapse, y = G3G4.matched.pData$resection, data=G3G4.matched.pData), silent= T)


try(logreg.G3G4.LCA <- logisticRegression(x= G3G4.matched.pData$OS.cat, y = G3G4.matched.pData$LCA, data=G3G4.matched.pData), silent= T)

try(logreg.G3G4.LCA.PFS <- logisticRegression(x= G3G4.matched.pData$relapse, y = G3G4.matched.pData$LCA, data=G3G4.matched.pData), silent= T)


try(logreg.G3G4.MYC <- logisticRegression (x= G3G4.matched.pData$OS.cat, y= G3G4.matched.pData$MYC.cat,  data=G3G4.matched.pData), silent= T)

try(logreg.G3G4.MYC.PFS <- logisticRegression (x= G3G4.matched.pData$relapse, y= G3G4.matched.pData$MYC.cat,  data=G3G4.matched.pData), silent= T)



try(logreg.G3G4.MYCN <- logisticRegression (x=G3G4.matched.pData$OS.cat , y=G3G4.matched.pData$MYCN.cat, data=G3G4.matched.pData), silent= T)

try(logreg.G3G4.MYCN.PFS <- logisticRegression (x=G3G4.matched.pData$relapse , y=G3G4.matched.pData$MYCN.cat, data=G3G4.matched.pData), silent= T)


try(logreg.G3G4.MYCMYCN <- logisticRegression (x=G3G4.matched.pData$OS.cat, y= G3G4.matched.pData$MYCMYCN.cat, data=G3G4.matched.pData), silent= T)

try(logreg.G3G4.MYCMYCN.PFS <- logisticRegression (x=G3G4.matched.pData$relapse, y= G3G4.matched.pData$MYCMYCN.cat, data=G3G4.matched.pData), silent= T)


try(logreg.G3G4.meth <- logisticRegression (x= G3G4.matched.pData$OS.cat, y= G3G4.matched.pData$meth, data = G3G4.matched.pData), silent= T)

try(logreg.G3G4.meth.PFS <- logisticRegression (x= G3G4.matched.pData$relapse, y= G3G4.matched.pData$meth, data = G3G4.matched.pData), silent= T)


try(logreg.G3G4.meth7 <- logisticRegression (x = G3G4.matched.pData$OS.cat, y= G3G4.matched.pData$meth7,  data = G3G4.matched.pData), silent= T)

try(logreg.G3G4.meth7.PFS <- logisticRegression (x = G3G4.matched.pData$relapse, y= G3G4.matched.pData$meth7,  data = G3G4.matched.pData), silent= T)


try(logreg.G3G4.TP53 <- logisticRegression (x = G3G4.matched.pData$OS.cat , y= G3G4.matched.pData$TP53.cat, data=G3G4.matched.pData), silent= T)

try(logreg.G3G4.TP53.PFS <- logisticRegression (x = G3G4.matched.pData$relapse , y=G3G4.matched.pData$TP53.cat, data=G3G4.matched.pData), silent= T)


try(logreg.G3G4.TERT <- logisticRegression (x= G3G4.matched.pData$OS.cat, y= G3G4.matched.pData$TERT.cat, data = G3G4.matched.pData), silent= T)

try(logreg.G3G4.TERT.PFS <- logisticRegression (x= G3G4.matched.pData$relapse, y= G3G4.matched.pData$TERT.cat, data = G3G4.matched.pData), silent= T)


logreg.G3G4.q13loss <- logisticRegression(x = G3G4.matched.pData$OS.cat, y = G3G4.matched.pData$q13loss, data = G3G4.matched.pData)


logreg.G3G4.q13loss.PFS <- logisticRegression(x = G3G4.matched.pData$relapse, y = G3G4.matched.pData$q13loss, data = G3G4.matched.pData)


logreg.G3G4.G3G4HR <- logisticRegression(x = G3G4.matched.pData$OS.cat, y = G3G4.matched.pData$G3G4_HR, data = G3G4.matched.pData)

logreg.G3G4.G3G4HR.PFS <- logisticRegression(x = G3G4.matched.pData$relapse, y = G3G4.matched.pData$G3G4_HR, data = G3G4.matched.pData)


### then list of results logreg.G3G4

logistic.reg.G3G4 <- as.list(mget(ls(pattern="logreg.G3G4"))) 

#################################################################################################
#################################################################################################
### analyse SHH factors

SHH.matched.pData <- matched.test.pData[(matched.test.pData$meth=="SHH"), ]

SHH.matched.pData <- SHH.matched.pData[!is.na(SHH.matched.pData$meth), ]

SHH.old.matched.pData <- SHH.matched.pData[SHH.matched.pData$meth7.cat=="SHH_Old", ]

### checking mstatus, resection, LCA, MYCN, TP53
# x <- SHH.matched.pData$OS.cat ### OS variable
# y <- SHH.matched.pData$sex ###"variable of interest" 
# z  <- SHH.matched.pData$relapse ### PFS variable
# data <- SHH.matched.pData

logreg.assoc <- function (x,y,z,data){
  log.reg.temp <- logisticRegression (x = x, y = y, data= data)
  log.reg.temp.PFS <- logisticRegression (x = z, y = y, data = data)
  log.reg.list <- list(OS = log.reg.temp, PFS = log.reg.temp.PFS)
  return(log.reg.list)
}

log.reg.SHH.sex <- logreg.assoc(x =  SHH.matched.pData$OS.cat, y = SHH.matched.pData$sex, z = SHH.matched.pData$relapse, data = SHH.matched.pData)

log.reg.SHH.mstatus <- logreg.assoc (x = SHH.matched.pData$OS.cat, y = SHH.matched.pData$mstatus, z = SHH.matched.pData$relapse, data = SHH.matched.pData)

log.reg.SHH.resection <- logreg.assoc (x = SHH.matched.pData$OS.cat, y = SHH.matched.pData$resection, z = SHH.matched.pData$relapse, data = SHH.matched.pData)

log.reg.SHH.LCA <- logreg.assoc (x = SHH.matched.pData$OS.cat, y = SHH.matched.pData$LCA, z = SHH.matched.pData$relapse, data = SHH.matched.pData)

### in SHH_old (SHH_child) group

log.reg.SHHold.mstatus <- logreg.assoc (x = SHH.old.matched.pData$OS.cat, y = SHH.old.matched.pData$mstatus, z = SHH.old.matched.pData$relapse, data = SHH.old.matched.pData)

log.reg.SHHold.resection <- logreg.assoc (x = SHH.old.matched.pData$OS.cat, y = SHH.old.matched.pData$resection, z = SHH.old.matched.pData$relapse, data = SHH.old.matched.pData)

log.reg.SHHold.LCA <- logreg.assoc (x = SHH.old.matched.pData$OS.cat, y = SHH.old.matched.pData$LCA, z = SHH.old.matched.pData$relapse, data = SHH.old.matched.pData)

log.reg.SHHold.MYCN <- logreg.assoc (x = SHH.old.matched.pData$OS.cat, y = SHH.old.matched.pData$MYCN.cat, z = SHH.old.matched.pData$relapse, data = SHH.old.matched.pData)

log.reg.SHH.TP53 <- logreg.assoc (x = SHH.old.matched.pData$OS.cat, y = SHH.old.matched.pData$TP53.cat, z = SHH.old.matched.pData$relapse, data = SHH.old.matched.pData)


### only need to run up to here 22/1/19

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################

### New dataframes for survival analyses to see effect of clinical variables in survival cohort 22/1/19 n=145. 
### THIS WAS NOT USED FOR THE FINAL ANALYSIS

### restrict analysis to age 3-16 yo, treated with curative intent including CSI

#cat ("first create dataframe for age 3-16 years, age.incl.df", sep = "\n") 

Age.incl <- matched.test.pData$age.filter== "TRUE"
Age.incl.df <- matched.test.pData [Age.incl,]
#summary(Age.incl.df)

#cat ("comparing with previous data frames for accuracy", sep = "\n")
#summary(test.pData)
#summary (matched.test.pData)

####################################

#cat ("processing summary stats for 3-16 yo for 4 molecular subgroups", sep = "\n")

### G3, aged 3-16 years
G3.incl <- Age.incl.df$meth =="G3"
G3.group.incl <- Age.incl.df [G3.incl, ]
#summary(G3.group.incl)
#nrow(G3.group.incl)

### G4, aged 3-16 years
G4.incl <- Age.incl.df$meth =="G4"
G4.group.incl <- Age.incl.df [G4.incl, ]
#summary(G4.group.incl)
#nrow(G4.group.incl)

### SHH, aged 3-16 years
SHH.incl <- Age.incl.df$meth == "SHH"
SHH.group.incl <- Age.incl.df [SHH.incl, ]
#summary (SHH.group.incl)
#nrow(SHH.group.incl)

### WNT, aged 3-16 years
WNT.incl <- Age.incl.df$meth == "WNT"
WNT.group.incl <- Age.incl.df [WNT.incl, ]
#summary (WNT.group.incl)
#nrow (WNT.group.incl)


### defining features of 7 molecular groups, group 3 and 4 subgroups

#cat ("processing summary stats for 7 molecular subgroup data, age 3-16 years", sep = "\n")

G3.high <- Age.incl.df$meth7 == "Grp3_HighRisk"
G3.high.incl <- Age.incl.df [G3.high, ]
#summary(G3.high.incl)

G3.low <- Age.incl.df$meth7 == "Grp3_LowRisk"
G3.low.incl <- Age.incl.df [G3.low, ]
#summary(G3.low.incl)


G4.high <- Age.incl.df$meth7 =="Grp4_HighRisk"
G4.high.incl <- Age.incl.df [G4.high, ]
#summary (G4.high.incl)

G4.low <- Age.incl.df$meth7 =="Grp4_LowRisk"
G4.low.incl <- Age.incl.df [G4.low, ]
#summary (G4.low.incl)

SHH.inf <- Age.incl.df$meth7 == "SHH_Inf"
SHH.inf.incl <- Age.incl.df [SHH.inf, ]
#summary (SHH.inf.incl)

SHH.old <- Age.incl.df$meth7 =="SHH_Old"
SHH.old.incl <- Age.incl.df [SHH.old, ]
#summary (SHH.old.incl)

##########################

### survival analysis using functions from this source file "clinical_data_functions_master.R" (i.e this file)

##########################

### creating dataframe for survival analysis
#cat ("restrict survival analysis for age 3-16 years, curative intent", sep = "\n")

###  matched.test.incl.pData is the dataframe that contains survival group for further analysis

index.incl <- match(names(goi.vsd), rownames(Age.incl.df)) 
matched.test.incl.pData.prelim <- Age.incl.df[index.incl[!is.na(index.incl)],] 
matched.goi.vsd.incl.prelim <- goi.vsd[!is.na(index.incl)] ### is this an error 31/10/18

###  include curative Yes/No as a filter 30/1/18 
### then create same size matched dataframes, then rerun the script on the server 6/2/18

matched.test.incl.pData <- matched.test.incl.pData.prelim[which(matched.test.incl.pData.prelim$curative =="Yes"),] ### this yields n=145 (primary database n=166, 16/1/18), this account for RNA drop off and is good to continue
index.final.incl <- match(names(matched.goi.vsd.incl.prelim), rownames(matched.test.incl.pData))
matched.goi.vsd.incl <- matched.goi.vsd.incl.prelim[!is.na(index.final.incl)]
matched.goi.vsd.cat.incl<- ifelse(matched.goi.vsd.incl>median(goi.vsd, na.rm = T), "high","low")

### G3 and G4 combined dataframe

#cat ("creating combined dataframe to assess biomarker in G3 G4 combined group, for survival cohort, aged 3-16 years, curative intent", sep = "\n")

G3.match <- matched.test.incl.pData$meth=="G3" 
G4.match <- matched.test.incl.pData$meth=="G4"

G3.match.df <- matched.test.incl.pData [G3.match, ]
G4.match.df <- matched.test.incl.pData [G4.match, ]

G3G4.match.df <- rbind(G3.match.df, G4.match.df)

nrow(G3G4.match.df)

index.incl <- match(names(goi.vsd), rownames(G3G4.match.df)) 
matched.G3G4.incl.pData <- G3G4.match.df[index.incl[!is.na(index.incl)],] 
is.vector(matched.G3G4.incl.pData) # unhashed 25/10/18
matched.goi.vsd.G3G4.incl <- goi.vsd[!is.na(index.incl)] 
matched.goi.vsd.cat.G3G4.incl <- ifelse(matched.goi.vsd.G3G4.incl>median(goi.vsd, na.rm = T), "high","low")

#summary(matched.G3G4.incl.pData$meth7)

### creating G3 and G4 high risk vs G3/G4 low risk for the subsequent multivariate Schwalbe modelling

matched.G3G4.incl.pData$G3G4.high.incl <- ifelse ((matched.G3G4.incl.pData$meth7.cat== "Grp3_HighRisk"| matched.G3G4.incl.pData$meth7.cat== "Grp4_HighRisk"), "G3G4_HighRisk", "G3G4_LowRisk")  
#######################################################################################################
### creating matched G4 dataframe (added 25/10/18) to allow KM survival analysis on selected biomarkers 

G4.index.incl <- match(names(goi.vsd), rownames(G4.match.df)) 
matched.G4.incl.pData <- G4.match.df[G4.index.incl[!is.na(G4.index.incl)],] 
matched.goi.vsd.G4.incl <- goi.vsd[!is.na(G4.index.incl)] 
matched.goi.vsd.cat.G4.incl <- ifelse(matched.goi.vsd.G4.incl>median(goi.vsd, na.rm = T), "high","low")
# length(matched.goi.vsd.cat.G4.incl) ### for numbers with adequate expression data for goi.vsd

#######################################################################################################
### creating SHH group in matched survival cohort

SHH.group.match <- matched.test.incl.pData$meth7.cat == "SHH_Old" | matched.test.incl.pData$meth7.cat == "SHH_Inf"
SHH.group.match.df <- matched.test.incl.pData[SHH.group.match, ]
SHH.index.incl <- match(names(goi.vsd), rownames(SHH.group.match.df))
matched.SHH.incl.pData <- SHH.group.match.df[SHH.index.incl[!is.na(SHH.index.incl)], ]

### creating matched goi data

matched.goi.vsd.SHH.incl <- goi.vsd[!is.na(SHH.index.incl)] 
matched.goi.vsd.cat.SHH.incl <- ifelse(matched.goi.vsd.SHH.incl>median(goi.vsd, na.rm = T), "high","low")

### creating SHH_Old group as continuous and categorical

SHH.old.match <- matched.test.incl.pData$meth7.cat == "SHH_Old" ### correlates to SHH_child cohort only (n=20) (does not include SHH_infant, n=8)
SHH.old.match.df <- matched.test.incl.pData [SHH.old.match, ]
SHH.old.index.incl <- match (names(goi.vsd), rownames(SHH.old.match.df))
matched.SHH.old.incl.pData <- SHH.old.match.df[SHH.old.index.incl[!is.na(SHH.old.index.incl)], ]

### creating matched goi data as continuous and categorical
matched.goi.vsd.SHH.old.incl <- goi.vsd[!is.na(SHH.old.index.incl)] 
matched.goi.vsd.cat.SHH.old.incl <- ifelse(matched.goi.vsd.SHH.old.incl>median(goi.vsd, na.rm = T), "high","low")

### WNT (not included in original analysis, for descriptive purposes only 19/6/18)
# WNT.match <- matched.test.incl.pData$meth=="WNT" 
# WNT.match.df <- matched.test.incl.pData [WNT.match, ]
# WNT.incl <- matched.test.incl.pData$meth =="WNT"
# WNT.incl.df <- matched.test.incl.pData[WNT.incl, ]

### Determining numbers for the dataset 19/6/18

# levels(matched.test.incl.pData$meth)
# matched.data.G3 <- matched.test.incl.pData [matched.test.incl.pData$meth =="G3", ]
# matched.data.G3.na <- na.omit(matched.data.G3$meth) 
# length(matched.data.G3.na)


#######################################################################################################

### Generating log regression results for survival cohort only 21/1/19

logreg.incl.sex <- logisticRegression(x=matched.test.incl.pData$OS.cat, y = matched.test.incl.pData$sex, data = matched.test.incl.pData)

logreg.incl.sex.PFS <- logisticRegression(x =matched.test.incl.pData$relapse, y= matched.test.incl.pData$sex, data = matched.test.incl.pData)


try(logreg.incl.mstatus <- logisticRegression(x= matched.test.incl.pData$OS.cat, y= matched.test.incl.pData$mstatus, data=matched.test.incl.pData), silent= T)

try(logreg.incl.mstatus.PFS <- logisticRegression(x= matched.test.pData$relapse, y= matched.test.pData$mstatus, data=matched.test.pData), silent= T)


try(logreg.incl.resection <- logisticRegression(x =matched.test.incl.pData$OS.cat, y = matched.test.incl.pData$resection, data=matched.test.incl.pData), silent= T)

try(logreg.incl.resection.PFS <- logisticRegression(x =matched.test.incl.pData$relapse, y = matched.test.incl.pData$resection, data=matched.test.incl.pData), silent= T)


try(logreg.incl.LCA <- logisticRegression(x= matched.test.incl.pData$OS.cat, y = matched.test.incl.pData$LCA, data=matched.test.incl.pData), silent= T)

try(logreg.incl.LCA.PFS <- logisticRegression(x= matched.test.incl.pData$relapse, y = matched.test.incl.pData$LCA, data=matched.test.incl.pData), silent= T)


try(logreg.incl.MYC <- logisticRegression (x= matched.test.incl.pData$OS.cat, y= matched.test.incl.pData$MYC.cat,  data=matched.test.incl.pData), silent= T)

try(logreg.incl.MYC.PFS <- logisticRegression (x= matched.test.incl.pData$relapse, y= matched.test.incl.pData$MYC.cat,  data=matched.test.incl.pData), silent= T)



try(logreg.incl.MYCN <- logisticRegression (x=matched.test.incl.pData$OS.cat , y=matched.test.incl.pData$MYCN.cat, data=matched.test.incl.pData), silent= T)

try(logreg.incl.MYCN.PFS <- logisticRegression (x=matched.test.incl.pData$relapse , y=matched.test.incl.pData$MYCN.cat, data=matched.test.incl.pData), silent= T)


try(logreg.incl.MYCMYCN <- logisticRegression (x=matched.test.incl.pData$OS.cat, y= matched.test.incl.pData$MYCMYCN.cat, data=matched.test.incl.pData), silent= T)

try(logreg.incl.MYCMYCN.PFS <- logisticRegression (x=matched.test.incl.pData$relapse, y= matched.test.incl.pData$MYCMYCN.cat, data=matched.test.incl.pData), silent= T)


try(logreg.incl.meth <- logisticRegression (x= matched.test.incl.pData$OS.cat, y= matched.test.incl.pData$meth, data = matched.test.incl.pData), silent= T)

try(logreg.incl.meth.PFS <- logisticRegression (x= matched.test.incl.pData$relapse, y= matched.test.incl.pData$meth, data = matched.test.incl.pData), silent= T)


try(logreg.incl.meth7 <- logisticRegression (x = matched.test.incl.pData$OS.cat, y= matched.test.incl.pData$meth7,  data = matched.test.incl.pData), silent= T)

try(logreg.incl.meth7.PFS <- logisticRegression (x = matched.test.incl.pData$relapse, y= matched.test.incl.pData$meth7,  data = matched.test.incl.pData), silent= T)


logreg.incl.q13loss <- logisticRegression(x=matched.test.incl.pData$OS.cat, y= matched.test.incl.pData$q13loss, data = matched.test.incl.pData)

logreg.incl.q13loss.PFS <- logisticRegression(x= matched.test.incl.pData$relapse, y= matched.test.incl.pData$q13loss, data = matched.test.pData)


try(logreg.incl.TP53 <- logisticRegression (x = matched.test.incl.pData$OS.cat , y= matched.test.incl.pData$TP53.cat, data=matched.test.incl.pData), silent= T)

try(logreg.incl.TP53.PFS <- logisticRegression (x = matched.test.incl.pData$relapse , y= matched.test.incl.pData$TP53.cat, data=matched.test.incl.pData), silent= T)


try(logreg.incl.TERT <- logisticRegression (x= matched.test.incl.pData$OS.cat, y= matched.test.incl.pData$TERT.cat, data = matched.test.incl.pData), silent= T)

try(logreg.incl.TERT.PFS <- logisticRegression (x= matched.test.incl.pData$relapse, y= matched.test.incl.pData$TERT.cat, data = matched.test.incl.pData), silent= T)


logreg.incl.WNT <- logisticRegression(x = matched.test.incl.pData$OS.cat, y = matched.test.incl.pData$WNT, data = matched.test.incl.pData)

logreg.incl.WNT.PFS <- logisticRegression(x = matched.test.incl.pData$relapse, y = matched.test.incl.pData$WNT, data = matched.test.incl.pData)
