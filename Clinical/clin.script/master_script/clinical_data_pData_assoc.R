### Run univariable logistic and cox regression in test.pData to determine which factors are significant in univariate analysis
### input: test.pData
### output: list of factors and p values and adjusted p values

### Author: Dr Marion Mateos
### Date: October 26 2018

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_master.R")

### loading in clinical data object = test.pData
load("/home/nmm199/R/MB_RNAseq/Clinical/test.pData")

### run the goi (MYC) and clinPath Assess function to generate matched.test.pData or can use the script below:

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

write.csv (summary.217, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/summary.217.csv")

#################################################################################################
### can read in the matched.test.pData object using MYC as the goi if have not generated it above

# load (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/summary.217.csv")

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


try(logreg.pData.TP53 <- logisticRegression (x = matched.test.pData$OS.cat , y= matched.test.pData$TP53.cat, data=matched.test.pData), silent= T)

try(logreg.pData.TP53.PFS <- logisticRegression (x = matched.test.pData$relapse , y= matched.test.pData$TP53.cat, data=matched.test.pData), silent= T)


try(logreg.pData.TERT <- logisticRegression (x= matched.test.pData$OS.cat, y= matched.test.pData$TERT.cat, data = matched.test.pData), silent= T)

try(logreg.pData.TERT.PFS <- logisticRegression (x= matched.test.pData$relapse, y= matched.test.pData$TERT.cat, data = matched.test.pData), silent= T)


### generate list of results

logistic.reg.results <- as.list(mget(ls(pattern="logreg.pData"))) 

logistic.reg.results




#################################################################################################

### Determine which results are significant in G3G4
### make dataframe

summary(matched.test.pData$meth)
G3G4.matched.pData <- matched.test.pData[(matched.test.pData$meth=="G3"| matched.test.pData$meth =="G4"), ]


try(logreg.G3G4.mstatus <- logisticRegression(x= matched.test.pData$OS.cat, y= matched.test.pData$mstatus, data=matched.test.pData), silent= T)

try(logreg.G3G4.mstatus.PFS <- logisticRegression(x= matched.test.pData$relapse, y= matched.test.pData$mstatus, data=matched.test.pData), silent= T)


try(logreg.G3G4.resection <- logisticRegression(x =matched.test.pData$OS.cat, y = matched.test.pData$resection, data=matched.test.pData), silent= T)

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


try(logreg.pData.TP53 <- logisticRegression (x = matched.test.pData$OS.cat , y= matched.test.pData$TP53.cat, data=matched.test.pData), silent= T)

try(logreg.pData.TP53.PFS <- logisticRegression (x = matched.test.pData$relapse , y= matched.test.pData$TP53.cat, data=matched.test.pData), silent= T)


try(logreg.pData.TERT <- logisticRegression (x= matched.test.pData$OS.cat, y= matched.test.pData$TERT.cat, data = matched.test.pData), silent= T)

try(logreg.pData.TERT.PFS <- logisticRegression (x= matched.test.pData$relapse, y= matched.test.pData$TERT.cat, data = matched.test.pData), silent= T)


