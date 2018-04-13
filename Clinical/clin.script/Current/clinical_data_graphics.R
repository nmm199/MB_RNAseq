### GRAPHICAL DEPICTION OF RESULTS
### created 6 March 2018
### author: Dr Marion Mateos

### Description: graphical depiction of p values against adjusted p values, with abline cutoffs for RNA seq expression data (transcript files) for childhood medulloblastoma cohort

### library(density)

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_master.R")

### loading in clinical data object = test.pData
load("/home/nmm199/R/MB_RNAseq/Clinical/test.pData")

### need results.master file to generate input file

results.master <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.filt.genefilter.20180220.rds")
# results.master <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.filt.genefilter.random.20180314.rds")
### results.filt.genefilter.20180220.rds: this is a gene filtered file for transcripts and relationship to survival and clinical characteristics within childhood medulloblastoma cohort)

### need to extract clinical data results

### run clinical_data_extract_DW.R on results.master (such as results.filt.genefilter.20180220.rds)

### may wish to alter graphics later


#########################################################################################################################################
### set output files

pdf ("/home/nmm199/R/MB_RNAseq/Clinical/clin.results/stats.report.pdf") ### adjust to which data are being used
# par(mfrow = c(2,2)) ### will output all graphs in 2 x 2 format
par (mfrow = c(1,1))

### set file for log output
log.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/pDatalog.report.data.txt"

### write out tables as required. 
# log.report.transpose



#########################################################################################################################################
### Functions in this file:
### plotEcdf function
### plotEcdf.double (for both p and adj p)
### plotHist



#########################################################################################################################################
#########################################################################################################################################
### FUNCTIONS
### MOVE TO separate functions file later
#########################################################################################################################################
### FUNCTION 1:
### plotEcdf function
### use plotEcdf function for p value/ adjusted p value and then create the relevant abline, colours. 
### input: variables (x, test.name ["name of statistics test" ], cutoff  is p value cut off)
### example inputs for testing function are:

# x <- cox.PFS.cat.G3G4.df[,1]
# test.name <- "cox PFS p values for G3G4"
# xlab <- "p value"
# cutoff=0.05 ### used to be c(-2,2)
# ylab = "Fn(p value)"


plotEcdf <- function(x, y = NULL, test.name, xlab = "p value", ylab = "Fn(p value)", cutoff=0.05){
  cdf.x <- ecdf(x)
  plot(ecdf(x), xlab = xlab, ylab = ylab, main = paste("cumulative density plot of", test.name), col = "dodgerblue")
  abline(a = 0, b = 1, lty = 2) ### where a = intercept, b = slope, lty = line type
  abline (v = 0.05, col = "red", lty = 2 )
  min(x, na.rm = T) -> min.x
  temp.no.dn.x <- length(which(x<cutoff[1])) ### define object cutoff
  max(x, na.rm = T) -> max.x
  # temp.no.up.x <- length(which(x>cutoff[2]))
  text(min.x-(0.1*min.x), 0.9, paste("Number transcripts p <", cutoff[1],"=", temp.no.dn.x), pos = 4) 
  # text(max.x-(0.1*max.x), 0.2, paste("Number transcripts z >", cutoff[2],"=", temp.no.up.x), pos = 2)
  if(!is.null(y)){
    cdf.y <- ecdf(y)
  }  
  
}


#########################################################################################################################################
### FUNCTION 2: plotEcdf.double function
### allows two statistical results to be plotted on the same graph for the same dataset
### example of inputs:

# x <- p.km.EFS.all  ### a statistical result for p value
# z <- adjusted.p.km.EFS.all  ### a statistical result for adjusted p value
# test.name <- "kaplan meier EFS"  ### name of statistical test
# xlab <- "p value" ### x axis label
# ylab <- "Fn(p value)" ### y axis label

plotEcdf.double <- function(x, y = NULL, z, test.name, xlab = "p value", ylab = "Fn(p value)", cutoff=0.05){
  cdf.x <- ecdf(x)
  plot(ecdf(x), xlab = xlab, ylab = ylab, main = paste("cumulative density plot of", test.name), col = "dodgerblue")
  lines((ecdf(z)), col = "red")
  abline(a = 0, b = 1, lty = 2) ### where a = intercept, b = slope, lty = line type
  abline (v = 0.05, col = "red", lty = 2 )
  min(x, na.rm = T) -> min.x
  temp.no.dn.x <- length(which(x<cutoff[1])) ### define object cutoff
  max(x, na.rm = T) -> max.x
  temp.no.dn.z <- length(which(z<cutoff[1]))
  text(min.x-(0.1*min.x), 0.9, paste("Number transcripts p <", cutoff[1],"=", temp.no.dn.x, "\nNumber transcripts adj p <", cutoff[1], "=", temp.no.dn.z), pos = 4) ### invalid position
  legend("bottomright", legend = c("p value", "adjusted p value"), fill = c("dodgerblue", "red"), border = FALSE) 
  if(!is.null(y)){
    cdf.y <- ecdf(y)
  } 
  return(as.data.frame(cbind("number(p<0.05)" = temp.no.dn.x, "number (adj p<0.05)" = temp.no.dn.z), row.names = test.name)) ### added 13/4/18
}  

# text(3.5, 150, paste("Mean =", round(MyMean, 1), "\n Median =", 
# round(MyMedian, 1), "\n Std.Dev =", round(MySd, 1)))

#########################################################################################################################################
### FUNCTION 3: plotEcdf.Zscore function
### developed separately because Z score labelling is different from p value labelling
### same concept as base "plotEcdf" function


plotEcdf.Zscore <- function(x, y = NULL, test.name, xlab = "z-score", cutoff=c(-2,2)){
  cdf.x <- ecdf(x)
  plot(ecdf(x), xlab = xlab, main = paste("cumulative density of", test.name), col = "red") ### cumulative density plot of XYZ
  abline(h = 0.5, v = 0)
  abline(v = cutoff, lty = 2) 
  min(x, na.rm = T) -> min.x
  temp.no.dn.x <- length(which(x<cutoff[1]))
  max(x, na.rm = T) -> max.x
  temp.no.up.x <- length(which(x>cutoff[2]))
  text(min.x-(0.1*min.x), 0.9, paste("Number transcripts z <", cutoff[1],"=", temp.no.dn.x), pos = 4)
  text(max.x-(0.1*max.x), 0.2, paste("Number transcripts z >", cutoff[2],"=", temp.no.up.x), pos = 2)
  if(!is.null(y)){
    cdf.y <- ecdf(y)
  }  
  
}

#########################################################################################################################################

### FUNCTION 4: plotHist function
### use plotHist function to plot histograms of number of transcripts below threshold (with abline, coloured for <0.05) 
### example for hardcoding: ### hist(p.km.EFS.all) ### hist(adjusted.p.km.EFS.all)
### automatically have included breaks = 100
### input 
### x <- "survival statistic result (p value or adjusted p value)"
### test.name <- "name of statistical test performed"


### input variable examples for plotHist function
# x <- cox.PFS.cat.G3G4.df[,1]
# test.name <- "Cox PFS categorical G3/G4"
# breaks = 100

plotHist <- function(x, test.name, breaks = 100, xlab = "p-value", cutoff = 0.05, text.pos = 0.9){
  hist.res <- hist(x, breaks = breaks, plot = F)
  max(hist.res$counts) -> temp.height
  length(which(x<0.05)) -> temp.no.sig
  
  if(length(cutoff)==1){
    ifelse(hist.res$breaks<cutoff,"red", "grey") -> hist.cols
  }else{
    ifelse(hist.res$breaks<cutoff[1]|hist.res$breaks>cutoff[2],"red", "grey") -> hist.cols
  }
  
  hist(x, breaks = breaks,  xlab = xlab, main = paste("Histogram of", test.name), col = hist.cols)
  
  if(length(cutoff)==1){
    abline(v= cutoff, lty = 2 , col = "red")
    text(text.pos, temp.height-(temp.height*0.1), paste("Number Genes p <", cutoff, "=", temp.no.sig), pos = 2)
  }else{
    abline(v= cutoff, lty = 2 , col = "red")
    text(text.pos, temp.height-(temp.height*0.1), paste("Number Genes p <", cutoff[1],"or p >", cutoff[2],"=", temp.no.sig), pos = 2) 
  }
  
  # return(as.data.frame(c("total number significant", temp.no.sig))) ### added new line 13/4/18
  return(as.data.frame(temp.no.sig, row.names = test.name)) ### added new line 13/4/18
}

#########################################################################################################################################
#########################################################################################################################################
### here is the start of the script to execute to generate graphics
#########################################################################################################################################
### TO CREATE WRAPPER FUNCTION:
### define the variable based on pulling out the subsetted value
### then do lapply 
### return a list or a dataframe
#########################################################################################################################################

### INDIVIDUAL FUNCTIONS BELOW WITH APPLICATION TO EACH SURVIVAL TEST, CLINICAL VARIABLES, TRANSCRIPT DATA (CATEGORICAL OR CONTINUOUS) 
### UNIVARIATE AND MULTIVARIATE ANALYSES
#########################################################################################################################################
### UNIVARIATE
### use plotEcdf function

### histogram, use plotHist function

### categorical expression data for OS, PFS, EFS for both overall (all data) and G3G4

plotHist(cox.OS.cat.all.df[,1],  "Cox OS categorical overall", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(cox.OS.cat.all.df[,2],  "Cox OS categorical overall", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)
plotHist(cox.OS.cat.G3G4.df[,1], "Cox OS categorical G3/G4", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(cox.OS.cat.G3G4.df[,2], "Cox OS categorical G3/G4", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(cox.PFS.cat.all.df[,1],  "Cox PFS categorical overall", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(cox.PFS.cat.all.df[,2],  "Cox PFS categorical overall", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)
plotHist(cox.PFS.cat.G3G4.df[,1], "Cox PFS categorical G3/G4", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(cox.PFS.cat.G3G4.df[,2], "Cox PFS categorical G3/G4", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(cox.EFS.cat.all.df[,1],  "Cox EFS categorical overall", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(cox.EFS.cat.all.df[,2],  "Cox EFS categorical overall", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)
plotHist(cox.EFS.cat.G3G4.df[,1], "Cox EFS categorical G3/G4", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(cox.EFS.cat.G3G4.df[,2], "Cox EFS categorical G3/G4", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

### continuous expression data for OS, PFS for both overall (all data) and G3G4. Cox EFS data does not exist for continuous transcript expression

plotHist(cox.OS.cont.all.df[,1],  "Cox OS continuous overall", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(cox.OS.cont.all.df[,2],  "Cox OS continuous overall", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)
plotHist(cox.OS.cont.G3G4.df[,1], "Cox OS continuous G3/G4", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(cox.OS.cont.G3G4.df[,2], "Cox OS continuous G3/G4", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(cox.PFS.cont.all.df[,1],  "Cox PFS continuous overall", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(cox.PFS.cont.all.df[,2],  "Cox PFS continuous overall", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)
plotHist(cox.PFS.cont.G3G4.df[,1], "Cox PFS continuous G3/G4", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(cox.PFS.cont.G3G4.df[,2], "Cox PFS continuous G3/G4", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

#########################################################################################################################################
### MULTIVARIATE
### for multivariate cox regression models 

### combined model

plotHist(multivar.cox.OS.combined.cat.df[,1],  "Cox OS combined multivar categorical", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(multivar.cox.OS.combined.cat.df[,2],  "Cox OS combined multivar categorical", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(multivar.cox.OS.combined.cont.df[,1],  "Cox OS combined multivar continuous", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(multivar.cox.OS.combined.cont.df[,2],  "Cox OS combined multivar continuous", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(multivar.cox.PFS.combined.cat.df[,1],  "Cox PFS combined multivar categorical", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(multivar.cox.PFS.combined.cat.df[,2],  "Cox PFS combined multivar categorical", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(multivar.cox.PFS.combined.cont.df[,1],  "Cox PFS combined multivar continuous", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(multivar.cox.PFS.combined.cont.df[,2],  "Cox PFS combined multivar continuous", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

### PNET5 model

plotHist(multivar.cox.OS.PNET5.cat.df[,1],  "Cox OS PNET5 multivar categorical", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(multivar.cox.OS.PNET5.cat.df[,2],  "Cox OS PNET5 multivar categorical", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(multivar.cox.OS.PNET5.cont.df[,1],  "Cox OS PNET5 multivar continuous", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(multivar.cox.OS.PNET5.cont.df[,2],  "Cox OS PNET5 multivar continuous", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(multivar.cox.PFS.PNET5.cat.df[,1],  "Cox PFS PNET5 multivar categorical", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(multivar.cox.PFS.PNET5.cat.df[,2],  "Cox PFS PNET5 multivar categorical", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(multivar.cox.PFS.PNET5.cont.df[,1],  "Cox PFS PNET5 multivar continuous", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(multivar.cox.PFS.PNET5.cont.df[,2],  "Cox PFS PNET5 multivar continuous", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)


### lancet G3G4 model

plotHist(multivar.cox.OS.lancetG3G4.cat.df[,1],  "Cox OS lancet G3G4 multivar categorical", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(multivar.cox.OS.lancetG3G4.cat.df[,2],  "Cox OS lancet G3G4 multivar categorical", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(multivar.cox.OS.lancetG3G4.cont.df[,1],  "Cox OS lancet G3G4 multivar continuous", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(multivar.cox.OS.lancetG3G4.cont.df[,2],  "Cox OS lancet G3G4 multivar continuous", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(multivar.cox.PFS.lancetG3G4.cat.df[,1],  "Cox OS lancet G3G4 multivar categorical", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(multivar.cox.PFS.lancetG3G4.cat.df[,2],  "Cox OS lancet G3G4 multivar categorical", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(multivar.cox.PFS.lancetG3G4.cont.df[,1],  "Cox OS lancet G3G4 multivar categorical", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(multivar.cox.PFS.lancetG3G4.cont.df[,2],  "Cox OS lancet G3G4 multivar categorical", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

#########################################################################################################################################

### wrapper function would contain similar to
### plotHist (x, test.name = test.name, xlab = "p value", ")
### might need to separate for p value and adjusted p value
### loop through different "test.name" options
### plotHist(x [,2], "test.name", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

#########################################################################################################################################
#########################################################################################################################################
### COX MULTIVARIATE using plotEcdf

### combined model
plotEcdf.double(x=multivar.cox.OS.combined.cat.df[,"cox.pval"], z= multivar.cox.OS.combined.cat.df[,"cox.adj.pval"], test.name = "multivar cox OS combined (cat)")
plotEcdf.double(x=multivar.cox.PFS.combined.cat.df[,"cox.pval"], z= multivar.cox.PFS.combined.cat.df[,"cox.adj.pval"], test.name = "multivar cox PFS combined (cat)")

plotEcdf.double(x=multivar.cox.OS.combined.cont.df[,"cox.pval"], z= multivar.cox.OS.combined.cont.df[,"cox.adj.pval"], test.name = "multivar cox OS combined (cont)")
plotEcdf.double(x=multivar.cox.PFS.combined.cont.df[,"cox.pval"], z= multivar.cox.PFS.combined.cont.df[,"cox.adj.pval"], test.name = "multivar cox PFS combined (cont)")

### PNET5 
plotEcdf.double(x=multivar.cox.OS.PNET5.cat.df[,"cox.pval"], z= multivar.cox.OS.PNET5.cat.df[,"cox.adj.pval"], test.name = "multivar cox OS PNET5 (cat)")
plotEcdf.double(x=multivar.cox.PFS.PNET5.cat.df[,"cox.pval"], z= multivar.cox.PFS.PNET5.cat.df[,"cox.adj.pval"], test.name = "multivar cox PFS PNET5 (cat)")

plotEcdf.double(x=multivar.cox.OS.PNET5.cont.df[,"cox.pval"], z= multivar.cox.OS.PNET5.cont.df[,"cox.adj.pval"], test.name = "multivar cox OS PNET5 (cont)")
plotEcdf.double(x=multivar.cox.PFS.PNET5.cont.df[,"cox.pval"], z= multivar.cox.PFS.PNET5.cont.df[,"cox.adj.pval"], test.name = "multivar cox PFS PNET5 (cont)")

### G3G4

plotEcdf.double(x=multivar.cox.OS.lancetG3G4.cat.df[,"cox.pval"], z= multivar.cox.OS.lancetG3G4.cat.df[,"cox.adj.pval"], test.name = "multivar cox OS G3G4(cat)")
plotEcdf.double(x=multivar.cox.PFS.lancetG3G4.cat.df[,"cox.pval"], z= multivar.cox.PFS.lancetG3G4.cat.df[,"cox.adj.pval"], test.name = "multivar cox PFS G3G4 (cat)")

plotEcdf.double(x=multivar.cox.OS.lancetG3G4.cont.df[,"cox.pval"], z= multivar.cox.OS.lancetG3G4.cont.df[,"cox.adj.pval"], test.name = "multivar cox OS G3G4 (cont)")
plotEcdf.double(x=multivar.cox.PFS.lancetG3G4.cont.df[,"cox.pval"], z= multivar.cox.PFS.lancetG3G4.cont.df[,"cox.adj.pval"], test.name = "multivar cox PFS G3G4 (cont)")


#########################################################################################################################################
### MULTIVARIATE COX Z SCORES

### combined model OS/PFS Z score

plotEcdf.Zscore(x=multivar.cox.OS.combined.cat.df[, "cox.Zscore"], test.name = "multivar cox OS comb Z score (cat)")
plotEcdf.Zscore(x=multivar.cox.OS.combined.cont.df[, "cox.Zscore"], test.name = "multivar cox OS comb Z score (cont)")

plotEcdf.Zscore(x=multivar.cox.PFS.combined.cat.df[, "cox.Zscore"], test.name = "multivar cox PFS comb Z score (cat)")
plotEcdf.Zscore(x=multivar.cox.PFS.combined.cont.df[, "cox.Zscore"], test.name = "multivar cox PFS comb Z score (cont)")


### OS Z scores for other models (PNET4, lancetG3G4, SHHold)
plotEcdf.Zscore(x=multivar.cox.OS.PNET5.cat.df[, "cox.Zscore"], test.name = "multivar cox OS PNET5 Z score  (cat)")
plotEcdf.Zscore(x=multivar.cox.OS.PNET5.cont.df[, "cox.Zscore"], test.name = "multivar cox OS PNET5 Z score (cont)")

plotEcdf.Zscore(x=multivar.cox.OS.lancetG3G4.cat.df[, "cox.Zscore"], test.name = "multivar cox OS G3G4 Z score  (cat)")
plotEcdf.Zscore(x=multivar.cox.OS.lancetG3G4.cont.df[, "cox.Zscore"], test.name = "multivar cox OS G3G4 Z score (cont)")

plotEcdf.Zscore(x=multivar.cox.OS.SHHold.cat.df[, "cox.Zscore"], test.name = "multivar cox OS SHH older Z score  (cat)")
plotEcdf.Zscore(x=multivar.cox.OS.SHHold.cont.df[, "cox.Zscore"], test.name = "multivar cox OS SHH older Z score (cont)")


### PFS Z scores for other models (PNET4, lancetG3G4, SHHold)

plotEcdf.Zscore(x=multivar.cox.PFS.PNET5.cat.df[, "cox.Zscore"], test.name = "multivar cox PFS PNET5 Z score  (cat)")
plotEcdf.Zscore(x=multivar.cox.PFS.PNET5.cont.df[, "cox.Zscore"], test.name = "multivar cox PFS PNET5 Z score (cont)")

plotEcdf.Zscore(x=multivar.cox.PFS.lancetG3G4.cat.df[, "cox.Zscore"], test.name = "multivar cox PFS G3G4 Z score  (cat)")
plotEcdf.Zscore(x=multivar.cox.PFS.lancetG3G4.cont.df[, "cox.Zscore"], test.name = "multivar cox PFS G3G4 Z score (cont)")

plotEcdf.Zscore(x=multivar.cox.PFS.SHHold.cat.df[, "cox.Zscore"], test.name = "multivar cox PFS SHH older Z score  (cat)")
plotEcdf.Zscore(x=multivar.cox.PFS.SHHold.cont.df[, "cox.Zscore"], test.name = "multivar cox PFS SHH older Z score (cont)")

#########################################################################################################################################
#########################################################################################################################################

### UNIVARIATE CHI

### chi.report.function for Chi squared test
### e.g output is chi.mstatus.result, names(chi.mstatus.result = "chi.p.value", "adjusted.pval")

chi.report.function <-  function (chi.test, results.master){
  chi.cohort <- nrow(chi.test)
  range.chi.cohort <- range(chi.test)
  grep.infinite <- grep ("inf", chi.test[, "chi.p.value"])
  p.sig.chi.test <- which(chi.test[, "chi.p.value"]<0.05)
  adj.p.sig.chi.test <- which (chi.test[, "adjusted.pval"]<0.05)

chi.list <- list ("transcript.number" = chi.cohort,
       "p.val.range" = range.chi.cohort, 
       "total.sig.pval.number" = length(p.sig.chi.test), 
       "total.sig.adj.pval.number" = length(adj.p.sig.chi.test)
       )
return (chi.list)
}      

chi.report.mstatus <- chi.report.function(chi.mstatus.result, results.master)
chi.report.relapse <- chi.report.function(chi.relapse.result, results.master)
chi.report.resection <- chi.report.function(chi.resection.result, results.master)
chi.report.MYCN <- chi.report.function (chi.MYCN.result, results.master)
chi.report.MYC <- chi.report.function (chi.MYC.result, result.master)
chi.report.MYCMYCN <- chi.report.function(chi.MYCMYCN.result, results.master)

chi.report.table <- as.data.frame(cbind(chi.report.mstatus, chi.report.MYC, chi.report.MYCN, chi.report.MYCMYCN, chi.report.relapse, chi.report.resection)) ### solution: as.data.frame to return a number not just "numeric" for p.val.range

print(chi.report.table)

# write.csv(chi.report.table, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/Complete_transcripts/chireport.gp.filt.mb.vsd.20180413.csv") ### returns an error "unimplemented type 'list' in EncodeElement'

### aim to run a function ( run lapply within the function and return (list)), and return object is then written out as a dataframe


#########################################################################################################################################
### LOGISTIC REGRESSION REPORTING FUNCTION


Log.report.function <-  function (log.reg.test, results.master){
  log.cohort <- nrow(log.reg.test)
  range.logreg.pval <- range(log.reg.test[,"logreg.pval"])
  range.logreg.adjpval <- range(log.reg.test[,"logreg.adj.pval"])
  range.logreg.OR <- range(log.reg.test [, "logreg.OR"])
  
  p.sig.log.reg.test <- which(log.reg.test[, "logreg.pval"]<0.05)
  adj.p.sig.log.reg.test <- which (log.reg.test[, "logreg.adj.pval"]<0.05)
  
 log.reg.list <- list ("transcript.number" = log.cohort,
                    "p.val.range" = range.logreg.pval,
                    "adj.p.val.range" = range.logreg.adjpval,
                    "total.sig.pval.number" = length(p.sig.log.reg.test), 
                    "total.sig.adj.pval.number" = length(adj.p.sig.log.reg.test)
  )
  return (log.reg.list)
}



log.report.list <- lapply(logistic.reg.results, Log.report.function)

log.report.table <- as.data.frame (log.report.list)
log.report.transpose <- t(log.report.table) ### transpose #
print(log.report.transpose)

write.csv(log.report.transpose, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/Complete_transcripts/logreport.gp.filt.mb.vsd.20180413.csv")
#########################################################################################################################################

### COX REGRESSION UNIVARIATE
### example: cox.PFS.cat.G3G4.df
# colnames(cox.PFS.cat.G3G4.df)

plotEcdf.double(x = cox.OS.cat.all.df[, "cox.pval"], z = cox.OS.cat.all.df[, "cox.adj.pval"], test.name = "cox OS in all (categorical)", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = cox.PFS.cat.all.df[, "cox.pval"], z = cox.PFS.cat.all.df[, "cox.adj.pval"], test.name = "cox PFS in all (categorical)", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = cox.EFS.cat.all.df[, "cox.pval"], z = cox.EFS.cat.all.df[, "cox.adj.pval"], test.name = "cox EFS in all (categorical)", xlab = "p value", ylab = "Fn(p value)")

plotEcdf.double(x = cox.OS.cat.G3G4.df[, "cox.pval"], z = cox.OS.cat.G3G4.df[, "cox.adj.pval"], test.name = "cox OS in G3G4 (categorical)", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = cox.PFS.cat.G3G4.df[, "cox.pval"], z = cox.PFS.cat.G3G4.df[, "cox.adj.pval"], test.name = "cox PFS in G3G4 (categorical)", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = cox.EFS.cat.G3G4.df[, "cox.pval"], z = cox.EFS.cat.G3G4.df[, "cox.adj.pval"], test.name = "cox EFS in G3G4 (categorical)", xlab = "p value", ylab = "Fn(p value)")

plotEcdf.double(x = cox.OS.cont.all.df[, "cox.pval"], z = cox.OS.cont.all.df[, "cox.adj.pval"], test.name = "cox OS in all (continuous)", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = cox.PFS.cont.all.df[, "cox.pval"], z = cox.PFS.cont.all.df[, "cox.adj.pval"], test.name = "cox PFS in all (continuous)", xlab = "p value", ylab = "Fn(p value)")

plotEcdf.double(x = cox.OS.cont.G3G4.df[, "cox.pval"], z = cox.OS.cont.G3G4.df[, "cox.adj.pval"], test.name = "cox OS in G3G4 (continuous)", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = cox.PFS.cont.G3G4.df[, "cox.pval"], z = cox.PFS.cont.G3G4.df[, "cox.adj.pval"], test.name = "cox PFS in G3G4 (continuous)", xlab = "p value", ylab = "Fn(p value)")


### SHH results (Some of these are significant for continuous variable for 9693 transcript file 13/4/18)

plotEcdf.double(x = cox.OS.cat.SHH.df[, "cox.pval"], z = cox.OS.cat.SHH.df[, "cox.adj.pval"], test.name = "cox OS in SHH (categorical)", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = cox.PFS.cat.SHH.df[, "cox.pval"], z = cox.PFS.cat.SHH.df[, "cox.adj.pval"], test.name = "cox PFS in SHH (categorical)", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = cox.OS.cont.SHH.df[, "cox.pval"], z = cox.OS.cont.SHH.df[, "cox.adj.pval"], test.name = "cox OS in SHH (continuous)", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = cox.PFS.cont.SHH.df[, "cox.pval"], z = cox.PFS.cont.SHH.df[, "cox.adj.pval"], test.name = "cox PFS in SHH (continuous)", xlab = "p value", ylab = "Fn(p value)")


plotEcdf.double(x = cox.OS.cat.SHH.old.df[, "cox.pval"], z = cox.OS.cat.SHH.old.df[, "cox.adj.pval"], test.name = "cox OS in SHH older (categorical)", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = cox.PFS.cat.SHH.old.df[, "cox.pval"], z = cox.PFS.cat.SHH.old.df[, "cox.adj.pval"], test.name = "cox PFS in SHH older (categorical)", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = cox.OS.cont.SHH.old.df[, "cox.pval"], z = cox.OS.cont.SHH.old.df[, "cox.adj.pval"], test.name = "cox OS in SHH older (continuous)", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = cox.PFS.cont.SHH.old.df[, "cox.pval"], z = cox.PFS.cont.SHH.old.df[, "cox.adj.pval"], test.name = "cox PFS in SHH older (continuous)", xlab = "p value", ylab = "Fn(p value)")


### Z scores

plotEcdf.Zscore(x = cox.OS.cat.all.df[,3], test.name = "cox OS Z scores all (cat)", xlab = "Z-score")
plotEcdf.Zscore(x = cox.PFS.cat.all.df[,3], test.name = "cox PFS Z scores all (cat)", xlab = "Z-score")
plotEcdf.Zscore(x = cox.EFS.cat.all.df[,3], test.name = "cox EFS Z scores all (cat)", xlab = "Z-score")

plotEcdf.Zscore(x = cox.OS.cat.G3G4.df[,3], test.name = "cox OS Z scores G3G4 (cat)", xlab = "Z-score")
plotEcdf.Zscore(x = cox.PFS.cat.G3G4.df[,3], test.name = "cox PFS Z scores G3G4 (cat)", xlab = "Z-score") ### Z score specific, cox uses categorical expression data
plotEcdf.Zscore(x = cox.EFS.cat.G3G4.df[,3], test.name = "cox EFS Z scores G3G4 (cat)", xlab = "Z-score")

plotEcdf.Zscore(x = cox.OS.cont.all.df[,3], test.name = "cox OS Z scores all (cont)", xlab = "Z-score")
plotEcdf.Zscore(x = cox.PFS.cont.all.df[,3], test.name = "cox PFS Z scores all (cont)", xlab = "Z-score")

plotEcdf.Zscore(x = cox.OS.cont.G3G4.df[,3], test.name = "cox OS Z scores G3G4 (cont)", xlab = "Z-score")
plotEcdf.Zscore(x = cox.PFS.cont.G3G4.df[,3], test.name = "cox PFS Z scores G3G4 (cont)", xlab = "Z-score")


### Z scores SHH specific
plotEcdf.Zscore(x = cox.OS.cat.SHH.df[,3], test.name = "cox OS Z scores SHH (cat)", xlab = "Z-score")
plotEcdf.Zscore(x = cox.PFS.cat.SHH.df[,3], test.name = "cox PFS Z scores SHH (cat)", xlab = "Z-score")

plotEcdf.Zscore(x = cox.OS.cat.SHH.old.df[,3], test.name = "cox OS Z scores SHH older (cat)", xlab = "Z-score")
plotEcdf.Zscore(x = cox.PFS.cat.SHH.old.df[,3], test.name = "cox PFS Z scores SHH older (cat)", xlab = "Z-score") ### Z score specific, cox uses categorical expression data

plotEcdf.Zscore(x = cox.OS.cont.SHH.df[,3], test.name = "cox OS Z scores SHH (cont)", xlab = "Z-score")
plotEcdf.Zscore(x = cox.PFS.cont.SHH.df[,3], test.name = "cox PFS Z scores SHH (cont)", xlab = "Z-score")

plotEcdf.Zscore(x = cox.OS.cont.SHH.old.df[,3], test.name = "cox OS Z scores SHH older (cont)", xlab = "Z-score")
plotEcdf.Zscore(x = cox.PFS.cont.SHH.old.df[,3], test.name = "cox PFS Z scores SHH older (cont)", xlab = "Z-score")


#########################################################################################################################################
### KAPLAN MEIER UNIVARIATE
### use plotEcdf function so that both p value and p adjust are on the same graph
### do for OS/PFS for km all and G3G4 (note that this is for continuous transcripts against a dichotomous survival outcome)

plotEcdf.double(x = km.OS.all.results[,1], z = km.OS.all.results[,2],  test.name = "kaplan meier OS", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = km.PFS.all.results[,1], z = km.PFS.all.results[,2],  test.name = "kaplan meier PFS", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = km.EFS.all.results[,1], z = km.EFS.all.results[,2],  test.name = "kaplan meier EFS", xlab = "p value", ylab = "Fn(p value)")


plotEcdf.double(x = km.OS.G3G4.results[,1], z = km.OS.G3G4.results [,2],  test.name = "kaplan meier OS G3G4", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = km.PFS.G3G4.results[,1], z = km.PFS.G3G4.results[,2],  test.name = "kaplan meier PFS G3G4", xlab = "p value", ylab = "Fn(p value)")
plotEcdf.double(x = km.EFS.G3G4.results[,1], z = km.EFS.G3G4.results[,2],  test.name = "kaplan meier EFS G3G4", xlab = "p value", ylab = "Fn(p value)")


### Z scores do not apply to the km survival (eg. columns are OS.p.value and OS.adjusted.pval) # colnames(km.OS.all.results)

#########################################################################################################################################
#########################################################################################################################################
### RESULTS THAT ARE GENERATED HOWEVER SMALL NUMBERS AND GRAPHICS PROBLEMATIC

### lancet SHH older model ### not all these graphs are useful

plotHist(multivar.cox.OS.SHHold.cat.df[,1],  "Cox OS SHH older multivar categorical", breaks = 100, xlab = "p-value", cutoff = 0.05)
# plotHist(multivar.cox.OS.SHHold.cat.df[,2],  "Cox OS SHH older multivar categorical", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(multivar.cox.OS.SHHold.cont.df[,1],  "Cox OS SHH older multivar continuous", breaks = 100, xlab = "p-value", cutoff = 0.05)
# plotHist(multivar.cox.OS.SHHold.cont.df[,2],  "Cox OS SHH older multivar continuous", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(multivar.cox.PFS.SHHold.cat.df[,1],  "Cox PFS SHH older multivar categorical", breaks = 100, xlab = "p-value", cutoff = 0.05)
# plotHist(multivar.cox.PFS.SHHold.cat.df[,2],  "Cox PFS SHH older multivar categorical", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)

plotHist(multivar.cox.PFS.SHHold.cont.df[,1],  "Cox PFS SHH older multivar continuous", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(multivar.cox.PFS.SHHold.cont.df[,2],  "Cox PFS SHH older multivar continuous", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)



### SHH older ### unusual results due to small numbers of patients? (cox.adj.pval = 1) ?exclude
plotEcdf.double(x=multivar.cox.OS.SHHold.cat.df[,"cox.pval"], z= multivar.cox.OS.SHHold.cat.df[,"cox.adj.pval"], test.name = "multivar cox OS SHH older (cat)")
plotEcdf.double(x=multivar.cox.PFS.SHHold.cat.df[,"cox.pval"], z= multivar.cox.PFS.SHHold.cat.df[,"cox.adj.pval"], test.name = "multivar cox PFS SHH older (cat)")

plotEcdf.double(x=multivar.cox.OS.SHHold.cont.df[,"cox.pval"], z= multivar.cox.OS.SHHold.cont.df[,"cox.adj.pval"], test.name = "multivar cox OS SHH older (cont)")
plotEcdf.double(x=multivar.cox.PFS.SHHold.cont.df[,"cox.pval"], z= multivar.cox.PFS.SHHold.cont.df[,"cox.adj.pval"], test.name = "multivar cox PFS SHH older (cont)")

#########################################################################################################################################
#########################################################################################################################################

### when wish to view PDF then need dev.off, sink()

dev.off()
# sink()
#########################################################################################################################################
#########################################################################################################################################
### GETTING STARTED WITH WRITING FUNCTIONS
### DETERMINE DATASET, RANGES, STRUCTURE, NAMES (names/colnames), length/nrow

#names(logistic.reg.results)
#head(logistic.reg.results)
#colnames(logistic.reg.results[[1]])
#"logreg.pval" 
#"logreg.adj.pval" 

#nrow(extract.logreg.LCA.df) ### 9693
#length(extract.logreg.LCA.df) ### 48465

#range(extract.logreg.LCA.df[,"logreg.pval"])
#range(extract.logreg.LCA.df[, "logreg.adj.pval"])
#range(extract.logreg.LCA.df[,"logreg.OR"])

# rm(log.reg.test)
# log.reg.test <- extract.logreg.LCA.df
# results.master <- results.master

#########################################################################################################################################

### REMOVE THIS HARDCODING WHEN FILE IS COMPLETE
### Hardcoding for graphs (examples)

### Ecdf plots that describe the empirical cumulative distribution frequency

# plot(ecdf(p.km.EFS.all)) 
# plot(ecdf(adjusted.p.km.EFS.all))  ### km.EFS.all.results[,"EFS.adjusted.pval"]

### Histograms

# hist(cox.PFS.cat.G3G4.df[,1])
# hist(cox.PFS.cat.G3G4.df[,2])


### redundant script removed above, for example:
# plotEcdf(x = cox.PFS.cat.G3G4.df[,1], test.name = "cox PFS p values for G3G4", xlab = "p value")
# plotEcdf (x = cox.PFS.cat.G3G4.df[,2], test.name = "cox PFS adj p values for G3G4", xlab = "adjusted p value")
# p.cox.PFS.cat.G3G4 <- cox.PFS.cat.G3G4.df[,1]
# p.adj.cox.PFS.cat.G3G4 <- cox.PFS.cat.G3G4.df[,2]
# plotEcdf.double(x = p.cox.PFS.cat.G3G4, z = p.adj.cox.PFS.cat.G3G4, test.name = "cox PFS in G3G4 (categorical)", xlab = "p value", ylab = "Fn(p value)")
### has been replaced by:
# plotEcdf.double(x = cox.EFS.cat.G3G4.df[, "cox.pval"], z = cox.EFS.cat.G3G4.df[, "cox.adj.pval"], test.name = "cox EFS in G3G4 (categorical)", xlab = "p value", ylab = "Fn(p value)")


### another example
# plotEcdf(x = p.km.EFS.all, test.name = "kaplan meier EFS", xlab = "p value") ### entire cohort fpr just p value
# plotEcdf(x = adjusted.p.km.EFS.all, test.name = "kaplan meier EFS", xlab = "adjusted p value", ylab = "Fn(adjusted p value)") ### entire cohort for adjusted p value
### has been replaced by:
# plotEcdf.double(x = p.km.EFS.all, z = adjusted.p.km.EFS.all,  test.name = "kaplan meier EFS", xlab = "p value", ylab = "Fn(p value)")


### another example of simplified coding
#p.km.EFS.all <- km.EFS.all.results[, "EFS.p.value"]
#adjusted.p.km.EFS.all <- km.EFS.all.results[,"EFS.adjusted.pval"]
#hist(p.km.EFS.all)
#hist(adjusted.p.km.EFS.all)

### replaced by
# plotHist(cox.PFS.cat.G3G4.df[,1], "Cox PFS categorical G3/G4", breaks = 100, xlab = "p-value", cutoff = 0.05)
# plotHist(cox.PFS.cat.G3G4.df[,2], "Cox PFS categorical G3/G4", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05) 
#########################################################################################################################################


# plot(density(x, na.rm = "T"))

### remove the density plot as is not as useful as ecdf

# plot(density(adjusted.p.km.EFS.all))
# lines(density(adjusted.p.km.EFS.all), col = "red")
# lines(density(p.km.EFS.all), col = "dodgerblue") ### this will overlay the unadjusted p value against the adjusted p value


#########################################################################################################################################
### INTERROGATING SPECIFIC CLINICAL DATA AND SPECIFIC TRANSCRIPTS
#########################################################################################################################################
### km survival curve graphing using significant transcripts e.g ENSG00000168772.10_1 for gp.filt.mb.vsd, results.sig.multivar.cox.OS.combined.cat
### reload relevant sections from clinical_data_master.R

library(NMF)

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

### To depict KM curves for defined transcript (goi) 

RNA.data <- "/home/dan/mygit/rna_seq_mb/paper/MB.vsd.txt" 

mb.vsd.novel <- read.delim(file="/home/dan/mygit/rna_seq_mb/paper/vsd.novel.txt") 

mb.vsd <- read.delim(RNA.data)

# mb.vsd.random <- randomize(mb.vsd) ### generate this first then run the clinPathAssess function on this

rownames(mb.vsd.random) <- rownames(mb.vsd)

### then assess within clinPathAssess 

goi <- "ENSG00000168772"
# goi <- "ENSG00000173818.16"  ### sig in G3G4 lancet model above current factors, in PFS and OS
goi.vsd <- as.numeric(mb.vsd[goi,]) ### 9/1/18 hashed when running the clinical_data_master.R; unhashed when interrogating clinPathAssess function ie. ### unhash ** when running individual goi

### the output would be a vector with a continuous variable, names equal NMB numbers
names(mb.vsd) -> names(goi.vsd)  ### added 16/1/18, unhash ** when running individual goi

test.pData = test.pData

## results for run of one
results.master.goi <- clinPathAssess(test.pData,
 goi.vsd,
 pdf.file = pdf.file,
 log.file = log.file)


# names(results.master.goi)<- row.names(goi.vsd)

saveRDS (results.master.goi, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.ENSG00000168772.rds")
readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/results.master.ENSG00000168772.rds")



#########################################################################################################################################

### evaluating the specifics of extracted data 30/1/18, which requires the clinical_data_master to be loaded 

# plot(ecdf(mb.vsd["ENSG00000188314",])) ### this showed a flat density curve, fn(x)=1
# mb.vsd["ENSG00000188314", ]  ### this showed that all expression values were the SAME, which explains why p value was 0 and cox z score infinite



#########################################################################################################################################
### Examples from external sources for plotting 2 or more sets of data on same graph
#plot(x, y1, ylim=range(c(y1,y2)))
# second plot  EDIT: needs to have same ylim
###par(new = TRUE)
#plot(x, y2, ylim=range(c(y1,y2)), axes = FALSE, xlab = "", ylab = "")

#matplot(x, cbind(y1,y2))
#matplot(x, cbind(y1,y2), pch=1)

# dev.off()



### This page aims to explain how to add a legend to R plot made in base R. It is done using the legend() function. The main arguments are:

# topright : where do you want to add the legend ? You can put : “bottomright”, “bottom”, “bottomleft”, “left”, “topleft”, “top”, “topright”, “right”, “center”).
#legend = c(“name1”, “name2”) : names you want to show.
#col = c(“red”, “blue”) : colors of the symbols
#pch = 15 : type of symbols (see graph # to know what symbol number you need
                        #    bty = “n” : If you don’t want a box around the legend. Write “o” if you want one
                         #   pt.cex = 2 : Size of the symbols
                         #   cex = 0.8 : Size of the text
                         #   text.col = “black” : color of the text
                         #   horiz = TRUE : legend in column or in row ?
                         #   inset = c(0.1, 0.1) : % (from 0 to 1) to draw the legend away from x and y axis
                          ###  You can also give the X and Y coordinate of the legend: legend(3, 5, legend = c(“A”, “B”))
                          ### Note that an equivalent page exist concerning legends with ggplot2.
                          ### see example of this in the pdf_graphics_example.R file 
                            
                            
                            
                            
                            
                            
                            
                         
                            