### GRAPHICAL DEPICTION OF RESULTS
### created 6 March 2018
### author: Dr Marion Mateos

### Description: graphical depiction of p values against adjusted p values, with abline cutoffs for RNA seq expression data (transcript files) for childhood medulloblastoma cohort
### Need input file: 
### run clinical_data_extract_DW.R on a defined file (such as results.filt.genefilter.20180220.rds, which is gene filtered file for transcripts and relationship to survival and clinical characteristics within childhood medulloblastoma cohort)
### may wish to alter graphics later

### library(density)

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_master.R")

### loading in clinical data object = test.pData
load("/home/nmm199/R/MB_RNAseq/Clinical/test.pData")

#########################################################################################################################################
### set output files

pdf ("/home/nmm199/R/MB_RNAseq/Clinical/clin.results/stats.report.file") ### adjust to which data are being used
par(mfrow = c(2,2))

### set file for log output
log.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/pDatalog.report.data.txt"

#########################################################################################################################################
# histo.p.adj.km.EFS.all <- hist(adjusted.p.km.EFS.all)


p.km.EFS.all <- km.EFS.all.results[, "EFS.p.value"]
adjusted.p.km.EFS.all <- km.EFS.all.results[,"EFS.adjusted.pval"]

plot(ecdf(p.km.EFS.all))
plot(ecdf(adjusted.p.km.EFS.all))  ### km.EFS.all.results[,"EFS.adjusted.pval"]

hist(p.km.EFS.all)
hist(adjusted.p.km.EFS.all)
### lines(density(km.EFS.p.extract.assembled.all), col = "red")

plot(density(adjusted.p.km.EFS.all))
lines(density(adjusted.p.km.EFS.all), col = "red")
lines(density(p.km.EFS.all), col = "blue") ### this will overlay the unadjusted p value against the adjusted p value

#########################################################################################################################################
### Chi squared test
### e.g output is chi.mstatus.result, names(chi.mstatus.result = "chi.p.value", "adjusted.pval")

chi.report.function <-  function (chi.test, results.master){
  chi.cohort <- nrow(chi.test)
  range.chi.cohort <- range(chi.test)
  grep.infinite <- grep ("inf", chi.test[, "chi.p.value"])
  p.sig.chi.test <- which(chi.test[, "chi.p.value"]<0.05)
  adj.p.sig.chi.test <- which (chi.test[, "adjusted.pval"]<0.05)

chi.list <- list ("cohort.number" = chi.cohort,
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

chi.report.table <- as.data.frame(cbind(chi.report.mstatus, chi.report.MYC, chi.report.MYC, chi.report.MYCMYCN, chi.report.relapse, chi.report.resection)) ### solution: as.data.frame to return a number not just "numeric" for p.val.range

print(chi.report.table)

# plot(chi.report.table[3,])


#########################################################################################################################################

### Univariate survival 
### example: cox.PFS.cat.G3G4.df

head(cox.PFS.cat.G3G4.df)

hist(cox.PFS.cat.G3G4.df[,1])
hist(cox.PFS.cat.G3G4.df[,2])

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
}

# x <-cox.PFS.cat.G3G4.df[,3]

plotEcdf <- function(x, y = NULL, test.name, xlab = "z-score", cutoff=c(-2,2)){
  cdf.x <- ecdf(x)
  plot(ecdf(x), xlab = xlab, main = paste("cumulative density plot of", test.name), col = "red")
  abline(h = 0.5, v = 0)
  abline(v = cutoff, lty = 2)
  min(x, na.rm = T) -> min.x
  temp.no.dn.x <- length(which(x<cutoff[1]))
  max(x, na.rm = T) -> max.x
  temp.no.up.x <- length(which(x>cutoff[2]))
  text(min.x-(0.1*min.x), 0.9, paste("Number Genes z <", cutoff[1],temp.no.dn.x), pos = 4)
  text(max.x-(0.1*max.x), 0.2, paste("Number Genes z >", cutoff[2],temp.no.up.x), pos = 2)
  if(!is.null(y)){
    cdf.y <- ecdf(y)
  }  
  
}



plotHist(cox.PFS.cat.G3G4.df[,1], "Cox PFS categorical G3/G4", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(cox.PFS.cat.G3G4.df[,2], "Cox PFS categorical G3/G4", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05)
plotHist(cox.PFS.cat.G3G4.df[,3], "Cox PFS categorical G3/G4", breaks = 100, xlab = "Z-score", cutoff = c(-2, 2))



plot(ecdf(cox.PFS.cat.G3G4.df[,3]))
# plot(density(x, na.rm = "T"))


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

dev.off()