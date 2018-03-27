### GRAPHICAL DEPICTION OF RESULTS
### created 6 March 2018
### author: Dr Marion Mateos

### Description: graphical depiction of p values against adjusted p values, with abline cutoffs for RNA seq expression data (transcript files) for childhood medulloblastoma cohort

### library(density)

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_master.R")

### loading in clinical data object = test.pData
load("/home/nmm199/R/MB_RNAseq/Clinical/test.pData")

### need results.master file to generate input file

# results.master <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.filt.genefilter.20180220.rds")
# results.master <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.filt.genefilter.random.20180314.rds")
### results.filt.genefilter.20180220.rds: this is a gene filtered file for transcripts and relationship to survival and clinical characteristics within childhood medulloblastoma cohort)

### run clinical_data_extract_DW.R on results.master (such as results.filt.genefilter.20180220.rds)

### may wish to alter graphics later


#########################################################################################################################################
### set output files

pdf ("/home/nmm199/R/MB_RNAseq/Clinical/clin.results/stats.report.pdf") ### adjust to which data are being used
par(mfrow = c(2,2)) ### will output all graphs in 2 x 2 format

### set file for log output
log.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/pDatalog.report.data.txt"

#########################################################################################################################################

### define the variable based on pulling out the subsetted value

p.km.EFS.all <- km.EFS.all.results[, "EFS.p.value"]
adjusted.p.km.EFS.all <- km.EFS.all.results[,"EFS.adjusted.pval"]

### use plotEcdf function



### histogram, use function if possible

hist(p.km.EFS.all)
hist(adjusted.p.km.EFS.all)




#########################################################################################################################################

### Chi squared test
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

chi.report.table <- as.data.frame(cbind(chi.report.mstatus, chi.report.MYC, chi.report.MYC, chi.report.MYCMYCN, chi.report.relapse, chi.report.resection)) ### solution: as.data.frame to return a number not just "numeric" for p.val.range

print(chi.report.table)

### aim to run a function ( run lapply within the function and return (list)), and return object is then written out as a dataframe

# plot(chi.report.table[3,])


#########################################################################################################################################

### Univariate survival 
### example: cox.PFS.cat.G3G4.df

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
}



# hist(p.km.EFS.all)
# hist(adjusted.p.km.EFS.all)


# x <-cox.PFS.cat.G3G4.df[,3]

plotEcdf.Zscore <- function(x, y = NULL, test.name, xlab = "z-score", cutoff=c(-2,2)){
  cdf.x <- ecdf(x)
  plot(ecdf(x), xlab = xlab, main = paste("cumulative density plot of", test.name), col = "red")
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


plotEcdf.Zscore(x = cox.PFS.cat.G3G4.df[,3], test.name = "cox PFS Z scores for G3G4", xlab = "Z-score") ### Z score specific, cox uses categorical expression data


### use plotEcdf function for p value/ adjusted p value and then create the relevant abline, colours. 

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

plotEcdf(x = cox.PFS.cat.G3G4.df[,1], test.name = "cox PFS p values for G3G4", xlab = "p value")

### try to do function so that both p value and p adjust are on the same graph
plotEcdf(x = p.km.EFS.all, test.name = "kaplan meier EFS", xlab = "p value") ### entire cohort
plotEcdf(x = adjusted.p.km.EFS.all, test.name = "kaplan meier EFS", xlab = "adjusted p value", ylab = "Fn(adjusted p value)")
plotEcdf.double(x = p.km.EFS.all, z = adjusted.p.km.EFS.all,  test.name = "kaplan meier EFS", xlab = "p value", ylab = "Fn(p value)")


### interrogate function
x <- p.km.EFS.all
z <- adjusted.p.km.EFS.all
test.name <- "kaplan meier EFS"
xlab <- "p value"
ylab <- "Fn(p value)"

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
}  


# text(3.5, 150, paste("Mean =", round(MyMean, 1), "\n Median =", 
# round(MyMedian, 1), "\n Std.Dev =", round(MySd, 1)))


#########################################################################################################################################

plotHist(cox.PFS.cat.G3G4.df[,1], "Cox PFS categorical G3/G4", breaks = 100, xlab = "p-value", cutoff = 0.05)
plotHist(cox.PFS.cat.G3G4.df[,2], "Cox PFS categorical G3/G4", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05) ### why does the number of adjusted p val with p<0.05 = 0
plotHist(cox.PFS.cat.G3G4.df[,3], "Cox PFS categorical G3/G4", breaks = 100, xlab = "Z-score", cutoff = c(-2, 2)) ### Z score


### when wish to view PDF then need dev.off, sink()

dev.off()
# sink()
#########################################################################################################################################

### Hardcoding for graphs (examples)

### Ecdf plots that describe the empirical cumulative distribution frequency

# plot(ecdf(p.km.EFS.all)) 
# plot(ecdf(adjusted.p.km.EFS.all))  ### km.EFS.all.results[,"EFS.adjusted.pval"]

### Histograms

# hist(cox.PFS.cat.G3G4.df[,1])
# hist(cox.PFS.cat.G3G4.df[,2])





# plot(density(x, na.rm = "T"))

### remove the density plot as is not as useful as ecdf

# plot(density(adjusted.p.km.EFS.all))
# lines(density(adjusted.p.km.EFS.all), col = "red")
# lines(density(p.km.EFS.all), col = "dodgerblue") ### this will overlay the unadjusted p value against the adjusted p value

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
                            
                            
                            
                            
                            
                            
                            
                         
                            