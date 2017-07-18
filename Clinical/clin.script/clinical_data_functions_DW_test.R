### File information:

### This is a list of functions that were written for the "clinical_data_4.R" script analysis
### Many of these functions are generic for survival analysis
### Aims of this document are to provide tools to perform univariate analysis in a group of children treated for medulloblastoma. 
### All kaplan-meier survival analyses provide graphical output with labelled axes and p value
### The survival analysis is designed in a cohort of children who received therapy with curative intent (including cranio-spinal irradiation), age 3-16 years

### Author: Dr Marion Mateos
### Date: July 14 2017


##########################################################################################

### Function Number 1
### Function entitled "chi.sq"
### provide 
## input for 2x2 table: 
## x <- variable 1
## y <- variable 2

## output: 
## percentage affected by variables in 2x2 table
## p value for pearson chi-squared analysis
## chi squared residuals 

chi.sq <- function(x,y){
  table.temp <- table(x, y) 
  table.temp.perc <- prop.table(table.temp)*100
  summary.table(table.temp)
  chi.test.temp <- chisq.test(table.temp) 
  chi.test.temp.stat <- c(stat=chi.test.temp$statistic, p.value=chi.test.temp$p.value) 
  chi.test.temp.res <- chi.test.temp$residuals
  aheatmap(chi.test.temp$residuals, Rowv=NA, Colv = NA)
  list.temp <- list  (table.temp, 
                      table.temp.perc,
                      chi.test.temp,
                      chi.test.temp.res
  )
  
  return(list.temp)
}



##########################################################################################



### Function Number 2
### Function entitled "cor.result"
### Aim of function: evaluation correlation between 2 variables, x and y

### input:
## x <- variable 1
## y <- variable 2

### output:
## list of correlation

cor.result <- function(x,y){
  cor.temp <- cor.test(x, y)
  cor.temp.summary <- summary(cor.temp)
  list.cor <- list(cor.temp, 
                   cor.temp.summary
  )
  return(list.cor)
}


##########################################################################################


### Function Number 3
### Function entitled "lin.reg"
### input
##  x <- variable 1
##  y <- variable 2

### output
## linear regression p value


#lin.reg <- function(x,y){
#  temp.reg <- lm (x ~y)
 # summary.temp <- summary(temp.reg)
 # summary.temp$coefficients -> out.stats
 # plot(x,y)
 # abline(lm (x ~y))
 # list.temp <- list (summary.temp,out.stats,p.val=out.stats[2,4])
 # return(list.temp)
#}


lin.reg <- function(x,y){
  temp.reg <- lm (x ~y)
  summary.temp <- summary(temp.reg)
  summary.temp$coefficients -> out.stats
  plot(x,y)
  abline(lm (x ~y), main = "Biomarker correlation with continuous variable",xlab = "Continuous variable", ylab = "Biomarker expression")
  list.temp <- list (summary.temp,out.stats,p.val=out.stats[2,4])
  return(list.temp)
}


##########################################################################################


### Function Number 4
### Function entitled "km.log.test" to create kaplan meier survival curves for age 3=16 year old children treated with curative intent, MB
### input:
## time (PFS)
## event (Progression/relapse)
## marker

### output:
## km survival curve
## p values plotted on graph
## y axis with values as %
## legend and p value for survival analysis 


km.log.test <- function(time, event, marker, out.file = "none"){
  if(out.file!="none"){
    pdf(out.file)
  }
  km.PFS.incl <- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
  plot(km.PFS.incl, yaxt="n", col = c("red", "blue"),xlab = "time to progression/relapse (years)", ylab = "PFS (%)", xlim = c(0,10), main = "Biomarker expression and progression-free survival (PFS)",  lty = 1:2)
  PFS.names <- c("biomarker - high", "biomarker - low")
  legend (x="topright", PFS.names,  lty= 1:2, col = c("red","blue"))
  axis(2, at=pretty(event), lab=pretty(event) * 100, las=TRUE)
  PFS.incl.logrank <- survdiff(Surv(time, event) ~ marker)
  1 - pchisq(PFS.incl.logrank$chisq, length(PFS.incl.logrank$obs)-1) -> surv.p.val
  text(4,0.1,paste("p =",round(surv.p.val, 3)), pos = 4, cex = 1)
  if(out.file!="none"){
    dev.off()
  }
}


##########################################################################################


### Function Number 5
### Function entitled "km.log.test.sub" to create kaplan meier survival curves for age 3=16 year old children treated with curative intent, for chosen MB subgroup(s)
### input:
## time (PFS)
## event (Progression/relapse)
## marker

### output:
## km survival curve
## p values plotted on graph
## y axis with values as %
## legend and p value for survival analysis 
## title to be adjusted according to subgroup, currently for G3/G4 combined


km.log.test.sub <- function(time, event, marker, out.file = "none"){
  if(out.file!="none"){
    pdf(out.file)
  }
  km.PFS.incl <- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
  plot(km.PFS.incl, yaxt="n", col = c("red", "blue"),xlab = "time to progression/relapse (years)", ylab = "PFS (%)", xlim = c(0,10), main = "Biomarker expression and progression-free survival (PFS) in G3/G4",  lty = 1:2)
  PFS.names <- c("biomarker - high", "biomarker - low")
  legend (x="topright", PFS.names,  lty= 1:2, col = c("red","blue"))
  axis(2, at=pretty(event), lab=pretty(event) * 100, las=TRUE)
  PFS.incl.logrank <- survdiff(Surv(time, event) ~ marker)
  1 - pchisq(PFS.incl.logrank$chisq, length(PFS.incl.logrank$obs)-1) -> surv.p.val
  text(4,0.1,paste("p =",round(surv.p.val, 3)), pos = 4, cex = 1)
  if(out.file!="none"){
    dev.off()
  }
}




##########################################################################################


### Function number 6
### Function entitled "km.log.test.OS" to create kaplan-meier survival curves for OS for children treated age 3-16years with curative intent
### input:
## time
## event
## marker

### output:
## km survival curve
## p values plotted on graph
## y axis with values as %
## legend and p value for survival analysis 


km.log.test.OS <- function(time, event, marker, out.file = "none"){
  if(out.file!="none"){
    pdf(out.file)
  }
  km.OS.incl <- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
  plot(km.OS.incl,yaxt="n", col = c("red", "blue"),xlab = "overall survival (years)", ylab = "OS (%)", xlim = c(0,10), main = "Biomarker expression and overall survival (OS)",  lty = 1:2)
  OS.names <- c("biomarker - high", "biomarker - low")
  legend (x="topright", OS.names,  lty= 1:2, col = c("red","blue"))
  axis(2, at=pretty(event), lab=pretty(event) * 100, las=TRUE)
  OS.incl.logrank <- survdiff(Surv(time, event) ~ marker)
  1 - pchisq(OS.incl.logrank$chisq, length(OS.incl.logrank$obs)-1) -> surv.p.val
  text(4,0.1,paste("p =",round(surv.p.val, 3)), pos = 4, cex = 1)
  if(out.file!="none"){
    dev.off()
  }
}


##########################################################################################


### Function number 7
### Function entitled "cox.result.OS"  to produce cox regression hazard ratio for overall survival in cohort of children
### input
## time
## event (dichotomous)
## marker of interest (univariate)
## strata = different levels e.g risk group or those with other categorical differences

### output
## "summary.cox" result which includes the coefficients, hazard ratio (exp(coef)), number of patients (n), number of events (nevents)


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

##################################################################################

### Function Number 8

### Function entitled "km.log.test.EFS" to create kaplan meier survival curves for EFS for age 3-16 year old children treated with curative intent, MB
### Event is defined as death due to any cause, relapse/progression, second malignancy 
### input:
## time
## event
## marker

### output:
## km survival curve
## p values plotted on graph
## y axis with values as %
## legend and p value for survival analysis 


km.log.test.EFS <- function(time, event, marker, out.file = "none"){
  if(out.file!="none"){
    pdf(out.file)
  }
  km.EFS.incl <- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
  plot(km.EFS.incl, yaxt="n", col = c("red", "blue"),xlab = "event-free survival (years)", ylab = "EFS (%)", xlim = c(0,10), main = "Biomarker expression and event-free survival (EFS)",  lty = 1:2)
  EFS.names <- c("biomarker - high", "biomarker - low")
  legend (x="topright", EFS.names,  lty= 1:2, col = c("red","blue"))
  axis(2, at=pretty(event), lab=pretty(event) * 100, las=TRUE)
  EFS.incl.logrank <- survdiff(Surv(time, event) ~ marker)
  1 - pchisq(EFS.incl.logrank$chisq, length(EFS.incl.logrank$obs)-1) -> surv.p.val
  text(4,0.1,paste("p =",round(surv.p.val, 3)), pos = 4, cex = 1)
  if(out.file!="none"){
    dev.off()
  }
}



#################################################################################

### Function number 9

### Function named "km.log.test.OS.sub" to examine potential biomarker in subgroup of childhood MB e.g G3G4 combined
### Function creates kaplan-meier survival curves for OS for children treated age 3-16years with curative intent, within a specified molecular subgroup(s)
### input:
## time
## event
## marker

### output:
## km survival curve
## p values plotted on graph
## y axis with values as %
## legend and p value for survival analysis 
## title to be adjusted for the specific subgroup, currently written for G3G4 combined


km.log.test.OS.sub <- function(time, event, marker, out.file = "none"){
  if(out.file!="none"){
    pdf(out.file)
  }
  km.OS.incl <- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
  plot(km.OS.incl,yaxt="n", col = c("red", "blue"),xlab = "overall survival (years)", ylab = "OS (%)", xlim = c(0,10), main = "Biomarker expression and overall survival (OS) in G3/G4",  lty = 1:2)
  OS.names <- c("biomarker - high", "biomarker - low")
  legend (x="topright", OS.names,  lty= 1:2, col = c("red","blue"))
  axis(2, at=pretty(event), lab=pretty(event) * 100, las=TRUE)
  OS.incl.logrank <- survdiff(Surv(time, event) ~ marker)
  1 - pchisq(OS.incl.logrank$chisq, length(OS.incl.logrank$obs)-1) -> surv.p.val
  text(4,0.1,paste("p =",round(surv.p.val, 3)), pos = 4, cex = 1)
  if(out.file!="none"){
    dev.off()
  }
}


######################################################################################

### Function number 10

### function name "km.log.test.EFS.sub" to test biomarker within subgroup, in this case it will be run in G3/G4 combined
### same inputs and outputs as parent function "km.log.test.EFS", see below

### Event is defined as death due to any cause, relapse/progression, second malignancy 
### input:
## time
## event
## marker

### output:
## km survival curve
## p values plotted on graph
## y axis with values as %
## legend and p value for survival analysis
## title to be adjusted for the specific subgroup, currently written for G3G4 combined


km.log.test.EFS.sub <- function(time, event, marker, out.file = "none"){
  if(out.file!="none"){
    pdf(out.file)
  }
  km.EFS.incl <- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
  plot(km.EFS.incl, yaxt="n", col = c("red", "blue"),xlab = "event-free survival (years)", ylab = "EFS (%)", xlim = c(0,10), main = "Biomarker expression and event-free survival (EFS) in G3/G4",  lty = 1:2)
  EFS.names <- c("biomarker - high", "biomarker - low")
  legend (x="topright", EFS.names,  lty= 1:2, col = c("red","blue"))
  axis(2, at=pretty(event), lab=pretty(event) * 100, las=TRUE)
  EFS.incl.logrank <- survdiff(Surv(time, event) ~ marker)
  1 - pchisq(EFS.incl.logrank$chisq, length(EFS.incl.logrank$obs)-1) -> surv.p.val
  text(4,0.1,paste("p =",round(surv.p.val, 3)), pos = 4, cex = 1)
  if(out.file!="none"){
    dev.off()
  }
}



##########################################################################################


### Function number 11

### A file named "updatepData" which is to create a new dataframe using most uptodate files, which is then used in Function number 9, to assess biomarker performance within the clinical group

### input 
## x.data, current database
## meth7, worksheet with 7 molecular group data
## cytogen, worksheet with cytogenetic arm calls

### output
## a dataframe called test.pData
## a log file called "temp.log.txt" (summary data for n=802 patients, minimal information which may not be used in final analysis)
## a pdf file called "temp.pdf"
## temp.pdf contains qq plots for age as a continuous variable, summary statistics for the entire cohort (n=802), boxplot for age (continuous variable) across 4 molecular subgroups )



updatepData <-  function(x.data, meth7, cytogen, pdf.file = NULL, log.file = NULL){
  #### convert pData object columns into categorical variables
  
  if(!is.null(pdf.file)){
    pdf(pdf.file)
  }
  if(!is.null(log.file)){
    sink(log.file)
  }
  
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
  RTX <- pData$RTX_R
  
  ### CSI
  CSI <-pData$RTXCSI_R
  
  ### molecular subgroups
  meth <-pData$X450K_R
  
  meth7.cat <- meth7[c(1)]
  
  cytogen.q13.cat <- ifelse(cytogen$q13 == "Loss", "q13 loss", "no q13 loss") 
  
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
  
  ### curative group
  
  curative <- pData$Curative_R
  
  
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
  
  cat ("reading out the summary statistics for the entire cohort n=802", sep ="\n")
  
  ### assessing the data, parametric etc
  qqnorm(age.cont)
  plot(age.cont, col='red')
  
  
  ### Getting summary characteristics
  print(summary(age.cont))
  mean(age.cont, na.rm = T)
  median(age.cont, na.rm = T)
  
  
  cat ("age summary statistics n=802", sep ="\n")
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
                           Relapsetodeath,
                           Event,
                           EFS
  )
  rownames(test.pData) <- rownames(pData)
  
  ### need to add the 7 group methylation status using dataframes
  test.meth7 <- data.frame(meth7.cat)
  test.meth7$Sample_Name <- rownames(test.meth7)
  test.pData$meth7 <- test.meth7[match(rownames(test.pData), rownames(test.meth7)), ]$all.calls
  dev.off()
  sink()
  return(test.pData)
}


###############################################################################

### Function number 12
### aim is to input a candidate marker (goi, gene of interest, and then to run univariate analysis for association with outcome of interest and as a survival biomarker in age 3-16yrs treated with curative intent)
### input files
## data frame called test.pData, goi.vsd (biomarker file)

### output files
## marker.results.pdf file for graphics including heatmap association, boxplots, survival curves
## marker.results.txt file for record of workings and output statistics including chi squared analysis, logistic regression, linear regression


clinPathAssess <- function(test.pData,
                           goi.vsd,
                           pdf.file = NULL,
                           log.file = NULL
){
  
  #############################################
  ### setting up output files for log, messages and PDF files
  if(!is.null(log.file)){
    sink(log.file)
  }
  # cat ()
  if(!is.null(pdf.file)){
    pdf(pdf.file, width = 10, height = 10) 
  }
  
  
  #############################################
  
  
  ### goi.vsd #### name of a putative biomarker
  ### matching the expression data (goi data) with the initially compiled clinical data with important variables
  
  index <- match(names(goi.vsd), rownames(test.pData)) 
  matched.test.pData <- test.pData[index[!is.na(index)],] 
  is.vector(matched.test.pData)
  matched.goi.vsd <- goi.vsd[!is.na(index)] 
  matched.goi.vsd.cat <- ifelse(matched.goi.vsd>median(goi.vsd, na.rm = T), "high","low") 
  
  
  ### add cytogenetic data to existing data frame, NMB650 to be resolved but will not be included in survival analysis therefore add cytogenetic data at this stage
  ### duplicate NMB650 data emailed Ed 
  
  cytogen.q13.cat <- cytogen [c("SampleID", "q13")]
  
  ### need to convert q13 loss into loss, and rest into "no loss"
  
  cytogen.q13 <- ifelse(cytogen.q13.cat$q13 =="Loss", "q13 loss", "No q13 loss")
  
  ### make q13 loss dataframe
  
  cytogen.q13.df <- data.frame(cytogen.q13.cat[,-1], 
                               row.names=cytogen.q13.cat [,1],
                               cytogen.q13)
  
  
  matched.test.pData$q13loss <- cytogen.q13.df[match(rownames(matched.test.pData), rownames(cytogen.q13.df)),]$cytogen.q13
  
  #############################################
  
  ### summary data 
  
  ### visualising means of the goi in the matched dataset
  
  age.mean.goi <- plotmeans(matched.test.pData$age.cont ~ matched.test.pData$meth, data=matched.test.pData)
  
  ### summary data for all variables within the age 3-16 yo group, to construct survival cohort
  age.df <- data.frame(NMB, age.filter, age.cont)
  print(summary(age.df))
  
  ##############################################
  
  ### Chi squared analysis
  
  matched.pData.RTX <- ifelse(matched.test.pData$RTX=="Yes", "RTX", "No RTX")
  matched.pData.CSI <-ifelse(matched.test.pData$CSI=="Yes", "CSI", "No CSI")
  matched.pData.cat.infant <- ifelse(matched.test.pData$age.cat.infant=="TRUE", "infant", "non infant")
  
  cat ("assessing age.cat.infant against biomarker", sep ="\n")
  # list.age.cat.infant <- chi.sq(x = matched.test.pData$age.cat.infant, y = matched.goi.vsd.cat)
  list.age.cat.infant <- chi.sq(x = matched.pData.cat.infant, y = matched.goi.vsd.cat)
  print(list.age.cat.infant)
  
  cat ("assessing sex against biomarker", sep ="\n")
  list.sex <- chi.sq (x = matched.test.pData$sex, y= matched.goi.vsd.cat)
  print(list.sex)
  
  cat ("assessing metastatic status against biomarker", sep ="\n")
  list.mstatus <- chi.sq(x = matched.test.pData$mstatus, y=matched.goi.vsd.cat) 
  print(list.mstatus)
  
  cat ("assessing relapse status against biomarker", sep ="\n")
  list.relapse <- chi.sq(x = matched.test.pData$relapse, y = matched.goi.vsd.cat) 
  print(list.relapse)
  
  cat ("assessing resection status against biomarker", sep ="\n")
  list.resection <- chi.sq (x = matched.test.pData$resection, y = matched.goi.vsd.cat)
  print(list.resection)
  
  cat ("assessing 4 molecular group status against biomarker", sep ="\n")
  list.meth.4 <- chi.sq (x= matched.test.pData$meth, y = matched.goi.vsd.cat)
  print(list.meth.4)
  
  cat ("assessing 7 molecular group status against biomarker", sep ="\n")
  list.meth.7 <- chi.sq (x = matched.test.pData$meth7, y = matched.goi.vsd.cat)
  print(list.meth.7)
  
  
  cat ("assessing histopath status against biomarker", sep ="\n")
  list.path <- chi.sq (x = matched.test.pData$histopath, y = matched.goi.vsd.cat)
  print(list.path)
  
  cat("assessing LCA vs non LCA status against biomarker", sep ="\n")
  list.LCA <- chi.sq (x = matched.test.pData$LCA, y= matched.goi.vsd.cat)
  print(list.LCA)
  
  cat("assessing MYC amplified status against biomarker", sep ="\n")
  list.MYC <- chi.sq (x = matched.test.pData$MYC.cat, y = matched.goi.vsd.cat)
  print(list.MYC)
  
  cat("assessing MYCN amplified status against biomarker", sep ="\n")
  list.MYCN <- chi.sq (x = matched.test.pData$MYCN.cat, y = matched.goi.vsd.cat)
  print (list.MYCN)
  
  cat("assessing MYC/MYCN amplified status against biomarker", sep ="\n")
  list.MYCMYCN <- chi.sq(x= matched.test.pData$MYCMYCN.cat, y= matched.goi.vsd.cat)
  print (list.MYCMYCN)
  
  cat("assessing TP53 mutational status against biomarker", sep ="\n")
  list.TP53 <- chi.sq (x = matched.test.pData$TP53.cat, y = matched.goi.vsd.cat)
  print(list.TP53)
  
  cat("assessing TERT mutational status against biomarker", sep ="\n")
  list.TERT <- chi.sq (x = matched.test.pData$TERT.cat, y = matched.goi.vsd.cat)
  print(list.TERT)
  
  cat ("assessing cytogenetic q13 loss against biomarker", sep = "\n")
  list.q13loss <- chi.sq (x = matched.test.pData$q13loss, y = matched.goi.vsd.cat)
  print(list.q13loss)
  
  
  ### is the biomarker overrepresented in poor prognostic groups or those who received different therapy
  cat("assessing RTX delivery status against biomarker, is the biomarker overrepresented in poor prognostic groups or those who received different therapy", sep ="\n")
  # list.RTX <- chi.sq (x = matched.test.pData$RTX, y = matched.goi.vsd.cat)
  list.RTX <- chi.sq (x = matched.pData.RTX, y = matched.goi.vsd.cat)
  print(list.RTX)
  
  cat("assessing CSI delivery status against biomarker, is the biomarker overrepresented in poor prognostic groups or those who received different therapy", sep ="\n")
  # list.CSI <- chi.sq (x = matched.test.pData$CSI, y = matched.goi.vsd.cat)
  list.CSI <- chi.sq (x = matched.pData.CSI, y = matched.goi.vsd.cat)
  print(list.CSI)
  
  
  ### run Fisher's exact test on those where count < 5 e.g pathology "other"
  histopath.table <- table(as.factor(matched.test.pData$histopath), as.factor(matched.goi.vsd.cat))
  histopath.result <- fisher.test(histopath.table)
  histopath.result
  
  
 ###################################
  
  ### Correlation coefficients
  
  cat ("correlation coefficients for association between variables", sep ="\n")
  
  list.cor.age <- cor.result(x = matched.test.pData$age.cont, y = matched.goi.vsd)
  list(list.cor.age)
  print(list.cor.age[[1]])
  
  
  ##################################
  
  ###  linear regression
  
  lin.reg.age <- lin.reg(x= matched.test.pData$age.cont, y= matched.goi.vsd)
  lin.reg.age
  
  ### question: how do get 95% confidence interval, tried predict(lin.reg.age, interval = "confidence"), object needs to be as a dataframe
  # predict(lm(matched.test.pData$age.cont ~ matched.goi.vsd, data = matched.test.pData, interval = "confidence"))
  
  
  
  #################################
  
  ### Mann-whitney U (aka wilcoxon rank sum test) for non-parametric data
  
  ### age continuous 
  age.cont.wilcox <- wilcox.test(matched.test.pData$age.cont ~ matched.goi.vsd.cat, exact = F, correct = F)
  age.cont.wilcox
  
  
  ##################################
  ### logistic regression
  
  cat ("processing logistic regression individually", sep ="\n")
  
  ### cannot make logistic regression function work, error message "y values must be 0 <= y <= 1"
  
  # x <- matched.test.pData$age.cat.infant
  # y <- matched.goi.vsd
  # data.source <-  matched.test.pData
  
  # log.reg <- function(x,y,data.source){
  # log.reg.temp <- glm(y ~ x, family = binomial (link = 'logit'), data = data.source)
  # summary.temp <-summary(log.reg.temp)
  # boxplot.temp <- boxplot (y ~ x, col = c("red", "blue"), xlab = "x", ylab = "Biomarker expression", main = "Correlation between biomarker and variable x")
  # }
  
  
  ### age categorical
  cat ("processing logistic regression for age.cat.infant & biomarker, OR, 95% CI result and boxplots", sep ="\n")
  
  log.reg.age.cat <- glm(age.cat.infant ~ matched.goi.vsd, family = binomial(link= 'logit'), data=matched.test.pData)
  print(summary(log.reg.age.cat))
  print(cbind(OR=exp(coef(log.reg.age.cat)), exp(confint(log.reg.age.cat))))
  age.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$age.cat.infant, col = c("red", "blue"), xlab = "Infant", ylab = "Biomarker expression", main = "Correlation between biomarker and age (infant vs non infant)")
  
  
  ### sex 
  cat ("processing logistic regression for sex & biomarker, OR, 95% CI result and boxplots", sep ="\n")
  
  log.reg.sex <- glm (matched.test.pData$sex ~ matched.goi.vsd, family = binomial (link = 'logit'), data= matched.test.pData)
  print(summary(log.reg.sex))
  print(cbind(OR=exp(coef(log.reg.sex)), exp(confint(log.reg.sex))))
  sex.boxplot <- boxplot (matched.goi.vsd ~ matched.test.pData$sex, col = c("red", "blue"), xlab = "Gender", ylab = "Expression of biomarker", main = "Biomarker expression and gender")
  
  ### metastatic status
  cat ("processing logistic regression for metastatic status & biomarker, OR, 95% CI result and boxplots", sep ="\n")
  
  log.reg.mstatus <- glm(matched.test.pData$mstatus ~ matched.goi.vsd,  family = binomial(link= 'logit'), data=matched.test.pData)
  print(summary(log.reg.mstatus))
  print(cbind(OR=exp(coef(log.reg.mstatus)), exp(confint(log.reg.mstatus))))
  mstatus.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$mstatus, col = c("red", "blue"), xlab = "M status", ylab = "Biomarker expression", main = "Correlation between biomarker and metastatic status")
  
  
  ### relapse 
  cat ("processing logistic regression for relapse & biomarker, OR, 95% CI result and boxplots", sep ="\n")
  
  log.reg.relapse <- glm(matched.test.pData$relapse ~ matched.goi.vsd, family = binomial(link= 'logit'), data=matched.test.pData)
  print(summary(log.reg.relapse))
  print(cbind(OR=exp(coef(log.reg.relapse)), exp(confint(log.reg.relapse))))
  relapse.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$relapse, col = c("red", "blue"), xlab = "Relapse status", ylab = "Biomarker expression",  main = "Correlation between biomarker and relapse")
  
  
  
  ### resection
  cat ("processing logistic regression for resection & biomarker, OR, 95% CI result and boxplots", sep ="\n")
  
  log.reg.resection <- glm (matched.test.pData$resection ~ matched.goi.vsd, family = binomial(link= 'logit'), data=matched.test.pData)
  print(summary(log.reg.resection))
  print(cbind(OR=exp(coef(log.reg.resection)), exp(confint(log.reg.resection))))
  resection.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$resection, col = c("red", "blue"), xlab = "Resection status", ylab = "Biomarker expression", main = "Correlation between biomarker and resection status")
  
 
  
  ### histopath
  cat ("processing logistic regression for histopath & biomarker, OR, 95% CI result and boxplots", sep ="\n")
  
  log.reg.histopath <- glm (matched.test.pData$histopath ~ matched.goi.vsd, family = binomial(link='logit'), data=matched.test.pData)
  print(summary(log.reg.histopath))
  print(cbind(OR=exp(coef(log.reg.histopath)), exp(confint(log.reg.histopath))))
  
  histopath.pw <- pairwise.t.test(matched.goi.vsd, matched.test.pData$histopath)
  histopath.pw
  histopath.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$histopath, col = c("red", "blue"), xlab = "Histopathology subtype", ylab = "Biomarker expression", main = "Correlation between biomarker and histopathology")
  
  
  ### visualise relationship between biomarker and LCA pathology
  cat ("processing logistic regression for LCA pathology & biomarker, OR, 95% CI result and boxplots", sep ="\n")
  
  log.reg.LCA <- glm (matched.test.pData$LCA ~ matched.goi.vsd, family = binomial(link='logit'), data=matched.test.pData)
  print(summary(log.reg.LCA))
  print(cbind(OR=exp(coef(log.reg.LCA)), exp(confint(log.reg.LCA))))
  
  LCA.pw <- pairwise.t.test(matched.goi.vsd, matched.test.pData$LCA, level = matched.test.pData$LCA == "Non LCA")
  LCA.pw
  LCA.boxplot <- boxplot (matched.goi.vsd~ matched.test.pData$LCA,col=c("red","blue"),  xlab = "LCA pathology", ylab = "Biomarker expression", main = "Correlation between biomarker and LCA pathology")
  
  
  ### MYC.cat
  cat ("processing logistic regression for MYC amplification & biomarker, OR, 95% CI result and boxplots", sep ="\n")
  
  log.reg.MYC <- glm (matched.test.pData$MYC.cat ~ matched.goi.vsd, family = binomial(link='logit'), data=matched.test.pData)
  print(summary(log.reg.MYC))
  print(cbind(OR=exp(coef(log.reg.MYC)), exp(confint(log.reg.MYC))))
  MYC.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$MYC.cat, col = c("Red", "Blue"), xlab = "MYC amplification", ylab = "Biomarker expression", main = "Correlation between biomarker and MYC amplification")
  
  ### MYCN.cat
  
  log.reg.MYCN <- glm (matched.test.pData$MYCN.cat ~ matched.goi.vsd, family = binomial (link='logit'), data=matched.test.pData)
  print(summary(log.reg.MYCN))
  MYCN.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$MYCN.cat, col = c("Red", "Blue"), xlab = "amplification", ylab = "Biomarker expression", main = "Correlation between biomarker and MYCN amplification")
  

  ### combined MYC / MYCN amplification group
  
  log.reg.MYCMYCN <- glm (matched.test.pData$MYCMYCN.cat ~ matched.goi.vsd, family = binomial (link='logit'), data=matched.test.pData)
  print(summary(log.reg.MYCMYCN))
  
  MYCMYCN.pw <- pairwise.t.test(matched.goi.vsd, matched.test.pData$MYCMYCN.cat)
  MYCMYCN.pw
  MYCMYCNN.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$MYCMYCN.cat, col = c("Red", "Blue"), xlab = "MYC or MYCN amplification", ylab = "Biomarker expression", main = "Correlation between biomarker and MYC/MYCN amplification")
  
  
  ### TP53
  
  log.reg.TP53 <- glm (matched.test.pData$TP53.cat ~ matched.goi.vsd,family = binomial (link='logit'), data=matched.test.pData)
  print(summary(log.reg.TP53))
  TP53.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$TP53.cat, col = c("Red", "Blue"), xlab = "TP53 mutational status", ylab = "Biomarker expression", main = "Correlation between biomarker and TP53 mutational status")
  
  
  ### q13 loss
  log.reg.q13loss <- glm (matched.test.pData$q13loss ~ matched.goi.vsd,family = binomial (link='logit'), data=matched.test.pData)
  print(summary(log.reg.TP53))
  q13.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$q13loss, col = c("Red", "Blue"), xlab = "Cytogenetic arm q13 status", ylab = "Biomarker expression", main = "Correlation between biomarker and q13 loss")
  
  
  
  ### additional subgroup specific tests
  
  ### TERT (may need to specific subgroup)
  
  log.reg.TERT <- glm (matched.test.pData$TERT.cat ~ matched.goi.vsd, family = binomial (link = 'logit'), data = matched.test.pData)
  print(summary(log.reg.TP53))
  TERT.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$TERT.cat, col = c("Red", "Blue"), xlab = "TERT status", ylab = "Biomarker expression", main = "Correlation between biomarker and TERT status")
  
  
  ##################################################
  
  ### logistic regression and molecular subgroups
  
  ### includes some more detailed assumption testing and post-hoc analysis
  cat ("logistic regression and association with molecular subgroups")
  
  log.reg.meth <- glm (matched.test.pData$meth ~ matched.goi.vsd, family = binomial(link='logit'), data = matched.test.pData)
  print(summary(log.reg.meth))                  
  
  
  log.reg.meth7 <- glm (matched.test.pData$meth7 ~ matched.goi.vsd, family = binomial(link='logit'), data = matched.test.pData)
  print(summary(log.reg.meth7)) 
  
  
  ### visualise distribution of biomarker in cohort
  
  # message ("visualisation of biomarker and relationship with methylation")
  
  # qqnorm(matched.goi.vsd)
  
  ### visualise relationship between biomarker and methylation groups
  
  meth.boxplot <- boxplot(matched.goi.vsd~matched.test.pData$meth, col=c("yellow","green","red","blue"), xlab = "Methylation subgroup", ylab = "Biomarker expression", main = "Correlation between biomarker and 4 molecular subgroups")
  meth7.boxplot <- boxplot(matched.goi.vsd~matched.test.pData$meth7, col=c("yellow","green","red","blue"),  xlab = "Methylation subgroup", ylab = "Biomarker expression", main = "Correlation between biomarker and 7 molecular subgroups")
  
  
  ### if biomarker is normally distributed, can use ANOVA (one-way ANOVA)
  
  # meth7.aov <- aov (matched.goi.vsd ~ matched.test.pData$meth7, data = matched.test.pData)
  # summary(meth7.aov)
  # plot(meth7.aov)
  
  
  ### post-hoc tests for methylation
  
  ### pairwise t-test to determine where the difference lies between the groups
  meth.pw <- pairwise.t.test(matched.goi.vsd, matched.test.pData$meth)
  meth.pw
  
  
  ###################################
  
  
  ### summary data for each methylation subgroup (n=4), no age restriction for first analyses
  
  cat ("processing summary stats for 4 molecular subgroups, no age restriction", sep ="\n")
  
  ### G3
  
  G3 <- matched.test.pData$meth =="G3" 
  
  G3.group <- matched.test.pData [G3, ]
  print(summary(G3.group))
  
  ## G4
  
  G4 <- matched.test.pData$meth =="G4"
  G4.group <- matched.test.pData [G4, ]
  print(summary(G4.group))
  
  
  ### SHH
  
  SHH <- matched.test.pData$meth =="SHH"
  SHH.group <- matched.test.pData [SHH, ]
  print(summary(SHH.group))
  
  ### WNT
  
  WNT <- matched.test.pData$meth == "WNT"
  WNT.group <- matched.test.pData [WNT, ]
  print(summary(WNT.group))
  
  cat ("processing summary stats for 7 molecular subgroups, no age restriction", sep = "\n")
  
  ### G3 subgroups
  
  G3.high.group <- matched.test.pData$meth7 =="Grp3_HighRisk"
  G3.high.data <- matched.test.pData [G3.high.group, ]
  print(summary(G3.high.data))
  
  G3.low.group <- matched.test.pData$meth7 =="Grp3_LowRisk"
  G3.low.data <- matched.test.pData [G3.low.group, ]
  print(summary(G3.low.data))
  
  
  ### G4 subgroups
  
  G4.high.group <- matched.test.pData$meth7 =="Grp4_HighRisk"
  G4.high.data <- matched.test.pData [G4.high.group, ]
  print(summary(G4.high.data))
  
  G4.low.group <- matched.test.pData$meth7 =="Grp4_LowRisk"
  G4.low.data <- matched.test.pData [G4.low.group, ]
  print(summary(G4.low.data))
  
  ### SHH
  SHH.inf.group <- matched.test.pData$meth7 == "SHH_Inf"
  SHH.inf.data <- matched.test.pData [SHH.inf.group, ]
  print (summary (SHH.inf.data))
  
  SHH.old.group <- matched.test.pData$meth7 =="SHH_Old"
  SHH.old.data <- matched.test.pData [SHH.old.group, ]
  print (summary (SHH.old.data))
  
  
  #############################################
  
  ### New dataframes for survival analyses
  
  ### restrict analysis to age 3-16 yo, treated with curative intent including CSI
  
  cat ("create dataframe for age 3-16 years, curative intent called age.incl.df", sep = "\n")
  
  Age.incl <- matched.test.pData$age.filter== "TRUE"
  Age.incl.df <- matched.test.pData [Age.incl,]
  print(summary(Age.incl.df))
  
  
  ### compare to prior dataframes to check accuracy of new dataframe
  
  cat ("comparing with previous data frames for accuracy", sep = "\n")
  summary(test.pData)
  summary (matched.test.pData)
  
  
  
  ####################################
  
  cat ("processing summary stats for 3-16 yo for 4 molecular subgroups", sep = "\n")
  
  ### G3, aged 3-16 years
  G3.incl <- Age.incl.df$meth =="G3"
  G3.group.incl <- Age.incl.df [G3.incl, ]
  print (summary(G3.group.incl))
  nrow(G3.group.incl)
  
  ### G4, aged 3-16 years
  G4.incl <- Age.incl.df$meth =="G4"
  G4.group.incl <- Age.incl.df [G4.incl, ]
  print(summary(G4.group.incl))
  nrow(G4.group.incl)
  
  ### SHH, aged 3-16 years
  SHH.incl <- Age.incl.df$meth == "SHH"
  SHH.group.incl <- Age.incl.df [SHH.incl, ]
  print(summary (SHH.group.incl))
  nrow(SHH.group.incl)
  
  ### WNT, aged 3-16 years
  WNT.incl <- Age.incl.df$meth == "WNT"
  WNT.group.incl <- Age.incl.df [WNT.incl, ]
  print(summary (WNT.group.incl))
  nrow (WNT.group.incl)
  
  
  ### defining features of 7 molecular groups, group 3 and 4 subgroups
  
  cat ("processing summary stats for 7 molecular subgroup data, age 3-16 years", sep = "\n")
  
  G3.high <- Age.incl.df$meth7 == "Grp3_HighRisk"
  G3.high.incl <- Age.incl.df [G3.high, ]
  print(summary(G3.high.incl))
  
  G3.low <- Age.incl.df$meth7 == "Grp3_LowRisk"
  G3.low.incl <- Age.incl.df [G3.low, ]
  print(summary(G3.low.incl))
  
  
  G4.high <- Age.incl.df$meth7 =="Grp4_HighRisk"
  G4.high.incl <- Age.incl.df [G4.high, ]
  print(summary (G4.high.incl))
  
  G4.low <- Age.incl.df$meth7 =="Grp4_LowRisk"
  G4.low.incl <- Age.incl.df [G4.low, ]
  print(summary (G4.low.incl))
  
  SHH.inf <- Age.incl.df$meth7 == "SHH_Inf"
  SHH.inf.incl <- Age.incl.df [SHH.inf, ]
  print(summary (SHH.inf.incl))
  
  SHH.old <- Age.incl.df$meth7 =="SHH_Old"
  SHH.old.incl <- Age.incl.df [SHH.old, ]
  print(summary (SHH.old.incl))
  
  
  
  ##########################
  
  ### survival analysis using functions from source file "clinical_data_functions.R"
  
  ##########################
  
  ### creating dataframe for survival analysis
  
  cat ("restrict survival analysis for age 3-16 years, curative intent", sep = "\n")
  cat ("creating matched data frame", sep = "\n")
  
  ### check whether to retain the two dataframes Age.incl.df and matched.test.incl.pData 
  
  index.incl <- match(names(goi.vsd), rownames(Age.incl.df)) 
  matched.test.incl.pData <- Age.incl.df[index.incl[!is.na(index.incl)],] 
  is.vector(matched.test.incl.pData)
  matched.goi.vsd.incl <- goi.vsd[!is.na(index.incl)] 
  matched.goi.vsd.cat.incl <- ifelse(matched.goi.vsd.incl>median(goi.vsd, na.rm = T), "high","low")
  
  ### G3 and G4 combined dataframe
  
  cat ("creating combined dataframe to assess biomarker in G3 G4 combined group, for survival cohort, aged 3-16 years, curative intent", sep = "\n")
  
  G3.match <- matched.test.incl.pData$meth=="G3" 
  G4.match <- matched.test.incl.pData$meth=="G4"
  
  G3.match.df <- matched.test.incl.pData [G3.match, ]
  G4.match.df <- matched.test.incl.pData [G4.match, ]
  
  G3G4.match.df <- rbind(G3.match.df, G4.match.df)
  
  # nrow(G3G4.match.df)
 
  index.incl <- match(names(goi.vsd), rownames(G3G4.match.df)) 
  matched.G3G4.incl.pData <- G3G4.match.df[index.incl[!is.na(index.incl)],] 
  is.vector(matched.G3G4.incl.pData)
  matched.goi.vsd.G3G4.incl <- goi.vsd[!is.na(index.incl)] 
  matched.goi.vsd.cat.G3G4.incl <- ifelse(matched.goi.vsd.G3G4.incl>median(goi.vsd, na.rm = T), "high","low")
  
  print(summary(matched.G3G4.incl.pData$meth7))
  
  
  ### creating binary relapse variables labelled 0,1 for event analysis
  
  relapse.bin <- ifelse(matched.test.pData$relapse == "relapse", 1, 0)
  relapse.bin.incl <- ifelse(matched.test.incl.pData$relapse == "relapse", 1,0)
  relapse.G3G4.bin.incl <- ifelse(matched.G3G4.incl.pData$relapse == "relapse", 1,0)
  
  OS.cat.bin <- ifelse(matched.test.pData$OS.cat == "Dead", 1, 0) 
  OS.cat.bin.incl <- ifelse(matched.test.incl.pData$OS.cat == "Dead", 1,0)
  OS.G3G4.bin.incl <- ifelse(matched.G3G4.incl.pData$OS.cat == "Dead", 1,0)
  
  EFS.cat.bin <- ifelse(matched.test.pData$Event == "Event", 1, 0)
  EFS.cat.bin.incl <- ifelse(matched.test.incl.pData$Event == "Event", 1,0)
  EFS.G3G4.bin.incl <- ifelse(matched.G3G4.incl.pData$Event == "Event", 1,0)
  
  
  #### PFS
  
  cat("running km survival curve for PFS and biomarker, graphical output to PDF", sep = "\n")
  
  km.log.test.all <- km.log.test(time = matched.test.incl.pData$PFS, event = relapse.bin.incl, marker = matched.goi.vsd.cat.incl)
  
  km.log.test.G3G4 <- km.log.test.sub(time = matched.G3G4.incl.pData$PFS, event = relapse.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl)
  
  ### cox regression analysis
  
  cat("cox regression analysis, PFS, biomarker, age 3-16 years treated with curative therapy", sep ="\n")
  
  cox.relapse.incl <- coxph (Surv(matched.test.incl.pData$PFS, relapse.bin.incl) ~ matched.goi.vsd.cat.incl)
  summary(cox.relapse.incl)$logtest
  summary(cox.relapse.incl)
  
  cox.relapse.incl.G3G4 <- coxph (Surv(matched.G3G4.incl.pData$PFS, relapse.G3G4.bin.incl) ~ matched.goi.vsd.cat.G3G4.incl)
  print(summary(cox.relapse.incl.G3G4)$logtest)
  print(summary(cox.relapse.incl.G3G4))
 
  ##########################################
  
  ### OS 
  
  cat ("km analysis for OS for children aged 3-16 years, treated with curative intent", sep = "\n")
  
  km.log.test.OS.all <- km.log.test.OS(time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.cat.incl )
  
  km.log.test.OS.G3G4 <- km.log.test.OS.sub(time = matched.G3G4.incl.pData$Followup, event = OS.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl)
  
  
  ### cox regression analysis
  
  cat("cox regression on age 3-16 years, curative", sep = "\n")
  
  
  cox.result.OS.all <- cox.result.OS (time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.cat.incl)
  print(summary(cox.result.OS.all))
  
  cox.result.OS.G3G4 <- cox.result.OS (time = matched.G3G4.incl.pData$Followup, event = OS.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl)
  print(summary(cox.result.OS.G3G4))
  
  
  
  #####################################
  
  ### EFS: 
  
  cat ("km analysis for EFS for children aged 3-16 years, treated with curative intent", sep = "\n")
  
  km.log.test.EFS.all <- km.log.test.EFS(time = matched.test.incl.pData$EFS, event = EFS.cat.bin.incl, marker = matched.goi.vsd.cat.incl)
  km.log.test.EFS.G3G4 <- km.log.test.EFS.sub(time = matched.G3G4.incl.pData$EFS, event = EFS.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl) 
  
  
  ### cox regression analysis
  
  cat ("cox regression for EFS on age 3-16 years", sep = "\n")
  
  cox.EFS.incl <- coxph (Surv(matched.test.incl.pData$EFS, EFS.cat.bin.incl) ~ matched.goi.vsd.cat.incl)
  summary(cox.EFS.incl)$logtest
  print(summary(cox.EFS.incl))
  
  
  cox.EFS.incl.G3G4 <- coxph (Surv(matched.G3G4.incl.pData$EFS, EFS.G3G4.bin.incl) ~ matched.goi.vsd.cat.G3G4.incl)
  summary(cox.EFS.incl.G3G4)$logtest
  print(summary(cox.EFS.incl.G3G4))
  
  
  ###################################
  
  sink()
  dev.off()
  
  ###### need to return some objects with results
  list() -> res
  return(res)
}





