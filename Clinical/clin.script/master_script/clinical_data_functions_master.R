library(NMF)
library(gplots)
library(car)
library(stats)
library(survival)

### File information:

### This is a list of functions that were written for the "clinical_data_4.R" script analysis and then updated for the clinical_data_master.R script
### Many of these functions are generic for survival analysis
### Aims of this document are to provide tools to perform univariate analysis in a group of children treated for medulloblastoma. 
### All kaplan-meier survival analyses provide graphical output with labelled axes and p value
### The survival analysis is designed in a cohort of children who received therapy with curative intent (including cranio-spinal irradiation), age 3-16 years

### Author: Dr Marion Mateos
### Date: July 14 2017 - September 20 2017

### names of functions:
### chi.sq
### cor.result
### lin.reg
### km.log.test
### km.log.test.OS
### cox.result.OS
### cox.result.surv
### cox.multivar.surv_5
### cox.multivar.surv.PNET5_4
### cox.multivar.surv_3
### km.log.test.EFS
### logisticRegression
### updatepData
### clinPathAssess
### cox.dataframe
### annotate.HTseq.IDs
### log.reg.dataframe
### gp.style.filter

##############################################################################################

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
  # x = matched.test.pData$resection 
  # y = matched.goi.vsd.cat
  table.temp <- table(x, y) ### check how to label x and y so outputted in list.temp
  table.temp.perc <- prop.table(table.temp)*100
  summary.table(table.temp) ### note that the pvalue for independence of all factors is not the p value for the Pearson's chi-squared test, verified 140917
  # chi.test.temp <- chisq.test(table.temp) 
  chi.test.temp <- try(chisq.test(table.temp), silent = T) 
  chi.test.temp.stat <- try(c(stat=chi.test.temp$statistic, p.value=chi.test.temp$p.value), silent = T) ### added 26/9/17
  chi.test.temp.res <- try(chi.test.temp$residuals, silent = T)                                         ### added 26/9/17
  try(aheatmap(chi.test.temp$residuals, Rowv=NA, Colv = NA), silent = T)
  list.temp <- list  (p.value = chi.test.temp$p.value,                                                   ### subset chi.test.temp$p.value may not exist
                      chi.test = chi.test.temp.stat, 
                      table.temp = table.temp, ## does not work if put in variable x or y, here
                      table.temp.perc = table.temp.perc,
                      chi.test.temp = chi.test.temp,
                      chi.test.temp.res = chi.test.temp.res
  )
  
  return(list.temp)
}


###########################################################################################

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


#####################################################################################

### Function Number 3
### Function entitled "lin.reg"
### input
##  x <- variable 1
##  y <- variable 2

### output
## linear regression p value


lin.reg <- function(x,y){
  temp.reg <- lm (x ~y)
  summary.temp <- summary(temp.reg)
  summary.temp$coefficients -> out.stats
  plot(x,y)
  abline(lm (x ~y))
  list.temp <- list (summary.temp,out.stats,p.val=out.stats[2,4])
  return(list.temp)
}

###########################################################################################


### Function Number 4
### Function entitled "km.log.test" to create kaplan meier survival curves for age 3=16 year old children treated with curative intent, MB
### input:
## time
## marker
## event

### output:
## km survival curve
## p values plotted on graph
## y axis with values as %
## legend and p value for survival analysis 

### example inputs
# time <- matched.test.incl.pData$PFS
# event = relapse.bin.incl
# marker = matched.goi.vsd.cat.incl

# time <- matched.G3G4.incl.pData$PFS  
# event <- relapse.G3G4.bin.incl  
#marker <- matched.goi.vsd.cat.G3G4.incl

km.log.test <- function(time, event, marker, out.file = "none"){
  if(out.file!="none"){
    pdf(out.file)
  }
  km.PFS.incl <- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
  # plot(km.PFS.incl, yaxt="n", col = c("red", "blue"),xlab = "time to progression/relapse (years)", ylab = "PFS (%)", xlim = c(0,10), main = "Biomarker expression and progression-free survival (PFS)",  lty = 1:2)
   plot(km.PFS.incl, yaxt="n", col = c("red", "blue"),xlab = "time to progression/relapse (years)", ylab = "PFS (%)", xlim = c(0,10), main = "CDK6 and PFS",  lty = 1:2) ### "Marker and progression-free survival(PFS)"
  PFS.names <- c("biomarker - high", "biomarker - low")
  legend (x="topright", PFS.names,  lty= 1:2, col = c("red","blue"))
  axis(2, at=pretty(event), lab=pretty(event) * 100, las=TRUE)
  PFS.incl.logrank <- survdiff(Surv(time, event) ~ marker)
  1 - pchisq(PFS.incl.logrank$chisq, length(PFS.incl.logrank$obs)-1) -> surv.p.val.PFS
  text(4,0.1,paste("p =",round(surv.p.val.PFS, 3)), pos = 4, cex = 1) ### if the p value does not paste, may be due to rounding, e.g was 3, set at 5 dp (25/10/18)
  #assign(paste0("Survival_pval_", marker, surv.p.val.PFS)) ### added 140917
  PFS.surv.table <- summary(km.PFS.incl)
  return (list(PFS.p.val = surv.p.val.PFS,
               PFS.surv.table)) ### added 140917
  if(out.file!="none"){
    dev.off()
  }
}


############################################################################################


### Function number 5
### Function entitled "km.log.test.OS"to create kaplan-meier survival curves for OS for children treated age 3-16years with curative intent
### input:
## time
## event
## marker

### output:
## km survival curve
## p values plotted on graph
## y axis with values as %
## legend and p value for survival analysis 


### example inputs
# time <- matched.test.incl.pData$Followup
# event = OS.cat.bin.incl
# marker = matched.goi.vsd.cat.incl
km.log.test.OS <- function(time, event, marker, out.file = "none"){
  if(out.file!="none"){
    pdf(out.file)
  }
  km.OS.incl <- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
  # plot(km.OS.incl,yaxt="n", col = c("red", "blue"),xlab = "overall survival (years)", ylab = "OS (%)", xlim = c(0,10), main = "Biomarker expression and overall survival (OS)",  lty = 1:2)
   plot(km.OS.incl,yaxt="n", col = c("red", "blue"),xlab = "overall survival (years)", ylab = "OS (%)", xlim = c(0,10), main = "CDK6 and overall survival (OS)",  lty = 1:2)
  OS.names <- c("biomarker - high", "biomarker - low")
  legend (x="topright", OS.names,  lty= 1:2, col = c("red","blue"))
  axis(2, at=pretty(event), lab=pretty(event) * 100, las=TRUE)
  OS.incl.logrank <- survdiff(Surv(time, event) ~ marker)
  1 - pchisq(OS.incl.logrank$chisq, length(OS.incl.logrank$obs)-1) -> surv.p.val.OS
  text(4,0.1,paste("p =",round(surv.p.val.OS, 3)), pos = 4, cex = 1)
  OS.surv.table <- summary(km.OS.incl)
  return (list(OS.p.val = surv.p.val.OS, 
               OS.surv.table))
  if(out.file!="none"){
    dev.off()
  }
}



#######################################################################################

### Function number 6
### Function entitled "cox.result.OS"  to produce cox regression hazard ratio for overall survival in cohort of children
### input
## time
## event (dichotomous)
## marker of interest (univariate)
## strata = different levels e.g risk group or those with other categorical differences

### output
## "summary.cox" result which includes the coefficients, hazard ratio (exp(coef)), number of patients (n), number of events (nevents)

cox.result.OS <- function (time, event, marker, strata = NULL, data)  
{
  if(is.null(strata)){
    cox.temp <- coxph (Surv(time, event)~marker, data= data)
  }else{
    cox.temp <- coxph (Surv(time, event)~marker, data= data)
  }
  #summary.cox <- c(rownames(summary(cox.temp)$coefficients),summary(cox.temp)$coefficients,
                  # summary(cox.temp)$n,
                  # summary(cox.temp)$nevent)
 # names(summary.cox) <- c("marker_name",colnames(summary(cox.temp)$coefficients), "n","nevents")
  cox.p.val <- summary(cox.temp)$logtest[3] ### p value can also be called within cox.temp$coefficients
  cox.HR <- summary(cox.temp)$conf.int[1] ### called within cox.temp$coefficients
  cox.lower.95CI <- summary(cox.temp)$conf.int[3] 
  cox.upper.95CI <- summary(cox.temp)$conf.int[4]
  cox.Zscore <- summary(cox.temp)$coefficients[4]
  cox.n <-summary(cox.temp)$n
  cox.nevent <-summary(cox.temp)$nevent
  summary.OS.cox <- list(cox.OS.p.val = cox.p.val,cox.OS.HR = cox.HR,  cox.OS.lower.95CI = cox.lower.95CI, cox.OS.upper.95CI =cox.upper.95CI, cox.OS.Zscore = cox.Zscore, n.OS = cox.n, n.event.OS = cox.nevent)
  
  return (summary.OS.cox)
}



####################################################################
### Function number 6a: universal version of function entitled "cox.result.OS" 

### to produce cox regression hazard ratio for survival in cohort of children
### input
## time
## event (dichotomous)
## marker of interest (univariate)
## strata = different levels e.g risk group or those with other categorical differences

### output
## "summary.cox" result which includes the coefficients, hazard ratio (exp(coef)), number of patients (n), number of events (nevents)

### example inputs to examine outputs
#  time <- matched.test.incl.pData$PFS 
# event <- relapse.bin.incl
# marker = matched.goi.vsd.cat.incl
# data <- matched.test.incl.pData


cox.result.surv <- function (time, event, marker, strata = NULL, data)  
{
  if(is.null(strata)){
    cox.temp <- coxph (Surv(time, event)~marker, data = data) ### removed data = data, still did not work
  }else{
    cox.temp <- coxph (Surv(time, event)~marker, data = data) ### removed data = data, still did not work
  }
  cox.p.val <- summary(cox.temp)$logtest[3] ### p value also in cox.temp$coefficients, this current logtest is the likelihood ratio p value which is correct
  cox.HR <- summary(cox.temp)$conf.int[1] ### called within cox.temp$coefficients
  cox.lower.95CI <- summary(cox.temp)$conf.int[3] 
  cox.upper.95CI <- summary(cox.temp)$conf.int[4]
  cox.Zscore <- summary(cox.temp)$coefficients[4] ### added this in to access Z score
  cox.n <-summary(cox.temp)$n
  cox.nevent <-summary(cox.temp)$nevent
  summary.cox <- list(cox.pval = cox.p.val,cox.HR = cox.HR, cox.lower.95CI = cox.lower.95CI, cox.upper.95CI =cox.upper.95CI, cox.Zscore = cox.Zscore, n = cox.n, n.event = cox.nevent)
  #names(summary.cox) <- c("marker_name",cox.p.val, cox.HR, cox.lower.95CI, cox.upper.95CI, "n","nevents")
    
  ### alternative output that generates everything except for 95%CI
  #summary.cox <- c(rownames(summary(cox.temp)$coefficients),summary(cox.temp)$coefficients,
                   #summary(cox.temp)$n,
                   #summary(cox.temp)$nevent)
  #names(summary.cox) <- c("marker_name","coef", "cox.HR", "se(coef)", "z", "p.val", "n","nevents")
  return (summary.cox)
}



###############################################################################################

### Function number 6b: multivariable cox model

### input variables according to a classic cox regression model Surv(time,event)~ marker, data

# time <- matched.test.incl.pData$PFS
# event <- relapse.bin.incl
# marker = matched.goi.vsd.cat.incl ### must use this and not an individual ENSG identity number. (### specify ENSG within clinPathAssess function first)

### Factors A to G to include in the multivariable cox regression
# FacA <- matched.test.incl.pData$LCA
# FacB <- matched.test.incl.pData$MYCMYCN.cat ### changed to MYCMYCN.cat rather than individual MYC.cat and MYCN.cat as per DW 4/10/17 & 25/10/18
# FacC <- matched.test.incl.pData$mstatus
### FacD <- matched.test.incl.pData$resection ### 31/10/18 resection removed
# FacD <- matched.test.incl.pData$q13loss
# FacE <- matched.test.incl.pData$meth7.cat
### FacF <- matched.test.incl.pData$TP53.cat
# data <- matched.test.incl.pData

###################################
### updated function below 14/11/17

 # cox.multivar.surv_7 <- function (time, event, marker, FacA, FacB, FacC, FacD, FacE, FacF, FacG, strata = NULL, data) {
 # if(is.null(strata)){
  #  cox.temp <- coxph(Surv(time, event)~marker + FacA + FacB + FacC + FacD + FacE + FacF + FacG, data=data)
 # }else {
 #   cox.temp <- coxph(Surv(time, event)~marker + FacA + FacB +FacC +FacD + FacE + FacF + FacG, data=data)
 # }  
 # cox.p.val <- summary(cox.temp)$coefficients[1,5] ### updated 14/11
 # cox.HR <- summary(cox.temp)$coefficients[1,2] ### updated 14/11
 # cox.lower.95CI <- summary(cox.temp)$conf.int[1,3] ### as now multivariate, therefore need to access 1st row results
 # cox.upper.95CI <- summary(cox.temp)$conf.int[1,4]
 # cox.Zscore <- summary(cox.temp)$coefficients[1,4] ### added this in to access Z score
 # cox.n <-summary(cox.temp)$n
 # cox.nevent <-summary(cox.temp)$nevent
#  summary.cox <- list(cox.pval = cox.p.val,cox.HR = cox.HR, cox.lower.95CI = cox.lower.95CI, cox.upper.95CI =cox.upper.95CI, cox.Zscore = cox.Zscore, n = cox.n, n.event = cox.nevent)
#  return (summary.cox)
# }
 
### updated function 31/10/18

 cox.multivar.surv_5 <- function (time, event, marker, FacA, FacB, FacC, FacD, FacE, strata = NULL, data) {
   if(is.null(strata)){
     cox.temp <- coxph(Surv(time, event)~marker + FacA + FacB + FacC + FacD + FacE, data=data)
   }else {
     cox.temp <- coxph(Surv(time, event)~marker + FacA + FacB +FacC +FacD + FacE, data=data)
   }  
   cox.p.val <- summary(cox.temp)$coefficients[1,5] 
   cox.HR <- summary(cox.temp)$coefficients[1,2] 
   cox.lower.95CI <- summary(cox.temp)$conf.int[1,3] 
   cox.upper.95CI <- summary(cox.temp)$conf.int[1,4]
   cox.Zscore <- summary(cox.temp)$coefficients[1,4] 
   cox.n <-summary(cox.temp)$n
   cox.nevent <-summary(cox.temp)$nevent
   summary.cox <- list(cox.pval = cox.p.val,cox.HR = cox.HR, cox.lower.95CI = cox.lower.95CI, cox.upper.95CI =cox.upper.95CI, cox.Zscore = cox.Zscore, n = cox.n, n.event = cox.nevent)
   return (summary.cox)
 }
 
 

###############################################################################################
### Function number 6c for PNET5 survival markers  (25/10/18 changes - gender removed, 31/10/18 WNT vs non WNT instead of 4 methylation groups for overall analysis and only including those that are significant in our cohort)
### Factors Fac A - FacD

# FacA <- matched.test.incl.pData$LCA
# FacB <- matched.test.incl.pData$MYCMYCN.cat ### keep this in algorithm because of flow diagram includes MYC.cat
### FacC <- matched.test.incl.pData$MYCN.cat ### do not include this, as D/W Dan 25/10/18
# FacC <- matched.test.incl.pData$mstatus
# FacD<- matched.test.incl.pData$WNT_PNET5

### FacD <- matched.test.incl.pData$resection
### FacE <- matched.test.incl.pData$TP53.cat ### included in PNET5 algorithm update, updated 5/10/17
### FacF <- matched.test.incl.pData$meth
# data <- matched.test.incl.pData
# time <- matched.test.incl.pData$PFS
# event <- relapse.bin.incl
# marker <- matched.goi.vsd.cat.incl

### age is not included as we are including only age 3-16 years (and age < 3 years is listed in PNET5)
### included factors :
### LCA, mstatus, MYC or MYCN amplification (MYCMYCN.cat as d/w DW and confirmed 25/10/18), R+,  methylation subgroup (n=4, rather than WNT vs non WNT). TP53 included (will be broadly stratified in SHH)
### not including isochromosome 17q and how can this data be derived from matched.test.incl.pData$q17   (levels Gain, Neutral, Loss)


# cox.multivar.surv.PNET5_6 <- function (time, event, marker, FacA, FacB, FacC, FacD, FacE, FacF, strata = NULL, data) {
# if(is.null(strata)){
#  cox.temp <- coxph(Surv(time, event)~marker + FacA +FacB +FacC +FacD +FacE +FacF, data=data)
# }else {
#   cox.temp <- coxph(Surv(time, event)~marker + FacA +FacB +FacC +FacD +FacE + FacF, data=data)
#  }
# cox.p.val <- summary(cox.temp)$coefficients[1,5] ### updated 21/11/17. 
#  cox.HR <- summary(cox.temp)$coefficients[1,2] ### updated 21/11/17. 
#  cox.lower.95CI <- summary(cox.temp)$conf.int[1,3] ### multivariate, access 1st row results
#  cox.upper.95CI <- summary(cox.temp)$conf.int[1,4]
#  cox.Zscore <- summary(cox.temp)$coefficients[1,4] ### added this in to access Z score
#  cox.n <-summary(cox.temp)$n
#  cox.nevent <-summary(cox.temp)$nevent
#  summary.cox <- list(cox.pval = cox.p.val,cox.HR = cox.HR, cox.lower.95CI = cox.lower.95CI, cox.upper.95CI =cox.upper.95CI, cox.Zscore = cox.Zscore, n = cox.n, n.event = cox.nevent)
#  return (summary.cox)
# }


cox.multivar.surv.PNET5_4 <- function (time, event, marker, FacA, FacB, FacC, FacD, strata = NULL, data) {
  if(is.null(strata)){
    cox.temp <- coxph(Surv(time, event)~marker + FacA +FacB +FacC +FacD, data=data)
  }else {
    cox.temp <- coxph(Surv(time, event)~marker + FacA +FacB +FacC +FacD, data=data)
  }
  cox.p.val <- summary(cox.temp)$coefficients[1,5] ### updated 21/11/17. 
  cox.HR <- summary(cox.temp)$coefficients[1,2] ### updated 21/11/17. 
  cox.lower.95CI <- summary(cox.temp)$conf.int[1,3] ### multivariate, access 1st row results
  cox.upper.95CI <- summary(cox.temp)$conf.int[1,4]
  cox.Zscore <- summary(cox.temp)$coefficients[1,4] ### added this in to access Z score
  cox.n <-summary(cox.temp)$n
  cox.nevent <-summary(cox.temp)$nevent
  summary.cox <- list(cox.pval = cox.p.val,cox.HR = cox.HR, cox.lower.95CI = cox.lower.95CI, cox.upper.95CI =cox.upper.95CI, cox.Zscore = cox.Zscore, n = cox.n, n.event = cox.nevent)
  return (summary.cox)
}

################################################################################################

### Function number 6d for G3G4 and SHH subgroup multivariate analsis
### Factors Fac A - FacC
### input factors for G3G4 Lancet subanalysis
# MYC, G3G4_HighRisk, q13loss, no adjustment for gender
### Function name changed 31/10/18 to reflect survival analysis with 3 input factors

### input factors for SHH.old (SHH_Child) subanalysis
# MYCN, TP53, mstatus (although the MYCN amplification p value was 0.084 on the multivariate, keep in as it is accepted RF in SHH-mb) 


cox.multivar.surv_3 <- function (time, event, marker, FacA, FacB, FacC, strata = NULL, data) {    ### updated name 31/10/18
  if(is.null(strata)){
    cox.temp <- coxph(Surv(time, event)~marker + FacA + FacB + FacC, data=data)
  }else {
    cox.temp <- coxph(Surv(time, event)~marker + FacA + FacB +FacC, data=data)
  }
  # cox.p.val <- summary(cox.temp)$logtest[3] ### p value can also be called within cox.temp$coefficients
  cox.p.val <- summary(cox.temp)$coefficients[1,5] ### updated 21/11/17
  # cox.HR <- summary(cox.temp)$conf.int[1] ### called within cox.temp$coefficients
  cox.HR <- summary(cox.temp)$coefficients[1,2] ### updated 21/11/17
  cox.lower.95CI <- summary(cox.temp)$conf.int[1,3] ### as now multivariate, therefore need to access 1st row results
  cox.upper.95CI <- summary(cox.temp)$conf.int[1,4]
  cox.Zscore <- summary(cox.temp)$coefficients[1,4] ### added this in to access Z score
  cox.n <-summary(cox.temp)$n
  cox.nevent <-summary(cox.temp)$nevent
  summary.cox <- list(cox.pval = cox.p.val,cox.HR = cox.HR, cox.lower.95CI = cox.lower.95CI, cox.upper.95CI =cox.upper.95CI, cox.Zscore = cox.Zscore, n = cox.n, n.event = cox.nevent)
  return (summary.cox)
}

################################################################################################
#################################################################################################

### Function Number 7

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

### example input
# time <- matched.test.incl.pData$EFS
# event <- EFS.cat.bin.incl
# marker <- matched.goi.vsd.cat.incl


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
  1 - pchisq(EFS.incl.logrank$chisq, length(EFS.incl.logrank$obs)-1) -> surv.p.val.EFS
  text(4,0.1,paste("p =",round(surv.p.val.EFS, 3)), pos = 4, cex = 1)
  EFS.surv.table <- summary(km.EFS.incl)
  return(list(EFS.p.val = surv.p.val.EFS, 
              EFS.surv.table))
  if(out.file!="none"){
    dev.off()
  }
}



###########################################################################

### Function number 8
### logisticRegression
### replaced the hard coded script within previous clinicalPathAssess function 14/9/17, to check if function working

### input
## x <- dependent variable (must be binary, 0,1)
## y <- biomarker (independent variable)
## data <- matched dataset

### output is a dataframe containing: 
## OR
## 95% confidence intervals
## p value for logistic regression


#x <- matched.test.pData$age.cat.infant
#y <- matched.goi.vsd
#data <- matched.test.pData

logisticRegression <- function(x,y, data) {
  temp.glm <- glm (x ~ y, family = binomial (link = 'logit'), data= data)
  summary.temp <- summary(temp.glm)$coefficients  ### contains p value, z value and standard error
  LR.OR = exp(coef(temp.glm))
  LR.OR.95CI <- exp(confint(temp.glm))
  interim.LR <- cbind (LR.OR, LR.OR.95CI, summary.temp)
  return(list (LR.pval= summary.temp[2,4],
               interim.LR = interim.LR,  
               LR.OR.val=interim.LR[2,1],
               lower.95CI=LR.OR.95CI[2,1], 
               upper.95CI=LR.OR.95CI [2,2]
               )
  )
}

### note previously was outputting the interim.LR (ie cbind dataframe, and the relevant outputs were:   
### pattern: OR =log.reg.age.cat[2,1], LCI 95 = log.reg.age.cat [2,2], UCI 95 = log.reg.age.cat[2,3], p value = log.reg.age.cat [2,7])
### 140917: have now updated logisticRegression function to output the required values, listed below: 
### including dataframe with raw values
### OR, pvalue and 95CI:  LR.OR.pval, LR.pval, lower.95CI, upper.95CI respectively


#logisticRegression <- function(x,y, data) {
# temp.glm <- glm (x ~ y, family = binomial (link = 'logit'), data= data)
# stats.temp <- cbind(OR=exp(coef(temp.glm)), exp(confint(temp.glm)))
# summary.temp <- list (stats.temp, summary(temp.glm)$coefficients, summary(temp.glm)$data$Event)
# return(summary.temp)
#}

############################################################################


### Function number 9
### updatepData
### called within the pData_input_2017_09_11.R file, only needed when clinical database file is updated. 

updatepData <-  function(pData, meth7, cytogen, pdf.file = NULL, log.file = NULL){
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
  RTX <- ifelse(RTX=="Yes", "RTX", "No RTX")  ## changed 12/9/17
  # RTX <- pData$RTX_R 
  
  ### CSI
  CSI <-pData$RTXCSI_R 
  CSI <- ifelse(CSI=="Yes", "CSI", "No CSI") ## changed 14/9/17
  #CSI <-pData$RTXCSI_R
  
  ### molecular subgroups
  meth <-pData$X450K_R
  
 
  meth7.cat <- meth7[match(rownames(pData), rownames(meth7)),c(1)]
  
  
##########################
 cytogen[match(rownames(pData), rownames(cytogen)), ] -> cytogen.ordered
  
  
  
  #cytogen.ordered <- cytogen[match(rownames(pData), rownames(cytogen)),]
  cytogen.q13.cat <- ifelse(cytogen.ordered$q13 == "Loss", "q13 Loss", "no q13 loss") ## changed 12/9/17
  
  
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
  # PFS <-pData$PFS_yr_R
  PFS <-as.numeric(as.character(pData$PFS_yr_R))
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
                           EFS,
                           meth7.cat,
                           q13loss= cytogen.q13.cat,
                           cytogen.ordered
  )
  
  dev.off()
  sink()
  return(test.pData)
}


###############################################################################

### clinPathAssess function is designed to input a goi.vsd and run that against clinical correlation tests
### output is a list of list of results
### important elements within the list of lists are then called out within the clinical_data_110917.R script

### goi.vsd #### name of a group of putative biomarkers or the converted numerical expression data within the RNA seq input file

clinPathAssess <- function(test.pData,
                           goi.vsd,
                           pdf.file = NULL,
                           log.file = NULL
)
  {
  
  ### attempt with Dan
     # test.pData = test.pData
     # x -> goi.vsd   
     # pdf.file = NULL
     # log.file = NULL
  
  ### will need to run clinical_data_master generate gp.filt.mb.vsd which is required for goi.vsd below. Can generate in clinical_data_master.R first, or below
  ### input required
      ### 1. a goi line
      ### 2. goi.vsd line 717 in this file that links to gp.filt.mb.vsd
      ### 3. names(goi.vsd) line 723 
  ### if then interrogate functions, will need to generate matched.test.pData
  ### unhash ** when running individual goi
  
  # goi <- "ENSG00000124588" ### NQO2 ** or equivalent goi of interest
  # goi <- "ENSG00000136997" ### MYC

  
  ### unhash here for goi.vsd
  # goi.vsd <- as.numeric(gp.filt.mb.vsd[goi,]) ### ** 9/1/18 hashed when running the clinical_data_master.R; unhashed if running individual goi
                                          
  ###  12/6/18 use filtered file( gp.filt.mb.vsd gives same results as filt.mb.vsd)
  ### the output would be a vector with a continuous variable, names equal NMB numbers
  
  # names(goi.vsd) <- names(gp.filt.mb.vsd) ### ** added 16/1/18, unhash ** when running individual goi

  # test.pData = test.pData ###**
  # pdf.file = NULL  ###**
  # log.file = NULL  ###**

  
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
  
  ### matching the expression data (goi data) with the initially compiled clinical data with important variables
  ### will need to exclude duplicates (two RNA samples separately listed, defined as "NMBnumberT") 16/1/18
  
  ##### looking to combine dataset to plot expression data for alive vs dead patients
  # matchdata <- test.pData[names(goi.vsd),]
  # dim(matchdata)
  # test <- cbind(matchdata,goi.vsd)
  # plot(test$goi.vsd, col = ifelse(levels(test$OS.cat) =="Alive","red", "blue"), xlab = "matched patient sample", ylab = "Expression", main = "Expression within the cohort compared to survival")
  # legend (x= "topright", lty=1:2, paste(text = c("alive","dead")), col = c("red", "blue"))
  
  #####
  
  index <- match(names(goi.vsd), rownames(test.pData)) 
  matched.test.pData <- test.pData[index[!is.na(index)],] 
  is.vector(matched.test.pData)
  #as.data.frame(matched.test.pData) ### added 070917 in attempt to avoid downstream" $ atomic in a vector" error
  matched.goi.vsd <- goi.vsd[!is.na(index)] 
  matched.goi.vsd.cat <- ifelse(matched.goi.vsd>median(goi.vsd, na.rm = T), "high","low") 
  
  #############################################
  
  ### summary data 
  
  ### visualising means of the goi in the matched dataset
  
  age.mean.goi <- plotmeans(matched.test.pData$age.cont ~ matched.test.pData$meth, data=matched.test.pData)
  
  ### summary data for all variables within the age 3-16 yo group, to construct survival cohort
  
  age.df <- data.frame(test.pData$NMB, test.pData$age.filter, test.pData$age.cont)
  summary(age.df)
  
  ##############################################
  
  ### Chi squared analysis
  
  # matched.test.pData$RTX <- ifelse(matched.test.pData$RTX=="Yes", "RTX", "No RTX") ### have renamed these objects within the main test.pData file therefore do not need these next 3 lines of code (140917)
  # matched.test.pData$CSI <-ifelse(matched.test.pData$CSI=="Yes", "CSI", "No CSI")
  # matched.test.pData$age.cat.infant <- ifelse(matched.test.pData$age.cat.infant=="TRUE", "infant", "non infant") 
  
  ### create WNT vs non WNT variable
  matched.test.pData$WNT_PNET5<- ifelse(matched.test.pData$meth == "WNT", "WNT", "non-WNT") ### added 31/10/18
  matched.test.pData$WNT1_0 <- ifelse(matched.test.pData$WNT_PNET5 =="WNT", 1, 0)
  
  ### run chi-squared analysis
  
  list.age.cat.infant <- chi.sq(x = matched.test.pData$age.cat.infant, y = matched.goi.vsd.cat)
  
  list.sex <- chi.sq (x = matched.test.pData$sex, y= matched.goi.vsd.cat)
  
  list.mstatus <- chi.sq(x = matched.test.pData$mstatus, y=matched.goi.vsd.cat) 
  
  list.relapse <- chi.sq(x = matched.test.pData$relapse, y = matched.goi.vsd.cat) 
  
  list.resection <- chi.sq (x = matched.test.pData$resection, y = matched.goi.vsd.cat)
  
  list.meth.4 <- chi.sq (x= matched.test.pData$meth, y = matched.goi.vsd.cat)
  
  list.meth.7 <- chi.sq (x = matched.test.pData$meth7, y = matched.goi.vsd.cat)
  
  list.WNTnon <- chi.sq (x = matched.test.pData$WNT_PNET5, y = matched.goi.vsd.cat) ### added 31/10/18
  
  list.path <- chi.sq (x = matched.test.pData$histopath, y = matched.goi.vsd.cat)
  
  list.LCA <- chi.sq (x = matched.test.pData$LCA, y = matched.goi.vsd.cat)
  
  list.q13loss <- chi.sq (x = matched.test.pData$q13loss, y = matched.goi.vsd.cat)
  
  list.MYC <- chi.sq (x = matched.test.pData$MYC.cat, y = matched.goi.vsd.cat)
  
  list.MYCN <- chi.sq (x = matched.test.pData$MYCN.cat, y = matched.goi.vsd.cat)
  
  list.MYCMYCN <- chi.sq(x= matched.test.pData$MYCMYCN.cat, y= matched.goi.vsd.cat)
  
  list.TP53 <- chi.sq (x = matched.test.pData$TP53.cat, y = matched.goi.vsd.cat)
  
  list.TERT <- chi.sq (x = matched.test.pData$TERT.cat, y = matched.goi.vsd.cat)
  
  list.RTX <- chi.sq (x = matched.test.pData$RTX, y = matched.goi.vsd.cat)
  
  list.CSI <- chi.sq (x = matched.test.pData$CSI, y = matched.goi.vsd.cat)
  
  #chi.sq.results <- list(list.age.cat.infant,
  #                       list.sex,
   #                      list.mstatus,
   #                      list.relapse, 
   #                      list.resection, 
   #                      list.meth.4,
   #                      list.meth.7, 
   #                      list.path,
    #                     list.LCA, 
     #                    list.MYC, 
      #                   list.MYCN, 
       #                  list.MYCMYCN, 
        #                 list.TP53,
         #                list.TERT, 
          #               list.q13loss,
           #              list.RTX,
            #             list.CSI
  #)
  
  chi.sq.list <- as.list(mget(ls(pattern="list."))) 
  
  ### run Fisher's exact test on those where count < 5 e.g pathology "other"
  
  try(histopath.table <- table(as.factor(matched.test.pData$LCA), as.factor(matched.goi.vsd.cat)), silent = T)
  try(histopath.result <- fisher.test(histopath.table), silent = T)
 
  ###################################
  
  ### Correlation coefficients
  
  #cat ("correlation coefficients for association between variables", sep ="\n")
  
  try(list.cor.age <- cor.result(x = matched.test.pData$age.cont, y = matched.goi.vsd), silent = T)
  #list(list.cor.age)
 
  
  ##################################
  
  ###  linear regression
  
  try(lin.reg.age <- lin.reg(x= matched.test.pData$age.cont, y= matched.goi.vsd), silent = T)
  
  
  #################################
  
  ### Mann-whitney U (aka wilcoxon rank sum test) for non-parametric data
  
  ### age continuous 
  try(age.cont.wilcox <- wilcox.test(matched.test.pData$age.cont ~ matched.goi.vsd.cat, exact = F, correct = F), silent = T)
  
  
  ##################################
  ### logistic regression
  
  #cat ("processing logistic regression within defined logisticRegression function for each variable", sep ="\n")
  ### have now updated logisticRegression function to output the required values, listed below: 
  ### pvalue [[1]][[1]], raw dataframe with intercept and y values, OR and 95CI:  LR.pval, interim.LR, LR.OR.pval,  lower.95CI, upper.95CI respectively
  
  try(log.reg.age.cat <- logisticRegression(x = matched.test.pData$age.cat.infant, y= matched.goi.vsd, data = matched.test.pData), silent= T)
  try(log.reg.sex <- logisticRegression (x = matched.test.pData$sex, y = matched.goi.vsd, data= matched.test.pData), silent= T)
  try(log.reg.mstatus <- logisticRegression(x= matched.test.pData$mstatus, y= matched.goi.vsd, data=matched.test.pData), silent= T)
  try(log.reg.relapse <- logisticRegression(x= matched.test.pData$relapse, y = matched.goi.vsd, data=matched.test.pData), silent= T)
  try(log.reg.resection <- logisticRegression(x =matched.test.pData$resection, y = matched.goi.vsd, data=matched.test.pData), silent= T)
  try(log.reg.histopath <- logisticRegression (x= matched.test.pData$histopath, y= matched.goi.vsd, data=matched.test.pData), silent= T)
  try(log.reg.LCA <- logisticRegression(x=matched.test.pData$LCA, y= matched.goi.vsd, data=matched.test.pData), silent= T)
  try(log.reg.MYC <- logisticRegression (x= matched.test.pData$MYC.cat, y= matched.goi.vsd,  data=matched.test.pData), silent= T)
  try(log.reg.MYCN <- logisticRegression (x= matched.test.pData$MYCN.cat, y=matched.goi.vsd, data=matched.test.pData), silent= T)
  try(log.reg.MYCMYCN <- logisticRegression (x=matched.test.pData$MYCMYCN.cat, y= matched.goi.vsd, data=matched.test.pData), silent= T)
  try(log.reg.meth <- logisticRegression (x=matched.test.pData$meth, y= matched.goi.vsd, data = matched.test.pData), silent= T)
  try(log.reg.meth7 <- logisticRegression (x = matched.test.pData$meth7, y= matched.goi.vsd, data = matched.test.pData), silent= T)
  try(log.reg.WNT1_0 <-  logisticRegression (x = matched.test.pData$WNT1_0, y =  matched.goi.vsd, data = matched.test.pData), silent = T) ### added 31/10/18 ### WNT = 1, nonWNT =0
  try(log.reg.TP53 <- logisticRegression (x= matched.test.pData$TP53.cat, y= matched.goi.vsd, data=matched.test.pData), silent= T)
  try(log.reg.TERT <- logisticRegression (x= matched.test.pData$TERT.cat, y= matched.goi.vsd, data = matched.test.pData), silent= T)
 
   ######################################################################
  
  reg.log.list <- as.list(mget(ls(pattern="log.reg"))) 
  
  ######################################################################

  ### visualise distribution of biomarker in cohort
  #cat ("visualisation of biomarker and relationship with methylation", sep = "\n")
  qqnorm(matched.goi.vsd)
  
  ### posthoc for multiple categories e.g histopath, meth, meth7
  ### pairwise t-test to determine where the difference lies between the groups
  try(histopath.pw <- pairwise.t.test(matched.goi.vsd, matched.test.pData$histopath), silent = T)
  try(LCA.pw <- pairwise.t.test(matched.goi.vsd, matched.test.pData$LCA, level = matched.test.pData$LCA == "Non LCA"), silent = T)
  MYCMYCN.pw <- pairwise.t.test(matched.goi.vsd, matched.test.pData$MYCMYCN.cat)
  meth.pw <- pairwise.t.test(matched.goi.vsd, matched.test.pData$meth)
  meth7.pw <- pairwise.t.test(matched.goi.vsd, matched.test.pData$meth7.cat)
  
  ### if biomarker is normally distributed, can use ANOVA (one-way ANOVA)
  meth7.aov <- aov (matched.goi.vsd ~ matched.test.pData$meth7, data = matched.test.pData)
  #summary(meth7.aov)
  plot(meth7.aov)
  
  #######################################################
  
  ### hardcoded boxplots, separated now from logistic regression coding 140917
  
  age.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$age.cat.infant, col = c("red", "blue"), xlab = "Infant", ylab = "Biomarker expression", main = "CDK6 and age (infant vs non infant)")
  sex.boxplot <- boxplot (matched.goi.vsd ~ matched.test.pData$sex, col = c("red", "blue"), xlab = "Gender", ylab = "Expression of biomarker", main = "CDK6 expression and gender")
  mstatus.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$mstatus, col = c("red", "blue"), xlab = "M status", ylab = "Biomarker expression", main = "CDK6 expression and metastatic status")
  relapse.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$relapse, col = c("red", "blue"), xlab = "Relapse status", ylab = "Biomarker expression",  main = "Correlation between CDK6 and relapse")
  resection.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$resection, col = c("red", "blue"), xlab = "Resection status", ylab = "Biomarker expression", main = "Correlation between CDK6 and resection status")
  histopath.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$histopath, col = c("red", "blue"), xlab = "Histopathology subtype", ylab = "Biomarker expression", main = "Correlation between biomarker and histopathology")
  LCA.boxplot <- boxplot (matched.goi.vsd~ matched.test.pData$LCA,col=c("red","blue"),  xlab = "LCA pathology", ylab = "Biomarker expression", main = "Correlation between CDK6 and LCA")
  MYC.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$MYC.cat, col = c("Red", "Blue"), xlab = "MYC amplification", ylab = "Biomarker expression", main = "CDK6 and MYC amplification")
  MYCN.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$MYCN.cat, col = c("Red", "Blue"), xlab = "amplification", ylab = "Biomarker expression", main = "CDK6 and MYCN amplification")
  MYCMYCN.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$MYCMYCN.cat, col = c("Red", "Blue"), xlab = "MYC or MYCN amplification", ylab = "Biomarker expression", main = "CDK6 and MYC/MYCN amplification")
  ### determine legend 
  
  meth.boxplot <- boxplot(matched.goi.vsd~matched.test.pData$meth, col=c("yellow","green","red","blue"), xlab = "Methylation subgroup", ylab = "Biomarker expression", main = "CDK6 and 4 molecular subgroups")
  meth7.boxplot <- boxplot(matched.goi.vsd~matched.test.pData$meth7, col=c("pink", "yellow","dark green","light green", "red","dark red", "blue"),  xlab = "Methylation subgroup", ylab = "Biomarker expression", main = "CDK6 and 7 molecular subgroups", names = c("G3_HR","G3_LR", "G4_HR", "G4_LR", "SHH-Inf", "SHH-Child", "WNT"))
  TP53.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$TP53.cat, col = c("Red", "Blue"), xlab = "TP53 mutational status", ylab = "Biomarker expression", main = "CDK6 and TP53 status")
  TERT.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$TERT.cat, col = c("Red", "Blue"), xlab = "TERT amplification (pos vs neg)", ylab = "Biomarker expression", main = "CDK6 and TERT amplifiation")
  WNT.boxplot <- boxplot (matched.goi.vsd ~ matched.test.pData$WNT1_0, col = c("Red", "Blue"), xlab = "WNT (1) vs non WNT(0)", ylab = "Biomarker expression", main = "CDK6 and WNT status")
  
  #######################################################################
  
  ### summary data for each methylation subgroup (n=4), no age restriction for first analyses
  
  #cat ("processing summary stats for 4 molecular subgroups, no age restriction", sep ="\n")
  
  ### G3
  
  G3 <- matched.test.pData$meth =="G3" 
  
  G3.group <- matched.test.pData [G3, ]
  # summary (G3.group)
  
  ## G4
  
  G4 <- matched.test.pData$meth =="G4"
  G4.group <- matched.test.pData [G4, ]
  #summary(G4.group)
  
  ### SHH
  
  SHH <- matched.test.pData$meth =="SHH"
  SHH.group <- matched.test.pData [SHH, ]
  #summary(SHH.group)
  
  ### WNT
  
  WNT <- matched.test.pData$meth == "WNT"
  WNT.group <- matched.test.pData [WNT, ]
  #summary(WNT.group)
  
  #cat ("processing summary stats for 7 molecular subgroups, no age restriction", sep = "\n")
  
  ### G3 subgroups
  
  G3.high.group <- matched.test.pData$meth7 =="Grp3_HighRisk"
  G3.high.data <- matched.test.pData [G3.high.group, ]
  #summary(G3.high.data)
  
  G3.low.group <- matched.test.pData$meth7 =="Grp3_LowRisk"
  G3.low.data <- matched.test.pData [G3.low.group, ]
  #summary(G3.low.data)
  
  ### G4 subgroups
  
  G4.high.group <- matched.test.pData$meth7 =="Grp4_HighRisk"
  G4.high.data <- matched.test.pData [G4.high.group, ]
  #summary(G4.high.data)
  
  G4.low.group <- matched.test.pData$meth7 =="Grp4_LowRisk"
  G4.low.data <- matched.test.pData [G4.low.group, ]
  #summary(G4.low.data)
  
  ### SHH
  SHH.inf.group <- matched.test.pData$meth7 == "SHH_Inf"
  SHH.inf.data <- matched.test.pData [SHH.inf.group, ]
  #summary (SHH.inf.data)
  
  SHH.old.group <- matched.test.pData$meth7 =="SHH_Old"
  SHH.old.data <- matched.test.pData [SHH.old.group, ]
  #summary (SHH.old.data)
  
  #############################################
  
  ### New dataframes for survival analyses
  
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
  
  ### have created WNT vs non WNT group earlier in matched.test.pData ### added 31/10/18
  
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
  ### previous script lines 1073-1076 
  
  # index.incl <- match(names(goi.vsd), rownames(Age.incl.df)) 
  # matched.test.incl.pData <- Age.incl.df[index.incl[!is.na(index.incl)],] 
  # matched.goi.vsd.incl <- goi.vsd[!is.na(index.incl)] 
  # matched.goi.vsd.cat.incl <- ifelse(matched.goi.vsd.incl>median(goi.vsd, na.rm = T), "high","low")
  
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
  
  ### creating binary relapse variables labelled 0,1 for event analysis. ### call dimensions of the variable by length(var)
  
  relapse.bin <- ifelse(matched.test.pData$relapse == "relapse", 1, 0)
  relapse.bin.incl <- ifelse(matched.test.incl.pData$relapse == "relapse", 1,0)
  relapse.G3G4.bin.incl <- ifelse(matched.G3G4.incl.pData$relapse == "relapse", 1,0)
  relapse.G4.bin.incl <- ifelse(matched.G4.incl.pData$relapse=="relapse", 1, 0)
  relapse.SHH.bin.incl <- ifelse(matched.SHH.incl.pData$relapse == "relapse", 1,0)
  relapse.SHH.old.bin.incl <- ifelse(matched.SHH.old.incl.pData$relapse == "relapse", 1,0)
  
  
  OS.cat.bin <- ifelse(matched.test.pData$OS.cat == "Dead", 1, 0) 
  OS.cat.bin.incl <- ifelse(matched.test.incl.pData$OS.cat == "Dead", 1,0)
  OS.G3G4.bin.incl <- ifelse(matched.G3G4.incl.pData$OS.cat == "Dead", 1,0)
  OS.G4.bin.incl <- ifelse (matched.G4.incl.pData$OS.cat == "Dead", 1, 0)
  OS.SHH.bin.incl <- ifelse(matched.SHH.incl.pData$OS.cat == "Dead", 1,0)
  OS.SHH.old.bin.incl <- ifelse(matched.SHH.old.incl.pData$OS.cat == "Dead", 1,0)
  
  EFS.cat.bin <- ifelse(matched.test.pData$Event == "Event", 1, 0)
  EFS.cat.bin.incl <- ifelse(matched.test.incl.pData$Event == "Event", 1,0)
  EFS.G3G4.bin.incl <- ifelse(matched.G3G4.incl.pData$Event == "Event", 1,0)
  EFS.G4.bin.incl <- ifelse(matched.G4.incl.pData$Event == "Event", 1,0)
  EFS.SHH.bin.incl <- ifelse(matched.SHH.incl.pData$Event == "Event", 1,0)
  EFS.SHH.old.bin.incl <- ifelse(matched.SHH.old.incl.pData$Event == "Event", 1,0)
  
  #### PFS ### this is where the graphs are generated again
  
  #cat("running km survival curve for PFS and biomarker, graphical output to PDF", sep = "\n")
  ### incldes overall group, G3G4 combined, SHH overall and SHH_Child (SHH.old)
  
  # km.log.test.all <- km.log.test(time = as.numeric(as.character(matched.test.incl.pData$PFS)), event = relapse.bin.incl, marker = matched.goi.vsd.cat.incl)
  try(surv.km.PFS.all <- km.log.test(time = matched.test.incl.pData$PFS, event = relapse.bin.incl, marker = matched.goi.vsd.cat.incl), silent = T)
  
  # km.log.test.G3G4 <- km.log.test(time = as.numeric(as.character(matched.G3G4.incl.pData$PFS)), event = relapse.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl)
  
  try(surv.km.PFS.G3G4 <- km.log.test(time = matched.G3G4.incl.pData$PFS, event = relapse.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl), silent = T)
  
  try(surv.km.PFS.G4 <- km.log.test(time = matched.G4.incl.pData$PFS, event = relapse.G4.bin.incl, marker = matched.goi.vsd.cat.G4.incl), silent = T)
  
  try(surv.km.PFS.SHH <- km.log.test (time = matched.SHH.incl.pData$PFS, event = relapse.SHH.bin.incl, marker = matched.goi.vsd.cat.SHH.incl), silent = T)
  
  try(surv.km.PFS.SHH.old <- km.log.test (time = matched.SHH.old.incl.pData$PFS, event = relapse.SHH.old.bin.incl, marker = matched.goi.vsd.cat.SHH.old.incl), silent = T)
  
  ################################################################################
  ### cox regression analysis
  
  #cat("cox regression analysis, PFS, biomarker, age 3-16 years treated with curative therapy", sep ="\n") including overall, G3G4 combined, SHH and SHH.old
  ### added in continuous variables for expression data 4/10/17
  
  try(surv.cox.relapse.incl.cat <- cox.result.surv (time=matched.test.incl.pData$PFS, event =  relapse.bin.incl, marker = matched.goi.vsd.cat.incl, data = matched.test.incl.pData), silent = T)
  try(surv.cox.relapse.incl.contin <- cox.result.surv (time=matched.test.incl.pData$PFS, event =  relapse.bin.incl, marker = matched.goi.vsd.incl, data = matched.test.incl.pData), silent = T)
  
  try(surv.cox.relapse.incl.G3G4.cat <- cox.result.surv (time=matched.G3G4.incl.pData$PFS, event =  relapse.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl, data = matched.test.incl.pData), silent = T)
  try(surv.cox.relapse.incl.G3G4.contin <- cox.result.surv (time=matched.G3G4.incl.pData$PFS, event =  relapse.G3G4.bin.incl, marker = matched.goi.vsd.G3G4.incl, data = matched.test.incl.pData), silent = T)
  
  try(surv.cox.relapse.incl.G4.cat <- cox.result.surv(time = matched.G4.incl.pData$PFS, event = relapse.G4.bin.incl, marker = matched.goi.vsd.cat.G4.incl, data = matched.test.incl.pData), silent = T)
  try(surv.cox.relapse.incl.G4.contin <- cox.result.surv(time = matched.G4.incl.pData$PFS, event = relapse.G4.bin.incl, marker = matched.goi.vsd.G4.incl, data = matched.test.incl.pData), silent = T)
  
  try (surv.cox.relapse.incl.SHH.cat <- cox.result.surv (time = matched.SHH.incl.pData$PFS, event = relapse.SHH.bin.incl, marker = matched.goi.vsd.cat.SHH.incl, data = matched.test.incl.pData), silent = T)
  try (surv.cox.relapse.incl.SHH.contin <- cox.result.surv (time = matched.SHH.incl.pData$PFS, event = relapse.SHH.bin.incl, marker = matched.goi.vsd.SHH.incl, data = matched.test.incl.pData), silent = T)
  
  try (surv.cox.relapse.incl.SHH.old.cat <- cox.result.surv (time = matched.SHH.old.incl.pData$PFS, event = relapse.SHH.old.bin.incl, marker = matched.goi.vsd.cat.SHH.old.incl, data = matched.test.incl.pData), silent = T)
  try (surv.cox.relapse.incl.SHH.old.contin <- cox.result.surv (time = matched.SHH.old.incl.pData$PFS, event = relapse.SHH.old.bin.incl, marker = matched.goi.vsd.SHH.old.incl, data = matched.test.incl.pData), silent = T)
   
  #####################################################################
  ### cox multivariate with established risk factors per PNET5 and Schwalbe combined (which adds MYC, TP53 and q13loss and 7 molecular groups)
  
  # cox.multivar.PFS.incl <- coxph (Surv(matched.test.incl.pData$PFS, relapse.bin.incl)~ matched.goi.vsd.cat.incl + matched.test.incl.pData$LCA + matched.test.incl.pData$MYC.cat + matched.test.incl.pData$MYCN.cat + matched.test.incl.pData$mstatus + matched.test.incl.pData$resection + matched.test.incl.pData$q13loss + matched.test.incl.pData$TP53.cat + matched.test.incl.pData$sex, data = matched.test.incl.pData)
  # summary(cox.multivar.PFS.incl)
  
  ############
  #  cox.p.val.marker <- summary(cox.multivar.PFS.incl)$coefficients[1,5] 
  #  cox.p.val.model <- summary(cox.multivar.PFS.incl)$logtest[3] ### this is the p likelihood ratio
  # cox.HR <- summary(cox.multivar.PFS.incl)$conf.int[1] ### this is the hazard ratio for the matched.goi.vsd.cat.incl
  # cox.lower.95CI <- summary(cox.multivar.PFS.incl)$conf.int[1,3] ### as now multivariate, therefore need to access 1st row details for 
  # cox.upper.95CI <- summary(cox.multivar.PFS.incl)$conf.int[1,4]
  # cox.Zscore <- summary(cox.multivar.PFS.incl)$coefficients[1,4] ### added this in to access Z score
  # cox.n <-summary(cox.multivar.PFS.incl)$n
  # cox.nevent <-summary(cox.multivar.PFS.incl)$nevent
  # summary.cox <- list(cox.pval.marker = cox.p.val.marker, cox.p.val.modle = cox.p.val.model,cox.HR = cox.HR, cox.lower.95CI = cox.lower.95CI, cox.upper.95CI =cox.upper.95CI, cox.Zscore = cox.Zscore, n = cox.n, n.event = cox.nevent)
  # summary.cox
  
  ############################
  ### re-analysis 31/10/18 using factors that are significant in the RNA cohort in univariate regression for combined PNET5 and Lancet: mstatus, LCA, MYCMYCN, meth7, q13loss
  
  # try(multivar.cox.PFS.combined.cat <- cox.multivar.surv_7(time = matched.test.incl.pData$PFS, event = relapse.bin.incl, marker = matched.goi.vsd.cat.incl, FacA= matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat, FacC= matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$resection, FacE= matched.test.incl.pData$q13loss, FacF= matched.test.incl.pData$TP53.cat, FacG = matched.test.incl.pData$meth7.cat,  data = matched.test.incl.pData), silent =T)
  # try(multivar.cox.PFS.combined.contin <- cox.multivar.surv_7(time = matched.test.incl.pData$PFS, event = relapse.bin.incl, marker = matched.goi.vsd.incl, FacA= matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat, FacC= matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$resection, FacE= matched.test.incl.pData$q13loss, FacF= matched.test.incl.pData$TP53.cat,  FacG = matched.test.incl.pData$meth7.cat,  data = matched.test.incl.pData), silent =T)
  
  try(multivar.cox.PFS.combined.cat <- cox.multivar.surv_5(time = matched.test.incl.pData$PFS, event = relapse.bin.incl, marker = matched.goi.vsd.cat.incl, FacA= matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat, FacC= matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$meth7.cat, FacE= matched.test.incl.pData$q13loss,  data = matched.test.incl.pData), silent =T)
  try(multivar.cox.PFS.combined.contin <- cox.multivar.surv_5(time = matched.test.incl.pData$PFS, event = relapse.bin.incl, marker = matched.goi.vsd.incl, FacA= matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat, FacC= matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$meth7.cat, FacE= matched.test.incl.pData$q13loss,  data = matched.test.incl.pData), silent =T)
  
  ### reanalysis 31/10/18 after determined that in overall group, mstatus, LCA, MYCMYCN (vs non MYC/non MYC amplified), WNT vs non WNT in PNET5 schema is significant in univariate log regression
  
  try(multivar.cox.PFS.PNET5.cat <- cox.multivar.surv.PNET5_4 (time = matched.test.incl.pData$PFS, event = relapse.bin.incl, marker = matched.goi.vsd.cat.incl, FacA = matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat, FacC = matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$WNT_PNET5, data = matched.test.incl.pData), silent =T) ### updated here
  try(multivar.cox.PFS.PNET5.contin <- cox.multivar.surv.PNET5_4 (time = matched.test.incl.pData$PFS, event = relapse.bin.incl, marker = matched.goi.vsd.incl, FacA = matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat, FacC = matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$WNT_PNET5, data = matched.test.incl.pData), silent =T)
  
  
  ### cox multivariate with established risk factors per PNET5 only ### analysis prior to 31/10/18
  ### LCA, mstatus,  MYCMYCN.cat, resection, meth (after D/W Dan, reinterpretation of the paper), TP53.cat (added 5/10/17)
  
  # try(multivar.cox.PFS.PNET5.cat <- cox.multivar.surv.PNET5_6 (time = matched.test.incl.pData$PFS, event = relapse.bin.incl, marker = matched.goi.vsd.cat.incl, FacA = matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat, FacC = matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$resection,  FacE = matched.test.incl.pData$meth, FacF = matched.test.incl.pData$TP53.cat, data = matched.test.incl.pData), silent =T)
  # try(multivar.cox.PFS.PNET5.contin <- cox.multivar.surv.PNET5_6 (time = matched.test.incl.pData$PFS, event = relapse.bin.incl, marker = matched.goi.vsd.incl, FacA = matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat, FacC = matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$resection, FacE = matched.test.incl.pData$meth, FacF = matched.test.incl.pData$TP53.cat, data = matched.test.incl.pData), silent =T)
  
  ### PNET5 RF in G3G4 subgroup (updated 31/10/18 to include factors significant in univariate regression : mstatus, LCA, MYC)
  
  try(multivar.cox.PFS.PNET5.G3G4.cat <- cox.multivar.surv_3 (time =matched.G3G4.incl.pData$PFS,event = relapse.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl, FacA = matched.G3G4.incl.pData$LCA, FacB = matched.G3G4.incl.pData$MYC.cat, FacC = matched.G3G4.incl.pData$mstatus, data =  matched.G3G4.incl.pData), silent =T )
  try(multivar.cox.PFS.PNET5.G3G4.contin <- cox.multivar.surv_3 (time =matched.G3G4.incl.pData$PFS,event = relapse.G3G4.bin.incl, marker = matched.goi.vsd.G3G4.incl, FacA = matched.G3G4.incl.pData$LCA, FacB = matched.G3G4.incl.pData$MYC.cat, FacC = matched.G3G4.incl.pData$mstatus,  data =  matched.G3G4.incl.pData), silent =T )
  
  ### cox multivariate with established risk factors per Schwalbe only
  ### analysis 31/10/18 decision to include 3 significant factors from PNET5/Lancet paper for SHH_old (SHH_child) cohort: LCA, MYCN and TP53. Resection and metastatic status were not signicant in univariate logistic regression analysis
  
  try(multivar.cox.PFS.SHH.old.cat <- cox.multivar.surv_3(time = matched.SHH.old.incl.pData$PFS, event = relapse.SHH.old.bin.incl, marker = matched.goi.vsd.cat.SHH.old.incl, FacA = matched.SHH.old.incl.pData$MYCN.cat, FacB = matched.SHH.old.incl.pData$LCA, FacC = matched.SHH.old.incl.pData$TP53.cat,  data = matched.SHH.old.incl.pData ), silent =T)
  try(multivar.cox.PFS.SHH.old.contin <- cox.multivar.surv_3(time = matched.SHH.old.incl.pData$PFS, event = relapse.SHH.old.bin.incl, marker = matched.goi.vsd.SHH.old.incl, FacA = matched.SHH.old.incl.pData$MYCN.cat, FacB = matched.SHH.old.incl.pData$LCA, FacC = matched.SHH.old.incl.pData$TP53.cat, data = matched.SHH.old.incl.pData ), silent =T)
  
  ### subgroup specific for SHH_Child (which is SHH_Old in our dataframe), analysis prior to 31/10/18
  # try(multivar.cox.PFS.SHH.old.cat <- cox.multivar.surv.Schwalbe_3(time = matched.SHH.old.incl.pData$PFS, event = relapse.SHH.old.bin.incl, marker = matched.goi.vsd.cat.SHH.old.incl, FacA = matched.SHH.old.incl.pData$MYCN.cat, FacB = matched.SHH.old.incl.pData$mstatus, FacC = matched.SHH.old.incl.pData$TP53.cat,  data = matched.SHH.old.incl.pData ), silent =T)
  # try(multivar.cox.PFS.SHH.old.contin <- cox.multivar.surv.Schwalbe_3(time = matched.SHH.old.incl.pData$PFS, event = relapse.SHH.old.bin.incl, marker = matched.goi.vsd.SHH.old.incl, FacA = matched.SHH.old.incl.pData$MYCN.cat, FacB = matched.SHH.old.incl.pData$mstatus, FacC = matched.SHH.old.incl.pData$TP53.cat, data = matched.SHH.old.incl.pData ), silent =T)
  
  
  ### Group 3 and 4 combined: Group 3/4 HR, MYC amplification, q13 loss - analysis unchanged 31/10/18
  ### high risk methylation group (ie HR Group 3 and HR Group 4 combined) in matched.G3G4.incl.pData$G3G4.high.incl for categorical and continuous variable
  
  try(multivar.cox.PFS.Lancet.G3G4.cat <- cox.multivar.surv_3 (time = matched.G3G4.incl.pData$PFS, event = relapse.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl, FacA = matched.G3G4.incl.pData$MYC.cat, FacB = matched.G3G4.incl.pData$G3G4.high.incl, FacC = matched.G3G4.incl.pData$q13loss,  data = matched.G3G4.incl.pData), silent =T)
  try(multivar.cox.PFS.Lancet.G3G4.contin <- cox.multivar.surv_3 (time = matched.G3G4.incl.pData$PFS, event = relapse.G3G4.bin.incl, marker = matched.goi.vsd.G3G4.incl, FacA = matched.G3G4.incl.pData$MYC.cat, FacB = matched.G3G4.incl.pData$G3G4.high.incl, FacC = matched.G3G4.incl.pData$q13loss,  data = matched.G3G4.incl.pData), silent =T)
  
  
  #####################################################################
  
  ### OS 
  
  #cat ("km analysis for OS for children aged 3-16 years, treated with curative intent", sep = "\n")
  
  try(surv.km.OS.all <- km.log.test.OS(time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.cat.incl ), silent = T)
  try(surv.km.OS.G3G4 <- km.log.test.OS(time = matched.G3G4.incl.pData$Followup, event = OS.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl), silent = T)
  try(surv.km.OS.G4 <- km.log.test.OS (time =matched.G4.incl.pData$Followup, event = OS.G4.bin.incl, marker = matched.goi.vsd.cat.G4.incl), silent = T)
  try(surv.km.OS.SHH <- km.log.test.OS(time = matched.SHH.incl.pData$Followup, event = OS.SHH.bin.incl, marker = matched.goi.vsd.cat.SHH.incl), silent = T)
  try(surv.km.OS.SHH.old <- km.log.test.OS(time = matched.SHH.old.incl.pData$Followup, event = OS.SHH.old.bin.incl, marker = matched.goi.vsd.cat.SHH.old.incl), silent = T)

  ### cox regression analysis
  
  #cat("cox regression on age 3-16 years, curative", sep = "\n")
  
  try(surv.cox.result.OS.all.cat <- cox.result.OS (time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.cat.incl, data =  matched.test.incl.pData), silent = T)
  try(surv.cox.result.OS.all.contin <- cox.result.OS (time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.incl, data =  matched.test.incl.pData), silent = T)
  
  try(surv.cox.result.OS.G3G4.cat <- cox.result.OS (time = matched.G3G4.incl.pData$Followup, event = OS.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl, data =  matched.test.incl.pData), silent = T)
  try(surv.cox.result.OS.G3G4.contin <- cox.result.OS (time = matched.G3G4.incl.pData$Followup, event = OS.G3G4.bin.incl, marker = matched.goi.vsd.G3G4.incl, data =  matched.test.incl.pData), silent = T)
  
  try(surv.cox.results.OS.G4.cat <- cox.result.OS (time = matched.G4.incl.pData$Followup,event = OS.G4.bin.incl, marker = matched.goi.vsd.cat.G4.incl, data = matched.G4.incl.pData), silent = T)
  try(surv.cox.results.OS.G4.contin <- cox.result.OS (time = matched.G4.incl.pData$Followup,event = OS.G4.bin.incl, marker = matched.goi.vsd.G4.incl, data = matched.G4.incl.pData), silent = T)
  
  try (surv.cox.result.OS.SHH.cat <- cox.result.surv (time = matched.SHH.incl.pData$Followup, event = OS.SHH.bin.incl, marker = matched.goi.vsd.cat.SHH.incl, data = matched.test.incl.pData), silent = T)
  try (surv.cox.result.OS.SHH.contin <- cox.result.surv (time = matched.SHH.incl.pData$Followup, event = OS.SHH.bin.incl, marker = matched.goi.vsd.SHH.incl, data = matched.test.incl.pData), silent = T)
  
  try (surv.cox.result.OS.SHH.old.cat <- cox.result.surv (time = matched.SHH.old.incl.pData$Followup, event = OS.SHH.old.bin.incl, marker = matched.goi.vsd.cat.SHH.old.incl, data = matched.test.incl.pData), silent = T)
  try (surv.cox.result.OS.SHH.old.contin <- cox.result.surv (time = matched.SHH.old.incl.pData$Followup, event = OS.SHH.old.bin.incl, marker = matched.goi.vsd.SHH.old.incl, data = matched.test.incl.pData), silent = T)
 
  
  #####################################
  
  ### multivariate cox regression for OS, based on categories combined (PNET + Schwalbe all risk factors, PNET5 alone, Schwalbe (lancet) subgroups)
  
  # try(multivar.cox.OS.combined.cat <- cox.multivar.surv_7(time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.cat.incl, FacA= matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat, FacC= matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$resection, FacE= matched.test.incl.pData$q13loss, FacF= matched.test.incl.pData$TP53.cat, FacG = matched.test.incl.pData$meth7.cat,  data = matched.test.incl.pData), silent =T)
  # try(multivar.cox.OS.combined.contin <- cox.multivar.surv_7(time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.incl, FacA= matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat, FacC= matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$resection, FacE= matched.test.incl.pData$q13loss, FacF= matched.test.incl.pData$TP53.cat, FacG= matched.test.incl.pData$meth7.cat,  data = matched.test.incl.pData), silent =T)
  
  ### overall cohort: 31/10/18 - mstatus, LCA, MYCMYCN, meth7, q13loss
  try(multivar.cox.OS.combined.cat <- cox.multivar.surv_5(time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.cat.incl, FacA= matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat, FacC= matched.test.incl.pData$mstatus, FacD =  matched.test.incl.pData$q13loss, FacE= matched.test.incl.pData$meth7.cat,  data = matched.test.incl.pData), silent =T)
  try(multivar.cox.OS.combined.contin <- cox.multivar.surv_5(time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.incl, FacA= matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat, FacC= matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$q13loss, FacE= matched.test.incl.pData$meth7.cat,  data = matched.test.incl.pData), silent =T)
  
  
  ### cox multivariate with established risk factors per PNET5 that are significant in univariate analysis in the RNA cohort (n=217): LCA, MYCMYCN, mstatus, WNT vs non WNT
  ### 31/10/18: 
  
  try(multivar.cox.OS.PNET5.cat <- cox.multivar.surv.PNET5_4 (time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.cat.incl, FacA = matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat,  FacC = matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$WNT1_0, data = matched.test.incl.pData), silent =T)
  try(multivar.cox.OS.PNET5.contin <- cox.multivar.surv.PNET5_4 (time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.incl, FacA = matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat,  FacC = matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$WNT1_0, data = matched.test.incl.pData), silent =T)
  
  ### pre-31/10/18: LCA, mstatus,  MYCMYCN.cat, resection, methylation status (4 subgroup), TP53.cat added 5/10/17
  
  # try(multivar.cox.OS.PNET5.cat <- cox.multivar.surv.PNET5_6 (time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.cat.incl, FacA = matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat,  FacC = matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$resection, FacE = matched.test.incl.pData$meth, FacF = matched.test.incl.pData$TP53.cat, data = matched.test.incl.pData), silent =T)
  # try(multivar.cox.OS.PNET5.contin <- cox.multivar.surv.PNET5_6 (time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.incl, FacA = matched.test.incl.pData$LCA, FacB = matched.test.incl.pData$MYCMYCN.cat,  FacC = matched.test.incl.pData$mstatus, FacD = matched.test.incl.pData$resection, FacE = matched.test.incl.pData$meth, FacF = matched.test.incl.pData$TP53.cat, data = matched.test.incl.pData), silent =T)
  
  ### cox multivariable in PNET 5, G3G4 group
  ### 31/10/18: using mstatus, LCA, MYC (significant in univariate logistic regression)
  
  try(multivar.cox.OS.PNET5.G3G4.cat <- cox.multivar.surv_3 (time = matched.G3G4.incl.pData$Followup, event = OS.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl, FacA = matched.G3G4.incl.pData$MYC.cat, FacB = matched.G3G4.incl.pData$LCA, FacC = matched.G3G4.incl.pData$mstatus, data = matched.G3G4.incl.pData), silent =T)
  try(multivar.cox.OS.PNET5.G3G4.contin <- cox.multivar.surv_3 (time = matched.G3G4.incl.pData$Followup, event = OS.G3G4.bin.incl, marker = matched.goi.vsd.G3G4.incl, FacA = matched.G3G4.incl.pData$MYC.cat, FacB = matched.G3G4.incl.pData$LCA, FacC = matched.G3G4.incl.pData$mstatus, data = matched.G3G4.incl.pData), silent =T)
  
  # try(multivar.cox.OS.PNET5.G3G4.cat <- cox.multivar.surv.PNET5_6 (time = matched.G3G4.incl.pData$Followup, event = OS.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl, FacA = matched.G3G4.incl.pData$MYCMYCN.cat, FacB = matched.G3G4.incl.pData$LCA, FacC = matched.G3G4.incl.pData$mstatus, FacD = matched.G3G4.incl.pData$resection, FacE = matched.G3G4.incl.pData$meth, FacF = matched.G3G4.incl.pData$TP53.cat, data = matched.G3G4.incl.pData), silent =T)
  # try(multivar.cox.OS.PNET5.G3G4.contin <- cox.multivar.surv.PNET5_6 (time = matched.G3G4.incl.pData$Followup, event = OS.G3G4.bin.incl, marker = matched.goi.vsd.G3G4.incl, FacA = matched.G3G4.incl.pData$MYCMYCN.cat, FacB = matched.G3G4.incl.pData$LCA, FacC = matched.G3G4.incl.pData$mstatus,  FacD = matched.G3G4.incl.pData$resection, FacE = matched.G3G4.incl.pData$meth, FacF = matched.G3G4.incl.pData$TP53.cat, data = matched.G3G4.incl.pData), silent =T)
  
  ### cox multivariate with established risk factors per Schwalbe only
  
  ### subgroup specific for SHH_Child (which is SHH_Old in our dataframe)
  ### 31/10/18: SHH_Old significant factors in univariate logistic regression were LCA, MYCN, TP53 mutational status
  
  try(multivar.cox.OS.SHH.old.cat <- cox.multivar.surv_3(time = matched.SHH.old.incl.pData$Followup, event = OS.SHH.old.bin.incl, marker = matched.goi.vsd.cat.SHH.old.incl, FacA = matched.SHH.old.incl.pData$MYCN.cat, FacB = matched.SHH.old.incl.pData$LCA, FacC = matched.SHH.old.incl.pData$TP53.cat, data = matched.SHH.old.incl.pData ), silent =T)
  try(multivar.cox.OS.SHH.old.contin <- cox.multivar.surv_3(time = matched.SHH.old.incl.pData$Followup, event = OS.SHH.old.bin.incl, marker = matched.goi.vsd.SHH.old.incl, FacA = matched.SHH.old.incl.pData$MYCN.cat, FacB = matched.SHH.old.incl.pData$LCA, FacC = matched.SHH.old.incl.pData$TP53.cat, data = matched.SHH.old.incl.pData ), silent =T)
  
  ### Mstatus, TP53 for SHH-child, MYCN, matched.SHH.old.incl.pData ### analysis pre 31/10/18
  # try(multivar.cox.OS.SHH.old.cat <- cox.multivar.surv.Schwalbe_3(time = matched.SHH.old.incl.pData$Followup, event = OS.SHH.old.bin.incl, marker = matched.goi.vsd.cat.SHH.old.incl, FacA = matched.SHH.old.incl.pData$MYCN.cat, FacB = matched.SHH.old.incl.pData$mstatus, FacC = matched.SHH.old.incl.pData$TP53.cat, data = matched.SHH.old.incl.pData ), silent =T)
  # try(multivar.cox.OS.SHH.old.contin <- cox.multivar.surv.Schwalbe_3(time = matched.SHH.old.incl.pData$Followup, event = OS.SHH.old.bin.incl, marker = matched.goi.vsd.SHH.old.incl, FacA = matched.SHH.old.incl.pData$MYCN.cat, FacB = matched.SHH.old.incl.pData$mstatus, FacC = matched.SHH.old.incl.pData$TP53.cat, data = matched.SHH.old.incl.pData ), silent =T)
  

  ### Group 3 and 4 combined: Group 3/4 HR, MYC amplification, q13 loss ### analysis unchanged 31/10/18
  ### high risk methylation group (ie HR Group 3 and HR Group 4 combined) in matched.G3G4.incl.pData$G3G4.high.incl for categorical and continuous variable
  
  try(multivar.cox.OS.Lancet.G3G4.cat <- cox.multivar.surv_3 (time = matched.G3G4.incl.pData$Followup, event = OS.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl, FacA = matched.G3G4.incl.pData$MYC.cat, FacB = matched.G3G4.incl.pData$G3G4.high.incl, FacC = matched.G3G4.incl.pData$q13loss,   data = matched.G3G4.incl.pData), silent =T)
  try(multivar.cox.OS.Lancet.G3G4.contin <- cox.multivar.surv_3 (time = matched.G3G4.incl.pData$Followup, event = OS.G3G4.bin.incl, marker = matched.goi.vsd.G3G4.incl, FacA = matched.G3G4.incl.pData$MYC.cat, FacB = matched.G3G4.incl.pData$G3G4.high.incl, FacC = matched.G3G4.incl.pData$q13loss,   data = matched.G3G4.incl.pData), silent =T)
  
  
  
  #####################################
  
  ### EFS: 
  
  #cat ("km analysis for EFS for children aged 3-16 years, treated with curative intent", sep = "\n")
  
  try(surv.km.EFS.all <- km.log.test.EFS(time = matched.test.incl.pData$EFS, event = EFS.cat.bin.incl, marker = matched.goi.vsd.cat.incl), silent = T)
  try(surv.km.EFS.G3G4 <- km.log.test.EFS(time = matched.G3G4.incl.pData$EFS, event = EFS.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl) , silent = T)
  
  
  ### cox regression analysis, decision not to update EFS with categorical/continuous and multivariate 4/10/17

  #cat ("cox regression for EFS on age 3-16 years", sep = "\n")
  
  try(surv.cox.EFS.incl <-  cox.result.surv (time=matched.test.incl.pData$EFS, event =  EFS.cat.bin.incl, marker = matched.goi.vsd.cat.incl, data = matched.test.incl.pData), silent = T)

  try(surv.cox.EFS.incl.G3G4 <- cox.result.surv (time= matched.G3G4.incl.pData$EFS, event = EFS.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl, data = matched.test.incl.pData), silent = T)
  
  ####################################
  
  cox.p.values.list <- as.list (mget(ls(pattern = "surv.cox.")))
  cox.multivar.p.values.list <- as.list (mget(ls(pattern = "multivar.cox")))
  
  ###################################

  surv.p.values.list <- as.list(mget(ls(pattern="surv.km.")))
  
  ######################################################################
  
  ###### need to return some objects with results, need to generate matched dataframes to be able to name outputs correctly 070917
 
   result.goi <- list (surv.p.values.list=surv.p.values.list, 
                      chi.sq.list = chi.sq.list,
                      reg.log.list = reg.log.list, 
                      cox.p.values.list = cox.p.values.list,
                      cox.multivar.p.values.list = cox.multivar.p.values.list
                      
                      #age.correlation = list.cor.age
                      #linear.reg.age = lin.reg.age ### linear regression and age correlation hashed out as causing problems with novel transcript expression data
   )
  
  sink()
  dev.off()
    return(result.goi) ### used to be # return(res)
}


##########################################################
##########################################################

#### Function called cox.dataframe

### to extract the relevant elements from within a list to generate the cox dataframe for each goi
### input
### p value
### Z score
### Hazard ratio (HR)
### 95CI - upper (U95CI) and lower (L95CI)

### output
### cox dataframe with labelled columns


cox.dataframe <- function(pval, Zscore, HR, L95CI, U95CI){
  cox.pval.assembled <- do.call(rbind, pval)
  adj.cox.pval <- p.adjust(cox.pval.assembled, method = "BH")
  #cox.pval.combined.df <- cbind(cox.pval.assembled, adj.cox.pval)
  #colnames(cox.pval.combined.df)<-c("cox.pval", "cox.adj.pval")
  cox.Zscore.assembled <- do.call(rbind, Zscore)
  cox.HR.assembled <- do.call(rbind, HR)
  cox.L95CI.assembled <- do.call(rbind, L95CI)
  cox.U95CI.assembled <- do.call(rbind, U95CI)
  cox.allresults.df <- cbind(cox.pval.assembled, adj.cox.pval, cox.Zscore.assembled, cox.HR.assembled, cox.L95CI.assembled, cox.U95CI.assembled)
  colnames(cox.allresults.df)<- c("cox.pval", "cox.adj.pval", "cox.Zscore","cox.HR", "cox.HR.L95CI", "cox.HR.U95CI" )
  return (cox.allresults.df)
}



########################################################################################

### Function called annotate.HTseq.IDs
### input is a datafile (results.master) generated from clinical_data_master.R 
# listDatasets()
# listMarts()
library(biomaRt)

annotate.HTseq.IDs<-function(HTseq.IDs){
  ensemblID <- gsub("\\..*", '', HTseq.IDs)
  #fetch the annotation information from the Ensembl database
  #ensembl <- useMart('ensembl', dataset='hsapiens_gene_ensembl')### changed 6/11/18
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset='hsapiens_gene_ensembl') ### added 6/11/18
  symbols <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters='ensembl_gene_id', ensemblID, mart=ensembl)
  #combine the annotation with the RNA-Seq data
  annotatedix <- match(ensemblID,symbols$ensembl_gene_id)
  symbols[annotatedix,]->annotatedGenes
  return(cbind(HTseq.IDs,annotatedGenes))
}

### mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "jul2016.archive.ensembl.org");
###  If you are a user of biomaRt (a part of the Bioconductor library) change the host from 'www.biomart.org' to 'www.ensembl.org'. 24/4/18)

#listMarts(mart = NULL, host = "www.ensembl.org", path = "/biomart/martservice", port = 80, includeHosts= FALSE, archive = FALSE, ssl.verifypeer = TRUE, ensemblRedirect = TRUE, verbose = FALSE)
#ensembl <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
#listAttributes(ensembl) ### if wish to add further components to the annotate function that exist within the ensembl database. Adds computing power/time

########################################################################################

### Function called log.reg.dataframe

log.reg.dataframe <- function(pval, OR, L95CI, U95CI){
  logreg.pval.assembled <- do.call(rbind, pval) ### check if recurrent error here when more than one input goi
  logreg.adj.pval <- p.adjust(logreg.pval.assembled, method = "BH")
  logreg.OR.assembled <- do.call(rbind, OR)
  logreg.L95CI.assembled <- do.call(rbind, L95CI)
  logreg.U95CI.assembled <- do.call(rbind, U95CI)
  logreg.allresults.df <- cbind(logreg.pval.assembled, logreg.adj.pval, logreg.OR.assembled, logreg.L95CI.assembled, logreg.U95CI.assembled)
  colnames(logreg.allresults.df)<- c("logreg.pval", "logreg.adj.pval","logreg.OR", "logreg.OR.L95CI", "logreg.OR.U95CI" )
  return (logreg.allresults.df)
}

########################################################################################

### Function called gp.style.filter that helps created a filtered file, using preset filters related to transcript expression (minimum, variance, fold chance)

gp.style.filter<- function (x,fold.change, delta, prop, base, prop.base, na.rm = TRUE, neg.rm = TRUE) 
{
  if (na.rm == TRUE){x <- x[!is.na(x)]}
  if (neg.rm == TRUE){x <- x[x > 0]}
  lenx <- length(x)
  if (lenx < 4 || lenx < prop + 1){return(FALSE)}
  srtd <- sort(x)
  if (prop < 1) {
    bot <- lenx * prop
  }else {bot <- prop}
  top <- lenx - bot
  fold.filter<-FALSE
  delta.filter<-FALSE
  base.filter<-FALSE
  if (max(srtd[bot:top]) / min(srtd[bot:top]) > fold.change){fold.filter<-TRUE}
  if (max(srtd[bot:top]) - min(srtd[bot:top]) > delta){delta.filter<-TRUE}
  if (length(which(srtd > base))>(prop.base*lenx)){base.filter<-TRUE}
  if (fold.filter ==TRUE & delta.filter == TRUE & base.filter == TRUE){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
