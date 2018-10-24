### VALIDATION FUNCTION
### AIM: create function for validation of goi in the Cavalli dataset, with data generated from ClinPathAssess function for transcripts analysed from gtf file
### create it for a non MELK goi first
### Author: Dr Marion Mateos
### Date: October 17 2018
### "ENSG00000124588" ### NQO2
### "ENSG00000165304" ### MELK

### input files required
### named goi
### eset ### from validate_transcript R script, generate by running lines 1-153

### output
### list.goi ### list of cox multivariable p values (with and without adjustment for gender, HR, 95 confidence intervals (L95CI, U95CI))

## this function allows validation of a named ENSG ID transcript in the Cavalli dataset G3G4 cohort using Lancet oncology paper. Graphical depiction if required needs to be added. 
### Need to work out on what items to return then run script from start to finish after refresh R (17/10/18, well done)
### 2 different functions depending on whether using categorical or continuous expression data

#######################################################################################################################################
### deriving values from categorical expression data (>median, < median)

 goi <- "ENSG00000124588" ### NQO2
# goi <- "ENSG00000245322" ### Error in data.exprs[goi, ] : subscript out of bounds # LOC256880 ### does not exist in Cavalli dataset
# goi <- "ENSG00000165304"
# goi <- "ENSG00000173818" ### ENDOV

data <- eset 

goi.validate <- function(goi, data){
data.exprs <- exprs(data)
rownames(data.exprs) <- gsub("_at", "", rownames(data.exprs))
goi.exp <- data.exprs[goi, ]
# goi.exp <- data [goi, ]  
goi.cat <- ifelse(goi.exp>median(goi.exp, na.rm = T), "high", "low")
summary <- summary(goi.exp)
index <- match(names(goi.cat), rownames(data.exprs))
# index <- match(names(goi.cat), rownames(data))
matched.goi <- data[index[!is.na(index)],]
matched.goi$goi <- goi.cat
matched.goi$goi.exp <- goi.exp


### generate km surv data for overall cohort see below function, see if can incorporate ### confirm that the data is matched to goi and adjust the variables eset$OS and eset$Dead
### working here 24/10/18 
km.OS.all <- survfit(Surv(matched.goi$OS, matched.goi$Dead)~matched.goi$goi, type = "kaplan-meier", conf.type= "log") ###  matched.goi
summary.km.OS.all <- summary(km.OS.all)
plot(km.OS.all)

OS.all.logrank <- survdiff(Surv(matched.goi$OS, matched.goi$Dead) ~ matched.goi$goi) 
## check what the difference is between includinng categorical vs continuous data

### if wish to modulate Time, Event, then can include these as variables in the function
surv.pval.OS.all <- 1 - pchisq(OS.all.logrank$chisq, length(OS.all.logrank$obs)-1) 

### up to here 24/10/18

### add in the multivariable for the overall cohort prior to dividing into G3G4
### the following is dividing into G3G4

sub <- matched.goi@phenoData@data$meth7 ### an alternate way to access the expression data subcolumns
matched.G3G4 <- matched.goi[,which(sub=="Grp3_LowRisk"|sub=="Grp3_HighRisk"|sub=="Grp4_LowRisk"|sub=="Grp4_HighRisk")] ### changed comma position ### generates G3G4 expression set
matched.G3G4$G3G4_HR <- matched.G3G4$meth7=="Grp4_HighRisk"|matched.G3G4$meth7 =="Grp3_HighRisk"

cox.G3G4.goi <- coxph(Surv(matched.G3G4$OS, matched.G3G4$Dead)~matched.G3G4$goi + matched.G3G4$q13loss + matched.G3G4$G3G4_HR +  matched.G3G4$Gender + matched.G3G4$MYC, data = data)
cox.G3G4.goi_nogender <- coxph(Surv(matched.G3G4$OS, matched.G3G4$Dead)~matched.G3G4$goi + matched.G3G4$q13loss + matched.G3G4$G3G4_HR + matched.G3G4$MYC, data = data)
summary_cox <- list (summary_nogender = summary(cox.G3G4.goi_nogender)$conf.int, summary_genderincl = summary(cox.G3G4.goi)$conf.int)


### km results in G3G4
km.OS.G3G4 <- survfit(Surv(matched.G3G4$OS, matched.G3G4$Dead)~matched.G3G4$goi, type = "kaplan-meier", conf.type= "log")
summary.km.G3G4 <- summary(km.OS.G3G4)
plot(km.OS.G3G4)
OS.G3G4.logrank <- survdiff(Surv(matched.G3G4$OS, matched.G3G4$Dead) ~ matched.G3G4$goi) ### if wish to modulate Time, Event, then can include these as variables in teh function
surv.pval.OS.G3G4 <- 1 - pchisq(OS.G3G4.logrank$chisq, length(OS.G3G4.logrank$obs)-1)

# plot.goi <- plot(matched.G3G4$goi, xlab = "individual samples", ylab = "goi expression", main = "Expression of "goi" in G3G4 validation cohort")### how to insert goi name
list.goi <- list(n = summary(cox.G3G4.goi_nogender)$n,  
                 nevent = summary(cox.G3G4.goi_nogender)$nevent,  
                 pval_nogender = summary(cox.G3G4.goi_nogender)$coefficients, 
                 pval_gender = summary(cox.G3G4.goi)$coefficients, 
                 cox_summary = summary_cox, 
                 km.OS.all = OS.all.logrank, 
                 pval.km.all = surv.pval.OS.all,
                 km.OS.G3G4 = km.OS.G3G4, 
                 pval.km.G3G4 = surv.pval.OS.G3G4)  
return (list.goi)
}




##########################################################################
### continuous expression data

goi.validate.contin <- function(goi, data){
  data.exprs <- exprs(data)
  rownames(data.exprs) <- gsub("_at", "", rownames(data.exprs))
  goi.exp <- data.exprs[goi, ]
  # goi.exp <- data [goi, ]  
  goi.cat <- ifelse(goi.exp>median(goi.exp, na.rm = T), "high", "low")
  summary <- summary(goi.exp)
  index <- match(names(goi.cat), rownames(data.exprs))
  # index <- match(names(goi.cat), rownames(data))
  matched.goi <- data[index[!is.na(index)],]
  matched.goi$goi <- goi.cat     ### categorical expression data
  matched.goi$goi.exp <- goi.exp ### continuous expression data
  
  
  ### generate km surv data for overall cohort see below function, see if can incorporate ### confirm that the data is matched to goi and adjust the variables eset$OS and eset$Dead
  km.OS.all <- survfit(Surv(eset$OS, eset$Dead)~goi.cat, type = "kaplan-meier", conf.type= "log")
  summary.km.OS.all <- summary(km.OS.all)
  plot(km.OS.all)
  OS.all.logrank <- survdiff(Surv(eset$OS, eset$Dead) ~ matched.goi$goi) ### check eset$OS and eset$Dead variables still are for matched data
  ### if wish to modulate Time, Event, then can include these as variables in teh function
  surv.pval.OS.all <- 1 - pchisq(OS.all.logrank$chisq, length(OS.all.logrank$obs)-1) 
  
  ### add in the multivariable for the overall cohort prior to dividing into G3G4
  ### the following is dividing into G3G4
  
  sub <- matched.goi@phenoData@data$meth7 ### an alternate way to access the expression data subcolumns
  matched.G3G4 <- matched.goi[,which(sub=="Grp3_LowRisk"|sub=="Grp3_HighRisk"|sub=="Grp4_LowRisk"|sub=="Grp4_HighRisk")] ### changed comma position ### generates G3G4 expression set
  matched.G3G4$G3G4_HR <- matched.G3G4$meth7=="Grp4_HighRisk"|matched.G3G4$meth7 =="Grp3_HighRisk"
  
  cox.G3G4.goi <- coxph(Surv(matched.G3G4$OS, matched.G3G4$Dead)~matched.G3G4$goi + matched.G3G4$q13loss + matched.G3G4$G3G4_HR +  matched.G3G4$Gender + matched.G3G4$MYC, data = data)
  cox.G3G4.goi_nogender <- coxph(Surv(matched.G3G4$OS, matched.G3G4$Dead)~matched.G3G4$goi + matched.G3G4$q13loss + matched.G3G4$G3G4_HR + matched.G3G4$MYC, data = data)
  summary_cox <- list (summary_nogender = summary(cox.G3G4.goi_nogender)$conf.int, summary_genderincl = summary(cox.G3G4.goi)$conf.int)
  
  ### km results in G3G4
  km.OS.G3G4 <- survfit(Surv(matched.G3G4$OS, matched.G3G4$Dead)~matched.G3G4$goi, type = "kaplan-meier", conf.type= "log")
  summary.km.G3G4 <- summary(km.OS.G3G4)
  plot(km.OS.G3G4)
  OS.G3G4.logrank <- survdiff(Surv(matched.G3G4$OS, matched.G3G4$Dead) ~ matched.G3G4$goi) ### if wish to modulate Time, Event, then can include these as variables in teh function
  surv.pval.OS.G3G4 <- 1 - pchisq(OS.G3G4.logrank$chisq, length(OS.G3G4.logrank$obs)-1)
  
  # plot.goi <- plot(matched.G3G4$goi, xlab = "individual samples", ylab = "goi expression", main = "Expression of "goi" in G3G4 validation cohort")### how to insert goi name
  list.goi <- list(n = summary(cox.G3G4.goi_nogender)$n,  
                   nevent = summary(cox.G3G4.goi_nogender)$nevent,  
                   pval_nogender = summary(cox.G3G4.goi_nogender)$coefficients, 
                   pval_gender = summary(cox.G3G4.goi)$coefficients, 
                   cox_summary = summary_cox, 
                   km.OS.all = OS.all.logrank, 
                   pval.km.all = surv.pval.OS.all,
                   km.OS.G3G4 = km.OS.G3G4, 
                   pval.km.G3G4 = surv.pval.OS.G3G4)  
  return (list.goi)
}


##########################################################################

##########################################################################

### example from km survival curves from clinical_data_functions_master.R
### can you somehow include this in the validate_function.R above

### example inputs
# time <- matched.test.incl.pData$Followup
# event = OS.cat.bin.incl
# marker = matched.goi.vsd.cat.incl

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
  1 - pchisq(OS.incl.logrank$chisq, length(OS.incl.logrank$obs)-1) -> surv.p.val.OS
  text(4,0.1,paste("p =",round(surv.p.val.OS, 3)), pos = 4, cex = 1)
  OS.surv.table <- summary(km.OS.incl)
  return (list(OS.p.val = surv.p.val.OS, 
               OS.surv.table))
  if(out.file!="none"){
    dev.off()
  }
}

