### VALIDATION FUNCTION
### AIM: create function for validation of goi in the Cavalli dataset, with data generated from ClinPathAssess function for transcripts analysed from gtf file
### create it for a non MELK goi first
### Author: Dr Marion Mateos
### Date: October 17 2018
### "ENSG00000124588" ### NQO2
### "ENSG00000165304" ### MELK

## AIM: Validation of a named ENSG ID transcript in the Cavalli dataset overall and G3G4 cohort using Lancet oncology paper, primarily using multivariable cox regression 


### 2 different functions depending on whether using categorical or continuous expression data
### Function 1: goi.validate 
### Function 2: goi.validate.cont


### input files required
### named goi
### eset ### from validate_transcript R script, generate by running lines 1-153

### output
### list.goi ### list of cox multivariable p values (without adjustment for gender, including HR, 95 confidence intervals (L95CI, U95CI)) and Kaplan Meier p values for overall and G3G4
### list.goi.cont ### list of cox multivariable p values for continuous expression variable

#######################################################################################################################################
### deriving values from categorical expression data (>median, < median)

# goi <- "ENSG00000124588" ### NQO2
# goi <- "ENSG00000245322" ### Error in data.exprs[goi, ] : subscript out of bounds # LOC256880 ### does not exist in Cavalli dataset
# goi <- "ENSG00000165304"
# goi <- "ENSG00000173818" ### ENDOV
# goi <- "ENSG00000128626" ### MRPS12

data <- eset 

### this function uses categorical expression data as the main input

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

### add in the multivariable for the overall cohort prior to dividing into G3G4
### factors based on Approach November 2018 utilising factors that were significant in univariate analysis in the multivariable model
### overall: mstatus, LCA, MYC or MYCN amplified (one variable), 7 molecular groups, q13loss
### new variables MYCMYCN amplified, LCA pathology

matched.goi$MYCMYCN <- ifelse (matched.goi$MYC=="1"| matched.goi$MYCN =="1", 1, 0) 
matched.goi$LCA <- ifelse (matched.goi$histology == "LCA", 1, 0)

# table (matched.goi$MYC) ### confirms that the new variable is correct as MYC or MYCN amplified
# table (matched.goi$MYCN)
# table (matched.goi$MYCMYCN) ### confirms that the new variable is correct as MYC or MYCN amplified ### n=69 either MYC or MYCN amplified
matched.goi$mstatus <- matched.goi$Met.status_.1.Met._0.M0.

cox.overall.goi <- coxph(Surv(matched.goi$OS, matched.goi$Dead)~matched.goi$goi + matched.goi$LCA + matched.goi$mstatus + matched.goi$q13loss + matched.goi$meth7 + matched.goi$MYCMYCN, data = data)

###### updated 13/12/18:
cox.pval.goi <- summary(cox.overall.goi)[[7]][1,5] ### this accesses the p value
cox.HR.goi <- summary(cox.overall.goi)[[7]][1,2]
cox.lower.95CI.goi <- summary(cox.overall.goi)[[8]][1,3]
cox.upper.95CI.goi <- summary(cox.overall.goi)[[8]][1,4]
cox.n.goi <- summary(cox.overall.goi)[[4]]   ### n             ###cox.OS.MELK$n          ### cox.OS.MELK[[11]] 
cox.nevent.goi <- summary(cox.overall.goi)[[6]]
summary.cox.goi <- list(pval = cox.pval.goi, HR = cox.HR.goi, L95CI = cox.lower.95CI.goi, U95CI = cox.upper.95CI.goi, n = cox.n.goi  , nevent = cox.nevent.goi, table = summary(cox.overall.goi)[[7]], HR_table = summary(cox.overall.goi)$conf.int)  



### the following is dividing into G3G4

sub <- matched.goi@phenoData@data$meth7 ### an alternate way to access the expression data subcolumns
matched.G3G4 <- matched.goi[,which(sub=="Grp3_LowRisk"|sub=="Grp3_HighRisk"|sub=="Grp4_LowRisk"|sub=="Grp4_HighRisk")] ### changed comma position ### generates G3G4 expression set
matched.G3G4$G3G4_HR <- matched.G3G4$meth7=="Grp4_HighRisk"|matched.G3G4$meth7 =="Grp3_HighRisk"

# cox.G3G4.goi <- coxph(Surv(matched.G3G4$OS, matched.G3G4$Dead)~matched.G3G4$goi + matched.G3G4$q13loss + matched.G3G4$G3G4_HR +  matched.G3G4$Gender + matched.G3G4$MYC, data = data)
cox.G3G4.goi_nogender <- coxph(Surv(matched.G3G4$OS, matched.G3G4$Dead)~matched.G3G4$goi + matched.G3G4$q13loss + matched.G3G4$G3G4_HR + matched.G3G4$MYC, data = data)
# summary_cox <- list (summary_nogender = summary(cox.G3G4.goi_nogender)$conf.int, summary_genderincl = summary(cox.G3G4.goi)$conf.int)


### this section was added in Dec 11 2018 to improve the cox output

cox.n.G3G4.goi <- summary(cox.G3G4.goi_nogender)[[4]]   ### n             ###cox.OS.MELK$n          ### cox.OS.MELK[[11]] 
cox.nevent.G3G4.goi <- summary(cox.G3G4.goi_nogender)[[6]] ### nevent   ###cox.OS.MELK$nevent     ### cox.OS.MELK[[12]] 

# summary(cox.OS.G3G4.MELK.Lancet)[[7]] ### this is the table of relevance p value
cox.pval.G3G4.goi <- summary(cox.G3G4.goi_nogender)[[7]][1,5] ### this accesses the p value for MELK (row 1, position 5)
cox.HR.G3G4.goi <- summary(cox.G3G4.goi_nogender)[[7]][1,2]
cox.lower.95CI.G3G4.goi <- summary(cox.G3G4.goi_nogender)[[8]][1,3]
cox.upper.95CI.G3G4.goi <- summary(cox.G3G4.goi_nogender)[[8]][1,4]
summary.cox.G3G4.goi <- list(pval = cox.pval.G3G4.goi, HR = cox.HR.G3G4.goi, L95CI = cox.lower.95CI.G3G4.goi, U95CI =cox.upper.95CI.G3G4.goi, n = cox.n.G3G4.goi , nevent = cox.nevent.G3G4.goi, table = summary(cox.G3G4.goi_nogender)[[7]], HR_table = summary(cox.G3G4.goi_nogender)$conf.int)  


### km results in G3G4
km.OS.G3G4 <- survfit(Surv(matched.G3G4$OS, matched.G3G4$Dead)~matched.G3G4$goi, type = "kaplan-meier", conf.type= "log")
summary.km.G3G4 <- summary(km.OS.G3G4)
plot(km.OS.G3G4)
OS.G3G4.logrank <- survdiff(Surv(matched.G3G4$OS, matched.G3G4$Dead) ~ matched.G3G4$goi) ### if wish to modulate Time, Event, then can include these as variables in teh function
surv.pval.OS.G3G4 <- 1 - pchisq(OS.G3G4.logrank$chisq, length(OS.G3G4.logrank$obs)-1)

# plot.goi <- plot(matched.G3G4$goi, xlab = "individual samples", ylab = "goi expression", main = "Expression of "goi" in G3G4 validation cohort")### how to insert goi name
list.goi <- list(goi_overall_cox = summary.cox.goi,
                 km.OS.all = OS.all.logrank, 
                 pval.km.all = surv.pval.OS.all,
                 goi_G3G4_cox = summary.cox.G3G4.goi, 
                 n = summary(cox.G3G4.goi_nogender)$n,  
                 nevent = summary(cox.G3G4.goi_nogender)$nevent,  
                 pval_nogender = summary(cox.G3G4.goi_nogender)$coefficients,
                 # pval_gender = summary(cox.G3G4.goi)$coefficients, 
                 # cox_summary = summary_cox,
                 km.OS.G3G4 = km.OS.G3G4, 
                 pval.km.G3G4 = surv.pval.OS.G3G4
                 )  
return (list.goi)
}




##########################################################################
### continuous expression data

goi.validate.contin <- function(goi, data){
  data.exprs <- exprs(data)
  rownames(data.exprs) <- gsub("_at", "", rownames(data.exprs))
  goi.exp <- data.exprs[goi, ]
  # goi.exp <- data [goi, ]  
  summary <- summary(goi.exp)
  index <- match(names(goi.exp), rownames(data.exprs))
  # index <- match(names(goi.cat), rownames(data))
  matched.goi <- data[index[!is.na(index)],]
  matched.goi$goi.exp <- goi.exp ### continuous expression data
  
  
  
  matched.goi$MYCMYCN <- ifelse (matched.goi$MYC=="1"| matched.goi$MYCN =="1", 1, 0) 
  matched.goi$LCA <- ifelse (matched.goi$histology == "LCA", 1, 0)
  
  # table (matched.goi$MYC) ### confirms that the new variable is correct as MYC or MYCN amplified
  # table (matched.goi$MYCN)
  # table (matched.goi$MYCMYCN) ### confirms that the new variable is correct as MYC or MYCN amplified ### n=69 either MYC or MYCN amplified
  matched.goi$mstatus <- matched.goi$Met.status_.1.Met._0.M0.
  
  cox.overall.goi.cont <- coxph(Surv(matched.goi$OS, matched.goi$Dead)~matched.goi$goi.exp + matched.goi$LCA + matched.goi$mstatus + matched.goi$q13loss + matched.goi$meth7 + matched.goi$MYCMYCN, data = data)
  
  ###### updated 13/12/18:
  cox.pval.goi <- summary(cox.overall.goi.cont)[[7]][1,5] ### this accesses the p value
  cox.HR.goi <- summary(cox.overall.goi.cont)[[7]][1,2]
  cox.lower.95CI.goi <- summary(cox.overall.goi.cont)[[8]][1,3]
  cox.upper.95CI.goi <- summary(cox.overall.goi.cont)[[8]][1,4]
  cox.n.goi <- summary(cox.overall.goi.cont)[[4]]   ### n             ###cox.OS.MELK$n          ### cox.OS.MELK[[11]] 
  cox.nevent.goi <- summary(cox.overall.goi.cont)[[6]]
  summary.cox.goi.cont <- list(pval = cox.pval.goi, HR = cox.HR.goi, L95CI = cox.lower.95CI.goi, U95CI = cox.upper.95CI.goi, n = cox.n.goi  , nevent = cox.nevent.goi, table = summary(cox.overall.goi.cont)[[7]], HR_table = summary(cox.overall.goi.cont)$conf.int)  
  
  
   ### add in the multivariable for the overall cohort prior to dividing into G3G4
  ### the following is dividing into G3G4
  
  sub <- matched.goi@phenoData@data$meth7 ### an alternate way to access the expression data subcolumns
  matched.G3G4 <- matched.goi[,which(sub=="Grp3_LowRisk"|sub=="Grp3_HighRisk"|sub=="Grp4_LowRisk"|sub=="Grp4_HighRisk")] ### changed comma position ### generates G3G4 expression set
  matched.G3G4$G3G4_HR <- matched.G3G4$meth7=="Grp4_HighRisk"|matched.G3G4$meth7 =="Grp3_HighRisk"
  
  
  cox.G3G4.goi.cont <- coxph(Surv(matched.G3G4$OS, matched.G3G4$Dead)~matched.G3G4$goi.exp + matched.G3G4$q13loss + matched.G3G4$G3G4_HR + matched.G3G4$MYC, data = data)
  summary_cox <- list (summary_goi_G3G4_cont = summary(cox.G3G4.goi.cont)$conf.int, summary_genderincl = summary(cox.G3G4.goi.cont)$conf.int)
  
  ### km results in G3G4 and overall is not being included as this calculates a log rank p value for EACH level of expression, therefore not meaningful
  
  # plot.goi <- plot(matched.G3G4$goi, xlab = "individual samples", ylab = "goi expression", main = "Expression of "goi" in G3G4 validation cohort")### how to insert goi name
  list.goi.cont <- list(goi_overall_cox_cont = summary.cox.goi.cont, 
                   cox_summary = summary_cox, 
                   pval_G3G4_cox = summary(cox.G3G4.goi.cont)$coefficients, 
                   n = summary(cox.G3G4.goi.cont)$n,  
                   nevent = summary(cox.G3G4.goi.cont)$nevent 
                   )
  return (list.goi.cont)
}


##########################################################################

##########################################################################

### example from km survival curves from clinical_data_functions_master.R for graphics
