### This information is being stored in separate file to validate_transcript.R file
### The script was used to develop the validate_transcript.R and validate_function.R scripts and as such there may be useful tips in amongst the code
### Created December 13 2018
### Author: Dr Marion Mateos
### Originally devleoped the code using MELK as an example goi

###################################################################################################################################################################

### see if data below can be used to improve the validate_function.R

# km.OS<- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
### keep this until determine km curves within the goi.validate function

km.OS.MELK <- survfit(Surv(eset$OS, eset$Dead)~MELK.cat, type = "kaplan-meier", conf.type= "log")
km.OS.MELK <- survfit(Surv(matched.eset$OS, matched.eset$Dead)~matched.eset$MELK, type = "kaplan-meier", conf.type= "log") ### generates same curve as above code
summary(km.OS.MELK)
plot(km.OS.MELK)

plot(matched.eset.G3G4$MELKexp, xlab = "individual samples", ylab = "MELK expression", main = "Expression of MELK in G3G4 validation cohort")
abline(h=median(matched.eset.G3G4$MELKexp), lty = 1, col = "red") ### v for vertical ie x axis ### h, for horizontal

km.OS.G3G4.MELK <- survfit(Surv(matched.eset.G3G4.incl$OS, matched.eset.G3G4.incl$Dead)~matched.eset.G3G4.incl$MELK, type = "kaplan-meier", conf.type= "log") ### generates same curve as above code
summary(km.OS.G3G4.MELK)
plot(km.OS.G3G4.MELK)


### i) Lancet
### ii) PNET 5

cox.OS.G3G4.MELK.Lancet <- coxph(Surv(matched.eset.G3G4.incl$OS, matched.eset.G3G4.incl$Dead)~matched.eset.G3G4.incl$MELK + matched.eset.G3G4.incl$q13loss + matched.eset.G3G4.incl$meth7_HR +  matched.eset.G3G4.incl$Gender + matched.eset.G3G4.incl$MYC, data = matched.eset.G3G4.incl)

#str(summary(cox.OS.MELK))

cox.n.G3G4.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[4]] ### n             ###cox.OS.MELK$n          ### cox.OS.MELK[[11]] 
cox.nevent.G3G4.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[6]] ### nevent   ###cox.OS.MELK$nevent     ### cox.OS.MELK[[12]] 

summary(cox.OS.G3G4.MELK.Lancet)[[7]] ### this is the table of relevance p value
cox.pval.G3G4.MELK.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[7]][1,5] ### this accesses the p value for MELK (row 1, position 5)
cox.HR.G3G4.MELK.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[7]][1,2]
cox.lower.95CI.G3G4.MELK.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[8]][1,3]
cox.upper.95CI.G3G4.MELK.Lancet <- summary (cox.OS.G3G4.MELK.Lancet)[[8]][1,4]
summary.cox.G3G4.MELK.Lancet <- list(pval = cox.pval.G3G4.MELK.Lancet, HR = cox.HR.G3G4.MELK.Lancet, L95CI = cox.lower.95CI.G3G4.MELK.Lancet, U95CI =cox.upper.95CI.G3G4.MELK.Lancet, n = cox.n.G3G4.Lancet, nevent = cox.nevent.G3G4.Lancet, table = summary(cox.OS.G3G4.MELK.Lancet)[[7]], HR_table = summary(cox.OS.G3G4.MELK.Lancet)$conf.int)  

write.csv (summary.cox.G3G4.MELK.Lancet, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/MELK.multivar.G3G4.csv")

#####################################################################################################################

### ORIGINAL DATA FOR MELK, CAN RUN FOR VISUALISATION HOWEVER CAN GENERATE COX MODEL FROM VALIDATE FUNCTION ABOVE
### CAVEAT: the G3G4 subsetting appears to give different results as only 369 results were included in the manual function below

### match in data with MELK
### this worked below 17/7/18 - 8/8/18 for MELK

MELK <- eset.exprs["ENSG00000165304", ]
plot(MELK)
qqnorm(MELK) ### demonstrates that it is normally distributed
summary(MELK)
MELK.cat <- ifelse(MELK>median(MELK, na.rm = T), "high", "low")

index <- match(names(MELK.cat), rownames(eset))
matched.eset <- eset[index[!is.na(index)],]

summary_cavalli <- summary(pData(matched.eset))
# write.csv (summary_cavalli, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/summary.cavalli.csv")

### add in MELK categorical variable into the matched dataset directly to compare with OS outcomes

matched.eset$MELK <- MELK.cat
matched.eset$MELKexp <- MELK

# km.OS<- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
km.OS.MELK <- survfit(Surv(eset$OS, eset$Dead)~MELK.cat, type = "kaplan-meier", conf.type= "log")
km.OS.MELK <- survfit(Surv(matched.eset$OS, matched.eset$Dead)~matched.eset$MELK, type = "kaplan-meier", conf.type= "log") ### generates same curve as above code
summary(km.OS.MELK)
plot(km.OS.MELK)

### create relevant variables to then use to adjust in multivariable analysis

matched.eset$LCA  <- ifelse((matched.eset$histology=="LCA"), "LCA", "non-LCA")
matched.eset$mets <- ifelse(matched.eset$Met.status_.1.Met._0.M0.=="1","metastatic","non-metastatic")

### include in multivariable cox regression analysis adding in MYC or MYCN

cox.OS.MELK <- coxph(Surv(matched.eset$OS, matched.eset$Dead)~matched.eset$MELK + matched.eset$LCA + matched.eset$Gender + matched.eset$Subgroup + matched.eset$mets + matched.eset$MYC + matched.eset$MYCN, data = matched.eset)

#str(summary(cox.OS.MELK))

cox.n <- summary(cox.OS.MELK)[[4]] ### n             ###cox.OS.MELK$n          ### cox.OS.MELK[[11]] 
cox.nevent <- summary(cox.OS.MELK)[[6]] ### nevent   ###cox.OS.MELK$nevent     ### cox.OS.MELK[[12]] 

summary(cox.OS.MELK)[[7]] ### this is the table of relevance p value
cox.pval.MELK <- summary(cox.OS.MELK)[[7]][1,5] ### this accesses the p value for MELK (row 1, position 5)
cox.HR.MELK <- summary(cox.OS.MELK)[[7]][1,2]
cox.lower.95CI.MELK <- summary(cox.OS.MELK)[[8]][1,3]
cox.upper.95CI.MELK <- summary (cox.OS.MELK)[[8]][1,4]
summary.cox.MELK <- list(pval = cox.pval.MELK, HR = cox.HR.MELK, L95CI = cox.lower.95CI.MELK, U95CI =cox.upper.95CI.MELK, n = cox.n, nevent = cox.nevent, table = summary(cox.OS.MELK)[[7]])  

### putting in the Lancet Oncology factors including q13 loss (added 15/10/18)

cox.OS.MELK.Lancet <- coxph(Surv(matched.eset$OS, matched.eset$Dead)~matched.eset$MELK + matched.eset$meth7 + matched.eset$MYC + matched.eset$q13loss + matched.eset$Gender)

summary(cox.OS.MELK.Lancet)

### create function to output these p value and HR characteristics

### Plot transcript expression in Cavalli dataset, unfiltered, with median expression as cutoff

plot(MELK, xlab = "individual samples", ylab = "MELK expression", main = "Expression of MELK in validation cohort")
abline(h=median(MELK), lty = 1, col = "red") ### v for vertical ie x axis ### h, for horizontal

MELK.high <- which (MELK > median(MELK))
length(MELK.high)
MELK.low <- which (MELK < median (MELK))
length (MELK.low)

#########################################################################################################
### look at expression in Group 3/Group 4. Generate dataframe for all expression data (rather than just matched.eset for MELK)
### rename the matched.eset here
matched.eset.all <- matched.eset 

# sub <- matched.eset.all$Subgroup
sub <- matched.eset.all@phenoData@data$Subgroup ### an alternate way to access the expression data subcolumns

Group3.exp <- matched.eset.all[,which(sub=="Group3")]
fdata <- featureData(matched.eset.all)
Group3.exp@featureData <- fdata ### upload feature data so that it contains data
pData(Group3.exp) ### 144 samples, all expression features (21641)


### Group 3 and 4
G3G4.exp <- matched.eset.all[, which(sub=="Group3"|sub=="Group4")]
pData(G3G4.exp)    ### 470 samples, expression set retained

### matched.eset.all is an expression set that can be used to match any expression set in (goi), pData(eset) converts to dataframe
View(pData(matched.eset.all))

### create matched G3G4 dataframe with MELK expression data
### note that MELK is called in from earlier, therefore already in G3G4.exp dataset

View(pData(G3G4.exp))

matched.eset.G3G4 <- pData(G3G4.exp)
# class(matched.eset.G3G4)


###############################

### Exploring characteristics of the dataset, note that expression dataset needs the prefix(pData(eset)) to convert to dataframe
# dim(G3G4.exp)
# names(pData(G3G4.exp))
# summary(pData(G3G4.exp))
# str(pData(G3G4.exp))
# class(pData(G3G4.exp))
# rownames(pData(G3G4.exp))
# plot(G3G4.exp$MELKexp)
# G3G4.exp$MELK

### MELK-specific script #################################################################

plot(matched.eset.G3G4$MELKexp, xlab = "individual samples", ylab = "MELK expression", main = "Expression of MELK in G3G4 validation cohort")
abline(h=median(matched.eset.G3G4$MELKexp), lty = 1, col = "red") ### v for vertical ie x axis ### h, for horizontal

# MELK.cat.G3G4 <- ifelse(matched.eset.G3G4$MELK>median(matched.eset.G3G4$MELK, na.rm = T), "high", "low")
# index <- match(names(MELK.cat.G3G4), rownames(G3G4.exp))
# matched.eset.G3G4$MELK.cat <- MELK.cat.G3G4
# matched.eset.G3G4$MELKexp <- matched.eset.G3G4$MELK

### create new dataframe to remove the two samples which are coded as SHH_Inf and WNT

matched.eset.G3G4.incl <- matched.eset.G3G4[(matched.eset.G3G4$meth7=="Grp3_LowRisk"|matched.eset.G3G4$meth7=="Grp3_HighRisk"|matched.eset.G3G4$meth7=="Grp4_LowRisk"|matched.eset.G3G4$meth7=="Grp4_HighRisk"),  ]
# matched.eset.G3G4.excl <- matched.eset.G3G4[(matched.eset.G3G4$meth7=="SHH_Inf"|matched.eset.G3G4$meth7=="WNT"),  ] ### note can also use != nomenclature (does not equal)


matched.eset.G3G4.incl$meth7_HR <- matched.eset.G3G4.incl$meth7=="Grp4_HighRisk"|matched.eset.G3G4.incl$meth7 =="Grp3_HighRisk"

# km.OS<- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
km.OS.G3G4.MELK <- survfit(Surv(matched.eset.G3G4.incl$OS, matched.eset.G3G4.incl$Dead)~matched.eset.G3G4.incl$MELK, type = "kaplan-meier", conf.type= "log") ### generates same curve as above code
summary(km.OS.G3G4.MELK)
plot(km.OS.G3G4.MELK)


### multivariable cox regression analysis for 
### i) Lancet
### ii) PNET 5

cox.OS.G3G4.MELK.Lancet <- coxph(Surv(matched.eset.G3G4.incl$OS, matched.eset.G3G4.incl$Dead)~matched.eset.G3G4.incl$MELK + matched.eset.G3G4.incl$q13loss + matched.eset.G3G4.incl$meth7_HR  + matched.eset.G3G4.incl$MYC + matched.eset.G3G4.incl$Gender, data = matched.eset.G3G4.incl)
cox.OS.G3G4.Lancet <- coxph(Surv(matched.eset.G3G4.incl$OS, matched.eset.G3G4.incl$Dead)~ matched.eset.G3G4.incl$q13loss + matched.eset.G3G4.incl$meth7_HR  + matched.eset.G3G4.incl$MYC, data = matched.eset.G3G4.incl )

cox.OS.G3G4.Lancet.output <- summary(cox.OS.G3G4.Lancet)$conf.int

# write.csv(cox.OS.G3G4.Lancet.output, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/Lancet.G3G4.csv")

#str(summary(cox.OS.MELK))

#####################################################################################################################################

cox.n.G3G4.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[4]] ### n             ###cox.OS.MELK$n          ### cox.OS.MELK[[11]] 
cox.nevent.G3G4.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[6]] ### nevent   ###cox.OS.MELK$nevent     ### cox.OS.MELK[[12]] 

summary(cox.OS.G3G4.MELK.Lancet)[[7]] ### this is the table of relevance p value
cox.pval.G3G4.MELK.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[7]][1,5] ### this accesses the p value for MELK (row 1, position 5)
cox.HR.G3G4.MELK.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[7]][1,2]
cox.lower.95CI.G3G4.MELK.Lancet <- summary(cox.OS.G3G4.MELK.Lancet)[[8]][1,3]
cox.upper.95CI.G3G4.MELK.Lancet <- summary (cox.OS.G3G4.MELK.Lancet)[[8]][1,4]
summary.cox.G3G4.MELK.Lancet <- list(pval = cox.pval.G3G4.MELK.Lancet, HR = cox.HR.G3G4.MELK.Lancet, L95CI = cox.lower.95CI.G3G4.MELK.Lancet, U95CI =cox.upper.95CI.G3G4.MELK.Lancet, n = cox.n.G3G4.Lancet, nevent = cox.nevent.G3G4.Lancet, table = summary(cox.OS.G3G4.MELK.Lancet)[[7]], HR_table = summary(cox.OS.G3G4.MELK.Lancet)$conf.int)  

write.csv (summary.cox.G3G4.MELK.Lancet, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/MELK.multivar.G3G4.csv")

# cox.list.95CI <- summary(cox.OS.G3G4.MELK.Lancet)$conf.int  ### this brings up all the hazard ratios and 95CI range



