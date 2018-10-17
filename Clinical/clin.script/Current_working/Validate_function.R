### VALIDATION FUNCTION
### AIM: create function for validation of goi in the Cavalli dataset, with data generated from ClinPathAssess function for transcripts analysed from gtf file
### create it for a non MELK goi first
### Author: Dr Marion Mateos
### Date: October 17 2018
### "ENSG00000124588" ### NQO2

## this function allows validation of a named ENSG ID transcript in the Cavalli dataset G3G4 cohort using Lancet oncology paper. Graphical depiction if required needs to be added. 
### Need to work out on what items to return then run script from start to finish after refresh R (17/10/18, well done)

goi <- "ENSG00000124588"
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
# km.OS.goi <- 

sub <- matched.goi@phenoData@data$meth7 ### an alternate way to access the expression data subcolumns

matched.G3G4 <- matched.goi[,which(sub=="Grp3_LowRisk"|sub=="Grp3_HighRisk"|sub=="Grp4_LowRisk"|sub=="Grp4_HighRisk")] ### changed comma position ### generates G3G4 expression set
matched.G3G4$G3G4_HR <- matched.G3G4$meth7=="Grp4_HighRisk"|matched.G3G4$meth7 =="Grp3_HighRisk"

cox.G3G4.goi <- coxph(Surv(matched.G3G4$OS, matched.G3G4$Dead)~matched.G3G4$goi + matched.G3G4$q13loss + matched.G3G4$G3G4_HR +  matched.G3G4$Gender + matched.G3G4$MYC, data = data)
cox.G3G4.goi_nogender <- coxph(Surv(matched.G3G4$OS, matched.G3G4$Dead)~matched.G3G4$goi + matched.G3G4$q13loss + matched.G3G4$G3G4_HR + matched.G3G4$MYC, data = data)
summary_cox <- list (summary_nogender = summary(cox.G3G4.goi_nogender)$conf.int, summary_genderincl = summary(cox.G3G4.goi)$conf.int)


# plot.goi <- plot(matched.G3G4$goi, xlab = "individual samples", ylab = "goi expression", main = "Expression of "goi" in G3G4 validation cohort")### how to insert goi name
# return (summary_cox, )
}



========= based on script below====================

MELK.cat <- ifelse(MELK>median(MELK, na.rm = T), "high", "low")

index <- match(names(MELK.cat), rownames(eset))
matched.eset <- eset[index[!is.na(index)],]


### add in MELK categorical variable into the matched dataset directly to compare with OS outcomes

matched.eset$MELK <- MELK.cat
matched.eset$MELKexp <- MELK

# km.OS<- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log")
km.OS.MELK <- survfit(Surv(eset$OS, eset$Dead)~MELK.cat, type = "kaplan-meier", conf.type= "log")
km.OS.MELK <- survfit(Surv(matched.eset$OS, matched.eset$Dead)~matched.eset$MELK, type = "kaplan-meier", conf.type= "log") ### generates same curve as above code
summary(km.OS.MELK)
plot(km.OS.MELK)



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
