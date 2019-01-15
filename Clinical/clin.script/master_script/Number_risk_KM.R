### Number at risk graphs
### Author: Marion Mateos 
### January 10 2019

library(rms)
# install.packages("survminer")
library(survminer)

### to use this after running clinical_data_master and clinical_data_function_master.R. 
### need to specify goi
### then specify surv curve of interest

### need to find way of adding title

### updating the npsurv to include title and pvalue pasting

### contains 2 functions
### KM_risk_OS for OS labelling
### KM_risk_PFS for PFS labelling

### input
### time e.g OS or PFS
### event e.g death/relapse
### marker (goi of interest matched to data frame of interest e.g overall vs G3G4)

### output
### KM survival curves with labelled axes (PFS probability, range 0-1)
### n at risk and event tables



###############################################################################################################################
###############################################################################################################################


### read in functions file

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/KM_number_risk_functions.R")


###############################################################################################################################
###############################################################################################################################

### using function
OS.km.overall.risk <- KM_risk_OS(matched.test.incl.pData$Followup, OS.cat.bin.incl, matched.goi.vsd.cat.incl )
OS.km.G3G4.risk <- KM_risk_OS(time = matched.G3G4.incl.pData$Followup, event = OS.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl )
PFS.km.overall.risk <- KM_risk_PFS(time = matched.test.incl.pData$PFS, event = relapse.bin.incl, marker = matched.goi.vsd.cat.incl)
PFS.km.G3G4.risk <- KM_risk_PFS(time = matched.G3G4.incl.pData$PFS, event = relapse.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl) ### need to adjust d.p in function
PFS.km.G4.risk <- KM_risk_PFS (time = matched.G4.incl.pData$PFS, event = relapse.G4.bin.incl, marker = matched.goi.vsd.cat.G4.incl)



###############################################################################################################################
###############################################################################################################################

### original script for km n at risk plots
surv.km.marker <- npsurv(Surv(matched.test.incl.pData$Followup, OS.cat.bin.incl) ~ matched.goi.vsd.cat.incl)

survplot(fit  = surv.km.marker, 
         conf = c("none","bands","bars")[1],
         col=c("red","blue"),
         lty=c(1,2), ### 2 is dashed line
         # title(main = "Biomarker and survival"),
         xlab ="Numbers at risk / Time (Years)",
         ylab = "OS probability",
         # main="Biomarker status and survival",
         xlim=c(0,10),
         label.curves = list(keys = "lines"),     # legend instead of direct label
         levels.only  = TRUE,#FALSE,              # show only levels, no label
         abbrev.label = FALSE,                    # if label used, abbreviate
         ## fun = function(x) {1 - x},            # Cumulative probability plot         
         loglog   = FALSE,                        # log(-log Survival) plot
         logt     = FALSE,                        # log time
         #time.inc = 100,                         # time increment
         dots     = FALSE,#TRUE,                         # dot grid
         n.risk   = TRUE,  # show number at risk
         #y.n.risk = 'auto'
         srt.n.risk = 1
         ## sep.n.risk = 0.056, 
         ##adj.n.risk = 1
         ## y.n.risk = 0, 
         ##cex.n.risk = 0.6
)

# axis(2, at=pretty(OS.cat.bin.incl), lab=pretty(OS.cat.bin.incl) * 100, las=TRUE) ### this alters the y label to *100 however need to remove default scale 
OS.incl.logrank <- survdiff(Surv(matched.test.incl.pData$Followup, OS.cat.bin.incl) ~ matched.goi.vsd.cat.incl)
surv.p.val.OS <- 1 - pchisq(OS.incl.logrank$chisq, length(OS.incl.logrank$obs)-1)
text(4,0.2,paste("p =",round(surv.p.val.OS, 3)), pos = 4, cex = 1)


###############################################################################################################################
###############################################################################################################################

ggsurvplot(fit=surv.km.marker, data = matched.test.incl.pData, 
           title = "Biomarker and survival", ### subtitle = "Based on Kaplan-Meier estimates",
           # caption = "created with survminer",
           font.title = c(16, "bold", "darkblue"),
           # font.subtitle = c(15, "bold.italic", "purple"),
           # font.caption = c(14, "plain", "orange"),
           font.x = c(14, "bold.italic", "red"),
           font.y = c(14, "bold.italic", "darkred"),
           font.tickslab = c(12, "plain", "darkgreen"),
           ########## risk table #########,
           risk.table = TRUE,
           # risk.table.title = "Note the risk set sizes",
           # risk.table.subtitle = "and remember about censoring.",
           # risk.table.caption = "source code: website.com",
           risk.table.height = 0.25) ### was 0.45

