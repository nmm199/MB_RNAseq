
getwd()
sink()
# "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/
pdf ("plot.pdf") ### adjust to which data are being used
par(mfrow = c(2,2))
### set file for log output
# log.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/pDatalog.report.data.txt"


try(surv.km.OS.all <- km.log.test.OS(time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.cat.incl ), silent = T)
try(surv.km.OS.G3G4 <- km.log.test.OS(time = matched.G3G4.incl.pData$Followup, event = OS.G3G4.bin.incl, marker = matched.goi.vsd.cat.G3G4.incl), silent = T)
try(surv.km.OS.SHH <- km.log.test.OS(time = matched.SHH.incl.pData$Followup, event = OS.SHH.bin.incl, marker = matched.goi.vsd.cat.SHH.incl), silent = T)
try(surv.km.OS.SHH.old <- km.log.test.OS(time = matched.SHH.old.incl.pData$Followup, event = OS.SHH.old.bin.incl, marker = matched.goi.vsd.cat.SHH.old.incl), silent = T)

dev.off()


### cox regression analysis

