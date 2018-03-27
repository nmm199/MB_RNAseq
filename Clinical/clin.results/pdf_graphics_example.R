### this file was used to start looking at basic PDF graphic generation, using .pdf extension

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

### graphing example


# Create data:
a=c(1:5)
b=c(5,3,4,5,5)
c=c(4,5,4,3,1)

# Make a basic graph
plot( b~a , type="b" , bty="l" , xlab="value of a" , ylab="value of b" , col=rgb(0.2,0.4,0.1,0.7) , lwd=3 , pch=17 , ylim=c(1,5) )
lines(c ~a , col=rgb(0.8,0.4,0.1,0.7) , lwd=3 , pch=19 , type="b" )

# Add a legend
legend("bottomleft", 
       legend = c("Group 1", "Group 2"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.8,0.4,0.1,0.7)), 
       pch = c(17,19), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))