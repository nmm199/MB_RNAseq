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







#########################################################################################################################################
### REMOVED THIS HARDCODING FROM clinical_data_graphics.R 
#########################################################################################################################################

### Hardcoding for graphs (examples)

### Ecdf plots that describe the empirical cumulative distribution frequency

# plot(ecdf(p.km.EFS.all)) 
# plot(ecdf(adjusted.p.km.EFS.all))  ### km.EFS.all.results[,"EFS.adjusted.pval"]

### Histograms

# hist(cox.PFS.cat.G3G4.df[,1])
# hist(cox.PFS.cat.G3G4.df[,2])


### redundant script removed above, for example:
# plotEcdf(x = cox.PFS.cat.G3G4.df[,1], test.name = "cox PFS p values for G3G4", xlab = "p value")
# plotEcdf (x = cox.PFS.cat.G3G4.df[,2], test.name = "cox PFS adj p values for G3G4", xlab = "adjusted p value")
# p.cox.PFS.cat.G3G4 <- cox.PFS.cat.G3G4.df[,1]
# p.adj.cox.PFS.cat.G3G4 <- cox.PFS.cat.G3G4.df[,2]
# plotEcdf.double(x = p.cox.PFS.cat.G3G4, z = p.adj.cox.PFS.cat.G3G4, test.name = "cox PFS in G3G4 (categorical)", xlab = "p value", ylab = "Fn(p value)")
### has been replaced by:
# plotEcdf.double(x = cox.EFS.cat.G3G4.df[, "cox.pval"], z = cox.EFS.cat.G3G4.df[, "cox.adj.pval"], test.name = "cox EFS in G3G4 (categorical)", xlab = "p value", ylab = "Fn(p value)")


### another example
# plotEcdf(x = p.km.EFS.all, test.name = "kaplan meier EFS", xlab = "p value") ### entire cohort fpr just p value
# plotEcdf(x = adjusted.p.km.EFS.all, test.name = "kaplan meier EFS", xlab = "adjusted p value", ylab = "Fn(adjusted p value)") ### entire cohort for adjusted p value
### has been replaced by:
# plotEcdf.double(x = p.km.EFS.all, z = adjusted.p.km.EFS.all,  test.name = "kaplan meier EFS", xlab = "p value", ylab = "Fn(p value)")


### another example of simplified coding
#p.km.EFS.all <- km.EFS.all.results[, "EFS.p.value"]
#adjusted.p.km.EFS.all <- km.EFS.all.results[,"EFS.adjusted.pval"]
#hist(p.km.EFS.all)
#hist(adjusted.p.km.EFS.all)

### replaced by
# plotHist(cox.PFS.cat.G3G4.df[,1], "Cox PFS categorical G3/G4", breaks = 100, xlab = "p-value", cutoff = 0.05)
# plotHist(cox.PFS.cat.G3G4.df[,2], "Cox PFS categorical G3/G4", breaks = 100, xlab = "adjusted p-value", cutoff = 0.05) 
#########################################################################################################################################


# plot(density(x, na.rm = "T"))

### remove the density plot as is not as useful as ecdf

# plot(density(adjusted.p.km.EFS.all))
# lines(density(adjusted.p.km.EFS.all), col = "red")
# lines(density(p.km.EFS.all), col = "dodgerblue") ### this will overlay the unadjusted p value against the adjusted p value

