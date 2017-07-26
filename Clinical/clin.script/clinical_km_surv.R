
### KM survival curve update
### date: 25/7/17
### Aim: 
### 1. to output relevant n, nevents, 5-year OS (or equivalent survival variable) 
### 2. to understand the structure of km survival curves and inbuilt functions to output relevant varialbes in a list
### 3. to plot a graph that shows n, n at risk, number of events on the graph
### 4. to determine if this can be added to the current function in the main function file (currently clinical_data_functions_v4.R) or whether it needs to be created as a new function

### Input files
### in order for the current script to work, load in the following

source (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/clinical_data_7.R")
source (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/clinical_data_functions_v4.R")


library(rms)

###########################################################################
### this formula below outputs the relevant km survival data without plotting it
### where kaplan meier survival curve is generated using the following formula
### survfit (Surv(time, event) ~ variable, data, conf.type = "log-log")

km.OS.test <- survfit(Surv(matched.test.incl.pData$Followup, OS.cat.bin.incl)~ matched.goi.vsd.cat.incl, data = matched.test.incl.pData, conf.type =
                        "log-log")



### the following will provide simple statistics for n, n.events (which is useful) however the median, 95% LCL/UCL are not useful from this object
km.OS.test
# str(km.OS.test)

### to then examine the relevant fields such as time, n.risk, n.event, survival (ie OS), SE, and 95% confidence intervals, then output summary(km.os.test)
### colnames(object) is useful for other statistics outputs, to be able to name the item to output. This is not useful in kaplan-meier survival objects/output

summary(km.OS.test)
str(summary(km.OS.test))
#data_int <- summary(km.OS.test[which(summary(km.OS.test$time >= 5))])
# names(summary(km.OS.test))

### recalling the individual items within the summary(km.OS.test) however this does not seem to be that helpful, the summary(km.OS.test) is the most helpful, so I think it's best generated within the output list
n <- km.OS.test$n
n_events <- km.OS.test$n.event
n_risk <-km.OS.test$n.risk
surv_percent <- km.OS.test$surv
OS_time <- km.OS.test$time
std_error <- km.OS.test$std.err
censor <- km.OS.test$n.censor

### put this data into a dataframe which will allow us to extract relevant data, however there is an error in the standard error column

# km.OS.df <- cbind( OS_time, n_risk, n_events, censor, surv_percent, std_error)

### p value is plotted from the km.log.test.OS function 
### p value was derived within that function using, see below
# surv.p.val <- 1 - pchisq(OS.incl.logrank$chisq, length(OS.incl.logrank$obs)-1)
  

### plotting the output as simple command

plot(km.OS.test)
survplot(km.OS.test)

### adding additional data and updating the axes etc

OS.km.plot <- plot (km.OS.test, yaxt = "n", col = c(""))

### original script which uses the function within the attached function files, which plots the km survival curve, currently without the n, nevent data
# km.log.test.OS.all <- km.log.test.OS(time = matched.test.incl.pData$Followup, event = OS.cat.bin.incl, marker = matched.goi.vsd.cat.incl )


### Therefore summary for Louise:
km.OS.test <- survfit(Surv(matched.test.incl.pData$Followup, OS.cat.bin.incl)~ matched.goi.vsd.cat.incl, data = matched.test.incl.pData, conf.type =
                        "log-log")

summary(km.OS.test)


n <- km.OS.test$n
n_events <- km.OS.test$n.event
n_risk <-km.OS.test$n.risk
surv_percent <- km.OS.test$surv
OS_time <- km.OS.test$time
std_error <- km.OS.test$std.err
censor <- km.OS.test$n.censor

### dataframe with errors
### km.OS.df <- cbind( OS_time, n_risk, n_events, censor, surv_percent, std_error)

### p value derivation

# OS.incl.logrank <- survdiff(Surv(time, event) ~ marker)
OS.incl.logrank <- survdiff(Surv(matched.test.incl.pData$Followup, OS.cat.bin.incl)~ matched.goi.vsd.cat.incl)
surv.p.val <- 1 - pchisq(OS.incl.logrank$chisq, length(OS.incl.logrank$obs)-1)
surv.p.val

###############################################################
### comparison script
### creating a Kaplan-Meier survival curve for overall survival (OS) as an example, 


library(rms)

### to create input object for survplot, must use npsurv

km.OS.test.2 <- npsurv(Surv(matched.test.incl.pData$Followup, OS.cat.bin.incl) ~ matched.goi.vsd.cat.incl)

survplot(fit = km.OS.test.2, col = c("red", "blue"),     
        conf = c("none","bands","bars")[1],
        type = "kaplan meier",
        xlab = "overall survival (years)", ylab = "OS (%)", main = "Biomarker expression (high vs  low) and overall survival (%)", 
        xlim= c(0,10),
        #main = "Biomarker expression and overall survival (OS)"
        #label.curves = TRUE,                                                                                   # label curves directly
        #label.curves = list(keys = c("biomarker-high", "biomarker-low"),  col = c("red", "blue")),             # legend instead of direct label
        #surv.p.value
        label.curves = list(keys = "lines",  col = c("red", "blue")),
        levels.only  = TRUE,                                          # show only levels, no label
        abbrev.label = FALSE,                                         # if label used, abbreviate
        ## fun = function(x) {1 - x},                                 # Cumulative probability plot         
        loglog   = FALSE,                                             # log(-log Survival) plot
        logt     = FALSE,                                             # log time
        time.inc = 2,                                                 # time increment
        dots     = FALSE,                                             # dot grid
        n.risk   = TRUE,                                              # show number at risk
        srt.n.risk = 0,  #strt.n.risk = 1                             # moves the n.risk away from the y axis (default is 0) 
        #sep.n.risk = 0.056, adj.n.risk = 1,
        y.n.risk = 0                                                  # moves the n.risk table below the x axis the distance of 1/3 of ylim as 'auto' (default)
        # cex.n.risk = 0.6
)
      
# lines(km.OS.test.2, col=c("red","black"))#add censor lines: does not work

### script option trying to adjust the legend and axes

library(rms)

km.OS.test.2 <- npsurv(Surv(as.numeric(matched.test.incl.pData$Followup),as.numeric(OS.cat.bin.incl)) ~ matched.goi.vsd.cat.incl)

survplot(fit = km.OS.test.2, col = c("red", "blue"),     
         conf = c("none","bands","bars")[1],
         type = "kaplan meier",
         main = "Biomarker expression (high vs  low) and overall survival (%)", xlab = "overall survival (years)", ylab = "OS (%)",
         xlim= c(0,10),
         #main = "Biomarker expression and overall survival (OS)"
         #label.curves = TRUE,                                                                                   # label curves directly
         #label.curves = list(keys = c("biomarker-high", "biomarker-low"),  col = c("red", "blue")),              # legend instead of direct label
         #surv.p.value
         label.curves = list(keys = "lines",  col = c("red", "blue")),
         levels.only  = TRUE,                                          # show only levels, no label
         abbrev.label = FALSE,                                        # if label used, abbreviate
         ## fun = function(x) {1 - x},                                 # Cumulative probability plot         
         loglog   = FALSE,                                             # log(-log Survival) plot
         logt     = FALSE,                                             # log time
         time.inc = 2,                                                 # time increment
         dots     = FALSE,                                             # dot grid
         n.risk   = TRUE                                               # show number at risk
         ## srt.n.risk = 0, sep.n.risk = 0.056, adj.n.risk = 1,
         ## y.n.risk = 0, cex.n.risk = 0.6
)



# lines(km.OS.test.2, col=c("red","blue"))#add censor lines, did not work

### script from Alice:

library(rms)

nsurvG4<-npsurv(Surv(as.numeric(cohIN$PFS_time), as.numeric(cohIN$PFS)) ~ cohIN$Risk)

survplot(fit  = nsurvG4,
         conf = c("none","bands","bars")[1],
         col=c("red","black"),
         lty=c(1,1),
         main="Standard Risk All, PNET4 trial scheme",           # main title did not work
         xlab ="Numbers at risk / Time (Months)",
         ##title(main="Standard Risk Group 4, Shih scheme"),
         ylab = "Survival Probability",
         ## xlim(0,100),
         ##label.curves = TRUE,                                  # label curves directly
         label.curves = list(keys = "lines"),                   # legend instead of direct label
         levels.only  = TRUE,#FALSE,                            # show only levels, no label
         abbrev.label = FALSE,                                  # if label used, abbreviate
         ## fun = function(x) {1 - x},                          # Cumulative probability plot         
         loglog   = FALSE,                                      # log(-log Survival) plot
         logt     = FALSE,                                      # log time
         #time.inc = 100,                                       # time increment
         dots     = FALSE,#TRUE,                                # dot grid
         n.risk   = TRUE,  # show number at risk
         #y.n.risk = 'auto'
         srt.n.risk = 1
         ## sep.n.risk = 0.056, 
         ##adj.n.risk = 1
         ## y.n.risk = 0, 
         ##cex.n.risk = 0.6
)




#mod <- survfit(Surv(time,status)~ph.ecog, data=lung)
#plot(mod, yaxt="n", col=1:4, ylab="Survival Percent", xlab="Days", main="Lung Cancer Data\nEffect of ph.ecog on Survival")
#axis(2, at=seq(0,1,0.2), labels=paste(seq(0,100,20), "%", sep=""), las=1)



### the below is an attempt to add in relevant labels and move legend, does not work

survplot(fit = km.OS.test.2, col = c("red", "blue"),     
         conf = c("none","bands","bars")[1],
         type = "kaplan meier",
         xlab = "overall survival (years)", ylab = "OS (%)", main = "Biomarker expression (high vs  low) and overall survival (%)", 
         xlim= c(0,10),ylim=NULL,
         #main = "Biomarker expression and overall survival (OS)"
         #label.curves = TRUE,                                                                                   # label curves directly
         #label.curves = list(keys = c("biomarker-high", "biomarker-low"),  col = c("red", "blue")),              # legend instead of direct label
         #surv.p.value#
         plot.new ((OS.names <- c("biomarker - high", "biomarker - low")),
         legend (x="topright", OS.names,  lty= 1:2, col = c("red","blue")),
         axis(2, at=pretty(OS.cat.bin.incl), lab=pretty(OS.cat.bin.incl) * 100, las=TRUE),
         OS.incl.logrank <- survdiff(Surv(matched.test.incl.pData$Followup, OS.cat.bin.incl) ~ matched.goi.vsd.cat.incl),
         surv.p.val <- 1 - pchisq(OS.incl.logrank$chisq, length(OS.incl.logrank$obs)-1),
         text(4,0.1,paste("p =",round(surv.p.val, 3)), pos = 4, cex = 1),
         #label.curves = list(keys = "lines",  col = c("red", "blue")),
         levels.only  = FALSE,                                          # show only levels, no label
         abbrev.label = FALSE,                                        # if label used, abbreviate
         ## fun = function(x) {1 - x},                                 # Cumulative probability plot         
         loglog   = FALSE,                                             # log(-log Survival) plot
         logt     = FALSE,                                             # log time
         time.inc = 2,                                                 # time increment
         dots     = FALSE,                                             # dot grid
         n.risk   = TRUE                                               # show number at risk
         ## srt.n.risk = 0, sep.n.risk = 0.056, adj.n.risk = 1,
         ## y.n.risk = 0, cex.n.risk = 0.6
)


OS.names <- c("biomarker - high", "biomarker - low")
legend (x="topright", OS.names,  lty= 1:2, col = c("red","blue"))
axis(2, at=pretty(event), lab=pretty(event) * 100, las=TRUE)
OS.incl.logrank <- survdiff(Surv(time, event) ~ marker)


### legend.pos = "topright"

        
### comparison script taken from "Drawing survival curves in R.pdf" and from Alice's script
### file input for fit needs to be created using the npsurv function instead of survfit

# tmp.survfitSHH  <- npsurv(Surv(SHH.table$OS_Time, SHH.table$OS_Status) ~ SHH.table$riskGroup) 

survplot(fit  = tmp.survfitSHH,
         conf = c("none","bands","bars")[1],
         xlab = "", ylab = "Survival",
         xlim(0,100),
         label.curves = TRUE,                     # label curves directly
         ## label.curves = list(keys = "lines"),  # legend instead of direct label
         levels.only  = FALSE,                    # show only levels, no label
         abbrev.label = FALSE,                    # if label used, abbreviate
         ## fun = function(x) {1 - x},            # Cumulative probability plot         
         loglog   = FALSE,                        # log(-log Survival) plot
         logt     = FALSE,                        # log time
         time.inc = 100,                          # time increment
         dots     = TRUE,                         # dot grid
         n.risk   = TRUE                         # show number at risk
         ## srt.n.risk = 0, sep.n.risk = 0.056, adj.n.risk = 1,
         ## y.n.risk = 0, cex.n.risk = 0.6
)
### km survival function for comparison

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
