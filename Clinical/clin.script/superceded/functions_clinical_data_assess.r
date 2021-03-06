### File information:

### This is a list of functions that were written for the "clinical_data_4_update.R" script analysis
### Many of these functions are generic for survival analysis
### Aims of this document are to provide tools to perform univariate analysis in a group of children treated for medulloblastoma. 
### All kaplan-meier survival analyses provide graphical output with labelled axes and p value
### The survival analysis is designed in a cohort of children who received therapy with curative intent (including cranio-spinal irradiation), age 3-16 years

### Author: Dr Marion Mateos
### Date: 19 July 2017



###########################################################

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
  table.temp <- table(x, y) 
  table.temp.perc <- prop.table(table.temp)*100
  summary.table(table.temp)
  chi.test.temp <- chisq.test(table.temp) 
  chi.test.temp.stat <- c(stat=chi.test.temp$statistic, p.value=chi.test.temp$p.value) 
  chi.test.temp.res <- chi.test.temp$residuals
  aheatmap(chi.test.temp$residuals, Rowv=NA, Colv = NA)
  list.temp <- list  (chi.test.temp,
                      table.temp, 
                      table.temp.perc,
                      chi.test.temp.res
  )
  
  return(list.temp)
}



###########################################################

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




###########################################################

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
  out.stats <- summary.temp$coefficients
  temp.stats <- cbind(OR=exp(coef(temp.reg)), exp(confint(temp.reg)))
  plot(x,y)
  abline(lm (x ~y))
  list.temp <- list (summary.temp,out.stats,p.val=out.stats[2,4], temp.stats)
  return(list.temp )
}


###########################################################

### Function number 4
### function called "logisticRegression" which runs glm (logistic regression) for association between biomarker (independent) and outcome variable

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


#logisticRegression <- function(x,y, data) {
  
 # temp.glm <- glm (x ~ y, family = binomial (link = 'logit'), data= data)
 # stats.temp <- cbind(OR=exp(coef(temp.glm)), exp(confint(temp.glm)))
 # summary.temp <- list (stats.temp, summary(temp.glm)$coefficients, summary(temp.glm)$data$Event)
 # return(summary.temp)
#}



logisticRegression <- function(x,y, data) {
  
 temp.glm <- glm (x ~ y, family = binomial (link = 'logit'), data= data)
 return(cbind(OR=exp(coef(temp.glm)), exp(confint(temp.glm)), summary(temp.glm)$coefficients))
}

### Note: previously unable to make logistic regression function work, error message "y values must be 0 <= y <= 1", this working function above used glm(x ~y) rather than glm(y~x)


###########################################################

### Function Number 5
### Function entitled "km.log.test" to create kaplan meier survival curves for age 3=16 year old children treated with curative intent, MB
### input:
## time
## event
## marker

### output:
## km survival curve
## p values plotted on graph
## y axis with values as %
## legend and p value for survival analysis 


km.log.test <- function(time, event, marker, out.file = "none"){
  if(out.file!="none"){
    pdf(out.file)
  }
  km.PFS.incl <- survfit(Surv(time, event)~marker, type = "kaplan-meier", conf.type = "log-log")
  plot(km.PFS.incl, yaxt="n", col = c("red", "blue"),xlab = "time to progression/relapse (years)", ylab = "PFS (%)", xlim = c(0,10), main = "Biomarker expression and progression-free survival (PFS)",  lty = 1:2)
  summary(km.PFS.incl)
  PFS.names <- c("biomarker - high", "biomarker - low")
  legend (x="topright", PFS.names,  lty= 1:2, col = c("red","blue"))
  axis(2, at=pretty(event), lab=pretty(event) * 100, las=TRUE)
  PFS.incl.logrank <- survdiff(Surv(time, event) ~ marker)
  1 - pchisq(PFS.incl.logrank$chisq, length(PFS.incl.logrank$obs)-1) -> surv.p.val
  text(4,0.1,paste("p =",round(surv.p.val, 3)), pos = 4, cex = 1)
  if(out.file!="none"){
    dev.off()
  }
}


###

km.OS.test <- survfit(Surv(matched.test.incl.pData$Followup, OS.cat.bin.incl)~ matched.goi.vsd.cat.incl, data = matched.test.incl.pData, conf.type =
                        "log-log")

summary(km.OS.test)
#sum <- data.frame(summary(km.OS.test))

#rite.table(summary(km.OS.test), file="kaplan meier results for OS.txt", sep="\t")