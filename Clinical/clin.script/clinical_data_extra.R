#######################

### File information: this is unused script taken from clinical data script (clinical_data_3.R)
### File aim: 
### 1. Provide back up information and code in case this is required for further analysis
### 2. Files that were used to cross-reference RNA seq cohort (NMB numbers)
### 3. Files that were used to run the cytogenetic arm change analysis

### Author: Dr Marion Mateos
### Date: 3/7/17




### input files:
file = "/home/nmm199/R/MB_RNAseq/Input data/arm_calls_clean280617.txt"




### Functions:

### frequencies of cytogenetic arm changes

cytogen.freq <- function (x){
  table.temp <- table (x)
  return (table.temp)
}

# x <- chrosomomal arm
### returns frequency of gain, loss, neutral for each chromosomal arm included in Ed Schwalbe's arm analysis



########################################
### matching NMB numbers to determine which included in RNA seq for reporting purposes

# test.pData$RNAdata <- mb.vsd[match(rownames(test.pData), colnames(mb.vsd)), ]

### can save it as a vector instead
RNAseq <- colnames (mb.vsd)
RNAseq

### save as readable text file, which can then be uploaded to computer harddrive and matched to existing clinical database in excel
write.table(RNAseq, "/home/nmm199/R/MB_RNAseq/Clinical/RNAseq.txt", sep="\t")

#########################################


### Distribution of the data


### note that s.d. is only used with mean, with normally distributed symmetric data
# var(data)
# sd(age)
# range(age)
# head(age)
# min(age)
# max(age)

##################################################

### other ways of constructing chi squared table and chisq test

# table(as.factor(MYC.cat),as.factor(relapse))
# chisq.test(table(as.factor(MYC.cat),as.factor(relapse)))
# tab.MYC.relapse <- table(MYC.cat, relapse)



######################################################

### pearson's correlation coefficient
# cor.test(x, y)
# plot(x,y)
# abline (lm(x ~y)



### if relationship does not appear linear then use spearman correlation coefficient

# cor.test(x, y, method = 'spearman')


#########################################


### linear regression
# this may not be required, could use correlation above

age.reg <- lm(matched.test.pData$age.cont ~ matched.goi.vsd)
age.reg
summary(age.reg)



#########################################



### Mann-whitney U (aka wilcoxon rank sum test) for non-parametric data
### note that if the sample size is large enough (e.g > 30) then could theoretically use t-test, even if the data is non-parametric
### however we will use the Wilcoxon rank sum test


### note that for the output object("output.wilcox"), can derive p value by output.wilcox$p.value   e.g 
# age.cont.wilcox <- wilcox.test(matched.test.pData$age.cont ~ matched.goi.vsd.cat, exact = F, correct = F)
# age.cont.wilcox$p.value

### comparison to t-test for understanding data
# t.test(matched.test.pData$age.cont ~ matched.goi.vsd.cat)

#######################################

### AOV 
### AOV depends on Levene's test for homogeneity of variances, H0= variances are equal, HA= variances are not equal. Therefore wish for H0 to be true ie cannot reject H0 , i.e variances are equal (homogenous)
### check variances are equal (i.e p>0.05, to confirm assumption upon which one-way ANOVA is created

cat ("confirm that variances are equal, ie. leveneTest is not significant", sep ="\n")

leveneTest(matched.goi.vsd ~ matched.test.pData$meth7)

# plot(matched.goi.vsd, matched.test.pData$meth)

# qqplot(matched.goi.vsd, matched.test.pData$meth7)


###########################
### boxplots

meth.boxplot <- boxplot(matched.goi.vsd~matched.test.pData$meth, col=c("yellow","green","red","blue"), xlab = "Methylation subgroup", ylab = "Biomarker expression", main = "Correlation between biomarker and 4 molecular subgroups")

meth7.boxplot <- boxplot(matched.goi.vsd~matched.test.pData$meth7, col=c("yellow","green","red","blue"),  xlab = "Methylation subgroup", ylab = "Biomarker expression", main = "Correlation between biomarker and 7 molecular subgroups")

### alternative coding, creates same boxplots as above

# plot(matched.test.pData$meth, matched.goi.vsd, col = c("Yellow", "Green", "red","blue"),  xlab = "methylation groups", ylab= "Expression of biomarker", main = "Correlation with methylation and biomarker")



### post-hoc tests for methylation

### pairwise t-test to determine where the difference lies between the groups
meth.pw <- pairwise.t.test(matched.goi.vsd, matched.test.pData$meth)

meth.pw


# can specify bonferroni correction, which may be over-conservative
# meth.pw.bon <- pairwise.t.test(matched.goi.vsd, matched.test.pData$meth, p.adj='bonf')
# meth.pw.bon


###########################

### logistic regression alternate coding


log.reg.mstatus <- glm(matched.test.pData$mstatus ~ matched.goi.vsd,  family = binomial(link= 'logit'), data=matched.test.pData)
summary(log.reg.mstatus)

mstatus.boxplot <- boxplot(matched.goi.vsd ~ matched.test.pData$mstatus, col = c("red", "blue"), xlab = "M status", ylab = "Biomarker expression", main = "Correlation between biomarker and metastatic status")

### this also works for metastatic status
# log.reg.mstatus2 <- glm(mstatus ~ matched.goi.vsd,  family = binomial(link= 'logit'), data=matched.test.pData)
# summary(log.reg.mstatus2)



########

### assessing the frequency for each arm changes (e.g p1 is short arm on chromosome 1)

cytogen.p1 <- cytogen.freq(x = cytogen$p1)
cytogen.p2 <- cytogen.freq(x = cytogen$p2)
cytogen.p3 <- cytogen.freq(x = cytogen$p3)
cytogen.p4 <- cytogen.freq(x= cytogen$p4)
cytogen.p5 <- cytogen.freq(x= cytogen$p5)
cytogen.p6 <- cytogen.freq(x= cytogen$p6)
cytogen.p7 <- cytogen.freq(x= cytogen$p7)
cytogen.p8 <- cytogen.freq(x= cytogen$p8)
cytogen.p9 <- cytogen.freq(x= cytogen$p9)
cytogen.p10 <- cytogen.freq(x= cytogen$p10)
cytogen.p11 <- cytogen.freq(x= cytogen$p11)
cytogen.p12 <- cytogen.freq(x= cytogen$p12)
cytogen.p16 <- cytogen.freq(x= cytogen$p16)
cytogen.p17 <- cytogen.freq(x= cytogen$p17)
cytogen.p18 <- cytogen.freq(x= cytogen$p18)
cytogen.p19 <- cytogen.freq(x= cytogen$p19)
cytogen.p20 <- cytogen.freq(x= cytogen$p20)
cytogen.p21 <- cytogen.freq(x= cytogen$p21)


cytogen.q1 <- cytogen.freq(x = cytogen$q1)
cytogen.q2 <- cytogen.freq(x = cytogen$q2)
cytogen.q3 <- cytogen.freq(x = cytogen$q3)
cytogen.q4 <- cytogen.freq(x= cytogen$q4)
cytogen.q5 <- cytogen.freq(x= cytogen$q5)
cytogen.q6 <- cytogen.freq(x= cytogen$q6)
cytogen.q7 <- cytogen.freq(x= cytogen$q7)
cytogen.q8 <- cytogen.freq(x= cytogen$q8)
cytogen.q9 <- cytogen.freq(x= cytogen$q9)
cytogen.q10 <- cytogen.freq(x= cytogen$q10)
cytogen.q11 <- cytogen.freq(x= cytogen$q11)
cytogen.q12 <- cytogen.freq(x= cytogen$q12)
cytogen.q13 <- cytogen.freq(x= cytogen$q13)
cytogen.q14 <- cytogen.freq(x= cytogen$q14)
cytogen.q15 <- cytogen.freq(x= cytogen$q15)
cytogen.q16 <- cytogen.freq(x= cytogen$q16)
cytogen.q17 <- cytogen.freq(x= cytogen$q17)
cytogen.q18 <- cytogen.freq(x= cytogen$q18)
cytogen.q19 <- cytogen.freq(x= cytogen$q19)
cytogen.q20 <- cytogen.freq(x= cytogen$q20)
cytogen.q21 <- cytogen.freq(x= cytogen$q21)
cytogen.q22 <- cytogen.freq(x= cytogen$q22)


# cytogen.df <- as.data.frame(cytogen.p1, 
#                         cytogen.p2, 
#                         cytogen.p3,
#                         cytogen.p4, 
#                         cytogen.p5,
#                         cytogen.p6,
#                         cytogen.p7,
#                         cytogen.p8, 
#                         cytogen.p9, 
#                         cytogen.p10,
#                         cytogen.p11,
#                         cytogen.p12,
#                         cytogen.p16,
#                         cytogen.p17, 
#                         cytogen.p18,
#                         cytogen.p19, 
#                         cytogen.p20,
#                         cytogen.p21,
#                         cytogen.q1, 
#                         cytogen.q2, 
#                         cytogen.q3, 
#                         cytogen.q4, 
#                         cytogen.q5, 
#                         cytogen.q6, 
#                         cytogen.q7)

# rownames(cytogen.df) <- colnames(cytogen)                        

# nrow(cytogen)




##############################################


# G4 <- matched.test.pData$meth =="G4"
# G4.group <- matched.test.pData [G4, ]
# summary(G4.group)

### confirmed that the above is valid through subsetting the object and defining removal of na
# str(G4.group)
# median(G4.group$age.cont, na.rm=T)
# summary(G4.group$age.cont)
# View(G4.group)



###############################################

### old dataframes that did not incorporate curative therapy specifically therefore superceded

cat ("restrict analysis to age 3-16 years", sep = "\n")
cat ("create dataframe for age 3-16 years")

Age.incl <-matched.test.pData$age.filter =="TRUE"
Age.incl.df <- matched.test.pData [Age.incl, ]
summary(Age.incl.df)
View(Age.incl.df)




#####################################################

### exp(coef) is the hazard ratio, with 95% CI. p value is listed above against the variable
# help ("coxph")

##############################################


### checking that variables are binary before running km survival
### check alive status is binary
matched.test.pData$OS.cat

# OS.cat.bin <- ifelse(matched.test.pData$OS.cat == "Dead", 1, 0)
# OS.cat.bin

##############################################

### old code for survival analysis that was not performed on curative age 3-16 yo cohort therefore defunct




#cox.OS <- coxph (Surv(matched.test.pData$Followup, OS.cat.bin) ~ matched.goi.vsd.cat)
#summary(cox.OS)$logtest
#summary(cox.OS)




### old code prior to creating function

cox.OS.incl <- coxph (Surv(matched.test.incl.pData$Followup, OS.cat.bin.incl) ~ matched.goi.vsd.cat.incl)
summary(cox.OS.incl)$logtest
summary(cox.OS.incl)



#################################################


### old PFS code


cat ("evaluate effect of biomarker on PFS", sep = "\n")

### change relapse to binary
relapse.bin <- ifelse(matched.test.pData$relapse == "relapse", 1, 0)

# matched.test.incl.pData$PFS -> time
# relapse.bin.incl -> event
# matched.goi.vsd.cat.incl -> marker


cat("cox regression analysis, PFS and biomarker", sep ="\n")

cox.relapse <- coxph (Surv(matched.test.pData$PFS, relapse.bin) ~ matched.goi.vsd.cat)
summary(cox.relapse)$logtest
summary(cox.relapse)


### 




### old EFS code

### change event-free status to binary
matched.test.pData$Event

EFS.cat.bin <- ifelse(matched.test.pData$Event == "Event", 1, 0)
EFS.cat.bin

chi.sq(as.factor(EFS.cat.bin), as.factor(matched.goi.vsd.cat))
### creating a table as km curve suggests that all in the biomarker high group had event. Ask Dan

# table.EFS<- table(as.factor(EFS.cat.bin), as.factor(matched.goi.vsd.cat))
# table.EFS.perc <- prop.table(table.EFS)*100

# km.EFS <- survfit(Surv(matched.test.pData$EFS, EFS.cat.bin)~matched.goi.vsd.cat)

# km.EFS <- survfit(Surv(matched.test.pData$EFS, EFS.cat.bin==1)~matched.goi.vsd.cat) is the same as when not specifying status==1
# km.EFS

### troubleshooting, make sure that event is binary and that binary event is included as "Event" in the km.EFS graph

#km.EFS.incl <- survfit(Surv(matched.test.incl.pData$EFS, EFS.cat.bin.incl)~matched.goi.vsd.cat.incl)
#km.EFS.incl



plot(km.EFS, col = c("red", "blue"), xlab = "time to event (years)", ylab = "EFS", main = "Biomarker expression and event-free survival (EFS)",  lty = 1:2)
EFS.names <- c("biomarker - high", "biomarker - low")
legend (15,1, EFS.names,  lty= 1:2)
# temp.logrank.pval <- 0.262
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))

EFS.logrank <- survdiff(Surv(matched.test.pData$EFS, EFS.cat.bin) ~ matched.goi.vsd.cat)
EFS.logrank
summary(EFS.logrank)

message ("will need to manually adjust the graph based on the log-rank p value")

plot(km.EFS, col = c("red", "blue"), xlab = "time to event (years)", ylab = "PFS", main = "Biomarker expression and event-free survival (EFS)",  lty = 1:2)
EFS.names <- c("biomarker - high", "biomarker - low")
legend (15,1, EFS.names,  lty= 1:2)
temp.logrank.pval <- 0.262
text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))


message ("restrict survival analysis for age 3-16 years")



plot(km.EFS.incl, col = c("red", "blue"), xlab = "event-free survival (years)", ylab = "EFS", main = "Biomarker expression and event-free survival (EFS)",  lty = 1:2)
EFS.names.incl <- c("biomarker - high", "biomarker - low")
legend (15,1, EFS.names.incl,  lty= 1:2)
# temp.logrank.pval <- 0.67
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))

EFS.incl.logrank <- survdiff(Surv(matched.test.incl.pData$EFS, OS.cat.bin.incl) ~ matched.goi.vsd.cat.incl)
EFS.incl.logrank

### now go back in and add p value from EFS.incl.logrank into graph
message ("need to add p value from EFS.incl.logrank into graph")

plot(km.EFS.incl, col = c("red", "blue"), xlab = "event-free survival (years)", ylab = "EFS", main = "Biomarker expression and event-free survival (EFS)",  lty = 1:2)
EFS.names.incl <- c("biomarker - high", "biomarker - low")
legend (15,1, EFS.names.incl,  lty= 1:2)
temp.logrank.pval <- 0.543
text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))



### 
message("cox regression for EFS on entire cohort")

cox.EFS <- coxph (Surv(matched.test.pData$EFS, EFS.cat.bin) ~ matched.goi.vsd.cat)
summary(cox.EFS)$logtest
summary(cox.EFS)


################################################
### visualisation options to explore

# plot(ks.test(x[y],x[!y]))

# plot(ecdf(x[y]), col = "red")

# lines(ecdf(x[!y]), col = "blue")

# plot(density(x[y], na.rm = T), col = "red")
# lines(density(x[!y], na.rm = T), col = "blue")
# rug(x[y], col = "red")
# rug(x[!y], col = "blue")
# temp.pval <- 0.01
# text(0.06, 0.06, paste("P.val  = ", temp.pval))
# abline(v = median(x[y], na.rm = T), lty = 2, col = "red")
# abline(v = median(x[!y], na.rm = T), lty = 2, col = "blue")





### redundant or non functional code

##############################
### troubleshooting tips

### example codes below

# pData <- read.table(file= "/home/nmm199/test/data/Test_data.txt",  header=T)
# pData <- read.table(file= "/home/nmm199/test/data/Test2.txt",  header=T, sep = "\t")
# read.table(file="/home/nmm199/test/data/Test2.txt", header = TRUE, sep="", quote="'", dec=".")



### if unsure of how function works, can try specifying the variables, and explore the temporary output.

# x = matched.test.pData$age.cat.infant
# y = matched.goi.vsd.cat
# chi.sq <- function(x,y){
#  table.temp <- table(x, y) 
#  table.temp.perc <- prop.table(table.temp)*100
#  summary.table(table.temp)
#  chi.test.temp <- chisq.test(table.temp) 
#  chi.test.temp.stat <- c(stat=chi.test.temp$statistic, p.value=chi.test.temp$p.value) 
#  chi.test.temp.res <- chi.test.temp$residuals
#  aheatmap(chi.test.temp.res, Rowv=NA, Colv = NA)
# list.temp <- list  (table.temp, 
#                     table.temp.perc,
#                     chi.test.temp.stat,
#                     chi.test.temp.res
#                      )

#  return(list.temp)
#}

# list.temp



### exploring the commands within the chi square function

# table.age.cat <- table(matched.test.pData$age.cat.infant, matched.goi.vsd.cat) 
# table.age.cat.perc <- prop.table(table.age.cat)*100
# summary.table(table.age.cat)
# chi.test.age.cat <- chisq.test(table.age.cat) 
# chi.test.age.cat.res <- c(stat=chi.test.age.cat$statistic, p.value=chi.test.age.cat$p.value) 



### relabelling heatmap if needed for output

# aheatmap(chi.test.age.cat$residuals, Rowv=NA, Colv = NA, labRow = "Infant", labCol="Expression", main="Correlations")
# help(aheatmap)

# aheatmap(chi.test.age.cat$residuals, Rowv=NA, Colv = NA, main = "Expression correlated to Infant status") 


### generate list, without the function above
# list.age <- list(age= table.age.cat, 
#             perc=table.age.cat.perc,
#            residual=chi.test.age.cat.res,
#             chi=chi.test.age.cat
#             )
# list.age

### subsetting list results for equivalent slice of object

# list.age[[1]]
# list.age$age



### other examples for t.test
# t.test(x~y, data=data, var.equal=T)
# chisq.test(tab, correct=F)
# fisher.test(tab)
# glm(formula = meth ~ age, family = binomial(link = "logit"), 
#  data = data)


########################################

### record of code that does not work

### error in trying to use age.cont with categorical biomarker or is use biomarker is dependent and age (categorical) as independent factor
# log.reg.age.cat2<- glm(matched.goi.vsd ~ matched.test.pData$age.cat.infant, family = binomial(link='logit'), data=matched.test.pData)

# log.reg.age.cont <- glm(age.cont ~ matched.goi.vsd.cat, family = binomial(link= 'logit'), data=matched.test.pData)
# matched.goi.vsd.cat.bin <- ifelse(matched.goi.vsd.cat=="high", 1, 0)
# matched.goi.vsd.cat.bin.2 <- matched.goi.vsd.cat.bin, na.rm=T
# also na.action = omit as another option to incorporate

# log.reg.age.cont <- glm (formula = matched.test.pData$age.cont ~ matched.goi.vsd.cat.bin , family = binomial(link='logit'), data=matched.test.pData)
# matched.goi.vsd.cat.bin <- as.factor(matched.goi.vsd.cat.bin)

# is.factor(matched.goi.vsd.cat.bin)


# log.reg <- function(x,y){
# glm(formula = x~y, family = binomial(link = 'logit'), data = matched.test.pData)
# summary.log.temp <- summary(log.reg.temp)
# aov.temp <- aov(x~y)
# sum.aov.temp <- summary(aov.temp) 
# list.log.reg <- list (summary.log.temp, 
#                       anova.temp,
#                       aov.temp)
# return(list.log.reg)
#}

# x <- age.cat.infant.bin
# y <- matched.goi.vsd

# age.cat.infant.bin<- ifelse(matched.test.pData$age.cat.infant, TRUE==1, 0)
# log.reg.age <- log.reg(x = matched.test.pData$age.cat.infant, y = matched.goi.vsd)
# log.reg.age.bin <- log.reg (x = age.cat.infant.bin, y = matched.goi.vsd)

### attempt with Dan below for records
# x <- matched.test.pData$age.cont
# y <- matched.goi.vsd.cat=="high"


# cat.variable.cont <- function(x, y){
# data <- data.frame(x,y)
#log.reg <- function(x,y){
#glm(formula = x~y, family = binomial(link = "logit"), data = data)
#anova(model, test="Chisq")

#aov(x~y) -> aov.out
#str(summary(aov.out))
#summary(aov.out) -> sum.aov.out
#sum.aov.out


# subgroup.meth<- function (x,y){
# temp.sub <- x =="temp.group"
# temp.group <- matched.test.pData [temp.sub, ]
# summary(temp.group)
# return(summary(temp.group))
# }

# G3.subgroup <- subgroup.meth(x = matched.test.pData$meth )
# G3.subgroup

### cox ph

### i think this code below is redundant 

# cox.relapse.df <- data.frame(matched.goi.vsd.cat, matched.test.pData$relapse, matched.test.pData$PFS)
# cox.relapse.df <- data.frame(matched.goi.vsd.cat, relapse.bin, matched.test.pData$PFS)
# surv.cox.rel <- survfit (cox.relapse, newdata=cox.relapse.df)
# summary(surv.cox.rel)                

######################

# message ("adding in cytogenetic 13q data") - superceded because 13q loss data was added earlier on to test.pData

# index.incl.cytogen <- match(rownames (Age.incl.df),row.names(cytogen.q13.cat [,1]))
# matched.test.incl.pData <- Age.incl.df[index.incl.cytogen[!is.na(index.incl.cytogen)],]



#######################
### redundant PFS code 



# km.PFS <- survfit(Surv(matched.test.pData$PFS, relapse.bin)~matched.goi.vsd.cat)
# km.PFS

# plot(km.PFS, col = c("red", "blue"), xlab = "time to progression/relapse (years)", ylab = "PFS", main = "Biomarker expression and progression-free survival (PFS)",  lty = 1:2)
# PFS.names <- c("biomarker - high", "biomarker - low")
# legend (12,1, PFS.names,  lty= 1:2)
# temp.logrank.pval <- 0.494
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))

#PFS.logrank <- survdiff(Surv(matched.test.pData$PFS, relapse.bin) ~ matched.goi.vsd.cat)
# PFS.logrank

# message ("will need to manually adjust the graph based on the log-rank p value")

# plot(km.PFS, col = c("red", "blue"), xlab = "time to progression/relapse (years)", ylab = "PFS", main = "Biomarker expression and progression-free survival (PFS)",  lty = 1:2)
# PFS.names <- c("biomarker - high", "biomarker - low")
# legend (15,1, PFS.names,  lty= 1:2)
# temp.logrank.pval <- 0.494
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))
### now go back in and add p value from PFS.incl.logrank into graph
# message ("need to add p value from PFS.incl.logrank into graph")

#plot(km.PFS.incl, col = c("red", "blue"), xlab = "time to progression/relapse (years)", ylab = "PFS", main = "Biomarker expression and progression-free survival (PFS)",  lty = 1:2)
#PFS.names <- c("biomarker - high", "biomarker - low")
#legend (15,1, PFS.names,  lty= 1:2)
# temp.logrank.pval <- 0.016
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))



###############################
### redundant OS code

# km.OS <- survfit(Surv(matched.test.pData$Followup, OS.cat.bin)~matched.goi.vsd.cat)
# km.OS

# plot(km.OS, col = c("red", "blue"), xlab = "time to progression/relapse (years)", ylab = "OS", main = "Biomarker expression and overall survival (OS)",  lty = 1:2)
# OS.names <- c("biomarker - high", "biomarker - low")
# legend (15,1, OS.names,  lty= 1:2)
# temp.logrank.pval <- 0.494
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))


# km.log.test(matched.test.incl.pData$Followup,OS.cat.bin, matched.goi.vsd.cat.incl )

# OS.logrank <- survdiff(Surv(matched.test.pData$Followup, OS.cat.bin) ~ matched.goi.vsd.cat)
# OS.logrank
# summary(OS.logrank)

# message ("will need to manually adjust the graph based on the log-rank p value")

# plot(km.OS, col = c("red", "blue"), xlab = "time to progression/relapse (years)", ylab = "PFS", main = "Biomarker expression and progression-free survival (PFS)",  lty = 1:2)
# OS.names <- c("biomarker - high", "biomarker - low")
#legend (15,1, PFS.names,  lty= 1:2)
#temp.logrank.pval <- 0.331
#text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))


#message ("restrict survival analysis for age 3-16 years")



#plot(km.OS.incl, col = c("red", "blue"), xlab = "overall survival (years)", ylab = "OS", main = "Biomarker expression and overall survival (PFS)",  lty = 1:2)
#OS.names.incl <- c("biomarker - high", "biomarker - low")
#legend (15,1, OS.names.incl,  lty= 1:2)
# temp.logrank.pval <- 0.67
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))

#OS.incl.logrank <- survdiff(Surv(matched.test.incl.pData$Followup, OS.cat.bin.incl) ~ matched.goi.vsd.cat.incl)
#OS.incl.logrank

### now go back in and add p value from OS.incl.logrank into graph
# message ("need to add p value from OS.incl.logrank into graph")

# plot(km.OS.incl, col = c("red", "blue"), xlab = "overall survival (years)", ylab = "OS", main = "Biomarker expression and overall survival (PFS)",  lty = 1:2)
# OS.names.incl <- c("biomarker - high", "biomarker - low")
# legend (15,1, OS.names.incl,  lty= 1:2)
# temp.logrank.pval <- 0.625
# text(11, 0.8, paste("Log-rank p value  = ", temp.logrank.pval))




### code from Dan below

matched.test.pData$PFS
matched.test.pData$relapse

View(matched.test.pData)
y <- matched.goi.vsd
y.cat <- matched.goi.vsd.cat
time <- matched.test.pData$PFS
event <-  matched.test.pData$relapse=="relapse"


km <- survfit(Surv(time, event)~y.cat)
temp.survdiff <- survdiff(Surv(time, event)~y.cat)
str(temp.survdiff)
plot(km, col = c("red","blue"), xlab = "time (years)", ylab = "PFS")
legend
summary(km)
summary(temp.survdiff)
cox.result <- coxph (Surv(time, event)~y.cat)
summary(cox.result)$logtest
summary(cox.result)


