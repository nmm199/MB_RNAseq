# library(readr)


### confirm wd
# setwd("/home/nmm199/R/MB_Data/Input data/")
library(readr)
data <- read.csv("~/R/MB_Data/Input data/database020617.csv")
# data <-read.txt('/home/nmm199/..../.txt')

### assessing the data, parametric etc
# qqnorm(data$Age_R)
plot(data$Age_R, col='red')

age<-data$Age_R
meth<-data$X450K_R

glm(formula = meth ~ age, family = binomial(link = "logit"), 
    data = data)

### Getting summary characteristics
summary(age)
# mean(age)
# median(age)

### note that s.d. is only used with mean, with normally distributed symmetric data
### question 1: why do the following not work? 
# var(data)
# sd(age)
# range(age)
# head(age)
# min(age)
# max(age)

t.test(age, conf.level=0.95)
boxplot(age~meth, col="orange")

### other examples to t.test
# t.test(x~y, data=data, var.equal=T)

### install.packages('gplots')
library('gplots')
plotmeans(age~meth, data=data)


# data <- read.table(file= "/home/nmm199/test/data/Test_data.txt",  header=T)
# data <- read.table(file= "/home/nmm199/test/data/Test2.txt",  header=T)
### read.table coding - optimise this
# read.table(file="/home/nmm199/test/data/Test2.txt", header = TRUE, sep="", quote="'", dec=".")


# MYC <- data$MYC_R
# PFS<- data$PFS_R
# chisq.test(MYC~PFS)

# tab<- as.table(cbind(MYC, PFS))

# chisq.test(tab, correct=F)
# chisq.test(tab, correct=F)
# fisher.test(tab)

install.packages('survival')
library('survival')

# km<- survfit(Surv(sweeks, event)~treat)

# MB <- coxph (Surv(sweeks, event)~treat)
# summary (MB)




