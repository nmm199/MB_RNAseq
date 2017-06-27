# library(readr)
library(NMF)

#### basic fundtion providing a putative biomarker
mb.vsd <- read.delim(file = "/home/dan/mygit/rna_seq_mb/paper/MB.vsd.txt")
goi <- "ENSG00000136997"
as.numeric(mb.vsd[goi,]) -> goi.vsd
gsub("T","",names(mb.vsd)) -> names(goi.vsd)
### the output would be a vector with a continuous variable names equal NMB numbers

###################################################

#### read in database
pData <- read.csv("~/R/MB_Data/Input data/database020617.csv", row.names = 1)

#### convert pdata object into sig.test variables
### age 
pData$Age_R -> age.cont
pData$Age_R<3 -> age.cat.infant
pData$Age_R>16 -> age.cat.adult.16
pData$Age_R>21 -> age.cat.adult.21
!age.cat.infant&!age.cat.adult.16 -> age.filter

### metastatic
pData$Mstatus_R -> mstatus
#### I wonder if there are more room for increased categories

### relapse

pData$PFS_R -> relapse
ifelse(relapse==1, "relapse","non-relapse") -> relapse

#### more of the same please continue....


test.pData <- data.frame(age.cont,
                    age.cat.infant,
                    age.cat.adult.16,
                    age.cat.adult.21,
                    age.filter)
rownames(test.pData) <- rownames(pData)
#.............. please continue)
save(test.pData, file = ".....")                    

#############################################

### goi.vsd #### name of a putative biomarker

match(names(goi.vsd), rownames(test.pData)) -> index
test.pData[index[!is.na(index)],] -> matched.test.pData
goi.vsd[!is.na(index)] -> matched.goi.vsd
ifelse(matched.goi.vsd>median(goi.vsd, na.rm = T), "high","low") -> matched.goi.vsd.cat
#### test that you want to do

### name each test
sink(file "text file")
table(matched.test.pData$age.cat.infant, matched.goi.vsd.cat) -> table.age.cat
prop.table(table.age)*100 -> table.age.cat.perc
summary.table(table.age.cat)
chisq.test(table.age) -> chi.test.age.cat
c(stat=chi.test.age.cat$statistic, p.value=chi.test.age.cat$p.value) -> chi.test.age.cat.res
aheatmap(chi.test.age.cat$residuals, Rowv=NA, Colv = NA, xlab = "Expression")
unsink()

#### continu constructng test




### assessing the data, parametric etc
# qqnorm(data$Age_R)
plot(pdata$Age_R, col='red')

age<-data$Age_R
meth<-data$X450K_R

glm(formula = meth ~ age, family = binomial(link = "logit"), 
    data = data)

### Getting summary characteristics
summary(age)
mean(age, na.rm = T)
median(age, na.rm = T)

### note that s.d. is only used with mean, with normally distributed symmetric data
### question 1: why do the following not work? 
# var(data)
# sd(age)
# range(age)
# head(age)
# min(age)
# max(age)

t.test(age, conf.level=0.95)
boxplot(age~meth, col=c("yellow","green","red","blue"))

### other examples to t.test
# t.test(x~y, data=data, var.equal=T)

### install.packages('gplots')
library('gplots')
plotmeans(age~meth, data=data)


# data <- read.table(file= "/home/nmm199/test/data/Test_data.txt",  header=T)
# data <- read.table(file= "/home/nmm199/test/data/Test2.txt",  header=T, sep = "\t")
### read.table coding - optimise this
# read.table(file="/home/nmm199/test/data/Test2.txt", header = TRUE, sep="", quote="'", dec=".")


MYC <- data$MYC_R
PFS<- data$PFS_R
table(as.factor(MYC),as.factor(PFS))
chisq.test(table(as.factor(MYC),as.factor(PFS)))

tab <- table(MYC, PFS)

chisq.test(tab, correct=F)
# chisq.test(tab, correct=F)
fisher.test(tab)

#install.packages('survival')
library('survival')

km<- survfit(Surv(sweeks, event)~treat)
MB <- coxph (Surv(sweeks, event)~treat)
# summary (MB)


