load("pData_saved.RData")
#gene_list <- read.csv("")

gene_list <- as.list("ENSG00000136997")
for (gene in 1:length(gene_list)){
#### create an ordered list of variables measured in the matched data 
matched.test.pData <- test.pData[index[!is.na(index)],] 
  
variable_list <- list()
variable_list <- append(variable_list, names(matched.test.pData))
#### Exlude the nmb column from the list
#variable_list <- variable_list[!variable_list =="NMB"]

#### itterate through the variables list generated and complete a chi squared test for each item in the list 
for(v in 2:length(variable_list)){
  pdf(paste("Heatmap of chi.sq ",variable_list[v],".pdf"))
  #### assign the results of the chi squared test to res, [, v] defines a column number 1:
  res <- chi.sq(matched.test.pData[, v], matched.goi.vsd.cat)
  #### now rename res according to the name of the item in the variable list being assessed
  assign(paste0("chi_res_",variable_list[[v]]), res)
  dev.off()
}
#### create a list containing all of the results of the chi.squared test
chi.sq.results <- as.list(mget(ls(pattern="chi_res_")))
#### determine the results in the list that are significant 
significant_chi_results <- list()
stats_chi_results <- list()
for (c in 1:length(chi.sq.results)){
  #### first create a individual data frame (indf) for each of the results in the list ( makes it easier to work with)
  indf <- chi.sq.results[[c]]
  #assign(paste(names(chi.sq.results[c])), indf)
  #### Now write a logical for the 3rd item in the list ([[3]] is table.temp which contains the p.vale) and specifcy we want the p.value less tha 0.05 from it
  sig <- indf[which(indf[[1]]$p.value < 0.05)]
  #stats <- indf[[1]]
  #names(sig) <- paste(names(chi.sq.results[c]))
  #### now add all the significant results identified into a significant results list
  significant_chi_results <- append(significant_chi_results, sig)
  #stats_chi_results <- append(stats_chi_results, stats)
}

######
#### Because we have a list of significant results each item within contains thet table data, we can now run fishers test on each table in the list 
#### First create a list to hold the results of the anlyses 
fishers_results <- list()
#### Now itterate through each of the items in the significant  results list 
for (s in 1:length(significant_chi_results)){
  #### create the individual tables
  tab <- significant_chi_results[[s]]
  #### run the fishers test on the extravcted table
  result <- fisher.test(tab)
  ####
  fishers_results <- append(fishers_results, result)
}

### Correlation coefficients
#### run the tests you want to run on the continuous age data 
x <- matched.test.pData$age.cont
y <- matched.goi.vsd
age_results.cor <- cor.result(x,y)
age_results.lin.reg <- lin.reg(x,y)
age_results.wilcox <- wilcox.test(x,y)


##################################
### logistic regression

cat ("processing logistic regression for each variable", sep ="\n")
#### create a vector of test factors
test_factors <- as.vector(names(matched.test.pData))
#### generate an empty list to hold the factors 
fac <- list()
#### itterate through the list of test factorss and generate the heading names add them into the fac list
for(i in 1:length(test_factors)){
  fac <- append(fac, as.name(paste0("matched.test.pData$",test_factors[[i]], sep="")))
}
#### itterate through each of the factors in the list s and run a lofgistic regression on each 
for(t in 3:length(fac)){
  regression <- logisticRegression(matched.test.pData[, t], matched.goi.vsd, matched.test.pData)
  #### renmae the output variable from each run through the loop so that the factor tested is attached to the object 
  assign(paste0("log_reg_",test_factors[[t]]), regression)
}
#### create a list of all the logistic regression outputs based on the prefix assigned to the objects
reg.log.list <- as.list(mget(ls(pattern="log_reg_")))
#print(reg.log.list)

for(t in 3:length(fac)){
  pairwise <- pairwise.t.test(matched.goi.vsd, matched.test.pData[, t])
  #### renmae the output variable from each run through the loop so that the factor tested is attached to the object 
  assign(paste0("pairwise_t_test",test_factors[[t]]), regression)
}

pairwise_t_tests <- as.list(mget(ls(pattern="pairwise_t_test")))

#### stratify samples by age group, gender and subgroup
strat_infants <- test.pData[which(test.pData$agefac == "Infant"),]
strat_Juniors <- test.pData[which(test.pData$agefac == "Junior"),]
strat_Teenagers <- test.pData[which(test.pData$agefac == "Teenager"),]
strat_Males <- test.pData[which(test.pData$sex == "male"),]
strat_Females <- test.pData[which(test.pData$sex == "female"),]
strat_Male_Juniors <- test.pData[which(test.pData$sex == "male" & test.pData$agefac == "Junior"),]
strat_Female_Juniors <- test.pData[which(test.pData$sex == "female" & test.pData$agefac == "Junior"),]
strat_4subgroup3 <- test.pData[which(test.pData$subgroup4fac == "G3"),]
strat_4subgroup4 <- test.pData[which(test.pData$subgroup4fac == "G4"),]
strat_4subgroupWNT <- test.pData[which(test.pData$subgroup4fac == "WNT"),]
strat_4subgroupSHH <- test.pData[which(test.pData$subgroup4fac == "SHH"),]
strat_7subgroup4HR <- test.pData[which(test.pData$subgroup7fac == "4HR"),]
strat_7subgroup4LR <- test.pData[which(test.pData$subgroup7fac == "4LR"),]
strat_7subgroup3HR <- test.pData[which(test.pData$subgroup7fac == "3HR"),]
strat_7subgroup3LR <- test.pData[which(test.pData$subgroup7fac == "3HR"),]
strat_7subgroupWNT <- test.pData[which(test.pData$subgroup7fac == "WNT"),]
strat_7subgroupSHH <- test.pData[which(test.pData$subgroup7fac == "SHH"),]


stratification_groups <- as.list(mget(ls(pattern="strat_")))

strats <- names(stratification_groups)

#for (c in 2:ncol(test.pData)){
#  index <- match(names(goi.vsd), rownames(test.pData)) 
#  for (s in 1:length(stratification_groups)){
#  summary(stratification_groups[[s]])
#  matched.test.pData <- stratification_groups[[s]][index[!is.na(index)],] 
#  assign(paste0("groups_stratified_test.pData_",strats[[s]]), matched.test.pData)
#  #matched.test.pData <- test.pData[index[!is.na(index)],] 
#  matched.goi.vsd <- goi.vsd[!is.na(index)] 
#  matched.goi.vsd.cat <- ifelse(matched.goi.vsd>median(goi.vsd, na.rm = T), "high","low") 
#}
#}

all_stratif_groups <- as.list(mget(ls(pattern="groups_stratified_test.pData_")))

#for (grp in 1:length(all_stratif_groups)){
#  index.incl <- match(names(goi.vsd), rownames(all_stratif_groups[[grp]])) 
#  matched.incl.pData <- all_stratif_groups[[s]][index.incl[!is.na(index.incl)],] 
  
#}

cat ("creating combined dataframe to assess biomarker in G3 G4 combined group, for survival cohort, aged 3-16 years, curative intent", sep = "\n")
Treated_all_groups <- test.pData[which(test.pData$curative == "curative" & test.pData$childfac == "Child.M" | test.pData$childfac == "Child.F"),]

Treated_G3G4 <- test.pData[which(test.pData$curative == "curative" & test.pData$childfac == "Child.M" | test.pData$childfac == "Child.F" & test.pData$subgroup4fac == "G3" | test.pData$subgroup4fac == "G4"),]
#### Creating matched data frames containing RNAseq expression data and curative data for samples in test.pData
treat_grps <- as.list(mget(ls(pattern="Treated_")))
names_treat_grps <- names(treat_grps)

for (i in 1:length(treat_grps)){
index.incl <- match(names(goi.vsd), rownames(treat_grps[[i]])) 
matched.test.incl.pData <- treat_grps[[i]][index.incl[!is.na(index.incl)],]
assign(paste0("matched_curative_",names_treat_grps[[i]]),matched.test.incl.pData)
matched.goi.vsd.incl <- goi.vsd[!is.na(index.incl)] 
assign(paste0("matched_goi_curative_",names_treat_grps[i]),matched.goi.vsd.incl)
matched.goi.vsd.cat.incl <- ifelse(matched.goi.vsd.incl>median(goi.vsd, na.rm = T), "high","low")
assign(paste0("matched_goi_cat_curative_",names_treat_grps[i]),matched.goi.vsd.cat.incl)
}

#### create the lists for downstream processing 
curatives <- as.list(mget(ls(pattern="matched_curative_")))
names_curatives <- names(curatives)
genesofinterest <- as.list(mget(ls(pattern="matched_goi_cat_curative")))

#### change the name list to use otherwise re-run will cause the list to be regenrated containing binary data 
### creating binary relapse variables labelled 0,1 for event analysis
for (c in 1:length(curatives)){
  relapse_binary <- ifelse(curatives[[c]]$relapse == "relapse", 1, 0)
  assign(paste0("binary_relapse_",names_curatives[[c]]),relapse_binary)
  Overall_survival_binary <- ifelse(curatives[[c]]$OS.cat == "Dead", 1, 0)
  assign(paste0("binary_OS_",names_curatives[[c]]), Overall_survival_binary)
  EventFreeSurvival_binary <- ifelse(curatives[[c]]$Event == "Event", 1, 0)
  assign(paste0("binary_EFS_",names_curatives[[c]]), EventFreeSurvival_binary)
}

#### create the lists for downstream processing 
binaries <- as.list(mget(ls(pattern="binary_")))
named_binaries <- names(binaries)
#### Run Kaplan-Meier estimates for curative data using the binary survical data for the gene of interest 
for (bin in 1:length(binaries)){
  for (cur in 1:length(curatives)){
    for (goi in 1:length(genesofinterest)){
      km_results_log.test<- km.log.test(time = curatives[[cur]]$PFS, event = binaries[[bin]], marker = genesofinterest[[goi]])
      assign(paste0("Kaplan_Meier",names_curatives[[cur]]), km_results_log.test)
    }
  }
}