for (cur in 1:length(curatives)){
  cox.relapse.incl <- coxph (Surv(curatives[[cur]]$EFS, EFS_binaries[[cur]]) ~ genesofinterest[[cur]])
  assign(paste0("cox_relapse_",named_EFS_binaries[[cur]]), cox.relapse.incl)
}

library(foreach)
library(doParallel)
registerDoParallel(10) #####
all.cox.relapse.incl <- foreach (cur = 1:length(curatives))%dopar%{
  cox.relapse.incl <- coxph (Surv(curatives[[cur]]$EFS, EFS_binaries[[cur]]) ~ genesofinterest[[cur]])
  assign(paste0("cox_relapse_",named_EFS_binaries[[cur]]), cox.relapse.incl)
  return(cox.relapse.incl)
}



head(mb.vsd)

goi <- "ENSG00000136997.15_1"

apply(mb.vsd,1,function(x){cor(as.numeric(x),as.numeric(mb.vsd[goi,]))}) -> res


library(foreach)
library(doParallel)
registerDoParallel(10)
res <- foreach(i = 1:nrow(mb.vsd), .combine = rbind)%dopar%{
#res <- foreach(i = 1:5, .combine =rbind)%dopar%{
  as.numeric(mb.vsd[goi,]) -> y
  as.numeric(mb.vsd[i,]) -> x
  cor.test(x,y) -> cor.res
  cor(x,y) -> coef
  return(c(coef, cor.res$p.value))
}
res <- as.data.frame(res)
cbind(res,p.adjust(res[,2], method = "BH")) -> res
c("coef", "p.value", "p.adj") -> colnames(res)
rownames(res) <- rownames(mb.vsd)
#rownames(res) <- rownames(mb.vsd)[1:5]

res[order(res$p.value),] -> ordered.res
hist(res$p.value)
hist(res$p.adj)
