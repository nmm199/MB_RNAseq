save(ex.no, "ex_no.RData")
load("~/ex_no.RData")

saveRDS(ex.no, "ex_no.rds")
ex.no <- readRDS("~/ex_no.rds")