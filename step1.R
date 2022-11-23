rm(list = ls())  
options(stringsAsFactors = F)
f='GSE2880_eSet.Rdata'
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42872
library(GEOquery)
getwd()
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)

options("download.file.method.GEOquery"="libcurl")
gset <- getGEO('GSE2880', destdir=".",
               AnnotGPL = T,    
               getGPL = T)       
save(gset,file=f)   



load('GSE2880_eSet.Rdata')  
class(gset)  
length(gset)  
class(gset[[1]])
gset
a=gset[[1]] 
dat=exprs(a)
dim(dat)

dat[1:4,1:4] 
boxplot(dat,las=2)
pd=pData(a) 
