
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


######
library(stringr)
pd$title
group_list1=str_split(pd$title,',',simplify = T)[,2]
group_list2=str_split(group_list1,"_",simplify = T)[,1]
group_list3=gsub(": ","_",group_list2)
length(group_list3)
group_data_frame=data.frame(group_list3,1:62)
substr(group_data_frame$group_list3,1,nchar(group_data_frame$group_list3)-1)
??substr
g3=group_data_frame[c(8:14,36:42),]
g4=substr(g3$group_list3,1,nchar(g3$group_list3)-1)
dat3=dat[,c(8:14,36:42)]
View(dat3)
table(g4)

###
library(GEOquery)
gpl <- getGEO('GPL1355', destdir=".")
colnames(Table(gpl))  
head(Table(gpl)[,c(1,11)]) ## you need to check this , which column do you need
probe2gene=Table(gpl)[,c(1,11)]
head(probe2gene)
save(probe2gene,file = "probe2gene.Rdata")
load(file='probe2gene.Rdata')
ids=probe2gene 



library(stringr)  
colnames(ids)=c('probe_id','symbol')  
ids=ids[ids$symbol != '',]
ids=ids[ids$probe_id %in%  rownames(dat3),]
dat3[1:4,1:4]   
dat3=dat3[ids$probe_id,] 

###save the high expression sample from duplicated terms
ids$median=apply(dat3,1,median) 
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]
dat3=dat3[ids$probe_id,] 
rownames(dat3)=ids$symbol
dat3[1:4,1:4]  

save(dat3,g4,file = 'step2-output.Rdata')
