rm(list = ls())  
options(stringsAsFactors = F)
load("deg4.Rdata")
head(deg)
up=deg[deg$logFC>0.5,]
down=deg[deg$logFC<(-0.5),]
updown=rbind(up,down)
tail(updown)
nrow(down)
head(deg)
##33
library(clusterProfiler)
gene=bitr(rownames(deg),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Rn.eg.db")
head(gene)
nrow(gene2)
gene2=dplyr::distinct(gene,SYMBOL,.keep_all = TRUE)
head(gene2)
head(deg)
deg$SYMBOL=rownames(deg)

library(tidyverse)
deg2=deg %>% inner_join(gene2,by="SYMBOL")
head(deg2)
deg2_sort=deg2 %>% arrange(desc(logFC))
head(deg2_sort)

genelist=deg2_sort$logFC
names(genelist)=deg2_sort$ENTREZID
head(genelist)
gsekegg=gseKEGG(geneList = genelist,
                organism = "rno",
                pvalueCutoff = 1,
                pAdjustMethod = "BH")
head(gsekegg)
View(data.frame(gsekegg))
library(ggplot2)
library(enrichplot)
gseaplot2(gsekegg,
          title = "Steroid hormone biosynthesis",  #set title
          "rno00140", #plot "hsa04740"
          color="red", #line color
          base_size = 20, #font size
          subplots = 1:2, 
          pvalue_table = T) 




