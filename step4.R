#Enrichment
load(file = "deg4.Rdata")
head(deg)
up=deg[deg$logFC>0.5,]
down=deg[deg$logFC<(-0.5),]
updown=rbind(up,down)
tail(updown)
library(clusterProfiler)
g=bitr(geneID = rownames(updown),
       fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = "org.Rn.eg.db")

go_all=enrichGO(gene = g$ENTREZID,
                OrgDb = "org.Rn.eg.db",
                keyType = "ENTREZID",ont = "ALL")
library(ggplot2)
barplot(go_all,split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale = "free")
write.csv(go_all,"go_all_终.csv")

ekegg=enrichKEGG(gene = g$ENTREZID,
                 organism = "rno",
                 pvalueCutoff = 0.1)
dotplot(ekegg,showCategory=20)
write.csv(ekegg,"kegg_终.csv")

ewiki = enrichWP(gene = g$ENTREZID,
                 organism = "Rattus norvegicus",
                 pvalueCutoff = 0.5,
                 pAdjustMethod="BH",
                 minGSSize=3,
                 maxGSSize=1000,
                 qvalueCutoff=1)
dotplot(ewiki,showCategory=20)
library(ReactomePA)
ereac = enrichPathway(gene = g$ENTREZID,
                      pvalueCutoff = 0.5,
                      qvalueCutoff = 1,
                      minGSSize = 5,
                      maxGSSize = 1000,
                      pAdjustMethod = "none",
                      organism = "rat")


dotplot(ereac,showCategory=20)
head(ereac)
x=data.frame(ereac)
head(x[order(x$BgRatio,decreasing = T),])
View()

write.csv(ereac,"reac_zhogn.csv")
