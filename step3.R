rm(list = ls())  
options(stringsAsFactors = F)
load(file = 'step2-output.Rdata')
dat4=log(dat3)

# Check the data every time

group_list=g4
dat4[1:4,1:4] 
table(group_list) 

boxplot(dat4[1,]~group_list) 

bp=function(g){         
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
bp(dat4[1,]) 
bp(dat4[2,])
bp(dat4[3,])
bp(dat4[10030,])

dim(dat4)
dat4=log(dat4)
head(dat4)
design
write.csv(dat4,file = "dat4.csv")
class(dat4)
head(dat4)
dat4 <- read.csv("~/PXR/dat4.csv", header=T)
head(dat4)
dat4[,1:14][c("Nr3c1"),]

dat4=as.dat4a.frame(dat4)
rownames(dat4)=dat4[,1]

library(limma)
design=model.matrix(~factor( group_list ))
fit=lmFit(dat4,design)
fit=eBayes(fit)
options(digits = 4) 
#topTable(fit,coef=2,adjust='BH') 
topTable(fit,coef=2,adjust='BH') 

design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
head(design)
exprSet=dat4
rownames(design)=colnames(dat4)
design
colnames(design)=c("hp","hc")
contrast.matrix<-makeContrasts("hp-hc",
                               levels = design)
contrast.matrix 

deg = function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  ##important step
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}

deg = deg(exprSet,design,contrast.matrix)

head(deg)

save(deg,file = 'deg4.Rdata')

write.csv(deg,"deg5.csv")
####enhancevolcano
library(EnhancedVolcano)

selectLab_activatess=c(
  "Kat2b","Per2","Per1","Foxo3",
  "Nr3c2","Ugt1a1","Nfkbia","Dusp1","G6pc1",
  "Cebpa ","Cav1", "Hdac1", "Klf9"," Klf5", "Ncoa1", "Nfil3", "Trim63",
  "Tsc22d3"," Klf15")
selectLab_inactivates=c("Atp1b1", "Jun", "Mapk8", "Ar", 
                        " Nr4a1", "Mapk3", "Mapk1", "Nr4a2", "Nr4a3s")
lab_italics=rownames(deg)
p <- EnhancedVolcano(
  deg,
  x = "logFC",
  y = "P.Value",
  lab = lab_italics,    
  selectLab = selectLab_inactivates, 
  boxedLabels = TRUE,
  #pCutoff = 0.01,      
  #FCcutoff = 2,         
  #cutoffLineWidth = 0.7,      
  #cutoffLineType = "twodash", 
  xlim = c(-2.1, 2.1),  
  ylim = c(0, 4), 
  pointSize = 2.5,        
  labSize = 4,        
  #xlab = bquote(~Log[2] ~ "fold change"),    
  #ylab = bquote(~-Log[10] ~ italic(p-value)), 
  #axisLabSize = 17,     # 坐标轴字体大小
  # title = "mRNA differential expression",    
  titleLabSize = 18,    
  subtitle = bquote(italic("Volcano plot")), 
  subtitleLabSize = 17, 
  legendLabSize = 13,  
  col = c("darkgreen", "red3", "royalblue", "red3"), 
  colAlpha = 0.3,      
  gridlines.major = FALSE,     
  gridlines.minor = FALSE,
  drawConnectors = T,
  widthConnectors = 1.6,
  colConnectors = "yellow"
)
p

View(deg)

deg2=deg
deg2$name=rownames(deg)
View(deg2)

rm(list = ls())
## for heatmap 

## for heatmap 

load(file = "deg4.Rdata")
dat=dat3
group_list=g4
if(T){ 
  load(file = 'step1-output.Rdata')
  dat[1:4,1:4]
  table(group_list)
  x=deg$logFC 
  names(x)=rownames(deg) 
  cg=c(names(head(sort(x),100)),
       names(tail(sort(x),100)))
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F)
  n=t(scale(t(dat[cg,])))
  
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n) 
  pheatmap(n,show_colnames =F,
           show_rownames = F,
           cluster_cols = F, 
           annotation_col=ac,filename = 'heatmap_top200_DEG.png') 
  
}

write.csv(deg,file = 'deg.csv')
###################3
