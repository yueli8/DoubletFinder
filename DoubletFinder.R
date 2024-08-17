#https://www.jianshu.com/p/550f62642443如何修改R包源代码及使用修改后的R包-以DoubletFinder为例
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(dplyr)
library(Matrix)
library(clustree)
library(devtools)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(reshape2)
library(tidyverse)
library(garnett)
library(data.table)
library(DoubletFinder)
library(harmony)
library(patchwork)

setwd("~/nanjing/HEV_scRNA-seq_nankeda")

AHE01<-readRDS("ahe01.rds")
#assay
AHE01@assays
dim(AHE01@meta.data)
#AHE01<-RunHarmony(AHE01,"orig.ident")
AHE01<- RunUMAP(AHE01, dims = 1:20)
AHE01<- FindNeighbors(AHE01, dims = 1:20)
AHE01<- FindClusters(AHE01, resolution = 0.5)
clusters <- AHE01@meta.data$seurat_clusters
Idents(AHE01)<-"RNA_snn_res.05"
AHE01$seurat_clusters<-AHE01@active.ident
DimPlot(AHE01,label=T)

#remove doublet
# 对"scRNA_harmony"这个单细胞对象进行pN-pK参数扫描，以生成人工双细胞并计算每个细胞的pANN值
sweep.res.list <- paramSweep(AHE01, PCs = 1:20, sct = T)
# 对参数扫描的结果进行汇总，计算每种pN和pK组合下的双细胞检测指标。这些指标可以用来评估不同参数下的双细胞检测效果，并选择最优的参数。参数GT表示是否提供了真实的双细胞标签，此处没有提供
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
# 对汇总结果按照真实双胞胎比例（BCreal）进行升序排序，并显示排序后的数据框。这可以帮助我们找到双胞胎检测效果最好的参数组合
sweep.stats[order(sweep.stats$BCreal),]
# 根据汇总结果找到最优的pK参数。
bcmvn <- find.pK(sweep.stats) 
# 提取出全局最优的pK值，储存于"pK_bcmvn"
pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)]) 

# 估计的同源双细胞（由两个或多个相同类型的细胞融合而成的假阳性细胞，它们通常比异源双细胞更难以检测和去除）的比例     
homotypic.prop <- modelHomotypic(AHE01$seurat_clusters) 
# 计算总的双细胞数量（假设双细胞形成率为 7.5%）
nExp_poi <- round(0.004 *nrow(AHE01@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # 计算异源双细胞数量

# 使用确定好的参数鉴定doublets
scRNA_harmony <- doubletFinder(AHE01, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                               nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
#pK_bcmvn=14
DimPlot(scRNA_harmony,group.by = "DF.classifications_0.25_14_0", raster = FALSE)
AHE01<- subset(scRNA_harmony, subset = DF.classifications_0.25_14_0== "Singlet")
dim(AHE01)

dim(AHE02@meta.data)
AHE02<- RunUMAP(AHE02, dims = 1:20)
AHE02<- FindNeighbors(AHE02, dims = 1:20)
AHE02<- FindClusters(AHE02, resolution = 0.5)
clusters <- AHE02@meta.data$seurat_clusters
Idents(AHE02)<-"RNA_snn_res.05"
AHE02$seurat_clusters<-AHE02@active.ident
DimPlot(AHE02,label=T)
#remove doublet
sweep.res.list <- paramSweep(AHE02, PCs = 1:20, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
sweep.stats[order(sweep.stats$BCreal),]
bcmvn <- find.pK(sweep.stats) 
pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)]) 
homotypic.prop <- modelHomotypic(AHE02$seurat_clusters) 
nExp_poi <- round(0.023 *nrow(AHE02@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # 计算异源双细胞数量
AHE02 <- doubletFinder(AHE02, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
#always 0.25    pK_bcmvn=3  nExp_poi.adj=0
DimPlot(AHE02,group.by = "DF.classifications_0.25_3_0", raster = FALSE)
AHE02<- subset(AHE02, subset = DF.classifications_0.25_3_0== "Singlet")
dim(AHE02)

dim(AHE02@meta.data)
AHE02<- RunUMAP(AHE02, dims = 1:20)
AHE02<- FindNeighbors(AHE02, dims = 1:20)
AHE02<- FindClusters(AHE02, resolution = 0.5)
clusters <- AHE02@meta.data$seurat_clusters
Idents(AHE02)<-"RNA_snn_res.05"
AHE02$seurat_clusters<-AHE02@active.ident
DimPlot(AHE02,label=T)
#remove doublet
sweep.res.list <- paramSweep(AHE02, PCs = 1:20, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
sweep.stats[order(sweep.stats$BCreal),]
bcmvn <- find.pK(sweep.stats) 
pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)]) 
homotypic.prop <- modelHomotypic(AHE02$seurat_clusters) 
nExp_poi <- round(0.023 *nrow(AHE02@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # 计算异源双细胞数量
AHE02 <- doubletFinder(AHE02, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
#always 0.25    pK_bcmvn=3  nExp_poi.adj=0
DimPlot(AHE02,group.by = "DF.classifications_0.25_3_0", raster = FALSE)
AHE02<- subset(AHE02, subset = DF.classifications_0.25_3_0== "Singlet")
dim(AHE02)


