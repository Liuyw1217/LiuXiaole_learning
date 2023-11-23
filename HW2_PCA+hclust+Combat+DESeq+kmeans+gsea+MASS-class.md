---
title: "STAT115/215 BIO/BST282 HW2"
Subtitle: "Differential Expression Analysis & Sample Classification"
author: "Your Name"
date: "Due Date: Sun 31/7/2022 midnight "
output:
  html_document: default
  pdf_document: default
  word_document: default

---

## Part I: Differential expression

In this HW, we will evaluate the differentially expressed genes and pathways between breast cancer and normal breast tissues. Our collaborator generated RNA-seq on ten pairs of breast tumors and matched adjacent normal tissues, located at /mnt/data/data/HW2/raw_data1. The experiments were run in two batches, each batch with 5 pairs of samples, and the batch information is provided in batch.csv. We have run Salmon for gene quantification which is provided in Cannon at /mnt/data/data/HW2/raw_data1/Salmon_results. Remember to convert Ensembl ID to gene symbol, mapping between Ensembl id and gene symbol using "org.Hs.eg.db" package in R.

### Problem I.1

Please install the following R/Bioconductor packages. Then try "library(package)" to make sure the package works. 
Note: sva package with Combat function is used for batch effect removal; 
DESeq2 package is used for differential gene expression analysis; 
tximport package is used for importing transcript-level abundance into gene-level matrices for downstream analysis
ggplot2 package is used for general plotting in R; 
pheatmap package is used for heatmap plotting;
dplyr package is used for data frame manipulations in R;
fgsea package is used for gene-set enrichment analysis.

```{r
# ```{r install, eval = FALSE}
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("sva")
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")
# install.packages(c("ggplot2", "dplyr",
#                    "fgsea","pheatmap"))
library(ggplot2)
library(sva)
library(DESeq2)
library(tximport)
library(dplyr)
library(fgsea)
library(pheatmap)
library(ComplexHeatmap)

library(clusterProfiler)
library(factoextra)
library(FactoMineR)
library(ape)
library(ggtree)
library(Cairo)
library(org.Hs.eg.db)
```

### Problem I.2

For RNA-seq analysis, visualization using principle component analysis (PCA) or hierarchical clustering is an essential step of exploratory analysis to identify potental batch effect. Please import transcript-level TPM values from Salmon output and convert to gene-level TPM values. Perform PCA of the samples using log2 transformed TPM values. Indicate different tissues (tumor or normal, using shape to show) and batch (1 or 2, using color to show) of each sample. Next try to use hierarchical clustering for these samples.
Do the results suggest the presence of a batch effect?
For this question, you will load Salmon output at /mnt/data/data/HW2/raw_data1/Salmon_results. You also need to read in batch information provided in /mnt/data/data/HW2/raw_data1/batch.csv. Remember to convert Ensembl ID to gene symbol, mapping between Ensembl id and gene symbol using "org.Hs.eg.db" package in R.

```r
setwd("E:/learnimg/JLab/xiaoleLiu/week2")
path = "E:/learnimg/JLab/xiaoleLiu/week2/HW2_sf/"

batch <- read.csv("batch.csv")
##读取数据
T1 <- read.table("HW2_sf/T1.sf",sep="\t",header=T,fill=T)
#convert Ensembl ID to gene symbol
tx2gene <- bitr(substr(T1$Name,1,15),fromType = "ENSEMBLTRANS",toType = c("ENSEMBLTRANS","SYMBOL"),OrgDb = org.Hs.eg.db)

##获取HW2_sf文件夹下面的所有文件的文件名
filelist <- list.files(path)  
files <- paste(path,filelist,sep="") ##files为所有的路径
#使用tximport一步完成导入salmon生成的quant.sf文件，并将数据从transcripts水平转换成gene水平，获得non-nromalized count矩阵
#将转录水平转换为基因水平  https://blog.csdn.net/weixin_46585008/article/details/109362501
txi <- tximport(files,type = "salmon",tx2gene = tx2gene,countsFromAbundance = "lengthScaledTPM",ignoreTxVersion = TRUE)

#TPM获取
tpm <- txi$abundance
expr.dat <- t(tpm)
expr.dat <- log2(expr.dat+1)

#---PCA of the samples using log2 transformed TPM values----------------------------------
# scale. = TRUE表示分析前对数据进行归一化；
# 去除方差为零的列：
expr.dat <- expr.dat[, apply(expr.dat, 2, var) > 0]
# 去除重复的列：
expr.dat <- expr.dat[, !duplicated(expr.dat)]
# PCA
pca_result <- prcomp(expr.dat, center = TRUE, scale. = TRUE)
summary(pca_result)
# Importance of components:
#   PC1     PC2     PC3      PC4      PC5      PC6      PC7      PC8      PC9      PC10
# Standard deviation     45.1708 36.6899 27.6209 26.09908 25.06347 21.97337 21.70589 19.41518 16.91960 3.246e-13
# Proportion of Variance  0.2884  0.1902  0.1078  0.09626  0.08878  0.06823  0.06658  0.05327  0.04046 0.000e+00
# Cumulative Proportion   0.2884  0.4786  0.5864  0.68268  0.77145  0.83969  0.90627  0.95954  1.00000 1.000e+00

# 查看主成分的结果
head(pca_result$x)
# PC1         PC2         PC3         PC4         PC5 
# [1,] -62.56691  25.7230357 -21.4822088  17.2258942  
# [2,] -48.16334  12.9585128  -8.4761015  -8.8709393  
# [3,]  -7.91455 -72.7614181   2.2716855  -5.8765608

# 绘制PC1 PC2二维图
plot(pca_result$x, main="after PCA")

# 绘制PCA图
colnames(pca_result$x) <- substr(filelist,1,2)
plot_data <- data.frame(Label=colnames(pca_result$x),Pc1=pca_result$x[,1],Pc2=pca_result$x[,2],Tissue=batch$tissue,Batch=as.factor(batch$batch))

ggplot(plot_data,aes(Pc1,Pc2,col=Batch,shape=Tissue)) + 
  geom_point(size=3) + 
  geom_text(aes(label=Label),vjust = "outward") + 
  geom_hline(yintercept = 0,lty=2,col="red") + 
  geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
  theme_bw() + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="PCA_1",y="PCA_2",title = "PCA analysis")

ggsave("PCA_p1p2.png",width=8, height=8) 
```

![](D:\software\MarkText\picture\2023-11-20-10-04-17-PCA_p1p2.png)

```r
#---hierarchical clustering for these samples. Do the results suggest the presence of a batch effect?--------------------
data_tpm <- as.data.frame(t(tpm))
rownames(data_tpm) <- paste(substr(filelist,1,2),batch$batch,sep="-")
data_tpm <-log2(data_tpm+1)##标准化数据
tree.dist <- dist(scale(data_tpm))#计算距离 默认euclidean.
heatmap(as.matrix(tree.dist))

tree.hc = hclust(tree.dist)#层次聚类  默认complete
dend <- as.dendrogram(tree.hc) #聚类树

pdf("hclust.pdf",width=7, height=7)  
plot(dend,main="hierarchical clustering  before remove batch effect ",xlab="samples",ylab="height")
dev.off()
```

![](D:\software\MarkText\picture\2023-11-20-10-04-52-image.png)

### Problem I.3

Run COMBAT on the samples to remove the batch effect. Visualize the results using a similar PCA and hierarchical clustering as Problem 2. Provide evidence that the batch effects are successfully adjusted. 

```r
##------------去除批次xiaoying------------------------
#一个数据框或矩阵。数据的行应该是特征，列应该是样本,log2-TPM
my_data <- expr.dat
#创建一个设计矩阵，其中包含关于每个样本的批次信息。这告诉 ComBat 哪些样本属于哪个批次。
design <- model.matrix(~batch$tissue)
combat_result <- ComBat(dat = my_data, batch = batch$batch,mod = design)
exp_combat <- as.data.frame(t(combat_result))
exp_combat <- exp_combat[, apply(exp_combat, 2, var) > 0]
# 去除重复的列：
exp_combat <- exp_combat[, !duplicated(exp_combat)]
# PCA
pca_result <- prcomp(exp_combat, center = TRUE, scale. = TRUE)
summary(pca_result)
# Importance of components:
#   PC1     PC2     PC3     PC4      PC5      PC6      PC7
# Standard deviation     47.1930 31.5364 29.6661 26.7444 25.24607 23.68910 21.22897
# Proportion of Variance  0.3147  0.1406  0.1244  0.1011  0.09007  0.07931  0.06369
# Cumulative Proportion   0.3147  0.4553  0.5797  0.6808  0.77084  0.85014  0.91383

# 绘制PCA图
colnames(pca_result$x) <- substr(filelist,1,2)
plot_data <- data.frame(Label=colnames(pca_result$x),Pc1=pca_result$x[,1],Pc2=pca_result$x[,2],Tissue=batch$tissue,Batch=as.factor(batch$batch))

ggplot(plot_data,aes(Pc1,Pc2,col=Batch,shape=Tissue)) + 
  geom_point(size=3) + 
  geom_text(aes(label=Label),vjust = "outward") + 
  geom_hline(yintercept = 0,lty=2,col="red") + 
  geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
  theme_bw() + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="PCA_1",y="PCA_2",title = "After Combat PCA analysis")

ggsave("PCA_p1p2_combat.png",width=8, height=8)
```

![](D:\software\MarkText\picture\2023-11-20-10-06-00-PCA_p1p2_combat.png)

```r
##---层次聚类----
data_tpm <- exp_combat
rownames(data_tpm) <- paste(substr(filelist,1,2),batch$batch,sep="-")
data_tpm <-log2(data_tpm+1)##标准化数据
tree.dist <- dist(scale(data_tpm))#计算距离 默认euclidean.
tree.hc = hclust(tree.dist)#层次聚类  默认complete
dend <- as.dendrogram(tree.hc) #聚类树

pdf("hclust_combat.pdf",width=7, height=7)  
plot(dend,main="hierarchical clustering remove batch effect ",xlab="samples",ylab="height")
dev.off()
```

![](D:\software\MarkText\picture\2023-11-20-10-06-48-image.png)

### Problem I.4

Run DESeq2 based on paired samples adjusting for the batch effect to identify differentially-expressed genes between tumor and normal tissues. How many genes are expressed higher in tumors than normal. Let's use 1) FDR < 0.01 and 2) Log2 fold change > 1 as the cutoff. 
Note: please use the raw_count columns of the Salmon result and convert these to integer values for DESeq2.
Identify the top 5 most (by Log2FC) over expressed genes (FDR < 0.01) in tumor and normal, respectively.  

```r
#----DESeq2---------------------------
data_count <- ceiling(txi$counts)
colData <- data.frame(batch=batch$batch, tissue=factor(batch$tissue))
# 去除批次效应
dds <- DESeqDataSetFromMatrix(countData = data_count, colData = colData, design = ~ batch+tissue)
#design = ~ batch+tissue ，这里面的有两个变量，一个变量代表了我们的处理组和对照组的信息，另外一个变量就是批量信息。

head(dds) #查看构建好的矩阵
dds1 <- DESeq(dds)
#将结果用result()函数来获取,使用results函数设置padj cutoff(alpha<0.01)
res <- results(dds1,alpha = 0.01, contrast = c("tissue","Tumor","Normal"))
summary(res) # res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# 依次按照log2FoldChange和pvalue进行排序
res1 <- res1[order(res1$log2FoldChange,res1$pvalue, decreasing = c(TRUE,FALSE)), ]

# 获取padj(FDR)小于0.01，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
res1_up<- res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.01),]      # 表达量显著上升的基因 632个
res1_down<- res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.01),]    # 表达量显著下降的基因 552个
res_diff <- rbind(res1_up,res1_down)  #1184个差异表达基因

top_5 <- rownames(res1_up)[1:5]
# "EEF1A1P44" "SOX1"      "NGB"       "SIX6"      "ZDHHC22" 
head(res1_up,5)
```

![](D:\software\MarkText\picture\2023-11-20-10-27-11-image.png)

### Problem I.5

Visualize the differential gene expression values by making a volcano and an MA plot to summarize the differences between tumor and normal. In the volcano plot, draw the statistically significant (FDR < 0.01, Log2FC > 1) tumor up-genes red and down-genes blue.
Note: Be sure to use the lfcShrink function to get more robust estimates of the fold-changes for genes.

```r
# 可视化
resultsNames(dds1)
#logFoldChange缩小有助于基因的可视化和比较。可指定apeglm方法将dds对象传递给lfcShrink来缩小。缩小LFC可消除low readcount基因中与LFC变化相关的噪声。
library(apeglm)
res_shrink <- lfcShrink(dds1,coef = "tissue_Tumor_vs_Normal",res=res,type="apeglm")

#DESeq2可视化-------------------------------------------
#--1--MA plot,图中x轴为baseMean，y轴为logFoldChange，若padj<0.01点标为红色。
pdf("MA plot.pdf",width=9, height=7)
plotMA(res_shrink,ylim=c(-15,15),alpha=0.01,xlab="Mean of Normalized Counts",main="Tumor vs Normal")
dev.off()

#--2--volcano
# 数据整合
data_shrink <- data.frame(res_shrink) %>% na.omit() %>% 
  mutate(logp=-log10(padj)) %>%
  mutate(group=ifelse(logp >= 2 & log2FoldChange >= 1,"Up",ifelse(log2FoldChange <= -1 & logp >= 2 ,'Down','No change')))

pdf("volcano_diff.pdf",width=9, height=7)
ggplot(data_shrink,aes(log2FoldChange, logp))+
  geom_hline(yintercept = 2, linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = group),size = 1) +
  scale_colour_manual(values=c("Up"="red","Down"="blue","No change"="grey"))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  xlab("log2FoldChange")+
  ylab("-Log10(FDR q-value)")
dev.off()
```

![](D:\software\MarkText\picture\2023-11-20-10-13-53-image.png)

![](D:\software\MarkText\picture\2023-11-20-10-14-37-image.png)### Problem I.6

Try kmeans (try k = 4 or 7) clustering to group differentially expressed genes into different clusters. How many genes are there in each cluster? Draw a heatmap showing gene clusters. 

```r
##----KNN----------------TPM数据---------------------
diff.expr <- as.data.frame(t(exp_combat[,rownames(res_diff)]))
colnames(diff.expr) <- batch$X

library(cluster)
library(fpc)
library(circlize)
# iter.max: 最大迭代次数
# nstart: 选择的随机集的数目
# centers: 最优类数目
set.seed(123)
km_7 <- kmeans(diff.expr, centers=7, iter.max=100, nstart=25)
table(km_7$cluster)
# 1    2    3    4    5    6    7 
# 1155    1    1    1   20    1    5 
km_4 <- kmeans(diff.expr, centers=4, iter.max=100, nstart=25)
table(km_4$cluster)
# 1    2    3    4 
# 1    1 1181    1 
pdf("kmean_4_7.pdf",width=6, height=9)
p7 <- Heatmap(diff.expr,name = "heatmap", 
             km = 7,
             column_names_side = "bottom",
             col = colorRamp2(c(0, 200, 400), c("blue", "white", "red")),
             cluster_columns = FALSE,
             row_dend_side = "left",
             show_row_names = FALSE
             # km_title = "%i"
)
p7
p4 <- Heatmap(diff.expr,name = "heatmap", 
              km = 4,
              column_names_side = "bottom",
              col = colorRamp2(c(0, 200, 400), c("blue", "white", "red")),
              cluster_columns = FALSE,
              row_dend_side = "left",
              show_row_names = FALSE
              # km_title = "%i"
)
p4
dev.off()
```

<img title="" src="file:///D:/software/MarkText/picture/2023-11-20-10-17-52-image.png" alt="" width="319"><img title="" src="file:///D:/software/MarkText/picture/2023-11-20-10-18-08-image.png" alt="" width="320">

### Problem I.7: For graduate students only

If you run DESeq2 without removing batch effect, how many differential genes do you get? How do their K-means clustering look? Does batch effect removal gives more interpretable results? 

```r
#####-----对未进行批次矫正的数据进行DESeq,并进行Kmeans,比较差别------------
colData <- data.frame(batch=batch$batch, tissue=factor(batch$tissue))
dds <- DESeqDataSetFromMatrix(countData = data_count,colData = colData,design = ~tissue)
dds1 <- DESeq(dds)
res <- results(dds1,alpha = 0.01, contrast = c("tissue","Tumor","Normal"))
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$log2FoldChange,res1$pvalue, decreasing = c(TRUE,FALSE)), ]

res1_up<- res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.01),]      # 表达量显著上升的基因 564个
res1_down<- res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.01),]    # 表达量显著下降的基因 498个
res_diff <- rbind(res1_up,res1_down)  #1062个差异表达基因

##Kmeans
diff.expr <- as.data.frame(t(exp_combat[,rownames(res_diff)]))
colnames(diff.expr) <- batch$X

set.seed(123)
km_7 <- kmeans(diff.expr, centers=7, iter.max=100, nstart=25)
table(km_7$cluster)
km_4 <- kmeans(diff.expr, centers=4, iter.max=100, nstart=25)
table(km_4$cluster)
```

<img src="file:///D:/software/MarkText/picture/2023-11-20-10-23-00-image.png" title="" alt="" width="286">

### Problem I.8

From the batch-removed DESeq2 run results, extract the top 200 tumor-upregulated genes (by Log2FC, FDR < 0.01). Run DAVID GO analysis (http://david.abcc.ncifcrf.gov/) to see whether these genes are enriched in specific biological process (BP), pathways, etc.

```r
##---From the batch-removed DESeq2 run results, extract the top 200 tumor-upregulated genes (by Log2FC, FDR < 0.01)-------
top_200 <- rownames(res1_up)[1:200]
write.table(res1_up[1:200,],"top_200_UP.txt",sep="\t",col.names = T,row.names = T,quote = F)

##GSEA---fgsea--------------------
library(fgsea)
data_gsea <- res_diff[,2]
#提取log2FC,名字对应基因名，进行GSEA富集分析
names(data_gsea) <- rownames(res_diff)
gmt<-gmtPathways("c2.cp.kegg.v7.4.symbols.gmt")
set.seed(123)
fgseaRes <- fgsea(pathways=gmt,stats=data_gsea, nperm=1000)

# 挑选显著性通路,
# 一般认为|NES|>1，p-value<0.05，FDR<0.25的通路是显著富集的。 |NES|值越大，FDR值就越小，说明分析的结果可信度越高。
# FDRVal_up <- subset(gsea.res, padj < 0.05 & NES >0) %>% arrange(.,padj)
# FDRVal_down <- subset(gsea.res, padj < 0.05 & NES <0) %>% arrange(.,padj)

top_5pathway <- fgseaRes[ES!=0][head(order(pval),n=5),pathway]

pdf("gsea_top_5pathway.pdf",width=10, height=5)
plotGseaTable(gmt[top_5pathway],data_gsea,fgseaRes = fgseaRes,gseaParam=0.5)
dev.off()

pdf("KEGG_enrichment_pathway.pdf",width=18, height=7)
p1 <- plotEnrichment(gmt[["KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS"]],data_gsea)+ 
  labs(title = "KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS")
p2 <- plotEnrichment(gmt[["KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"]],data_gsea)+ 
  labs(title = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION")
p1+p2
dev.off()
```

![](D:\software\MarkText\picture\2023-11-20-10-23-49-image.png)

![](D:\software\MarkText\picture\2023-11-20-10-29-02-image.png)

### Problem I.9: For graduate students only

Run Gene Set Enrichment analysis (http://www.broadinstitute.org/gsea/index.jsp) using the summary statistics from problem 4. Show the top five gene sets or pathways that best capture the differentially expressed genes between tumor than normal. Comment on the biological relevance of the results. Plot GSEA enrichment plot for an interesting pathway. 
Mapping between gene sets and pathways is provided in /mnt/data/data/HW2/raw_data1/c2.cp.kegg.v7.4.symbols.gmt file.

![](D:\software\MarkText\picture\2023-11-20-10-29-53-image.png)

![](D:\software\MarkText\picture\2023-11-20-10-30-22-image.png)

## Part II: Sample classification

We provide you z-score normalized expression data of 50 breast tumor samples, 50 normal breast samples (your training and cross-validation data), and 20 samples without diagnosis (your testing data). We want to use the 100 samples with known diagnosis to train machine learning models in order to predict the 20 unknown samples. 
You will need the following libraries in R: ggplot2 and ggfortify for plotting, MASS and caret for machine learning, and pROC is for evaluating testing performance. The YouTube video on caret (https://youtu.be/z8PRU46I3NY) and the package documentation (http://topepo.github.io/caret/index.html) might be helpful.
All data for Part II are provided at /mnt/data/data/HW2/raw_data2.

### Problem II.1

Run PCA for dimension reduction on the 100 samples with known labels, and draw these 100 samples in a 2D plot. Do cancer and normal separate from the first two PCs? Would this be sufficient to classify the unknown samples?
z-score normalized data are provided in BRCA_zscore_data.txt. Phenotype data is in BRCA_phenotype.txt.

```r
library(ggplot2)
library(ggfortify)
library(MASS)
library(caret)
library(pROC)

setwd("E:/learnimg/JLab/xiaoleLiu/week2")

BRCA_data <- read.table("BRCA_zscore_data.txt",sp="\t",header = T)
BRCA_phenotype <- read.table("BRCA_phenotype.txt",sep="\t",header=T)

##--PCA降维---------
PCA_res <- prcomp(BRCA_data,scale. = T)
summary(PCA_res)

data_plot_pca <- data.frame(sample=rownames(PCA_res$x),PC1 = PCA_res$x[,1],PC2=PCA_res$x[,2],phenotype=BRCA_phenotype$phenotype)
p1 <- ggplot(data_plot_pca,aes(x=PC1,y=PC2,col=phenotype))+
  geom_point(size=3) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="PC1",y="PC2",title = "BRCA zscore data PCA analysis")
p1
```

![](D:\software\MarkText\picture\2023-11-20-10-34-09-image.png)

### Problem II.2

Draw a plot showing the cumulative % variance captured from the top 100 PCs. How many PCs are needed to capture 90% of the variance? 

```r
# Draw a plot showing the cumulative % variance captured from the top 100 PCs.--------------------
cumulative_var <- summary(PCA_res)$importance[3,]
which(cumulative_var>0.9)[1]  #25
# Plot the cumulative variance
plot(1:100, cumulative_var, type = "b", xlab = "Number of PCs", ylab = "Cumulative % Variance",
     main = "Cumulative % Variance Explained by Top 100 PCs", grid = TRUE)
abline(h=0.9)
```

![](D:\software\MarkText\picture\2023-11-20-10-34-54-image.png)

### Problem II.3

Apply machine learning methods (KNN, logistic regression, Ridge regression, LASSO, ElasticNet, random forest, and support vector machines) on the top 25 PCs of the training data and 5-fold cross validation to classify the samples. caret and MASS already implemented all of the machine learning methods, including cross-validation, so calling each is only one command. In order to get consistent results from different runs, use set.seed(). 

![](D:\software\MarkText\picture\2023-11-20-10-40-37-image.png)

```r
####----Apply machine learning methods------5-fold cross validation------------
## top 25 PCs of the training data
TrainData <- as.data.frame(PCA_res$x[,1:25])
# 对训练集进行PCA降维，保留前25个PCs
# PCA_res <- prcomp(BRCA_data, center = TRUE, scale. = TRUE)
# TrainData <- predict(PCA_res, newdata = BRCA_data)[, 1:25]

# 创建控制参数
mycontrol <- trainControl(method = "cv",   # 交叉验证方法
                          number = 5,      # 折数
                          verboseIter = FALSE,   # 显示每次迭代的详细信息
                          classProbs = TRUE)    # 输出类别概率                          
set.seed(123)

#--1--# KNN
model_knn <- train(TrainData,
                 y = BRCA_phenotype$phenotype,
                 method = "knn",
                 trControl = mycontrol)

#--2--# logistic regression
model_logic <- train(TrainData,   # 训练数据
               y = BRCA_phenotype$phenotype,   # 目标变量
               method = "glm",   # 使用逻辑回归方法
               trControl = mycontrol)   # 控制参数对象
# family="binomial"指定了逻辑回归的二元分类问题

#--3--# Ridge regression
model_ridge <- train(TrainData,   
               y = BRCA_phenotype$phenotype,   
               method = "ridge",   # 使用岭回归方法
               trControl = mycontrol)   
# Error: wrong model type for classification
# 岭回归（Ridge Regression）是一种用于回归问题的方法，而不是分类问题。因此，在使用岭回归时，不能将其用于分类任务
model_ridge <- train(TrainData,
               y = BRCA_phenotype$phenotype,
               method = "glmnet",
               trControl = mycontrol,
               tuneGrid = expand.grid(alpha = 0, lambda = seq(0.1, 1, by = 0.1)))  # 设置alpha为0，表示使用Ridge回归

#--4--# LASSO
model_lasso <- train(x = TrainData,   
               y = BRCA_phenotype$phenotype,  
               method = "glmnet",   # 使用glmnet方法
               trControl = mycontrol,
               tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 1, length = 100))) #alpha设置为1表示只使用L1正则化（即Lasso回归）

#--5--# ElasticNet
#glmnet方法默认使用的是弹性网络（Elastic Net），它是一种同时使用L1（Lasso）和L2（Ridge）正则化的线性回归方法。
# ElasticNet是一种结合了L1正则化（LASSO）和L2正则化（岭回归）的线性模型方法，用于解决具有高维特征和相关性的数据集的问题。
model_elasticnet <- train(x = TrainData,   
               y = BRCA_phenotype$phenotype, 
               method = "glmnet",   # 使用glmnet方法
               trControl = mycontrol,   # 控制参数对象
               tuneGrid = expand.grid(alpha = 0.5, lambda = seq(0.1, 1, by = 0.1)))   # 超参数网格
#`method`参数设置为"glmnet"以指定使用ElasticNet方法，`tuneGrid`参数用于指定超参数的网格搜索范围，其中`alpha`是ElasticNet混合比例，取值范围为0到1，`lambda`是正则化项的惩罚参数。

#--6--# random forest
model_rf <- train(x = TrainData,   
               y = BRCA_phenotype$phenotype, 
               method = "rf",   # 使用随机森林方法
               trControl = mycontrol)   # 控制参数对象

#--7--# support vector machines
#SVM分类分析的基本思想是将样本映射到高维特征空间，并在该空间中找到一个最优超平面，使得不同类别的样本能够被最大间隔所分开。对于线性可分的情况，使用线性核函数（例如线性核或多项式核）进行分类。对于非线性可分的情况，使用非线性核函数（例如径向基函数核）进行分类。
model_svm <- train(x = TrainData,   # 训练数据
               y = BRCA_phenotype$phenotype, 
               method = "svmRadial",   # 使用径向基函数核的支持向量机
               trControl = mycontrol)   # 控制参数对象
```

### Problem II.4

Summarize the performance of each machine learning method, in terms of accuracy and kappa. 

```r
## 从Accuracy和Kappa两方面总结了每种机器学习方法的性能。
#通过train函数训练了多个模型，存储在models列表中
models <- list(model_knn,model_logic,model_ridge,model_lasso,model_elasticnet,model_rf,model_svm)
# 使用resamples函数计算模型的性能度量结果
results <- resamples(models)
# 查看每个模型的性能度量结果
summary(results)

# 特定性能度量指标
values_model <- results$values

# lambda: lambda是弹性网络（Elastic Net）中的正则化参数。对于Lasso回归（L1正则化），lambda控制了稀疏性和特征选择的程度。在train函数中，通过交叉验证或网格搜索过程，找到了最佳的lambda值。
# Accuracy: 准确率是分类模型中常用的性能度量指标之一，表示模型预测的正确分类样本数与总样本数之间的比例。它是通过对预测结果与实际类别进行比较计算得出的。
# Kappa: Kappa系数（Cohen's Kappa）是一种用于衡量分类模型性能的统计度量。它考虑了模型的准确率，并根据预期的随机一致性对准确率进行校正。Kappa系数的取值范围在[-1, 1]之间，值越接近1表示模型性能越好。
```

![](D:\software\MarkText\picture\2023-11-20-10-47-08-image.png)

![](D:\software\MarkText\picture\2023-11-20-10-44-25-image.png)

### Problem II.5: For Graduate students only

Compare the performance difference between logistic regression, Ridge, LASSO, and ElasticNet. In LASSO, how many PCs have non-zero coefficient? What is the lambda for Ridge and LASSO, respectively? 

```r
# 获取Lasso回归模型的系数
lasso_coef <- coef(model_lasso$finalModel, s = model_lasso$bestTune$lambda)

# 打印每个变量的系数
print(lasso_coef)
```

![](D:\software\MarkText\picture\2023-11-20-10-45-42-image.png)

### Problem II.6

Use the PCA projections in Q1 to obtain the first 25 PCs of the 20 unknown samples. Use one method that performs well in Q4 to make predictions. Caret already used the hyper-parameters learned from cross-validation to train the parameters of each method on the full 100 training data. You just need to call this method to make the predictions. 
Expression data for the 20 unknown samples are provided in unknown_samples.txt.

```r
## 预测未知样本类别-----------------
un_sample <- read.table("unknown_samples.txt",sep="\t",header=T)
diagnose <- read.table("diagnosis.txt",sep="\t",header = T)

# 使用之前训练的PCA模型将测试集数据映射到前25个PCs
testData <- predict(PCA_res, newdata = un_sample)[, 1:25]
# 使用训练得到的模型对测试集的PCA结果进行预测
predictions <- predict(model_svm, newdata = testData)
```

![](D:\software\MarkText\picture\2023-11-20-10-50-42-image.png)

### Problem II.7: For Graduate students only

Can you find out the top 3 genes that are most important in this prediction method in Q6? Do they have some known cancer relevance? 

```r
##找出重要性前3的基因------------------------
# 对于SVM模型，通常没有直接的特征重要性度量。然而，你可以考虑使用特征权重（feature weights）来近似表示特征的重要性。特征权重反映了SVM模型中每个特征对决策边界的贡献。
# 对于使用"svmRadial"方法训练的ksvm对象，无法直接提取特征权重或系数。使用varImp()函数来计算特征的相对重要性，该函数可以基于特征在模型中的重要程度进行排序。
# 提取特征相对重要性
feature_importance <- varImp(model_svm, scale = FALSE)
feature_scores <- feature_importance$importance

# 找出重要性前3的基因
top_PCs <- dplyr::arrange(feature_scores,desc(Normal))[1:3,]
# Normal     Tumor
# PC2 0.9807923 0.9807923
# PC8 0.7482993 0.7482993
# PC5 0.6810724 0.6810724

#提取主成分贡献率：
#主成分分析的结果中，可以通过pca_result$rotation来获取主成分的旋转（加载）矩阵
sort(PCA_res$rotation[,2],decreasing = TRUE)[1:3]
# TRIM11       PAK4       CDK5 
# 0.02007697 0.01971736 0.01965764 
sort(PCA_res$rotation[,8],decreasing = TRUE)[1:3]
# MANSC1      RPRD2     GALNT2 
# 0.03584126 0.03080776 0.03036478 
sort(PCA_res$rotation[,5],decreasing = TRUE)[1:3]
# RPS7      RPS19      RPL13 
# 0.02950072 0.02930869 0.02906633 
```

### Problem II.8

Suppose a pathologist later made diagnosis on the 20 unknown samples (load the diagnosis.txt file). Based on this gold standard, draw an ROC curve of your predictions in Q6. What is the prediction AUC? 

```r
## ROC---------------------------------------
library(pROC)

# Obtain raw scores from the SVM model
scores <- predict(model_svm, newdata = testData, type = "raw")
# 将预测结果转换为二元变量
predictions_score <- as.numeric(scores)  # 将scores转换为数值类型的预测结果
predictions_score <- ifelse(predictions_score > 1, 1, 0)  # 大于1的预测结果视为正类（Tumor），否则为负类（Normal）

# 计算ROC曲线和AUC
library(pROC)
roc_obj <- roc(diagnose$phenotype, predictions_score)
auc <- auc(roc_obj) #Area under the curve: 0.899

# 绘制ROC曲线
plot(roc_obj, main = "ROC Curve", xlab = "False Positive Rate", ylab = "True Positive Rate")
text(x=0.8,y=0.8,"Area under the curve: 0.899")
```

![](D:\software\MarkText\picture\2023-11-20-10-52-45-image.png)
