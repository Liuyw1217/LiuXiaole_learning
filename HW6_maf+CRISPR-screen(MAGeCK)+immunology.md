# Spring 2022 STAT115/215 BIO/BST282

#### Ya Han

## Part I: Data exploration on TCGA

The Cancer Genome Atlas (TCGA) is an NCI project to comprehensively 
profile over 10K tumors in 33 cancer types. In this homework, we are 
going to explore TCGA data analysis.

Q1. Go to TCGA GDC website ([https://portal.gdc.cancer.gov/](https://portal.gdc.cancer.gov/)) and explore the GDC data portal. How many glioblastoma (GBM) cases in TCGA meet ALL of the following requirements?

1. Male;

2. Diagnosed at the age above 45;

3. Still alive.

Answers:

Q2. TCGA GDC ([https://portal.gdc.cancer.gov/](https://portal.gdc.cancer.gov/)) and Broad Firehose ([http://firebrowse.org/](http://firebrowse.org/))
 both provide processed TCGA data for downloading and downstream 
analysis. Download clinical data of GBM. What’s the average diagnosed 
age of all GBM patients?

Answers:

## Part II – Tumor Subtypes

Q1. GBM is one of the earliest cancer types to be processed by TCGA, 
and the expression profiling was initially done with Affymetrix 
microarray. Also, with brain cancer, it is hard to get sufficient number
 of normal samples. We provide the pre-processed expression matrix in 
(GBM_expr.txt) where samples are columns and genes are rows. Do a 
K-means (k=3) clustering from all the genes and the most variable 2000 
genes. Do tumor and normal samples separate in different clusters? Do 
the tumors samples consistently separate into 2 clusters, regardless of 
whether you use all the genes or most variable genes?

Answers:

```r
exp <- read.table("GBM_expr.txt",header=T,sep="\t",row.names = 1)
cli <- read.table("GBM_clin.txt",header=T,sep="\t")
########################---------K-means------####################################################################
# 设置随机种子，让结果可以重现

set.seed(123)

# 调用kmeans聚类算法 k = 4

km_sample <- kmeans(t(exp), centers = 2)

#选择高变异的2000个基因重新进行聚类分析
gene_sd <- apply(exp, 1, sd) # 计算每个基因在样本中的标准差

# 排序基因

sorted_genes <- sort(gene_sd, decreasing = TRUE) # 按照标准差从高到低排序

# 选择前2000个基因

top_2000_exp <- exp[names(sorted_genes)[1:2000],] # 获取前2000个基因的表达谱

#进行看K-means
km_2000 <- kmeans(t(top_2000_exp),centers = 2)

#发现结果并不理想，没有用全部基因聚类的结果好，肿瘤和正常样本没有被分开
```

Q2. LIMMA is a BioConductor package that does differential expression
 between microarrays, RNA-seq, and can remove batch effects (especially 
if you have experimental design with cmplex batches). Use LIMMA to see 
how many genes are differentially expressed between the two GBM subtypes
 (with FDR < 0.05% and logFC > 1.5)?

Answers:

```r
# NO have ?
```

![](D:\software\MarkText\picture\2023-11-10-12-20-52-image.png)

```r
library(limma)
library(edgeR) #edgeR将同时引入limma

group <- c(rep("Tumor", 60), rep("Normal",10)) %>% factor(., levels = c("Tumor", "Normal"))
# 创建设计矩阵，指定组别信息
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(exp)# 创建 DGEList 对象
dge <- DGEList(counts = exp)

# 这里我们使用上面提到的 filterByExpr() 进行自动过滤，去除低表达基因,edgeR包函数 
keep <- filterByExpr(dge)
dge <- dge[keep, ,keep.lib.sizes = FALSE]

# 归一化，得到的归一化系数被用作文库大小的缩放系数
dge <- calcNormFactors(dge)

# 使用 voom 方法进行标准化
v <- voom(dge, design, plot = TRUE, normalize = "quantile")

# 如果是芯片数据、TPM数据或已标准化的数据，不需要再进行标准化，可直接从这里开始进行差异分析
# 使用线性模型进行拟合
fit <- lmFit(v, design)

# 和上面两个包一样，需要说明是谁比谁
con <- paste(rev(levels(group)), collapse = "-")
con
# [1] "tumor-normal"

# 创建对比矩阵
cont.matrix <- makeContrasts(contrasts = c(con), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 获取差异表达基因结果
tempOutput <- topTable(fit2, coef = con, n = Inf)
DEG_limma_voom <- na.omit(tempOutput)

# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 1
FDR = 0.05
k1 <- (DEG_limma_voom$adj.P.Val < FDR) & (DEG_limma_voom$logFC < -logFC)
k2 <- (DEG_limma_voom$adj.P.Val < FDR) & (DEG_limma_voom$logFC > logFC)
DEG_limma_voom <- mutate(DEG_limma_voom, change = ifelse(k1, "down", ifelse(k2, "up", "no-change")))
table(DEG_limma_voom$change)

# 火山图
p <- ggplot(data = DEG_limma_voom, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p
```

Q3. For Graduate Students: From the DNA methylation profiles 
(GBM_meth.txt), how many genes are significantly differentially 
methylated between the two subtypes? Are DNA methylation associated with
 higher or lower expression of these genes? How many differentially 
expressed genes have an epigenetic (DNA methylation) cause?

Answers:

#不知道怎么比对样本呢

![](D:\software\MarkText\picture\2023-11-10-16-34-17-image.png)

Q4. With the survival data of the GBM tumors (GBM_clin.txt), make a 
Kaplan-Meier Curve to compare the two subtypes of GBM patients. Is there
 a significant difference in patient outcome between the two subtypes?

Answers:

```r
library("survival")
library("survminer")
cli$time <- as.numeric(ifelse(cli$vital.status==0,cli$days.to.last.followup,cli$days.to.death))
t <- gsub("\\.","-",colnames(exp))[1:60]
cli$group <- c(rep("group1",30),rep("group2",30))

fit <- survfit(Surv(time, vital.status) ~ group, data = cli)
summary(fit)
# 查看看完整的生存表格
summary(fit)$table

#按分层更改图形颜色，线型等
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # 添加风险表
           risk.table.col = "strata", # 根据分层更改风险表颜色
           linetype = "strata", # 根据分层更改线型
           surv.median.line = "hv", # 同时显示垂直和水平参考线
           ggtheme = theme_bw(), # 更改ggplot2的主题
           palette = c("#E7B800", "#2E9FDF"))#定义颜色
```

#怎么分组？

Q5. For Graduate Students: Use the differential genes (say this is Y 
number of genes) between the two GBM subtypes as a gene signature to do a
 Cox regression of the tumor samples. Does it give significant 
predictive power of patient outcome?

Answer:

```r
library("survival")
data(lung)
head(lung)
# inst time status age sex ph.ecog ph.karno pat.karno meal.cal wt.loss
# 1    3  306      2  74   1       1       90       100     1175      NA
# 2    3  455      2  68   1       0       90        90     1225      15

#分析离散型变量
res.cox <- coxph(Surv(time, status) ~ sex, data = lung)
#分析连续性变量
res.cox <- coxph(Surv(time, status) ~ age, data = lung)
#分析多个变量
res.cox <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data = lung)

summary(res.cox)
#coef就是偏回归系数，exp(coef)就是HR
```

Q6. For Graduate Students: Many studies use gene signatures to predict prognosis of patients. Take a look at this paper: [Most Random Gene Expression Signatures Are Significantly Associated with Breast Cancer Outcome | PLOS Computational Biology](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002240).
 It turns out that most published gene signatures are not significantly 
more associated with outcome than random predictors. Write a script to 
randomly sample Y genes in this expression data as a gene signature and 
do Cox regression on the sampled signature to predict patient outcome. 
Automate the script and random sample followed by Cox regression 100 
times. How does your signature in Q5 compared to random signatures in 
predicting outcome?

Answer:

```r
# 加载所需的库
library(survival)
library(plyr)

# 读取基因表达数据
gene_expression_data <- as.data.frame(t(read.table("GBM_expr.txt",header=T,sep="\t",row.names = 1)))
gene_expression_data$X <- gsub("\\.","-",rownames(gene_expression_data))
# 假设您的数据有一个名为 'SurvivalTime' 的时间列和一个名为 'Status' 的事件状态列
# 如果列名不同，请用实际的列名替换这些列名
survival_data <- as.data.frame(dplyr::left_join(cli[,c(1,2,6)],gene_expression_data,by="X"))

# 要生成的随机签名的数量
num_random_signatures <- 100

# 在每个签名中抽样的基因数量Y
num_genes_in_signature <- 13

# 执行Cox回归的函数
cox_regression <- function( signature_genes, survival_data) {
  signature_data <- survival_data[,c("time","vital.status",signature_genes)]
  cox_model <- coxph(Surv(time, vital.status) ~ ., data = signature_data)
  return(summary(cox_model)$logtest["pvalue"])
}

# 执行100次随机抽样和Cox回归
random_results <- replicate(num_random_signatures, {
  random_genes <- sample(colnames(survival_data)[-c(1,2)], num_genes_in_signature)
  cox_regression(random_genes, survival_data)
})

# 对您的基因签名执行Cox回归
your_genes <- c("BRCA1","BRCA2","TP53", "ERBB2","ESR1","PGR","MKI67","CHEK2","ATM","PALB2","PTEN","STK11","NBN")  # 用实际的基因名替换
your_result <- cox_regression(your_genes, survival_data)

# 将您的签名与随机签名进行比较
p_value_comparison <- sum(random_results <= your_result) / num_random_signatures

# 打印结果
cat("您的signatures的P值：", your_result, "\n")
#您的signatures的P值： 0.00194615 
cat("具有相等或较低P值的random signatures的比例：", p_value_comparison, "\n")                         
#具有相等或较低P值的random signatures的比例： 0.28 
```

## Part III – Tumor mutation analyses and precision medicine

Q1. The MAF files contain the mutations of each tumor compared to the
 normal DNA in the patient blood. Write a script to parse out the 
mutations present in each tumor sample and write out a table. The table 
should rank the genes by how many times mutation happens in the tumor 
samples provided. Submit the table with the top 20 genes.

![](D:\software\MarkText\picture\2023-11-15-20-08-11-image.png)

![](D:\software\MarkText\picture\2023-11-15-20-10-24-image.png)

```r
#编写一个脚本，解析每个肿瘤样本中出现的突变，并写出一个表格。该表应根据所提供的肿瘤样本中发生突变的次数对基因进行排序。提交包含前 20 个基因的表格。
library(data.table)
library(tidyverse)

setwd("/home4/liuyw/test/Liuxiaole_Training_homework_data/HW6")
maf_files <- list.files(pattern = "maf.txt")

maf_data <- read_tsv(maf_files[1], comment = "#", col_types = cols())
for (i in 2:length(maf_files)) {
  # 读取MAF文件
  f = fread(maf_files[i],sep="\t",header = T,fill = T)
  maf_data <- rbind(maf_data,f)
}

# 选择需要的列
selected_columns <- c("Hugo_Symbol", "Tumor_Sample_Barcode")
filtered_data <- maf_data %>% select(all_of(selected_columns))

# 过滤出具有意义的突变信息（排除空值等）
filtered_data <- filtered_data %>% drop_na()

# 计算每个基因的突变次数
gene_mutations_count <- filtered_data %>% count(Hugo_Symbol)

# 创建包含前20个基因的表格
top_genes_table <- gene_mutations_count %>% top_n(20, wt = n) %>% arrange(desc(n))

# 输出表格
print(top_genes_table)
# Hugo_Symbol     n
#   <chr>       <int>
# 1 TP53           16
# 2 IDH1           12
# 3 ATRX            9
# 4 NF1             7
# 5 EGFR            6
# 6 TTN             6
# 7 Unknown         6
# 8 NBPF10          5
# 9 APOB            4
# 10 MUC17          4
```

Q2. Existing clinical genetic testing laboratories use information 
about the frequency of a mutation in cohorts, like from the GBM cohort 
in TCGA, to assess a mutation’s clinical significance (guidelines: https://www.ncbi.nlm.nih.gov/pubmed/27993330).
 Of the top 20 mutated genes in Q1, what is the most frequent protein 
change (hint: count mutations with the exact same amino acid change as 
the same)? In what gene does the most frequent protein change occur? Do 
you think this gene and protein mutation pair forms a genetic subtype of
 GBM (i.e. Is this gene-protein mutation pair seen more in one GBM 
subtype than the other)?

Answer:

```r
#20 个突变基因中，最常见的蛋白质变化是什么?最常见的蛋白质变化发生在哪个基因中？
#您认为这种基因和蛋白质突变对是否形成了 GBM 的遗传亚型（即这种基因-蛋白质突变对是否在一种 GBM 亚型中比在另一种亚型中出现得更多）？
top20_maf <- maf_data %>% filter(Hugo_Symbol %in% top_genes_table$Hugo_Symbol) 

var_type <- top20_maf %>% count(Variant_Type) %>% top_n(3,wt=n) %>% arrange(desc(n))
# Variant_Type       n
# 1 SNP            128
# 2 DEL             17
# 3 INS              2
```

Q3. CBioPortal has a comprehensive list of tumor profiling results for interactive visualization. Go to cBioPortal ([http://www.cbioportal.org](http://www.cbioportal.org)),
 and select either “Glioblastoma” under “CNS/Brian” (left) or select 
“TCGA PanCancer Atlas Studies” under “Quick Select” (middle). Input each
 gene in Q1 and click Submit. From the OncoPrint tab, you can see how 
often each gene is mutated in GBM or all TCGA cancer types. Based on 
this, which of the genes in Q1 is likely to be a cancer driver gene?

Answer:

![](D:\software\MarkText\picture\2023-11-16-10-07-08-image.png)

Q4. From the Mutation tab on the cBioPortal result page, is this 
mutation a gain or loss of function mutation on the gene you identified 
from Q2?

Answer: 

![](D:\software\MarkText\picture\2023-11-16-10-22-44-image.png)

![](D:\software\MarkText\picture\2023-11-16-10-23-15-image.png)

Q5. From cBioPortal, select Glioblastoma (TCGA provisional, 
which has the largest number of samples) and enter the driver mutation 
gene in Q2. From the Survival tab, do GBM patients with this mutation 
have better outcome in terms of progression free survival and overall 
survival?

Answer:

![](D:\software\MarkText\picture\2023-11-16-10-32-57-image.png)

Q6. You are working with an oncologist collaborator to decide the 
treatment option for a GBM patient. From exome-seq of the tumor, you 
identified the top mutation in Q2. To find out whether there are drugs 
that can target this mutation to treat the cancer, go to [https://www.clinicaltrials.gov](https://www.clinicaltrials.gov) to find clinical trials that target the gene in Q2. How many trials are
 related to glioblastoma? How many of these are actively recruiting 
patients which this patient could potentially join? Hint: Search by the 
disease name and gene name.

Answer:

![](D:\software\MarkText\picture\2023-11-16-10-40-08-image.png)

## Part IV- CRISPR screens

We will learn to analyze CRISPR screen data from this paper: [Genome-wide CRISPR-Cas9 Screens Reveal Loss of Redundancy between PKMYT1 and WEE1 in Glioblastoma Stem-like Cells - PubMed](https://www.ncbi.nlm.nih.gov/pubmed/?term=26673326).
 To identify therapeutic targets for glioblastoma (GBM), the author 
performed genome-wide CRISPR-Cas9 knockout (KO) screens in 
patient-derived GBM stem-like cell line (GSCs0131).

MAGeCK tutorial: https://sourceforge.net/p/mageck/wiki/Home/ [MAGeCK download | SourceForge.net](https://sourceforge.net/projects/mageck/)

Q1. Use MAGeCK to do a basic QC of the CRISPR screen data (e.g. read 
mapping, ribosomal gene selection, replicate consistency, etc). Comment 
on the quality of the data based on your results.

Directory: /mnt/data/data/HW6/crispr_data

Answer:

```shell
# 从fastq文件中收集sgRNA read count信息
mageck count -l library.txt -n res --sample-label "Day0,Day23" --pdf-report --fastq Day0_rep1.fastq,Day0_rep2.fastq Day23_rep1.fastq,Day23_rep2.fastq 
# 比较样本间sgRNA和基因的差异
mageck test -k res.count.txt -t Day23 -c Day0 -n res
```

<img src="file:///D:/software/MarkText/picture/2023-11-16-20-33-22-image.png" title="" alt="" width="121"><img src="file:///D:/software/MarkText/picture/2023-11-16-20-33-47-image.png" title="" alt="" width="120"><img title="" src="file:///D:/software/MarkText/picture/2023-11-16-20-32-53-image.png" alt="" width="179">

Q2.  Analyze CRISPR screen data with MAGeCK to identify positive and negative
 selection genes. How many genes are selected as positive or negative 
selection genes, respectively, and what are their respective enriched 
pathways?

Answer:

- MAGeCK will generate output files containing results for each gene, including a ranking of genes based on their enrichment or depletion scores.

- Positive selection indicates genes that are enriched in the treatment group, suggesting a potential functional advantage or necessity.

- Negative selection indicates genes that are depleted in the treatment group, suggesting a potential inhibitory effect on cell survival or growth.

- **正向选择：**
  
  - 正向选择的基因是在处理组中敲除或抑制会给细胞提供生长或存活优势的基因。
  - 在CRISPR筛选中，正向选择的基因通常与细胞必需功能或通路相关，当其被干扰时，细胞获得了选择性优势。

- **负向选择：**
  
  - 负向选择的基因是在处理组中敲除或抑制导致细胞生长劣势或死亡的基因。
  - 负向选择通常指对细胞生存或增殖至关重要的基因。敲除这些基因会妨碍细胞的生长或存活。

```shell
# xxx.gene_summary.txt 文件
# 阳性筛选排序
sort -k 11,11n res.gene_summary.txt | less
```

<img title="" src="file:///D:/software/MarkText/picture/2023-11-16-20-37-20-image.png" alt="" width="316"><img title="" src="file:///D:/software/MarkText/picture/2023-11-16-20-39-23-image.png" alt="" width="324">

Q3. For Graduate Students: Genes negatively selected in this CRISPR 
screen could be potential drug targets. However, if they are always 
negatively selected in many cells, targeting such genes might create too
 much toxicity to the normal cells. Go to depmap (DepMap.org) which has 
CRISPR / RNAi screens of over 500 human cell lines, Click “Tools”  Data
 Explorer. Pick the top 3 negatively selected genes to explore. Select 
Gene Dependency from CRISPR (Avana) on the X axis and Omics from 
Expression on the Y axis, to see the relationship between the expression
 level of the gene and dependency (CRISPR screen selection) of the gene 
across ~500 cell lines. Are the top 3 genes good drug targets?

(在 CRISPR 筛选中被负选择的基因可能是潜在的药物靶点。但是，如果这些基因在许多细胞中总是被负选择，那么靶向这些基因可能会对正常细胞产生过多毒性。访问 depmap (DepMap.org)，其中有 500 多个人类细胞系的 CRISPR / RNAi 筛选结果，点击 "工具"  Data Explorer。选择前 3 个负向选择基因进行探索。在 X 轴上选择 CRISPR (Avana) 的基因依赖性，在 Y 轴上选择表达的 Omics，以查看约 500 个细胞系中基因表达水平与基因依赖性（CRISPR 筛选选择）之间的关系。前 3 个基因是好的药物靶点吗？)

Answer:

![](D:\software\MarkText\picture\2023-11-16-20-55-06-image.png)Q4. For Graduate Students: Let’s filter out pan essential genes 
(PanEssential.txt) from the negatively selected genes in Q2. Take the 
remaining top 10 genes, and check whether those genes have drugs or are 
druggable from this website: [http://www.oasis-genomics.org/](http://www.oasis-genomics.org/).
 Go to Analysis -> Pan Cancer Report, enter the top 10 genes and 
check the table for druggability (more druggable for higher number on 
Dr). Which of these genes are druggable?

(让我们从 Q2 中的负向选择基因中筛选出泛必需基因（PanEssential.txt）。取剩下的前 10 个基因，并从这个网站 http://www.oasis-genomics.org/ 上查看这些基因是否有药物或可以药物治疗。进入 "分析"->"泛癌症报告"，输入前 10 个基因，然后查看表中的可药性（博士编号越高，可药性越高）。这些基因中哪些是可药用的？)

Note: If the above website cannot be opened, try [https://www.dgidb.org/](https://www.dgidb.org/). Go to Search Potential Druggability or Search Drug-Gene Interactions ann enter the top 10 genes.

![](D:\software\MarkText\picture\2023-11-16-21-02-05-image.png)

## Part V. Cancer immunology and immunotherapy

Immune checkpoint inhibitors, which primarily activate CD8 T cells, 
have shown remarkable efficacy in melanoma (SKCM), but haven’t worked as
 well in GBM patients. Let’s explore the tumor immune microenvironment 
from TCGA data. Although the cancer patients in TCGA were not treated 
with immunotherapy, their response to other drugs and clinical outcome 
might be influenced by pre-treatment tumor immune microenvironment

Q1. TIMER ([http://timer.cistrome.org/](http://timer.cistrome.org/))
 estimated the infiltration level of different immune cells of TCGA 
tumors using different immune deconvolution methods. CD8A and CD8B are 
two gene markers on CD8 T cells. On the Diff Exp tab, compare the 
expression level of either CD8A or CD8B between GBM and SKCM (Metastatic
 Melanoma). Based on this, which cancer type have more CD8 T cells?

Answer:

<img src="file:///D:/software/MarkText/picture/2023-11-17-09-33-05-image.png" title="" alt="" width="262"><img src="file:///D:/software/MarkText/picture/2023-11-17-09-33-29-image.png" title="" alt="" width="264">

Q2. On the Gene tab, select both GBM and SKCM (Metastatic Melanoma), 
include CD8 T cells as the cell infiltrate. Check the following genes, 
PDCD1(PD1), CD274(PDL1), CTLA4 which are the targets of immune 
checkpoint inhibitors, to see whether their expression level is 
associated with immune cell infiltration in the GBM and SKCM tumors. 
Their higher expression usually indicate that T cells are in a 
dysfunctional state, which immune checkpoint inhibitors aim to revive.

Answer:

<img title="" src="file:///D:/software/MarkText/picture/2023-11-17-09-37-49-image.png" alt="" width="252"><img title="" src="file:///D:/software/MarkText/picture/2023-11-17-09-37-22-image.png" alt="" width="262">

Q3. On the Survival tab, select both GBM and SKCM, include CD8 T cell
 as the cell infiltrate, add tumor stage and patient age as the clinical
 variables to conduct survival analyses. Based on the Cox PH model, what
 factors are the most significantly associated with patient survival in 
each cancer type? Plot the Kaplan-Meier curve to evaluate how the immune
 cell infiltrate is associated with survival. Is CD8 T cell associated 
with patient survival in each cancer type?

GBM data has not stage information

Answer:

![](D:\software\MarkText\picture\2023-11-17-09-40-36-image.png)

Q4. For Graduate Students: Based on the above observations, can you 
hypothesize why immune checkpoint inhibitors don’t work well for some 
GBM patients?

Answer:
