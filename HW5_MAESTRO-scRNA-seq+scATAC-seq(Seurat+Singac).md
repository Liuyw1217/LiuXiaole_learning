---
title: "2022 Homework5 Single-cell Multiome Data Analysis"
author: "your name"
date: "Date"
output:
  html_document: default
  pdf_document: default
---

## Introduction

The technology of gathering data from multi-modality within the same cell offers opportunities for gaining holistic views of cells individually. 10X [Chromium Single-cell Multiome ATAC + Gene Expression](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression) simultaneously profiling transcriptome and epigenome at single-cell level, enabled deeper characterization of cell types and states.
In addition to their power in addressing biological questions, single-cell multiome also provides good testing data for evaluating label transferring accuracy. In a single omic scATAC-seq analysis, people typically transfer cell type labels from existing scRNA-seq data to obtain the cell type annotation. Leveraging multiome ATAC + Gene Expression, we can link transcriptome and epigenome modalities cell by cell. Thereby, we know the "ground truth" cell type label for each cell in scATAC-seq data and compare the results from Seurat's label transferring methods.
在同一细胞内收集多模态数据的技术为获得细胞个体的整体视图提供了机会。10X [Chromium Single-cell Multiome ATAC + Gene Expression](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression)在单细胞水平上同时分析转录组和表观基因组，能够更深入地描述细胞类型和状态。
单细胞多组学除了能解决生物学问题外，还为评估标签转移的准确性提供了良好的测试数据。在单个 omic scATAC-seq 分析中，人们通常从已有的 scRNA-seq 数据中转移细胞类型标签，以获得细胞类型注释。利用多组 ATAC + 基因表达，我们可以将转录组和表观组模式逐个细胞联系起来。这样，我们就能知道 scATAC-seq 数据中每个细胞的 "地面实况 "细胞类型标签，并比较 Seurat 标签转移方法的结果。

In this homework, we will start from the 10X fastq raw reads, go through pre-processing pipeline and downstream analysis to get an overall idea of single-cell data analysis. For the data pre-processing, we will use the <mark>MAESTRO</mark> package installed on Huawei cloud. For downstream analysis, we will work on your local computer's Rstudio using<mark> [Seurat](https://satijalab.org/seurat/index.html) </mark>and<mark> [Signac]</mark>(https://satijalab.org/signac/index.html) R packages.

## Homework Questions

### Part I. Running MAESTRO

- Understand MAESTRO workflow

- Set up MAESTRO analysis pipelines
  
  #### Step0. Get 10X pbmc_granulocyte_sorted_3k multiome data
  
  Single-cell multiome data were downloaded from the 10X website. You can find raw data [here](https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_3k). Data were downloaded to the Huawei cloud. There are two data folders under `/mnt/data/data/HW5/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k/` directory. One is `gex/` for scRNA-seq data and the other one is `atac/` for scATAC-seq data. You will need to activate the conda environment to access all the required packages: `$ source activate MAESTRO`. If you want to learn more about conda environment management, you can read through the link [conda](https://conda.io/en/latest/index.html).
  
  #### Step1. Configure MAESTRO working directory.
  
  MAESTRO is a snakemake pipeline developed for streamlined pre-processing of single-cell data and downstream analysis. For more details, you can read through the documentation of [MAESTRO](https://github.com/liulab-dfci/MAESTRO). MAESTRO is written in Snakemake. Please learn more in the [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) documentation.
  This step will configure two working directories for scRNA-seq and scATAC-seq, respectively. A `Snakefile` and a `config.yaml` file will be configured within each directory. `Snakefile` includes all the snakemake rules that will be later executed. `config.yaml` file contains the detailed parameter settings you need to specify when initiating the directory. You can also manually change the contents of the config.yaml file to customize your run. Let's call the directory for processing scRNA-seq data as `multiome_scrna/` and scATAC-seq data as `multiome_scatac/`.
  此步骤将分别为 scRNA-seq 和 scATAC-seq 配置两个工作目录。每个目录中都将配置一个 Snakefile 和一个 config.yaml 文件。Snakefile 包含稍后将执行的所有 snakemake 规则。config.yml 文件包含启动目录时需要指定的详细参数设置。你也可以手动更改 config.yaml 文件的内容，自定义运行。我们将处理 scRNA-seq 数据的目录称为 multiome_scrna/，将处理 scATAC-seq 数据的目录称为 multiome_scatac/。
  
  **Please read the detailed tutorials in MAESTRO documentations to have a better idea of how to set each parameter:**
  [scRNA-seq](https://github.com/liulab-dfci/MAESTRO/blob/master/example/RNA_infrastructure_10x/RNA_infrastructure_10x.md)
  [scATAC-seq](https://github.com/liulab-dfci/MAESTRO/blob/master/example/ATAC_infrastructure_10x/ATAC_infrastructure_10x.md)
  **Before running, Make sure you already have the below files on
  your server.**
  
  For scRNA-seq:
  
  1. STAR mapping index file: This is the STAR genome reference file for
     mapping human single-cell RNA-seq data.
  - Path: `/mnt/data/data/HW5/references/Refdata_scRNA_MAESTRO_GRCh38_1.2.2`
  2. barcode whitelist: This is the complete list of multiome scRNA-seq
     cell barcodes from the 10X Cell Ranger ARC workflow. We use the barcode
     whitelist to correct the cell barcodes we get from reads.
  - Path: `/mnt/data/data/HW5/references/whitelist/rna/737K-arc-v1.txt`
  3. lisa TF annotation file: This file contains data used for running
     LISA2.
  - Path: `/mnt/data/data/HW5/references/lisa_data/hg38_1000_2.0.h5`
  
  For scATAC-seq:
  
  1. giggleannotation file: giggle annotation file is required for
     regulator identification.
  - Path: `/mnt/data/data/HW5/references/giggle.all`
  2. minimap2 reference file: This is the STAR genome reference file for
     mapping human single-cell ATAC-seq data.
  - Path: `/mnt/data/data/HW5/references/Refdata_scATAC_MAESTRO_GRCh38_1.1.0`
  3. barcode whitelist: This is the complete list of multiome scATAC-seq
     cell barcodes from the 10X Cell Ranger ARC workflow. We use the barcode
     whitelist to correct the cell barcodes we get from reads.
  - Path: `/mnt/data/data/HW5/references/whitelist/atac/737K-arc-v1.txt`
  
  **hints:**
  
  1. MAESTRO has <mark>two sub-commands </mark>to initiate the working
     directories.
  
  2. You will need to feed each parameter settings to each
     sub-commands for initialization.

#### Step2. Run the snakemake pipeline.

**hint**: Estimated running time and memory usage:
3k multiome scrna-seq: 120min; 60G; 16 cores
3k multiome scatac-seq: 180min; 60G; 16 cores

```shell
##安装 MAESTRO
conda config --add channels defaults
conda config --add channels liulab-dfci
conda config --add channels bioconda
conda config --add channels conda-forge

# To make the installation faster, we recommend using mamba
conda install mamba -c conda-forge
mamba create -n MAESTRO-2 maestro=1.5.0 -c liulab-dfci

# Activate the environment
conda activate MAESTRO
```

```shell
##--1-- 10X PBMC 8k scRNA-seq
# Step 1. Configure the MAESTRO workflow
# Initialize the MAESTRO scRNA-seq workflow using MAESTRO scrna-init command. This will install a Snakefile and a config file in this directory.

#(1)MAESTRO scrna-init
MAESTRO scrna-init --platform 10x-genomics --species GRCh38 \
--cores 16 --rseqc --directory multiome_scrna --outprefix scrna_pbmc_granulocyte_sorted_3k \
--mapindex /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/references/Refdata_scRNA_MAESTRO_GRCh38_1.2.2/GRCh38_STAR_2.7.6a \
--whitelist /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/references/whitelist/rna/737K-arc-v1.txt \
--umi-length 12 --lisadir /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/references/lisa_data/hg38_1000_2.0.h5 --signature human.immune.CIBERSORT

#(2)MAESTRO samples-init
cd multiome_scrna
MAESTRO samples-init --assay_type scrna --platform 10x-genomics --data_type fastq --data_dir /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k/gex

# Step 2. Run MAESTRO
# Before running the workflow, please check the config.yaml and see if it is configured correctly. Once config.yaml is configured, users can use snakemake to run the workflow.
#First, test with a dry run to see if the pipeline work
##Remember to add -np for a DRY run:
cd multiome_scrna
#First, test with a dry run to see if the pipeline work
##Remember to add -np for a DRY run:
snakemake -np --rerun-incomplete -j 1
#If a dry run is 100% complete and no red error messages reported, we will create a sbatch job script under the `multiome_scrna/` folder
#-j 16 is the total number of cores you will need to specify when creating the job.
nohup snakemake --rerun-incomplete -j 16 > 10X_PBMC_3k.out &


```

```shell
##--2-- 10x PBMC 10k scATAC-seq
chromap -i -r /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/references/Refdata_scATAC_MAESTRO_GRCh38_1.1.0/GRCh38_genome.fa -o index
# chromap buildindex -r /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/references/Refdata_scATAC_MAESTRO_GRCh38_1.1.0/GRCh38_genome.fa -o index

# Step 1. Configure the MAESTRO workflow
# Initialize the MAESTRO scATAC-seq workflow using MAESTRO scATAC-init command. 
# scATAC-seq initiation:
MAESTRO scatac-init --platform 10x-genomics --format fastq --species GRCh38 --input_path /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k/atac \
--cores 16 --directory multiome_scatac \
--peak_cutoff 100 --count_cutoff 1000 --frip_cutoff 0.2 --cell_cutoff 50 \
--giggleannotation /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/references/giggle.all \
--fasta /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/references/Refdata_scATAC_MAESTRO_GRCh38_1.1.0/GRCh38_genome.fa \
--whitelist /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/references/whitelist/atac/737K-arc-v1.txt --rpmodel Enhanced --annotation --method RP-based --signature human.immune.CIBERSORT \
--index index

#samples-init
cd multiome_scatac
MAESTRO samples-init --assay_type scatac --platform 10x-genomics --data_type fastq --data_dir /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/pbmc_granulocyte_sorted_3k//pbmc_granulocyte_sorted_3k/atac
#dry run
snakemake -np --rerun-incomplete -j 1
#run
nohup snakemake --rerun-incomplete -j 16 > 10X_PBMC_3k.out &
```

```shell
#用minimap2进行比对，添加参数--mapping minimap2
# scATAC-seq initiation:
MAESTRO scatac-init --platform 10x-genomics --format fastq --species GRCh38 --input_path /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k/atac \
--cores 16 --directory multiome_scatac \
--peak_cutoff 100 --count_cutoff 1000 --frip_cutoff 0.2 --cell_cutoff 50 \
--giggleannotation /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/references/giggle.all \
--fasta /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/references/Refdata_scATAC_MAESTRO_GRCh38_1.1.0/GRCh38_genome.fa \
--whitelist /home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data/HW5/references/whitelist/atac/737K-arc-v1.txt --rpmodel Enhanced --annotation --method RP-based --signature human.immune.CIBERSORT \
--mapping minimap2


```

**1 Reads mapping stats: After getting the `Result/` output folder, for scRNA-seq run, please copy the content of `/mnt/data/data/HW5/pbmc_granulocyte_sorted_3k/multiome_scrna/Result/QC/scrna_pbmc_granulocyte_sorted_3k_bam_stat.txt` below (1 pts;). For the scATAC-seq run, please copy the content of `/mnt/data/data/HW5/pbmc_granulocyte_sorted_3k/multiome_scatac/Result/QC/flagstat.txt` below (1 pts;).**

![](D:\software\MarkText\picture\2023-11-21-09-55-42-image.png)

![](D:\software\MarkText\picture\2023-11-21-09-57-30-image.png)

**2. Cell Filtering QC plot: Please attach the cell filtering QC plot for the scATAC-seq run, which is saved as `/mnt/data/data/HW5/pbmc_granulocyte_sorted_3k/multiome_scatac/Result/QC/scatac_pbmc_granulocyte_sorted_3k_scATAC_read_distr.png` (1 pts;). Please briefly describe what you observed from this figure and how the filtering was affected by the cutoff value you set in the `MAESTRO scatac-init` subcommand (1 pts;).**

师姐的结果，我的结果里没有。。。。

![](D:\software\MarkText\picture\2023-11-21-09-59-16-image.png)

Cell Filtering QC plot

####Answer

**3. Bonus Question: Suppose you have a list of 500 cell barcodes, can you sub-sample the scATAC-seq data by cell barcodes to get a raw fastq file with only 500 cells? Please provide the path to downsampled fastq file on Huawei cloud (3 pts;).**
**hints:**

1. After running MAESTRO, cells passed the QC filter will be stored in `/mnt/data/data/HW5/pbmc_granulocyte_sorted_3k/multiome_scatac/Result/QC/scatac_pbmc_granulocyte_sorted_3k_scATAC_validcells.txt`. You can sub-sample 500 cell barcodes using the command `shuf -n 500 scatac_pbmc_granulocyte_sorted_3k_scATAC_validcells.txt > 500_vallidcells.txt`.

2. Sinto has a function called [filterbarcodes](https://timoast.github.io/sinto/basic_usage.html#filter-cell-barcodes-from-bam-file) that can sub-sample bam file according to the cell barcodes.

3. There are several tools that can convert bam files back to fastq. 10X has a script called [bamtofastq](https://support.10xgenomics.com/docs/bamtofastq). Another tool [bam2fasta](https://github.com/czbiohub/bam2fasta) have similar functions.

4. You don't have to start from bam files. Any methods are appreciated.

```shell
# how you can sub-sample scATAC-seq data to get a raw fastq file with only 500 cells:
#师姐答案，这个没有用minimap2
#Write a random permutation of the input lines to standard output.
shuf -n 500 /mnt/data/data/HW5/pbmc_granulocyte_sorted_3k/multiome_scatac/Result/QC/scatac_pbmc_granulocyte_sorted_3k_scATAC_validcells.txt > /mnt/data/home/yuyuanyuan/traning_homework/week5/part1/sub_sample/500_vallidcells.txt

#给barcode添加分组
awk 'BEGIN{OFS="\t"} {print $0,"scatac_pbmc_granulocyte_sorted_3k.sortedByPos.rmdp.CBadded.bam"}' /mnt/data/home/yuyuanyuan/traning_homework/week5/part1/sub_sample/500_vallidcells.txt > /mnt/data/home/yuyuanyuan/traning_homework/week5/part1/sub_sample/500_vallidcells_tag.txt

ln -s /mnt/data/data/HW5/pbmc_granulocyte_sorted_3k/multiome_scatac/Result/minimap2/scatac_pbmc_granulocyte_sorted_3k.sortedByPos.rmdp.CBadded.bam ./

samtools index scatac_pbmc_granulocyte_sorted_3k.sortedByPos.rmdp.CBadded.bam

#sub-sample bam file according to the cell barcodes
sinto filterbarcodes -b /mnt/data/home/yuyuanyuan/traning_homework/week5/part1/sub_sample/scatac_pbmc_granulocyte_sorted_3k.sortedByPos.rmdp.CBadded.bam -c /mnt/data/home/yuyuanyuan/traning_homework/week5/part1/sub_sample/500_vallidcells_tag.txt -p 20

samtools index scatac_pbmc_granulocyte_sorted_3k.sortedByPos.rmdp.CBadded.bam.bam

#转成fastq文件
bedtools bamtofastq -i /mnt/data/home/yuyuanyuan/traning_homework/week5/part1/sub_sample/scatac_pbmc_granulocyte_sorted_3k.sortedByPos.rmdp.CBadded.bam.bam -fq /mnt/data/home/yuyuanyuan/traning_homework/week5/part1/sub_sample/scatac_pbmc_granulocyte_sorted_3k.sortedByPos.rmdp.CBadded.bam.fastqfastq


```

### Part II. Single-cell RNA-seq

- Create a scRNA-seq Seurat Object
- QC and Pre-processing
- Normalization and Dimension Reduction
- KNN and Clustering
- Finding Markers
- Cell Type Annotation

Now, we are done with all the work on the server. In this part, we’ll continue to work on single-cell Seurat objects generated by MAESTRO. After running MAESTRO, all the output files will be organized under a data folder named `Result/`. A Seurat object containing the count matrix and the metadata will be stored in a .rds data list. You
can always find them in the MAESTRO analysis results: `multiome_scrna/Result/Analysis` and `multiome_scatac/Result/Analysis`. In a real data analysis workflow, you will continue to work on the .rds file you got from
MAESTRO. But in this homework, let’s use the prepared Seurat objects to keep downstream analysis consistent.

Please go to your local computer and download two .rds file from `/mnt/data/data/HW5/pbmc_granulocyte_sorted_3k`. The single-cell RNA-seq object is `/mnt/data/data/HW5/pbmc_granulocyte_sorted_3k/multiome_scrna/Result/Analysis/scrna_pbmc_granulocyte_sorted_3k_scRNA_Object.rds` and the single-cell ATAC-seq object is stored as `/mnt/data/data/HW5/pbmc_granulocyte_sorted_3k/multiome_scatac/Result/Analysis/scatac_pbmc_granulocyte_sorted_3k_scATAC_Object.rds`.

**Please install and load the following packages:**

```r
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
```

**1. In Rstudio, you can use `readRDS()` to load the .rds data. You will find that the .rds is a data list containing an `RNA` Seurat object and a differentially expressed gene list. We will only need to extract the Seurat object for further analysis. Please Describe the raw dataset’s composition: what are the number of
genes and number of cells in your cell feature matrix (1 pts;) ?**

```r
setwd("E:/learnimg/JLab/xiaoleLiu/week5/data/scRNA")
##Read the .rds file for scRNA-seq data:
rna<- readRDS("pbmc_granulocyte_sorted_3k_scRNA_Object.rds")

#MAESTRO also attached differentially expressed genes as a table under the .rds file. We need to extract the Seurat object.
rna<- rna$RNA  #Seurat object
rna
```

![](D:\software\MarkText\picture\2023-11-21-10-23-52-image.png)

**2. Filtering cells with a high proportion of mitochondrial reads (potential dead cells) or outlier number of genes (possible low reactions or multiplets) are essential steps in single-cell analysis. Outlier cells with too high or low gene coverage should be removed. The cutoff depends on the scRNA-seq technology and the distribution of each dataset. MAESTRO has already filtered the cells based on the number of counts and genes. In this question, please calculate the percentage of UMIs mapped to the mitochondrial genes, save the results under `percent.mt` column in the Seurat `@metadata` slot (1 pts;). Please visualize the distribution of `nFeature_RNA`, `nCount_RNA` and `percent.mt` in a single violin plot (1 pts;). hints: You may want to use `Idents(rna) <- rna$orig.ident` before running violin plot for better visualization. Use a table to show how many of
the cells have mitochondrial rate > 20% (1 pts;).**

```r
set.seed(123)
#--1--# QC 
Idents(rna) <- rna$orig.ident
#Idents() 函数用于为 Seurat 对象中的细胞分配标识符（identifiers）。rna$orig.ident 是一个包含细胞原始标识符的向量。每个元素对应于 rna 数据集中的一个细胞，并包含该细胞的原始标识符。
#Idents(rna) <- rna$orig.ident 将细胞的标识符设置为 rna$orig.ident 中包含的值。这样，Idents() 函数将 Seurat 对象中的细胞标识符更新为原始标识符向量中的值。
#此操作可能有助于将外部信息（如样本来源、处理组、实验条件等）与 Seurat 对象中的细胞关联起来。

#计算线粒体reads比例
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")
#visualize the distribution of `nFeature_RNA`, `nCount_RNA` and `percent.mt` in a single violin plot
VlnPlot(rna,features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),pt.size = 0)
#cells have mitochondrial rate \> 20%
mt_cell <- subset(rna, subset = percent.mt > 20)
mt_cell
```

![](D:\software\MarkText\picture\2023-11-21-10-25-43-image.png)

**3. After removing unwanted cells from the dataset (already done by MAESTRO. No need to filter in this homework), the next step is to normalize the data. Please use the default settings in Seurat to do the normalization: `NormalizeData()` (1 pts;). We next calculate a subset of features that exhibit high cell-to-cell variation
in the data set. These features will be used for downstream analysis to reduce computing time. Please use `FindVariableFeatures(..., selection.method = "vst", nfeatures = 2000)` to return the top 2,000 variable features (1 pts;). Next, we will apply a linear transformation (also known as scaling), a standard pre-processing step prior to dimensional reduction techniques. Please perform scaling on the top 2,000 variable genes (default) using `ScaleData()` (1 pts;).**

```r
#--2--# normalization
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
#FindVariableFeatures
rna = FindVariableFeatures(rna , selection.method = "vst", nfeatures = 2000)
VariableFeatures(rna)[1:100]
#Scale top 2,000 variable genes
rna <- ScaleData(rna)# 默认是仅在高可变基因上运行标准化
```

![](D:\software\MarkText\picture\2023-11-21-10-28-54-image.png)

**4. Perform linear dimensional reduction. Next, we’ll perform PCA on the scaled data. Please use `RunPCA()` to do the dimension reduction on the 2,000 variable genes we detected in the previous question (1 pts;). You can take a look at PCA cell embeddingsat `rna[['pca']]@cell.embeddings`.**

```r
#--3--# reduction
rna <- RunPCA(rna,features=VariableFeatures(object=rna))#基于前面得到的high variable基因的scale矩阵
#You can take a look at PCA cell embeddings at rna[['pca']]@cell.embeddings
rna[['pca']]@cell.embeddings[1:20,1:5]

DimPlot(rna,reduction = "pca")
ElbowPlot(rna,ndims = 50)

```

![](D:\software\MarkText\picture\2023-11-21-10-29-02-image.png)

**5. Since not all the PCs we calculated will be used in downstream analysis. In this question, We will determine how many PCs we should use in downstream analysis:**

- **5.1 How much variability is explained in each of the first 50PCs? Please show a scree plot with each PCs on the x-axis and variation explained by each PC on the y-axis (1 pts;).**
- **hint** : [Is it possible to figure out PC variance explained in Seurat? · Issue #982 · satijalab/seurat · GitHub](https://github.com/satijalab/seurat/issues/982)

```
#Please show a scree plot with each PCs on the x-axis and variation explained by each PC on the y-axis
mat <- GetAssayData(rna, assay = "RNA", slot = "scale.data")  #大matrix，2000特征 2633细胞
pca <- rna[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat)) #matrixStats包提供的函数之一。该函数用于计算矩阵（或数据框）的每行方差。

eigValues = (pca@stdev)^2  ## EigenValues,特征值，主成分的特征值，特征值衡量了该主成分所解释的方差的大小
varExplained = eigValues / total_variance

qplot(c(1:50), varExplained) +
  geom_line() +
  geom_point(size=4)+
  xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle("Scree Plot") + theme_bw()
```

![](D:\software\MarkText\picture\2023-11-21-10-31-50-image.png)

- **5.2 How many PCs do you need to cover 20% of the total variance? Show a bar plot with PCs on the x-axis and the cumulative sum of variance explained by each PC on the y-axis (1 pts;). In bulk RNA-seq experiments, can you recall how many PCs you need to explain 20% of the variability? What do you think is the main reason that causes such difference between single-cell and bulk RNA-seq (Graduate level 2 pts;)?**

```r
varExplained_cusum=cumsum (varExplained)
varExplained_cusum
```

```
##  [1] 0.07265365 0.09341231 0.11181564 0.12654058 0.13800581 0.14585500
##  [7] 0.15040909 0.15459109 0.15802056 0.16094504 0.16358224 0.16603339
## [13] 0.16827453 0.17041856 0.17248729 0.17437511 0.17622072 0.17802903
## [19] 0.17979082 0.18151151 0.18321558 0.18489700 0.18656391 0.18822228
## [25] 0.18986865 0.19150576 0.19313067 0.19474103 0.19634057 0.19793408
## [31] 0.19952321 0.20110863 0.20268287 0.20425691 0.20582470 0.20738932
## [37] 0.20894354 0.21049373 0.21203397 0.21356981 0.21509881 0.21662339
## [43] 0.21814315 0.21966091 0.22117447 0.22267661 0.22417505 0.22567048
## [49] 0.22716386 0.22865281
```

```r
barplot(height=varExplained_cusum,names.arg =as.character(c(1:50)), 
        xlab = "Principal Component",
        ylab = "cumulative sum of variance Explained",
        ylim = c(0,0.3))+
  abline(h=0.2,lwd=2,col="red", lty = 3)
```

![](D:\software\MarkText\picture\2023-11-21-10-34-18-image.png)

####Answer

```
I need 32 PCs to cover 20% of the total variance.
In bulk RNA-seq experiments, maybe just 1 PC can explain 20% of the variability. I think that the number of cells is often much larger than the number of samples is the main reason that causes such difference between single-cell and bulk RNA-seq.
```

- **5.3 What are The top 5 genes with the most positive and negative coefficients in each of the first 10 PCs? Use a table to show the results (1 pts;).**

```r
# 然而，在bulk RNA-seq中，通常需要更少的主成分来捕捉与scRNA-seq相似数量的变异性。这种差异的主要原因是scRNA-seq数据中存在的噪音和异质性。
# 在scRNA-seq中，每个细胞被视为一个独立的样本，由于细胞之间的差异、技术噪声和其他因素，导致变异性增加。另一方面，批量RNA-seq对细胞群体中的表达值进行平均，减少了个体细胞变异的影响，产生了更平滑的信号。

#What are The top 5 genes with the most positive and negative coefficients in each of the first 10 PCs? Use a table to show the results (1 pts;).
print(rna[["pca"]], dims = 1:10, nfeatures = 5)
```

**6. Umap Visualization: Dimensionality reduction is a powerful tool for machine learning practitioners to visualize and understand large, high-dimensional datasets. Seurat offers several non-linear dimension reduction techniques, such as tSNE and UMAP (as opposed to PCA which is a linear dimensional reduction technique). UMAP is a new technique by McInnes et al. that offers many advantages over t-SNE, most notably increased speed and better preservation of the data’s global structure.**

- **6.1 Use the first 15 PCs, apply `RunUMAP()` to perform non-linear dimension reduction for visualization. The results will be stored in `rna[['umap]]` (1 pts;).**

```r
#--4--Umap Visualization
rna <- RunUMAP(object = rna, dims = 1:15)
rna[['umap']]
# DimPlot(rna,reduction = "umap")
```

```
## A dimensional reduction object with key UMAP_ 
##  Number of dimensions: 2 
##  Projected dimensional reduction calculated:  FALSE 
##  Jackstraw run: FALSE 
##  Computed using assay: RNA
```

- **6.2 Please Visualize the cells on the PCA and UMAP embeddings individually (2 pts;) and comment on the number of cell clusters that appear in each plot (1 pts;). hint: Use `DimPlot()`. Describe the difference between PCA and UMAP on 2D plots (2 pts;).**

```r
DimPlot(rna,reduction = "pca")
DimPlot(rna,reduction = "umap")
```

<img title="" src="file:///D:/software/MarkText/picture/2023-11-21-10-38-18-image.png" alt="" width="321"><img title="" src="file:///D:/software/MarkText/picture/2023-11-21-10-38-40-image.png" alt="" width="317">

####Answer

```
There are 2 cell clusters that appear in pca plot; and there are 3 big cell clusters that appear in umap plot. 
The PCA projection looks for the direction that maximizes the variance, and usually ignores the variance in the other directions. Therefore, PCA tends to find larger cell clusters, while ignoring smaller ones. However, the idea of UMAP is different from PCA, which will find more clusters, and try to separate these clusters.
```

**7. Clustering: Seurat v3 applies a graph-based clustering approach, building upon initial strategies from Macosko et al.. To cluster the cells, we first construct a KNN graph based on the euclidean distance in PCA space and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity).**

- **7.1 Use the `FindNeighbors()` function and take first 15 PCs as input to perform KNN (1 pts;).**

```r
#--5--Clustering
# Use the FindNeighbors() function and take first 15 PCs as input to perform KNN 
rna <- FindNeighbors(rna, dims = 1:15)
```

- **7.2 We next apply modularity optimization techniques to iteratively group cells together. Use `FindClusters()` function with different `resolution` to perform clustering and draw the resulting clusters in different colors on UMAP (`resolution` = 0.4, 0.6, 0.8) (3 pts;).**

```r
#resolution越大cluster越多。Use FindClusters() function with different resolution to perform clustering and draw the resulting clusters in different colors on UMAP (resolution = 0.4, 0.6, 0.8) 
rna <- FindClusters(rna, resolution = 0.4)
rna <- FindClusters(rna, resolution=0.6)
rna <- FindClusters(rna, resolution=0.8)
DimPlot(rna, reduction = "umap", group.by=c('RNA_snn_res.0.4','RNA_snn_res.0.6','RNA_snn_res.0.8'))
```

![](D:\software\MarkText\picture\2023-11-21-10-41-39-image.png)

- **7.3 How does resolution influence the number of clusters and the number of cells assigned to each cluster? Please provide a table to show the number of cells in each cluster (Graduate 1 pts;). Is there a correct number of clusters in a particular data set? why or why not (Graduate level 1pts;)?**

```r
#the number of cells in each cluster
table(rna$RNA_snn_res.0.4)
table(rna$RNA_snn_res.0.6)
table(rna$RNA_snn_res.0.8)
```

**8. Find Cluster Markers: For further analysis, please keep using resolution = 0.6 to cluster the cells. Seurat has a function called `FindMarkers()` that can perform several tests to identify differentially expressed genes for a single cluster.**

- **8.1 Use Wilcox Rank Sum test (default) in `FindMarkers(..., min.pct = 0.25)` to identify all markers of
  cluster 5 (1 pts;). Print the top5 markers in cluster 5. hints: You need
  to rerun `FindClusters()` with resolution = 0.6 to get
  correct number of cell clusters.**

```r
#--6-- Find Cluster Markers
#Use Wilcox Rank Sum test (default) in FindMarkers(..., min.pct = 0.25) to identify all markers of cluster 5 (1 pts;). Print the top5 markers in cluster 5. 
rna <- FindClusters(rna, resolution = 0.6)
Idents(rna) <- rna$RNA_snn_res.0.6

cluster5.markers <- FindMarkers(rna, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 5)
```

```
##               p_val avg_log2FC pct.1 pct.2     p_val_adj
## CCL5  5.225108e-219   3.770360 0.978 0.118 7.446301e-215
## GZMH  4.429950e-186   3.017419 0.617 0.035 6.313122e-182
## GZMA  5.284179e-185   2.528293 0.850 0.086 7.530484e-181
## NKG7  5.037707e-173   3.103419 0.950 0.145 7.179237e-169
## TRGC2 3.764142e-162   2.364719 0.567 0.035 5.364279e-158
```

- **8.2 `FindAllMarkers()` function can find markers for every cluster compared to all remaining cells (one vs. the rest). Apply `FindAllMarkers()` function with `min.pct = 0.25` and `logfc.threshold = 0.25` to report only positive genes in each cluster (1 pts;). Please print top 2 markers in each cluster weighted by log2FC (1pts;).**

```r
#FindAllMarkers() function can find markers for every cluster compared to all remaining cells (one vs. the rest). Apply FindAllMarkers() function with min.pct = 0.25 and logfc.threshold = 0.25 to report only positive genes in each cluster (1 pts;). 
#Please print top 2 markers in each cluster weighted by log2FC
rna.markers <- FindAllMarkers(rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
rna.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) #wt指定用于排序的权重列


```

- **8.3 Visualize the gene expression values of these potential markers (top2) on your UMAP plots using `FeaturePlot()` (2 pts;).**

```r
gene.markers <- rna.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

FeaturePlot(rna, features = gene.markers$gene)
```

![](D:\software\MarkText\picture\2023-11-21-10-50-03-image.png)

**9. Annotation: Based on markers in each cluster, MAESTRO referenced human.immune.CIBERSORT data to do the cell type annotations. The annotated cell types are stored in the `assign.ident` column in Seurat metadata. Please set `Idents(rna) <- rna$assign.ident` and plot the cell type annotation results in UMAP (1 pts;)**

```r
#--7--Annotation
# Based on markers in each cluster, MAESTRO referenced human.immune.CIBERSORT data to do the cell type annotations. The annotated cell types are stored in the assign.ident column in Seurat metadata. 
# Please set Idents(rna) <- rna$assign.ident and plot the cell type annotation results in UMAP
Idents(rna) <- rna$assign.ident
DimPlot(rna, reduction = "umap",label = TRUE, pt.size = 0.5)

```

![](D:\software\MarkText\picture\2023-11-21-10-50-51-image.png)

### Part III. Single-cell ATAC-seq

- Create a scATAC-seq Seurat Object
- Data Normalization and Dimension Reduction
- Clustering
- Cell Type Annotation
- Finding Differentially Accessible Peaks

**Please install and load following the packages:**

```
library(data.table)
library(dplyr)
library(Seurat)
library(Signac)
library(ggplot2)
```

**1. Read the multiome scATAC-seq Seurat object: Same as scRNA-seq .rds file, we need to extract Seurat object from .rds data list. How many assays are stored in this Seurat object and what is the current active assay (1 pts;)? Now switch the default assay to ‘ATAC’ using `DefaultAssay(atac) <- 'ATAC'`. How many cells and peaks are retained in the peak count matrix (1 pts;)? Important note:
Please use `Idents(atac) <- atac$orig.ident` to update the cell ID before moving to Question 2.**

```r
setwd("E:/learnimg/JLab/xiaoleLiu/week5/data/scATAC")
set.seed(123)

#--1-- Create a scATAC-seq Seurat Object
atac <- readRDS("pbmc_granulocyte_sorted_3k_scATAC_Object.rds")
str(atac)
atac <- atac$ATAC
atac
Idents(atac) <- atac$orig.ident  #update the cell ID
#将对象 atac 中的默认测量数据设置为名为 'ATAC' 的测量数据。这样，在后续的分析中，可以通过使用默认测量数据的名称来引用和访问 'ATAC' 数据。
DefaultAssay(atac) <- 'ATAC'   #switch the default assay to ‘ATAC’ 
# In .rds data, there are 2 list —— ATAC(SeuratObject) and peaks(data.frame). In Seurat object, there are 2 assays —— "ATAC" and "ACTIVITY","ACTIVITY" is the current active assay. 
# There are 44292 cells and 2664 peaks are retained in the peak count matrix.

```

**2. Normalization and Linear Dimensional Reduction: scATACseq data are very sparse. It is sparser than scRNAseq. Thus, It is necessary to do pre-processing steps before clustering scATACseq data. Signac performs a two-step normalization method called [TF-IDF](https://divingintogeneticsandgenomics.rbind.io/post/clustering-scatacseq-data-the-tf-idf-way/) (term frequency-inverse document frequency). To get rid of the low dynamic range of scATAC-seq data, we further use `FindTopFeatures()` to only choose the top n% of features (peaks) for dimension reduction. Next, we run SVD (singular value decomposition) on the normalized TD-IDF matrix for linear dimensional reduction. The combined steps of TF-IDF followed by SVD are also known as LSI (latent semantic indexing).**

- **2.1: Before performing normalization, let’s plot a UMAP on the first 15 PCs to see how the plot looks(1 pts;). hints: Same as we have done in scRNA-seq analysis, do `ScaleData()` and `RunPCA()` (Running PCA will take some time), then `RunUMAP(..., reduction = 'pca', dims =1:15)`. Plot the UMAP.
  Does this UMAP look good? Leave a brief comment on what you have observed (1 pts;). We will perform normalization later and you can compare their differences.**

```r
#--2--Normalization and Linear Dimensional Reduction
# Before performing normalization, let’s plot a UMAP on the first 15 PCs to see how the plot looks
atac <- ScaleData(atac)
atac <- RunPCA(atac)
atac <- RunUMAP(atac, reduction = 'pca', dims =1:15)
DimPlot(atac, reduction = "umap", pt.size = 0.5)
```

- **2.2: Let’s apply `RunTFIDF()`, `FindTopFeatures(..., min.cutoff = 'q0')`, and `RunSVD()` sequentially to do the normalization and linear dimension reduction (1 pts;). After running SVD, you will find the cell imbedding information under `atac[['lsi']]`**

```r
# 归一化及线性降维 LSI
#Let’s apply RunTFIDF(), FindTopFeatures(..., min.cutoff = 'q0'), and RunSVD() sequentially to do the normalization and linear dimension reduction
#After running SVD, you will find the cell imbedding(中间体) information under atac[['lsi']]
#ATAC-seq数据经过TF-IDF转换、特征选择和SVD降维
atac <- RunTFIDF(atac)  #Term Frequency-Inverse Document Frequency (TF-IDF)转换,用于衡量每个基因或每个ATAC峰对于每个细胞中的重要程度, #行、列归一化
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(atac)  #降维

atac[['lsi']]
```

- **2.3: Let’s use `DepthCor()` to measure the correlation of sequencing depth with each LSI component (1 pts;). Did you observe any strong correlation in this plot (1 pts;)?**

```r
DepthCor(atac)
# # 看降维之后的维度与测序深度相关性
```

![](D:\software\MarkText\picture\2023-11-21-10-59-05-image.png)



**3. Non-Linear Dimension Reduction and Clustering: Like scRNA-seq analysis, we will project the cell on a 2D plot using UMAP and do the clustering.**

- **3.1 Run UMAP dimension reduction on 2:30 LSI components (1 pts;) and plot the UMAP on a 2D graph (1 pts;). Compared to the UMAP you got from question 2.1 (without normalization), does this UMAP looks better? Leave your comment below (1 pts;).**

```r
#--3--Non-Linear Dimension Reduction and Clustering
#运行 UMAP用于降维和可视化数据时，通常会选择使用除去第一维度（PC1或LSI1等）之外的维度。
#这是因为第一维度通常捕获了数据中的最大方差，而其他维度则提供了更详细和特异的信息。
atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 2:30)
DimPlot(object = atac, reduction = "umap",pt.size = 0.5) 
```

- **3.2 Use the same dimension (reduction = 'lsi', dims = 2:30) to perform KNN `FindNeighbors()`, cluster the cells with `resolution = 0.6 and algorithm = 3` (1 pts;), and visualize the clustering results on a UMAP embedding (1 pts;). How many clusters did you get (1 pts;)? **
  
  ```r
  atac <- FindNeighbors(object = atac, reduction = 'lsi', dims = 2:30)
  atac <- FindClusters(object = atac, resolution = 0.6, algorithm = 3)
  DimPlot(object = atac, reduction = "umap",label=T,pt.size = 0.5) 
  ```

- **3.3 MAESTRO performed cell type annotation and stored the information as `assign.ident` in metadata. Please use `Idents()` to change the cell id as cell type labels and plot the cell annotation in UMAP (1 pts;). How many cell types did you get (1 pts)?**

```r
#--4--annotation
#MAESTRO performed cell type annotation and stored the information as assign.ident in metadata. Please use Idents() to change the cell id as cell type labels and plot the cell annotation in UMAP
Idents(atac) <- atac$assign.ident
DimPlot(atac, reduction = "umap",label = TRUE, pt.size = 0.5)

```

![](D:\software\MarkText\picture\2023-11-21-12-13-32-image.png)

**4. Differentially accessible peaks: To find differentially accessible regions between clusters of cells, we can perform a differential accessibility (DA) test.**

- **4.1 Let's compare CD8T cells and B cells to find differentially accessible peaks (1 pts;). hint: use `FindMarkers(..., min.pct = 0.2, test.use = 'LR')`. Print the top5 regions differentially expressed between CD8T and B cells (1 pts;).**
  
  ```r
  #--5--Differentially accessible peaks
  #To find differentially accessible regions between clusters of cells, we can perform a differential accessibility (DA) test.
  cd8t_b.diff <- FindMarkers(object = atac,
                             ident.1 = "CD8T",
                             ident.2 = "B", min.pct = 0.2, test.use = 'LR')
  
  head(cd8t_b.diff,n=5)
  ```

- **4.2 Use Violin plot to show the first DA peaks in the last question. Pleas show a violin plot with CD8T and B cell types on the x-axis, and expression level on the y-axis (1 pts;).**

```r
#show the first DA peaks
VlnPlot(object = atac, features = rownames(cd8t_b.diff)[1], pt.size = 0.1, idents = c("CD8T", "B"))

```

![](D:\software\MarkText\picture\2023-11-21-12-13-50-image.png)

**5. Based on the peak count matrix (assay 'ATAC'), MAESTRO created a gene activity matrix according to [regulatory potential model](https://github.com/liulab-dfci/MAESTRO/blob/master/example/Gene_activity_modelling/Gene_activity_modelling.md). This matrix is stored as 'ACTIVITY' assay under Seurat object. Please switch the default assay to 'ACTIVITY' using `DefaultAssay(atac) <- 'ACTIVITY'`, nNormalize and scale the gene activity matrix by using `NormalizeData()` (1 pts;). This gene activity matrix will be used for label transferring later.**

```r
#--6--gene activity
# 通过peak*cell 估计gene*cell ,然后通过gene注释cell类型
#Based on the peak count matrix (assay ‘ATAC’), MAESTRO created a gene activity matrix according to regulatory potential model. This matrix is stored as ‘ACTIVITY’ assay under Seurat object
#Please switch the default assay to ‘ACTIVITY’ using DefaultAssay(atac) <- 'ACTIVITY', Normalize and scale the gene activity matrix by using NormalizeData() (1 pts;). 
#This gene activity matrix will be used for label transferring later.
DefaultAssay(atac) <- 'ACTIVITY'
atac <- NormalizeData(atac,normalization.method = "LogNormalize",scale.factor = 10000)
```

### Part IV. Integrating scRNA-seq and scATAC-seq data

- Integrative analysis of multiome scRNA-seq and scATAC-seq
- Barcodes matching in multiome experiments data
- Label Transferring evaluation

**Please install and load following the packages:**

```r
library(circlize)
library(readr)
library(purrr)

##如果细胞在 scRNA 和 scATAC 数据中具有不同的cell barcode，我们该如何将它们配对呢？
# 可使用colnames()函数来探索scRNA-seq和scATAC-seq数据集中的列名（细胞条形码或其他标识符）
#通过检查列名，您可以确定在scRNA-seq和scATAC-seq数据集之间用于配对细胞的潜在公共标识符或元数据。这些标识符可能包括样本ID、细胞类型注释或任何两个数据集之间共有的其他信息。
#10X multiome data provides two separate lists of barcodes (barcode whitelist we used for barcode correction in MAESTRO):
#虽然列在两个不同的文件中，但两个文件中条形码的位置编码表示配对，这意味着它们具有相同的长度，您可以并排匹配从 RNA 到 ATAC 的细胞。

```

**1. When we pre-processed scRNA and scATAC separately by MAESTRO, some low-quality cells were filtered out during the analysis (Recall Part I. Question 2). In this Part, let's work on the cells common in single-cell gene expression and ATAC profiles. How are we going to pair the cells if they have different cell barcodes in scRNA and scATAC data (use `colnames()` to explore)? 10X multiome data provides two separate lists of barcodes (barcode whitelist we used for barcode correction in MAESTRO): one for gene expression and another for ATAC. Please download two barcodes lists from `/mnt/data/data/HW5/references/whitelist/atac/737K-arc-v1.txt` and `/mnt/data/data/HW5/references/whitelist/rna/737K-arc-v1.txt`to your local computer.**

- **1.1 Though listed in two separate files, The positional encoding of the barcodes in two files indicates the pairing, which means they have the same length, and you can match the cells from RNA to ATAC side by side. Now, Please match the scRNA-seq and scATAC-seq barcodes to find the list of common cells (2 pts;). How many cells are in common (1 pts;)?**
  
  ```r
  #Read multiome scatac-seq cell barcodes
  atac_whitelist<- read_tsv("scATAC/737K-arc-v1.txt", col_names = FALSE)
  #Read multiome scrna-seq cell barcodes
  rna_whitelist<- read_tsv("scRNA/737K-arc-v1.txt", col_names = FALSE)
  
  #Match scatac barcodes whitelists with scrna barcodes whitelist
  barcode_map <- set_names(atac_whitelist$X1, nm = rna_whitelist$X1) #为一个向量或列表的元素设置名称或标签
  atac_barcode <- as.character(barcode_map[Cells(rna)])
  # Attach atac cell barcodes to rna
  #将名为 rna 的对象中的细胞名称（cell names）重新命名为来自名为 atac_barcode 的对象的新名称。
  #new.names 的长度必须与 object 中的细胞数量相等，并且顺序应与要替换的细胞名称的顺序相对应。
  renamed_rna.obj <- RenameCells(object = rna,new.names = atac_barcode)
  #Find Overlapped cells between rna and atac
  barcodes.overlap <- intersect(Cells(renamed_rna.obj),Cells(atac))
  ```

- **1.2 Please use the list of common cells to subset scRNA and scATAC Seurat object (2 pts;). We will use these sub-sampled Seurat objects for the rest of the homework. hint: use `subset()`** 
  
  ```{r}
  #Taking common subset of scRNA-seq and scATAC-seq
  rna_filter <- subset(renamed_rna.obj, cells = barcodes.overlap)
  atac_filter <- subset(atac, cells = barcodes.overlap)
  ```

**2. Label Transferring: In Part III Question 3.3, we displayed a UMAP with cell type annotated by MAESTRO. You probably have noticed that the annotation result is not promising. In this section, let’s use the label transferring methods from Seurat to redo the cell annotation. You can find the detailed tutorial from [Seurat](https://satijalab.org/seurat/articles/atacseq_integration_vignette.html).**

- **2.1 First, let’s identify anchors between scATAC-seq data set and scRNA-seq data set (1 pts;). Set `DefaultAssay()` as ‘ACTIVITY’ for scATAC object. In `FindTransferAnchors()`, set ‘RNA’ as reference assay and ‘ACTIVITY’ as query assay.**

```
#--8--Label Transferring
# identify anchors between scATAC-seq data set and scRNA-seq data set
#在单细胞多组学数据集中，"anchors" 通常指的是两个不同单细胞数据集之间的共享细胞亚群或群簇。
#Set DefaultAssay() as ‘ACTIVITY’ for scATAC object. In FindTransferAnchors(), set ‘RNA’ as reference assay and ‘ACTIVITY’ as query assay.
#将来自 scATAC-seq 数据的基因活性评分与来自 scRNA-seq 的基因表达量化一起用作典型相关分析的输入
DefaultAssay(atac_filter) <- "ACTIVITY"
transfer.anchors <- FindTransferAnchors(reference = rna_filter, query = atac_filter, features = VariableFeatures(object = rna_filter), reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
```

- **2.2 Based on anchors between two modalities, please use 2:30 LSI as `weight.reduction` to transfer scRNA-seq cell types stored in `assign.ident` to scATAC-seq data (1 pts;). Attach predicted cell types to scATAC-seq metadata and Visualize the predicted cell types on UMAP (1 pts;).**

```
#将注释从 scRNA-seq 数据集中转移到 scATAC-seq 细胞上
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna_filter$assign.ident, weight.reduction = atac_filter[["lsi"]], dims = 2:30)
#添加到atac的meta.data里
atac_filter <- AddMetaData(atac_filter, metadata = celltype.predictions)

DimPlot(atac_filter, reduction = "umap",group.by = "predicted.id",label=T, pt.size = 0.5)
```

![](D:\software\MarkText\picture\2023-11-21-12-17-15-image.png)

**3. Though we can do label transferring to annotate scATAC-seq cells, back to the time without multiome data, it is hard to say if this label transferring method is accurate. The good thing in multiome data is that, by matching cell barcodes between scRNA-seq and scATAC-seq data, we actually know the “true label” for cells in scATAC-seq. We have a ground-truth annotation that can be used for evaluating the accuracy of label transferring.**

- **3.1 Since we subset the scRNA-seq and scATAC-seq objects with the list of common cells, the `assign.ident` column stored in scRNA-seq metadata is exactly the true cell type label for cells in scATAC-seq. So we can easily add `assign.ident` in scATAC-seq metadata and plot a UMAP. Please visualize the ground-truth cell type labels of scATAC-seq data on UMAP (1 pts;). Compare the results with predicted cell type UMAP from the last question. Virtually, does label transferring seem like an accurate method (1 pts;)?**

```r
#存储在 scRNA-seq 元数据中的 assign.ident 列正是 scATAC-seq 中细胞的真实细胞类型标签。
# 预测的细胞类型 UMAP 进行比较。实际上，标签转移是否是一种准确的方法
atac_filter <- AddMetaData(atac_filter, metadata = rna_filter$assign.ident,col.name = 'assign.ident.rna')
DimPlot(atac_filter, reduction = "umap",group.by = "assign.ident.rna",label=T, pt.size = 0.5)

```

![](D:\software\MarkText\picture\2023-11-21-12-18-25-image.png)

```{r}
I think label transferring is an accurate method. Because the cell type predicted from label transferring is similar to true cell type label.
```

- **3.2 What is the accuracy of the label transferring method based on the RP-enhanced model we used to calculate gene activity score? In Part IV Question 2.1, We used the 'ACTIVITY' slot to find anchors. This gene activity matrix is derived from the peak count matrix based on the RP-enhanced model we set in MAESTRO. Please show a table of what percent of the cells are correctly labeled (predicted.id == True_label) (1 pts;) ?**

```{r}
# What is the accuracy of the label transferring
table(atac_filter$assign.ident.rna,atac_filter$predicted.id)
table(atac_filter$assign.ident.rna==atac_filter$predicted.id)
```

- **3.3 For scATAC-seq data, create a heatmap with true cell type labels on x-axis and predicted cell type labels on y-axis. Fill with gradient colors to indicate the fraction of cells matched between corresponding cell types (Graduate level 2 pts;). Describe what clusters appear to map 1 to 1 between the two modalities and which clusters appear to split (Graduate level 2 pts;)?**
  
  ```r
  #For scATAC-seq data, create a heatmap with true cell type labels on x-axis and predicted cell type labels on y-axis. Fill with gradient colors to indicate the fraction of cells matched between corresponding cell types
  predictions <- table(atac_filter$assign.ident.rna,atac_filter$predicted.id)
  predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
  predictions <- as.data.frame(predictions)
  p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() +
    scale_fill_gradient(name = "Fraction of cells",
                        low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  p1
  #B, CD4Tconv, NK clusters appear to map 1 to 1 between the two modalities. pDC cluster appears to split. CD8T, Others appear to made the wrong prediction.
  
  ```
  
  ![](D:\software\MarkText\picture\2023-11-21-12-19-43-image.png)

- **3.4 Generate a density plot for prediction score generated in label transferring steps. Use two different colors for correctly annotated cells and incorrectly annotated cells (Graduate level 2 pts;).** 
  
  ```{r}
  #Generate a density plot for prediction score generated in label transferring steps. Use two different colors for correctly annotated cells and incorrectly annotated cells
  correct <- length(which(atac_filter$assign.ident.rna==atac_filter$predicted.id))
  incorrect <- length(which(atac_filter$assign.ident.rna!=atac_filter$predicted.id))
  annotation_correct <- atac_filter$assign.ident.rna==atac_filter$predicted.id
  atac_filter <- AddMetaData(atac_filter, metadata = annotation_correct,col.name = 'annotation_correct')
  data <- FetchData(atac_filter, vars = c("prediction.score.max", "annotation_correct"))
  p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
    geom_density(alpha = 0.5) + theme_bw() + scale_fill_discrete(name = "Annotation Correct",
                                                                 labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
                                                                                                                                                                               labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
  p2
  ```
  
  ![](D:\software\MarkText\picture\2023-11-21-12-20-20-image.png)
