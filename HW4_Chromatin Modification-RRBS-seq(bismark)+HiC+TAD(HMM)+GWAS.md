---
title: "STAT 115 Homework 4"
author: "(your name)"
date: "Due: Sunday 08/14/2022 by 11:59 pm"
output: html_document

---

# Part I. Chromatin Modification

DNA methylation patterns are altered in many diseases, including cancer, which makes this epigenetic mark an attractive target for various studies. Genome-wide detection of 5mC by bisulfite sequencing is regarded as the current gold standard for DNA methylation detection. In the HW3, we have performed the gene regulation analysis on prostate cancer using transcription factor ChIP-seq data. To better understand the effect of methylation on gene expression, we can utilize <mark>BS-seq</mark> to detect the DNA methylation positions in the genome.
**I.1 Reduced-representation bisulfite sequencing (<mark>RRBS-Seq</mark>) is a technique that uses one or multiple restriction enzymes on the genomic DNA to produce sequence-specific fragmentation. The fragmented genomic DNA is then treated with bisulfite and sequenced. RRBS-Seq is particularly effective on the sites with high methylation, such as in promoters and repeat regions. Given a subsampled RRBS-seq file(`/mnt/data/data/HW4/bs_subreads.fastq`) in a prostate cancer cell line (LNCaP),  perform the reads mapping with <mark>bismarker</mark>(https://github.com/FelixKrueger/Bismark/tree/master/Docs). We have generated a bismark hg38 index for you. One of the main jobs is to conduct the <mark>C-T or G-A conversion.</mark> Why are we doing this? Bismark have a function to deduplicate reads. Is it necessary for this sample? Why? After performing the mapping of the subsampled file, how many reads are unaligned and how many <mark>CpGs, CHGs, and CHHs </mark>are methylated?**
亚硫酸氢盐还原测序（RRBS-Seq）是一种在基因组DNA上使用一种或多种限制性酶来产生序列特异性片段的技术。然后用亚硫酸氢盐处理基因组 DNA 片段并进行测序。RRBS-Seq 对启动子和重复区等甲基化程度较高的位点特别有效。给定一个前列腺癌细胞系（LNCaP）的 RRBS-seq 子取样文件，用 bismarker进行读数映射。我们已为您生成了 bismark hg38 索引。其中一项主要工作是进行 C-T 或 G-A 转换。为什么要这样做？Bismark 有重复读数的功能。这个样本有必要这样做吗？为什么？对子采样文件进行映射后，有多少读数未对齐，有多少 CpGs、CHGs 和 CHHs 被甲基化？

**Hint:**  
Now that we prepared the genome for the question, follow the Bismark docs to run the main bismark functions and extract the methylations. Use directly the converted genome `/mnt/data/data/HW4/index/hg38/`. What do you see from `bismark2report`?

测序 DNA BS 转化后经过 PCR 扩增发生了 C->T 转化，那么反向互补链则是 G->A 转化。普通DNA文库经过PCR扩增会生成正负2种链，并且这两条链是反向互补的，而经过亚硫酸盐处理后，DNA链由于C->T的转换导致两条链并不互补，所以经过PCR扩增后会形成正负4种链，即所谓的OT、CTOB、OB、CTOB

![](D:\software\MarkText\picture\2023-11-20-18-57-45-image.png)

<mark>Bismark</mark> 是一款比较有名的甲基化测序比对软件，Bismark 比对前会将所有 reads 进行 C->T 和 G->A 转换，并且将参考基因组同样进行这 2 种转换。非链特异性建库reads的方向无法确定，则转换后的2种reads要分别与参考基因组的2种情况分别进行比对，相当于每条 reads 进行 4 次不同比对，从里面选择最佳比对作为结果。如果是链特异性建库由于已知reads的属于哪一条链，则会只比对两种情况，这样比对速度上要很多。由于Bismark会将reads和参考基因组做C->T或G->A转换，就某一种情况来说转换后只剩下3种碱基，故Bismark的工作原理又被俗称为“三碱基比对”。

![](D:\software\MarkText\picture\2023-11-20-18-59-10-image.png)

链接：https://www.jianshu.com/p/0ac601dcb48f

deduplicate_bismark --bam [options] <filenames>
对 Bismark 的比对结果**去重**，去除以相同方向比对到相同位置的 reads。**建议应用于 WGBS 数据，不适用与 RRBS ，amplicon，以及 target 富集文库的数据**。

```she
conda install -c bioconda bismark -y
# 准备比对基因组。需要指定一个目录，其中包含你想要比对reads的基因组(请注意，bismark_genome_prepare脚本期望该文件夹中包含FastA文件，扩展名为.fa或. FastA，每个文件有单个或多个序列条目)
bismark_genome_preparation --bowtie2 hg38/

nohup bismark --genome hg38/ bs_subreads.fastq --parallel 5  -o bismark_output/  &
# 该步骤会生成两个文件：
# bs_subreads_bismark_bt2.bam (包含比对结果和甲基化水平)
# bs_subreads_bismark_bt2_SE_report.txt (包含比对和甲基化水平的统计结果)


# deduplicate_bismark --bam [options] <filenames>
## 对 Bismark 的比对结果去重，去除以相同方向比对到相同位置的 reads。建议应用于 WGBS 数据，不适用与 RRBS ，amplicon，以及 target 富集文库的数据。

# 未被甲基化的C在BS的处理下会被为U，经历PCR之后，U会变为T，互补链上的G会变为A。因此在比对时，要将测序的reads进行C-T转换，参考基因组进行C-T、G-A转换，以此做比对时，reads才能比对上基因组。
# 不用deduplicate去重，因为RRBS-seq用的是切割同一位点的限制性内切酶，酶切位点是固定的，从同一位点切割，自然会产生较多的重复片段，所以reads总是从那些剪切位点开始测的，若去重会丢掉大量有用的reads，因此不需去重。
# 4173条reads没有比对上，1000条reads没有uniquely比对，4827条reads uniquely比对，比对率为48.3%
# 1343个CpGs, 21个CHGs, 99个CHHs被methylated
```

![](D:\software\MarkText\picture\2023-11-20-19-12-01-image.png)

**I.2 Methylation in cytosine at promoter regions normally suppresses the gene expression, while H3K4me3 and H3K27ac histone modifications at promoter regions imply higher gene expression. We have processed RRBS-seq data on the prostater cancer cell line dataset(https://www.encodeproject.org/experiments/ENCSR859PDD/) and report the high methylation signal sites in `/mnt/data/data/HW4/data/Methylation.bed`. Draw violin plots of the expression level of genes with methylation, H3K4me3, and H3K27ac in their promoters. Could you find that the higher methylation signals repress the gene expression?**  

启动子区域胞嘧啶的甲基化通常会抑制基因的表达，而启动子区域的 H3K4me3 和 H3K27ac 组蛋白修饰则意味着基因的高表达。我们处理了前列腺癌细胞系数据集的 RRBS-seq 数据，并在 `/mnt/data/data/HW4/data/Methylation.bed`中报告了高甲基化信号位点。绘制启动子中存在甲基化、H3K4me3 和 H3K27ac 的基因表达水平的小提琴图。您能否发现甲基化程度较高的信号会抑制基因的表达？

**Hint:**  

- The `/mnt/data/data/HW4/ProstateCancer_H3K4me3_peaks.bed` and `/mnt/data/data/HW4/ProstateCancer_H3K27ac_peaks.bed` are the H3K4me3 and H3K27ac peaks files of prostate cancer. Try to find the<mark> intersection of loops </mark>and histone modification signal intervals. In the `/mnt/data/data/HW4/Expr_loc.txt`, the first three columns are chromosome coordinates, the fourth column is the gene expression score of a prostate cancer tissue, and the rest columns are the gene symbols.  

- Use `>` to write your results into an output file. Use `uniq` command to remove the duplicates from your intersections. e.g. `bedtools intersect {YOUR FLAGS} | uniq > OUTPUT.TXT`  

- The methylation group might have a lower signal-noise ratio. So drop the genes with very low expression (< 0.0001) and use the log scale when you generate the violin plots for better visualization. You can use `ylim` and `scale_y_log10` from `ggplot`.
  
  ```shell
  # To find the intersection of loops and histone modification signal intervals
  # filter genes with expression scores greater than or equal to 0.0001, then pipes the result to bedtools intersect to find the intersection with the previous intersection.bed file. Finally, uniq is used to remove duplicate entries, and the result is saved in filtered_intersect.bed.
  # 有无methylation, H3K4me3, and H3K27ac的基因表达量
  awk '{if ($4 >= 0.0001) print $0}' Expr_loc.txt | bedtools intersect -a stdin -b ProstateCancer_H3K4me3_peaks.bed -wa| uniq > H3K4me3_intersect.txt
  awk '{if ($4 >= 0.0001) print $0}' Expr_loc.txt | bedtools intersect -a stdin -b ProstateCancer_H3K4me3_peaks.bed -v | uniq > H3K4me3_not_intersect.txt
  awk '{if ($4 >= 0.0001) print $0}' Expr_loc.txt | bedtools intersect -a stdin -b ProstateCancer_H3K27ac_peaks.bed -wa| uniq > H3K27ac_intersect.txt
  awk '{if ($4 >= 0.0001) print $0}' Expr_loc.txt | bedtools intersect -a stdin -b ProstateCancer_H3K27ac_peaks.bed -v | uniq > H3K27ac_not_intersect.txt
  awk '{if ($4 >= 0.0001) print $0}' Expr_loc.txt | bedtools intersect -a stdin -b Methylation.bed -wa| uniq > Methylation_intersect.txt
  awk '{if ($4 >= 0.0001) print $0}' Expr_loc.txt | bedtools intersect -a stdin -b Methylation.bed -v | uniq > Methylation_not_intersect.txt
  
  #-wa 选项表示输出a文件交集区域的完整行。-wb 则输出b文件交集区域的完整行。
  #-v 输出在-a参数文件中没有overlap的区域
  ```
  
  ```r
  # 
  setwd("/home4/liuyw/test/Liuxiaole_Training_homework_data/HW4")
  
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(ggsignif)
  
  # read files
  H3K4me3_expr <- read.table("H3K4me3_intersect.txt", sep = "\t", header = 0,col.names=c("chrom", "start", "end","score","Ensembl_ID","Symbol"))
  H3K27ac_expr <- read.table("H3K27ac_intersect.txt", sep = "\t", header = 0,col.names=c("chrom", "start", "end","score","Ensembl_ID","Symbol"))
  methy_expr <- read.table("Methylation_intersect.txt", sep = "\t", header = 0,col.names=c("chrom", "start", "end","score","Ensembl_ID","Symbol"))
  not_H3K4me3_expr <- read.table("H3K4me3_not_intersect.txt", sep = "\t", header = 0,col.names=c("chrom", "start", "end","score","Ensembl_ID","Symbol"))
  not_H3K27ac_expr <- read.table("H3K27ac_not_intersect.txt", sep = "\t", header = 0,col.names=c("chrom", "start", "end","score","Ensembl_ID","Symbol"))
  not_methy_expr <- read.table("Methylation_not_intersect.txt", sep = "\t", header = 0,col.names=c("chrom", "start", "end","score","Ensembl_ID","Symbol"))
  
  H3K4me3_expr["group"] = "H3K4me3"
  H3K27ac_expr["group"] = "H3K27ac"
  methy_expr["group"] = "DNA methylation"
  not_H3K4me3_expr["group"] = "H3K4me3"
  not_H3K27ac_expr["group"] = "H3K27ac"
  not_methy_expr["group"] = "DNA methylation"
  
  H3K4me3_expr["modification"] = "modified"
  H3K27ac_expr["modification"] = "modified"
  methy_expr["modification"] = "modified"
  not_H3K4me3_expr["modification"] = "unmodified"
  not_H3K27ac_expr["modification"] = "unmodified"
  not_methy_expr["modification"] = "unmodified"
  
  data <- rbind(H3K4me3_expr,H3K27ac_expr,methy_expr,not_H3K4me3_expr,not_H3K27ac_expr,not_methy_expr)
  ```
  
  ggplot(data, aes(x=group, y=score,fill = modification)) +
    geom_violin(scale="width") +
    scale_y_log10() +
    labs(x="Group", y="Expression Score (log10 scale)", title="Gene Expression Levels") +
    theme_bw()
  
  p <- ggplot(data, aes(group,score,fill=modification))+
    geom_violin(width = 0.8,alpha=0.2,position = position_dodge(0.8))+
    geom_boxplot(width=0.3,position = position_dodge(0.8))+
    scale_fill_lancet()+
    scale_y_log10() +
    labs(x="Group", y="Expression Score (log10 scale)", title="Gene Expression Levels") +
    theme_bw()+
    stat_compare_means(aes(group=modification),                       #按分组进行统计检验
  
                       method = "t.test",
                       paired = F,                             #非配对t检验
                       symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                          symbols=c("***","**","*","ns")),
                       label = "p.signif",
                       size=4.5)                              #显著性符号的大小 
  
  ggsave("violin.png",width = 6,height = 5)  ##ggplot 中直接保存

```
![](D:\software\MarkText\picture\2023-11-20-19-23-05-image.png) 

# Part II. HiC

Genome architecture plays a key role in nuclear functions. The spatial arrangement and proximity of genes have been linked to biological functions, such as gene replication, regulation, and transcription. The Hi-C technique allows researchers to extract the interaction frequency for all loci of a genome at high-throughput and at a genome-wide scale. In this part, we will learn how the HiC data integrates with other epigenetic data and genome architecture affects gene expression in prostate cancer.

基因组结构在核功能中起着关键作用。基因的空间排列和邻近程度与基因复制、调控和转录等生物功能有关。通过 Hi-C 技术，研究人员可以在全基因组范围内高通量提取基因组所有位点的相互作用频率。在本部分中，我们将了解 HiC 数据如何与其他表观遗传数据整合，以及基因组结构如何影响前列腺癌的基因表达。

**II.1 Given subsampled example fastq files generated by Hi-C technique(`/mnt/data/data/HW4/p2_data/HiC_fq`), perform reads <mark>mapping</mark>, <mark>filtering</mark>, and <mark>binning</mark> to this data with <mark>runHiC</mark>(http://xiaotaowang.github.io/HiC_pipeline/quickstart.html). How about the <mark>quality</mark> of this subsampled data? Could you explain the QC criteria of Hi-C data?**
**Hint:**  

- The runHiC script takes the bwa index as the input. You can use the one we used for HW3 `/mnt/data/data/HW3/bwa_index_hg38/bwa_index_hg38`.

- Note that the index file is **NOT*** the same as the one we used for bismark. What is the main difference between bowtie and bwa in practice? (Hint: https://www.biostars.org/p/134418/).

**Hint:**  

- Prepare your working directory with the index files and your metadata file `datasets.tsv`. See http://xiaotaowang.github.io/HiC_pipeline/quickstart.html#create-the-meta-data-file. In our case we would use these following setup configuration. Make sure you have more than one replicate sample in your metadata to avoid future ambiguity. You may need this file `/mnt/data/data/HW4/hg38.chrom.sizes` as well. For saving time and storage on the server, you can use soft link to link file to your work directory`ln -s /mnt/data/data/HW3/bwa_index_hg38/ path_to_your_work_directory/`.
```

  HiC_subreads GM06990 R1 HindIII
  HiC_subreads GM06990 R2 HindIII

```
Run the runHiC commands and generate the stats table. Put here the output `*.stats` and plots `*.png` for the summary group `allReps` **ONLY**. What is the message you can read from the results.

```shell
##首先 -p DATAFOLDER, --dataFolder DATAFOLDER: 根数据文件夹的路径，包含测序 reads 和参考基因组。
mkdir HiC
# https://xiaotaowang.github.io/HiC_pipeline/quickstart.html 按照这个教程重新建了bwa index
mkdir hg38
cd hg38
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/*
ln -s /home4/liuyw/test/Liuxiaole_Training_homework_data/HW4/HiC_fq/* /home4/liuyw/test/Liuxiaole_Training_homework_data/HW4/HiC/fastq
# -m METADATA, --metadata METADATA: 描述每个 SRA/FASTQ 文件的元数据文件。包含四列：SRA/FASTQ 文件名前缀，细胞系/样本名，生物复制标签和限制酶名称。
# 在运行 runHiC 之前，您需要做的另一件事是创建一个 TXT 文件 在工作区子文件夹下命名为“datasets.tsv”
# 在这里，“datasets.tsv” 是描述您的 Hi-C 数据的元数据文件，它应该包含 4 列。依次是：SRA 文件名的前缀（对于 FASTQ 读取格式，它应该是文件名除去“_1.fastq”或“_2.fastq”子字符串的前导部分），细胞系名称，生物复制标签和限制酶名称。
# 请注意，对于 Arima Hi-C，您可以将酶名设置为 Arima；对于不使用限制酶进行 DNA 切割的实验，您可以为记录任意设置酶名。例如，对于 Micro-C，您可以将其设置为 MNase；对于 ChIA-PET，您可以将其设置为sonication。
vim datasets.tsv
# HiC_subreads GM06990 R1 HindIII
# HiC_subreads GM06990 R2 HindIII
```

```shell
#--1--# mapping 
nohup runHiC mapping -m datasets.tsv -p input \
    -g hg38 -f fastq -F FASTQ -A bwa-mem  \
    -C ./hg38.chrom.sizes \
    --include-readid --drop-seq --include-sam --add-frag  \
    -t 10 --tmpdir tmp_dir \
    --logFile runHiC-mapping.log &

# 将在当前下创建名为 alignments-hg38 和 pairs-hg38 的两个子文件夹 工作目录
# reads对将用bwa-mem映射到hg38参考基因组，比对结果将以BAM格式报告，并置于pair-hg38下
# BAM文件将使用pairs-hg38下的pairtools解析为.pairs。
#--2--# filtering 
nohup runHiC filtering -m datasets.tsv --pairFolder pairs-hg38 --logFile runHiC-filtering.log --tmpdir tmp_dir --nproc 10 &
# 创建一个名为 filtered-hg38 的新子文件夹。有效的 .pairs.gz 文件
#--3--# binning
## an .mcool file will be produced under the coolers-hg38 sub-folder for each .pairs.gz file using cooler. 
nohup runHiC binning -f filtered-hg38 --logFile runHiC-binning.log --nproc 10 &
# 将在 coolers-hg38 子文件夹下为每个生成一个 .mcool 文件

##一步进行所有--pileup
runHiC pileup -p ./input -g hg38 -f fastq -F FASTQ -A bwa-mem -t 10 --include-readid --drop-seq --chunkSize 1500000 --logFile runHiC.log

#--4--# quliaty
#Then statistic tables about your data can be found in .stats files under the *filtered-hg38 sub-folder.
nohup runHiC quality -m datasets.tsv -L filtered-hg38 --logFile runHiC-quality.log &
```

  **II.2 In the II.1, you have learned how to generate a genomic interaction matrix on the example file. Here we provided a genomic interaction data on chr21 in prostate cancer(`/mnt/data/data/HW4/chr21.chr21_10000.cool`). Normalize this data at 10k resolution and perform loop calling with 10% CTCF model. How many loops with >0.9 confidence can you find? Then draw a genomic contact heatmap to show the interaction in chr21. Are there any highly interactive regions?**

  在第 II.1 节中，你已经学会了如何在示例文件中生成基因组相互作用矩阵。这里我们提供了前列腺癌中 chr21 的基因组相互作用数据。以 10k 分辨率对该数据进行归一化处理，然后使用 10% CTCF 模型执行循环调用。您能找到多少置信度大于 0.9 的loops？然后绘制基因组接触热图，显示 chr21 中的相互作用。是否存在高度交互区域？

  **Hint:** cooler(https://cooler.readthedocs.io/en/latest/) can perform normalization and peakachu(https://github.com/tariks/peakachu) can perform loop calling on Hi-C data. higlass(https://higlass.io/)(You may need to install this one on your local machine.) may help with visualization.

```shell
# # cooler标准化
cooler zoomify --resolutions 10000 --balance --nproc 10 -o chr21_10000.zoomify.cool chr21.chr21_10000.cool
# peakachu call loop
peakachu score_genome -r 10000 --balance \
    -m ./down10.ctcf.pkl \
    -p chr21_10000.zoomify.cool::/resolutions/10000 \
    -O chr21_10kb-scores.bedpe 

## ValueError: Number of features of the model must match the input. Model n_features is 243 and input n_features is 225 

# peakachu寻找显著loop
peakachu pool -r 10000 -i chr21_10kb-scores.bedpe -o chr21_10kb-scores.0.9.bedpe -t 0.9 &

# 置信区间90%内，有106个loop  # higlass包安装不成功，此处代码跑不出来，也可以用IGV看bedpe文件
knitr::include_graphics("p2-2-1.png")
```

```python
# 师姐代码
from higlass.client import View, Track
from higlass.tilesets import cooler
import higlass

ts1 = cooler('./bash/chr21_10000.zoomify.cool')
tr1 = Track('heatmap', tileset=ts1)
view1 = View([tr1])
display, server, viewconf = higlass.display([view1])

display
# higlass包安装不成功，此处代码跑不出来，也可以用IGV看bedpe文件
```

  **II.3 Transcription factors help construct the chromatin loops. With the provided prostate cancer ATAC-seq peaks file(`/mnt/data/data/HW4/tumor_ATAC_peaks.bed`), could you find the open regions in the loop anchors? What factors bind in the loop anchors? What potential roles may these factors play?**

```shell
# 师姐代码
#loop ATAC
cut -f1,2,3 ../p2-2/chr21_10kb-scores.0.9.bedpe > chr21_10kb_loop_anchor.bed
cut -f4,5,6 ../p2-2/chr21_10kb-scores.0.9.bedpe >> chr21_10kb_loop_anchor.bed
bedtools intersect -a /mnt/data/data/HW4/tumor_ATAC_peaks.bed -b chr21_10kb_loop_anchor.bed -wa |uniq > chr21_10kb_loop_anchor_ATAC.bed

#call motif
awk '{print $4"\t"$1"\t"$2"\t"$3"\t""."}' chr21_10kb_loop_anchor_ATAC.bed > chr21_10kb_loop_anchor_ATAC.homer.bed
nohup findMotifsGenome.pl chr21_10kb_loop_anchor_ATAC.homer.bed hg38 open_in_loop_motif -preparsedDir homer &
```

  **II.4 Based on the loop file, could you find the genes within loops on chr21? Do these target genes express higher than the genes without loops structure?**

```shell
# 师姐daima
awk ' BEGIN{OFS="\t"}{print $1,$2,$6}' ../p2-2/chr21_10kb-scores.0.9.bedpe > chr21_10kb_loop_anchor.bed
bedtools intersect -a /mnt/data/data/HW4/Expr_loc_chr21.txt -b ../p2-3/chr21_10kb_loop_anchor.bed -wa |uniq > gene_in_chr21_loop_anchor.txt
bedtools intersect -a /mnt/data/data/HW4/Expr_loc_chr21.txt -b ../p2-3/chr21_10kb_loop_anchor.bed -v |uniq > gene_NOT_in_chr21_loop_anchor.txt
```

```r
# 师姐代码
library(ggplot2)

setwd("/Users/Yihongyu/Desktop/codes/BioTrainning/week4/Homework4")

# read files
gene_chr21 <- read.table("./bash/p2-4/gene_in_chr21_loop.txt", sep = "\t", header = 0)
gene_not_chr21 <- read.table("./bash/p2-4/gene_NOT_in_chr21_loop.txt", sep = "\t", header = 0)
gene_chr21["group"] = "IN"
gene_not_chr21["group"] = "OUT"
dat <- rbind(gene_chr21, gene_not_chr21)
colnames(dat) = c("chrom", "start", "end", "expr_score", "ENSG", "symbol", "group")

ggplot(dat, aes(x = group, y = log(expr_score), fill = group)) +
  geom_boxplot() +
  ylab("Gene expression")
```

# Part III. Hidden Markov Model and TAD boundaries

  Topologically associating domains (TADs) define genomic intervals, where sequences within a TAD physically interact more frequently with each other than with sequences outside the TAD. TADs are often defined by HiC (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3149993/), an experimental technique designed to study the three-dimensional architecture of genomes. HiC generates PE sequenced data, where the two mate pairs indicate two genomic regions that are might be far apart in the genome, but physically interact with each other. If we look across the genome in bins (40kb in the early paper, but now can go down to 5-10kb with deeper sequencing), we could find reads that are mapped there and check whether their interacting mate pairs are mapped upstream or downstream. In each bin, we can calculate a directional index (DI) to quantify the degree of upstream or downstream bias of a given bin (for more details, see the supplement- `Supplement_10.1038_nature11082.pdf` ). For this HW, we ask you to implement a hidden Markov Model (Viterbi) to find regions with upstream bias (DI < 0) and those with downstream bias (DI > 0), even though the DI in individual bins might have some noise. This way, TAD boundaries could be discovered as clusters of bins from negative DIs to positive DIs (see Supplementary Figure 12b).
  For simplicity, we will only have two hidden states (upstream, downstream), and use the following HMM parameters (these do not necessarily capture the real data distribution, but just to help your implementation):

  拓扑关联域（TAD）定义了基因组区间，TAD 内的序列之间的物理相互作用比 TAD 外的序列更频繁。TAD通常由HiC定义，HiC是一种实验技术，旨在研究基因组的三维结构。HiC 生成 PE 测序数据，其中的两个配对表示两个基因组区域，这两个区域在基因组中可能相距甚远，但在物理上相互影响。如果我们在整个基因组中寻找一个个bin（早期论文中为 40kb，但现在随着测序的深入可以缩小到 5-10kb），我们就可以找到映射在那里的读数，并检查它们相互作用的配对是映射在上游还是下游。在每个bin中，我们可以计算一个方向指数（DI），以量化特定bin的上游或下游偏倚程度（更多详情)。对于本 HW，我们要求您实施一个隐藏马尔可夫模型（Viterbi），以找到具有上游偏差（DI < 0）和下游偏差（DI > 0）的区域，即使单个bin的 DI 可能有一些噪声。这样，就可以发现 TAD 边界是由负 DIs 到正 DIs 的分区群（见补充图 12b）。
  为简单起见，我们将只有两个隐藏状态（上游和下游），并使用以下 HMM 参数（这些参数不一定能捕捉到真实的数据分布，但可以帮助您实现）：

```
Initial probability: upstream = 0.5, downstream = 0.5
Transition probability: Pb(up to up) = Pb(dn to dn) = 0.9, Pb(up to dn) = Pb(dn to up) = 0.1
Emission probabilities:
P{<-1200, [-1200,-800), [-800,-500), [-500,0), [0,500), [-500,800), [800, 1200), >= 1200 | upstream} = (0.01, 0.01, 0.02, 0.04, 0.65, 0.15, 0.08, 0.04)
P{<-1200, [-1200,-800), [-800,-500), [-500,0), [0,500), [-500,800), [800, 1200), >= 1200 | downstream} = (0.04, 0.08, 0.15, 0.65, 0.04, 0.02, 0.01, 0.01)
```

  <mark>隐马尔可夫模型</mark>（Hidden Markov Model，HMM）是一种统计模型，广泛应用于时序数据建模和分析。它包含两个主要部分：状态序列和观测序列。

1. **状态序列（Hidden States）：**
   
   - 隐马尔可夫模型包含一个隐藏的状态序列，表示系统在不同时间点的内部状态。这些状态对外部观测是不可见的，因此称为"隐"状态。
   - 状态序列通常是一个马尔可夫链，即未来状态只依赖于当前状态，与过去状态无关。这就是“马尔可夫”模型的含义。

2. **观测序列（Observations）：**
   
   - 对应于每个时间点的隐状态，存在一个可观测的符号或数值，构成观测序列。这些观测是我们可以测量或观察到的数据。
   - 观测序列的生成是由隐藏状态序列驱动的，每个隐藏状态对应一个观测值的概率分布。

3. **转移概率（Transition Probabilities）：**
   
   - 描述状态序列中状态转移的概率。在离散时间步内，系统从一个状态转移到另一个状态的概率。
   - 转移概率通常用状态转移矩阵表示，其中矩阵元素(a_ij)表示从状态i到状态j的转移概率。

4. **发射概率（Emission Probabilities）：**
   
   - 描述在给定状态的情况下观测到特定观测值的概率。
   - 发射概率通常用发射矩阵表示，其中矩阵元素(b_ij)表示在状态i的情况下观测到观测值j的概率。

5. **初始概率分布（Initial Probabilities）：**
   
   - 描述模型在时间序列开始时处于各个状态的概率分布。
     HMM主要应用于两个问题：
- **学习问题（Learning Problem）：** 根据观测序列，估计模型的参数，包括转移概率、发射概率和初始概率。

- **预测问题（Inference Problem）：** 根据模型和观测序列，预测状态序列或未来的观测。
  
  **Given the DI file (`/mnt/data/data/HW4/ESC.Dixon_2015.DI.chr21.txt`), implement and utilize the Viterbi algorithm to predict the hidden states of the Hi-C data. Visualize your result with a graph utilizing the following: midpoint of genomic bin on the <mark>x-axis</mark>; DI score per bin on the <mark>y-axis</mark>; color: hidden state of the HMM. Please do not use the "HMM" package when you solving this question but recommended when you try to validate your results.**  
  
  给定 DI 文件，实施并利用Viterbi算法预测 Hi-C 数据的隐藏状态。用以下图表展示您的结果：X 轴为基因组bin的中点；Y 轴为每个bin的 DI 分数；颜色：HMM 的隐藏状态。
  
  **Hint1**: Examples HMM code can be found at:
  http://www.adeveloperdiary.com/data-science/machine-learning/implement-viterbi-algorithm-in-hidden-markov-model-using-python-and-r/
  **Hint2**: The observations are continuous or have too many discrete values. Try binning them into a few discrete regions. Use `cut` function built in R.
  
  ```r
  #官网的python的代码有具体解释
  forward = function(v, a, b, initial_distribution){
  T = length(v)
  M = nrow(a)
  alpha = matrix(0, T, M)
  
  alpha[1, ] = initial_distribution*b[, v[1]]
  
  for(t in 2:T){
    tmp = alpha[t-1, ] %*% a
    alpha[t, ] = tmp * b[, v[t]]
  }
  return(alpha)
  }
  
  backward = function(v, a, b){
  T = length(v)
  M = nrow(a)
  beta = matrix(1, T, M)
  
  for(t in (T-1):1){
    tmp = as.matrix(beta[t+1, ] * b[, v[t+1]])
    beta[t, ] = t(a %*% tmp)
  }
  return(beta)
  }
  
  BaumWelch = function(v, a, b, initial_distribution, n.iter = 100){
  for(i in 1:n.iter){
    T = length(v)
    M = nrow(a)
    K=ncol(b)
    alpha = forward(v, a, b, initial_distribution)
    beta = backward(v, a, b)
    xi = array(0, dim=c(M, M, T-1))
  
    for(t in 1:T-1){
      denominator = ((alpha[t,] %*% a) * b[,v[t+1]]) %*% matrix(beta[t+1,]) 
      for(s in 1:M){
        numerator = alpha[t,s] * a[s,] * b[,v[t+1]] * beta[t+1,]
        xi[s,,t]=numerator/as.vector(denominator)
      }
    }
  
    xi.all.t = rowSums(xi, dims = 2)
    a = xi.all.t/rowSums(xi.all.t)
  
    gamma = apply(xi, c(1, 3), sum)  
    gamma = cbind(gamma, colSums(xi[, , T-1]))
    for(l in 1:K){
      b[, l] = rowSums(gamma[, which(v==l)])
    }
    b = b/rowSums(b)
  
  }
  return(list(a = a, b = b, initial_distribution = initial_distribution))
  }
  
  Viterbi=function(v,a,b,initial_distribution) {  
  T = length(v)
  M = nrow(a)
  prev = matrix(0, T-1, M)
  omega = matrix(0, M, T)
  
  omega[, 1] = log(initial_distribution * b[, v[1]])
  for(t in 2:T){
    for(s in 1:M) {
      probs = omega[, t - 1] + log(a[, s]) + log(b[s, v[t]])
      prev[t - 1, s] = which.max(probs)
      omega[s, t] = max(probs)
    }
  }
  
  S = rep(0, T)
  last_state=which.max(omega[,ncol(omega)])
  S[1]=last_state
  
  j=2
  for(i in (T-1):1){
    S[j]=prev[i,last_state] 
    last_state=prev[i,last_state] 
    j=j+1
  }
  
  S[which(S==1)]='A'
  S[which(S==2)]='B'
  
  S=rev(S)
  
  return(S)  
  }
  
  setwd("/home4/liuyw/test/Liuxiaole_Training_homework_data/HW4/part3_HMM",)
  chr21_DI <- read.table("ESC.Dixon_2015.DI.chr21.txt", sep = "\t", header = 0, col.names = c("chrom","start","end","DI"))
  chr21_DI$Observations <- cut(chr21_DI$DI, breaks=c(-Inf,-1200,-800,-500,0,500,800,1200,Inf), right = FALSE)
  chr21_DI$midpoint <- (chr21_DI$start + chr21_DI$end)/2
  
  ## Initial Probabilities
  initial_distribution <- c(0.5, 0.5)
  ## Transition Probabilities
  A <- matrix(c(0.9,0.1,0.1,0.9), nrow = 2, ncol = 2,
            dimnames = list(c("upstream", "downstream"), c("upstream", "downstream")))
  ## Emission probabilities
  B <- t(matrix(c(0.01, 0.01, 0.02, 0.04, 0.65, 0.15, 0.08, 0.04,0.04, 0.08, 0.15, 0.65, 0.04, 0.02, 0.01, 0.01), nrow = 8, ncol = 2))
  rownames(B) <- c("upstream", "downstream")
  colnames(B) <- c("[-Inf,-1200)","[-1200,-800)","[-800,-500)","[-500,0)","[0,500)","[500,800)", "[800,1200)","[1200, Inf)")
  
  # hidden state
  myout.hidden=Viterbi(chr21_DI$Observations,A,B,initial_distribution)
  chr21_DI$HMM = myout.hidden
  
  # 画图
  library(ggplot2)
  
  ggplot(chr21_DI, aes(x = midpoint, y = DI, fill = HMM)) +
  geom_bar(stat = 'identity') +
  ylab("Directional Index (DI)")
  ```

![](D:\software\MarkText\picture\2023-11-21-09-03-45-image.png)

# Part IV: GWAS Followup

The NHGRI-EBI <mark>GWAS</mark> Catalog is a curated dataset of trait-associated genetic variants for human. While it provides association between single-nucleotide polymorphisms (SNPs) and trait (i.e. cancer), the genetic variants in GWAS catalog are not necessarily causative or functional for a trait, since SNPs can be highly correlated measured by <mark>linkage disequilibrium</mark> (LD). To learn the potential functional effect of a certain SNP, especially the non-coding variants, we can use <mark>RegulomeDB</mark> to explore the potential function of the SNP.

NHGRI-EBI GWAS 目录是一个人类性状相关基因变异的编辑数据集。虽然它提供了单核苷酸多态性（SNPs）与性状（如癌症）之间的关联，但 GWAS 目录中的遗传变异并不一定是性状的因果关系或功能性变异，因为 SNPs 可以通过连锁不平衡（LD）测量高度相关。要了解某个 SNP（尤其是非编码变异）的潜在功能效应，我们可以使用 RegulomeDB 来探索 SNP 的潜在功能。

You will explore the following online resources: 

The <mark>NHGRI-EBI GWAS</mark> catalog (https://www.ebi.ac.uk/gwas/), 

<mark>dbSNP</mark> (https://www.ncbi.nlm.nih.gov/snp/ ), 

<mark>LDLink</mark> (https://ldlink.nci.nih.gov/), and 

<mark>RegulomeDB</mark> (the beta version http://regulomedb.org or the more stable older version http://legacy.regulomedb.org/).
**IV.1 Explore whether there are genetic variants within the gene BRCA2 which are associated with any traits. What traits are associated with the BRCA2 variants? Which SNP has the smallest p-value related to breast cancer? What is the risk allele?**  

探究 BRCA2 基因中是否存在与任何性状相关的遗传变异。哪些性状与 BRCA2 变异相关？哪个 SNP 与乳腺癌相关的 p 值最小？风险等位基因是什么？

```
基因BRCA2有62个genetic vairants，25个traits
和BRCA2 variants相关的traits包括low density lipoprotein cholesterol measurement、low density lipoprotein cholesterol measurement、breast carcinoma、total cholesterol measurement和lung carcinoma

和breast cancer相关SNP中，rs11571833的p值最小，为6 x 10(-6)；risk allele为T
```

**IV.2 For the BRCA2 SNP with the most significant association with breast cancer, what consequence does the risk allele have on the BRCA2 protein sequence? Based on 1000 Genomes in LDLink, what is the allele frequency of the risk allele among the 5 ethnicities In the population with the highest risk in the resource, what is the expected number of people with heterozygous genotype at this SNP, assuming linkage disequilibrium?**  

对于与乳腺癌关系最密切的 BRCA2 SNP，风险等位基因对 BRCA2 蛋白序列有什么影响？根据 LDLink 中的《1000 个基因组》，风险等位基因在 5 个种族中的等位基因频率是多少？ 在资源中风险最高的人群中，假设存在连锁不平衡，该 SNP 的杂合基因型的预期人数是多少？

```
risk allele rs11571833-T与breast cancer最显著相关，对BRCA2蛋白的consequences是top gained(导致终止密码子提前，转录本更短)
在5个人种中(AFR非洲人、EUR欧洲人、AMR(Ad mixed American)、EAS东亚人、ASA南亚人)，SNP rs11571833位点的A频率为0.996，T的频率为0.004
EUR人种的risk最高，T的频率为1.09%，A为98.91%。该人种样本为503人，由此计算杂合子为503*0.9891*0.0109*2 = 10.84，约11人
```

**IV.3 Explore a certain SNP, rs4784227, that was reported to be associated with breast cancer. Is it an intergenic, exonic or intronic variant? What gene does it fall in?**  

探究据报道与乳腺癌有关的某个 SNP（rs4784227）。它是基因间变异、外显子变异还是内含子变异？它属于哪个基因？

```
rs4784227是intronic variant，位于CASC16基因内
```

**IV.4 Explore the SNP rs4784227 in RegulomeDB. What functional category does the rank score (or Regulome DB Score) implicate? What factors does RegulomeDB take into consideration while scoring the potential function of SNPs?**

在 RegulomeDB 中探索 SNP rs4784227。等级得分（或 Regulome DB Score）涉及哪种功能类别？RegulomeDB 在对 SNP 的潜在功能进行评分时考虑了哪些因素？

Regulome DB Score是2b，表示TF binding + any motif + DNase Footprint + DNase peak

SNP rs4784227可能会影响binding，TF ChIP-seq, motif, DNase-seq和组蛋白修饰等都用来评估SNP rs4784227的潜在功能

while scoring potential function of SNPs，RegulomeDB考虑了eQTL, TF binding, motif, DNase footprint, DNase peak等

![](D:\software\MarkText\picture\2023-11-20-20-24-38-image.png)

**IV.5 Describe the evidence that implicates the regulatory potential of rs4784227, for example, list several transcription factors with binding peaks overlapping this SNP; report the cell types with open chromatin regions overlapping this SNP.**

描述与 rs4784227 潜在调控作用有关的证据，例如，列出与该 SNP 结合峰重叠的几种转录因子；报告与该 SNP 重叠的染色质开放区域的细胞类型。

TF ChIP-seq栏可知，TCF12、TCFL12、TEAD4、FOXA1、ARID3A、MYBL2、NFIC等转录因子peak和SNP rs4784227有overlap

TF motif栏可知，rs4784227富集了lsgf3g motif，可能会有lsgf3g转录因子到rs4784227上

DNase-seq栏可知，Mcf7、Hepg2、A549、H7es、Nhbera和Panislets等细胞的开放染色质区域与rs4784227有overlap

组蛋白修饰栏可知，quiescent/low activity in Blood & T cell, Brain, Heart and other cells

**IV.6 Read the paper by Cowper-Sal et al. (PMID 23001124) and summarize the potential mechanisms of the above SNP's function in terms of affecting transcription factor-DNA interaction and regulating genes.**

阅读 Cowper-Sal 等人的论文（PMID 23001124），总结上述 SNP 在影响转录因子-DNA 相互作用和调控基因方面的潜在功能机制。

1) The risk-associated SNPs are enriched in the cistromes of some genes and the epigenome of histone modification.
2) The majority of the risk-associated SNPs modulate the affinity of chromatin for genes at distal regulatory elements, thereby resulting in allele-specific gene expression. 
3) SNPs contained within the gene interact multiplicatively with mutations increasing cancer risk.
