# LiuXiaole_learning
刘小乐课程学习笔记

网址：https://wukong.tongji.edu.cn/new_home/sjshen/learning/2022_harvard_lecture/ 
练习：/home/share/resource/lecture/lxl_wcf_train_lecture2022/Training_homework_data

## 主要学习了比对分析：STAR、bwa等
tophat调用bowtie，tophat2调用bowtie2
HISAT2既能做RNA-seq也能做DNA-seq
STAR是ENCODE官方推荐的RNA-seq比对工具
**DNA-seq：bowtie；bowtie2；BWA**
**RNA-seq：STAR；HISAT2；Tophat**
## RNA-seq数据分析：RMSE、salmon、RSeQC(
## 去除批次：Combat(sva)、DESeq2
## 降维聚类：PCA(prcomp)、hclust(dist)、kmeans
## 差异表达：DESeq2、limma
## 分类聚类回归等：MASS。
“MASS“包是一个十分强大的统计包，可以进行各种统计分析，我也将围绕它来介绍判别分析。”MASS“包既可以进行线性判别，也可以进行二次判别。
## ChIP-Seq:bwa比对、samtools查看比对结果、macs2（filterdup、callpeak、差异峰）、homer(TF motif)
Cistrome DB总共收录了30451人和26013小鼠的转录因子、组蛋白修饰和染色质可及性样本，可以说是目前最全面的研究ChIP-seq和DNase-seq的数据库
Cistrome Toolkit 可以检索哪个因子调控了我们感兴趣的基因，哪个因子结合在用户的基因组区域或者和用户的peak有显著重叠。
Cistrome-GO是用于对转录因子ChIP-seq峰进行功能富集分析的网络服务。它有两种工作模式。如果用户同时提供了TF的ChIP-seq峰文件和差异表达分析文件（基于TF），则Cistrome-GO将基于两种数据类型的整合执行集成模式分析。如果我们仅上传TF ChIP-seq峰文件，则Cistrome-GO将以单独模式执行分析。
## ATAC-seq:bedtools取交并补差集、
Lisa主要是利用来自人和小鼠的DNase-seq与H3K27ac和 ChIP-seq的数据，来确定导致差异表达的基因集的转录因子和调控因子。
## Chromatin Modification：bismark比对
Bismark是一款比较有名的甲基化测序比对软件。对 Bismark 的比对结果**去重**，去除以相同方向比对到相同位置的 reads。**建议应用于 WGBS 数据，不适用与 RRBS ，amplicon，以及 target 富集文库的数据**。
## HiC:runHiC(mapping、filtering、binning、quality)、cooler标准化、peakachu寻找显著loop、
## 隐马尔可夫模型：状态序列（Hidden States）、观测序列（Observations）、转移概率（Transition Probabilities）、发射概率（Emission Probabilities）、初始概率分布（Initial Probabilities）
## GWAS：NHGRI-EBI GWAS、dbSNP、LDLink、RegulomeDB
## MASESTRO:scrna-init、sample-init scatac-init 、snakemake
## Seurat、Signac
## maf文件、gain or loss of function mutation
CBioPortal has a comprehensive list of tumor profiling results for interactive visualization.
## CRISPR screens: MAGeCK（count、test）、
正向选择的基因是在处理组中敲除或抑制会给细胞提供生长或存活优势的基因。
负向选择的基因是在处理组中敲除或抑制导致细胞生长劣势或死亡的基因。
depmap (DepMap.org)，其中有 500 多个人类细胞系的 CRISPR / RNAi 筛选结果.查询基因在各种细胞的中的必要性 （Perturbation effects），换而言之，在这个细胞中去除这个基因是否影响细胞生存，再进一步在这个细胞中敲除或者敲低该基因是否会影响细胞活力
TIMER是一个综合数据库，主要功能是通过TIMER算法系统地分析不同癌症类型中的六种肿瘤浸润免疫细胞
