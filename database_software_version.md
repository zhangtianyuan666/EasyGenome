# Bioinformatics Analysis Resources

## 1. Software

本流程使用的软件及其主要用途如下： The software used in this workflow and their primary applications are summarized below

| module   | software      | 版本参数Version   | Function                                                                         | 功能                                                        | Github                                          |
|:---------|:--------------|:--------------|:---------------------------------------------------------------------------------|:----------------------------------------------------------|:------------------------------------------------|
| module1  | FastP         | v0.23.2       | High-performance data filtering and preprocessing for FASTQ data                 | 数据过滤与预处理                                                  | https://github.com/OpenGene/fastp               |
| module1  | kraken2       | v2.1.2        | Taxonomic classification for short- and long-read sequencing data                | 基于 k-mer 比对算法，对二代短读长、三代长读长测序数据进行高精度物种分类注释                 | https://github.com/DerrickWood/kraken2          |
| module1  | krona         | -             | An interactive visualization tool for species classification results.            | 用于物种分类结果的交互式可视化工具                                         | https://github.com/marbl/Krona                  |
| module1  | Porechop      | v0.2.4        | Adapter trimming for long-read sequencing data (ONT/PacBio)                      | 针对 ONT/PacBio 三代长读长测序数据，精准识别并去除测序接头序列                     | https://github.com/rrwick/Porechop              |
| module1  | NanoPlot      | v1.25.0       | Visualization of third-generation sequencing data quality metrics                | 可视化展示三代测序数据质控指标                                           | https://github.com/wdecoster/NanoPlot           |
| module2  | Jellyfish     | 2.2.10        | K-mer frequency statistics for genome/sequencing data                            | 高效计算基因组 / 测序数据的 k-mer 频率分布                                | https://github.com/gmarcais/Jellyfish           |
| module2  | Unicycler     | v0.5.1        | Hybrid genome assembly (short + long reads) for bacterial genomes                | 整合二代短读长和三代长读长数据，实现细菌基因组的混合组装，生成高质量连续基因组序列                 | https://github.com/rrwick/Unicycler             |
| module2  | GenomeScope   | v2            | Genome size estimation                                                           | 估算物种基因组大小、杂合度、重复序列比例等核心基因组特征                              | https://github.com/schatzlab/genomescope        |
| module2  | Flye          | v2.9.2        | De novo genome assembly for long-read sequencing data (ONT/PacBio)               | 针对 ONT/PacBio 三代长读长数据的从头组装工具，适配高重复、高杂合基因组，生成连续的组装序列       | https://github.com/mikolmogorov/Flye            |
| module2  | minimap2      | v2.24         | Efficient alignment of long-read to a reference genome                           | 三代测序深度统计                                                  | https://github.com/lh3/minimap2                 |
| module2  | racon         | v2.0          | Genome polishing for long-read sequencing assemblies                             | 基于三代长读长数据对基因组组装序列进行纠错优化，提升组装序列的碱基准确性                      | https://github.com/isovic/racon                 |
| module2  | pilon         | v1.24         | Genome polishing for short-read sequencing assemblies                            | 利用二代短读长数据对基因组组装序列进行精细纠错，修正单碱基错误和小片段插入缺失                   | https://github.com/broadinstitute/pilon/wiki    |
| module2  | bwa           | v0.7.12       | Efficient alignment of short-read data to a reference genome                     | 采用 BWT 算法将二代短读长测序 reads 高效比对到参考基因组，生成 SAM/BAM 格式比对结果      | https://github.com/lh3/bwa                      |
| module2  | samtools      | v1.17         | Sequencing depth statistics and BAM/CRAM file manipulation                       | 统计二代测序数据的基因组覆盖深度，支持 BAM/CRAM 格式文件的排序、索引、格式转换等操作           | https://github.com/samtools/samtools            |
| module2  | PlasFlow      | v1.1          | Genome polish for long-read data                                                 | 基于机器学习模型，从基因组组装序列中识别并区分质粒序列与染色体序列                         | https://github.com/smaegol/PlasFlow             |
| module2  | Quast         | v5.2.0        | Quality evaluation and comparison of genome assemblies                           | 全面评估基因组组装质量（N50、Contig 数量、基因组完整性等）                        | https://github.com/ablab/quast                  |
| module2  | checkm        | v1.1.3        | Assess completeness and contamination of microbial genomes                       | 基于单拷贝标记基因，评估细菌 / 古菌基因组组装的完整性和污染度                          | https://github.com/Ecogenomics/CheckM           |
| module2  | Busco         | v5.4.0        | Assess genome/transcriptome/protein completeness with conserved orthologs        | 基于进化保守的单拷贝直系同源基因集，评估基因组组装、转录组组装或蛋白序列的完整性                  | https://gitlab.com/ezlab/busco/-/releases#6.0.0 |
| module2  | Merqury       | v1.3          | Evaluate genome assembly quality                                                 | 评估基因组组装质量，检测组装错误、重复序列、杂合区域等                               | https://github.com/marbl/merqury                |
| module2  | bedtools      | v2.30.0       | Process and analyze genome interval data (e.g., BED, GFF, VCF formats)           | 处理和分析基因组区间数据（BED/GFF/VCF 格式）                              | https://github.com/arq5x/bedtools2              |
| module3  | Prokka        | v1.14.5       | Rapid prokaryotic genome annotation (gene prediction & functional annotation)    | 快速完成原核生物基因组的基因预测、功能注释，生成 GFF3/GBK 等标准化注释文件                | https://github.com/tseemann/prokka              |
| module3  | Bakta         | v1.9.3        | Comprehensive prokaryotic genome annotation with database integration            | 整合多个功能数据库，全面注释原核生物基因组（基因预测、功能注释、抗性基因、毒力因子等）               | https://github.com/oschwengers/bakta            |
| module3  | pseudofinder  | v1.1.0        | Pseudogene annotation and characterization in prokaryotic genomes                | 识别并注释原核生物基因组中的假基因，分析假基因的来源、突变类型及分布特征                      | https://github.com/filip-husnik/pseudofinder    |
| module3  | RepeatMasker  | v4.1.2        | Repeat sequence annotation and classification in genomes                         | 识别并分类基因组中的重复序列（串联重复、散在重复等），注释重复序列的位置、长度和类型                | https://github.com/Dfam-consortium/RepeatMasker |
| module3  | EggNOG-mapper | v2.1.9        | Functional annotation based on EggNOG orthology database                         | 基于 EggNOG 直系同源基因数据库，对蛋白序列进行功能注释                           | https://github.com/eggnogdb/eggnog-mapper       |
| module3  | Diamond       | v2.0.14       | High-speed protein alignment for functional annotation                           | 基于 BLASTP 算法优化的高速蛋白序列比对工具，用于功能注释的数据库比对                    | https://github.com/bbuchfink/diamond            |
| module3  | hmmscan       | v2.1.9        | HMM domain annotation using Pfam/TIGRFAMs databases                              | 利用 Pfam/TIGRFAMs 等 HMM 模型数据库，对蛋白序列进行结构域注释，识别功能保守结构域       | https://github.com/eggnogdb/eggnog-mapper       |
| module3  | antiSMASH     | v8.0.2        | Identification and annotation of secondary metabolite biosynthetic gene clusters | 识别并注释基因组中的次级代谢产物合成基因簇（如抗生素、毒素、色素等）                        | https://github.com/antismash/antismash          |
| module3  | BioMGCore     | -             | Statistical analysis of antiSMASH annotation results                             | 统计 antiSMASH 注释的次级代谢产物基因簇数量、类型、分布位置等                      | https://github.com/xielisos567/BioMGCore        |
| module3  | island        | v2.0          | Prediction of genomic islands in bacterial genomes                               | 识别细菌基因组中的基因组岛（Genomic Island），分析其位置、长度、GC 含量及功能注释         | https://github.com/oasisfeng/island             |
| module3  | PhiSpy        | v3.7.8        | Identification of lysogenic prophages in bacterial/archaeal genomes              | 识别细菌 / 古菌基因组中的溶源噬菌体序列，注释噬菌体整合位点、基因组成及潜在功能                 | https://github.com/linsalrob/PhiSpy             |
| module3  | minced        | v0.4.2        | CRISPR-Cas locus identification and annotation in microbial genomes              | 识别并注释微生物基因组中的 CRISPR-Cas 位点，确定重复序列（Repeat）和间隔序列（Spacer）特征 | https://github.com/ctSkennerton/minced          |
| module4  | hmmscan       | v2.1.9        | HMM domain annotation using Pfam/TIGRFAMs databases                              | 利用 Pfam/TIGRFAMs 等 HMM 模型数据库，对蛋白序列进行结构域注释，识别功能保守结构域       | https://github.com/eggnogdb/eggnog-mapper       |
| module4  | Diamond       | v2.0.14       | High-speed protein alignment for functional annotation                           | 基于 BLASTP 算法优化的高速蛋白序列比对工具，用于功能注释的数据库比对                    | https://github.com/bbuchfink/diamond            |
| module4  | Blastp        | v2.16.0       | Protein sequence alignment for PHI database annotation                           | 利用 BLASTP 算法将蛋白序列比对到 PHI 数据库                              | https://github.com/ncbi/blast_plus_docs         |
| module4  | rgi           | v5.2.1        | Functional annotation based on CARD database(antibiotic resistance)              | 基于 CARD 数据库注释抗生素抗性基因，预测抗性机制、抗生素类别及抗性基因型                   | https://github.com/arpcard/rgi                  |
| module4  | signalp6      | v6.0          | Prediction of signal peptides in bacterial/archaeal proteins                     | 预测细菌 / 古菌蛋白序列中的信号肽（Signal Peptide）                        | https://github.com/fteufel/signalp-6.0          |
| module4  | tmhmm         | v2.0          | Prediction of transmembrane helices in protein sequences                         | 预测蛋白序列中的跨膜螺旋结构（Transmembrane Helices）                     | https://github.com/richelbilderbeek/tmhmm       |
| module4  | EffectiveT3   | v1.0.1        | Annotation of Type III secretion system (T3SS) effector proteins                 | 注释细菌 Ⅲ 型分泌系统（T3SS）效应蛋白                                    | https://github.com/sm18lr88/OpenAI_TTS_GUI      |
| module5  | circlize      | v0.4.15       | Generate comprehensive circular genome visualization plots                       | 绘制全基因组综合圈图                                                | https://github.com/jokergoo/circlize            |
| module5  | cgview        | v2.0.3        | Generate circular genome plots for long-read assembly results                    | 针对三代长读长组装结果绘制基因组圈图                                        | https://github.com/paulstothard/cgview          |
| module6  | orthofinder   | v2.5.5        | Identification of orthologs and comparative genomics analysis                    | 识别多个物种间的直系同源基因 / 旁系同源基因                                   | https://github.com/davidemms/OrthoFinder        |
| module6  | jcvi          | v1.5.6        | Synteny analysis and visualization between multiple genomes                      | 开展多基因组共线性分析                                               | https://github.com/tanghaibao/jcvi              |
| module6  | seqkit        | v2.10.0       | Comprehensive sequence statistics and manipulation for FASTA/FASTQ               | 对 FASTA/FASTQ 格式序列进行多维度统计                                 | https://github.com/shenwei356/seqkit            |
| module6  | gtdbtk        | v2.3.0        | Taxonomic classification of bacterial/archaeal genomes                           | 基于通用单拷贝标记基因，对细菌 / 古菌基因组进行标准化物种分类注释                        | https://github.com/Ecogenomics/GTDBTk           |
| module6  | Pyani         | v0.2.13       | Calculate Average Nucleotide Identity (ANI) between microbial genomes            | 计算微生物基因组之间的平均核苷酸一致性（ANI）                                  | https://github.com/widdowquinn/pyani            |
| module6  | samtools      | v1.17         | Sequencing depth statistics and BAM/CRAM file manipulation                       | 统计二代测序数据的基因组覆盖深度，支持 BAM/CRAM 格式文件的排序、索引、格式转换等操作           | https://github.com/samtools/samtools            |
| module6  | minimap2      | v2.24         | Efficient alignment of long-read to a reference genome                           | 三代测序深度统计                                                  | https://github.com/lh3/minimap2                 |
| module6  | syri          | v1.7.1        | Identification of synteny and genome rearrangement between two assemblies        | 比较两个基因组组装结果，识别其共线性和基因组重排事件                                | https://github.com/schneebergerlab/syri         |
| module7  | roary         | v3.9.1        | Construction of bacterial pan-genome profiles                                    | 构建细菌泛基因组图谱                                                | https://sanger-pathogens.github.io/Roary/       |
| module8  | mlst          | v2.9          | Multi-locus sequence typing (MLST) for bacterial strain classification           | 基于多位点序列分型（MLST）对细菌菌株进行初步分类                                | https://github.com/tseemann/mlst                |
| module8  | chewBBACA     | v3.4          | Bacterial whole-genome multi-locus sequence typing (wgMLST)                      | 执行细菌全基因组多位点序列分型（wgMLST）                                   | https://github.com/tseemann/mlst                |


## 2. Databases

本流程使用的数据库及其主要用途如下：The databases used in this workflow are listed below

<div align="center">

<table width="80%" align="center">
  <thead>
    <tr style="text-align: right;">
      <th>module</th>
      <th>database</th>
      <th>版本参数Version</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>module1</td>
      <td>kraken2</td>
      <td>2025-10-15</td>
    </tr>
    <tr>
      <td>module2</td>
      <td>busco</td>
      <td>2024-11-14</td>
    </tr>
    <tr>
      <td>module2</td>
      <td>checkm</td>
      <td>2015-01-16</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>Bakta</td>
      <td>2024-01-19</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>nr_cluster</td>
      <td>2026-01-05</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>eggnog</td>
      <td>V5.0</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>Uniprot</td>
      <td>2025-10-15</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>tigerfam</td>
      <td>2025-08-06</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>antismash</td>
      <td>V4.0</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>cazy</td>
      <td>2025-08-26</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>tcdb</td>
      <td>2025-08-06</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>phi</td>
      <td>2025-05-01</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>card</td>
      <td>2025-05-29</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>cyped</td>
      <td>2025-10-15</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>VFDB</td>
      <td>2024-11-14</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>ICE</td>
      <td>2015-01-16</td>
    </tr>
    <tr>
      <td>module6</td>
      <td>gtdbtk</td>
      <td>2024-01-19</td>
    </tr>
  </tbody>
</table>

</div>



## 2. Databases

本流程使用的python包如下： The Python packages used in this workflow are listed below

<div align="center">

<table width="80%" align="center">
  <thead>
    <tr style="text-align: right;">
      <th>module</th>
      <th>python packages</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>module1</td>
      <td>sys</td>
    </tr>
    <tr>
      <td>module1</td>
      <td>argparse</td>
    </tr>
    <tr>
      <td>module1</td>
      <td>os</td>
    </tr>
    <tr>
      <td>module1</td>
      <td>subprocess</td>
    </tr>
    <tr>
      <td>module1</td>
      <td>json</td>
    </tr>
    <tr>
      <td>module1</td>
      <td>glob</td>
    </tr>
    <tr>
      <td>module1</td>
      <td>gzip</td>
    </tr>
    <tr>
      <td>module1</td>
      <td>statistics</td>
    </tr>
    <tr>
      <td>module1</td>
      <td>typing</td>
    </tr>
    <tr>
      <td>module1</td>
      <td>Bio</td>
    </tr>
    <tr>
      <td>module2</td>
      <td>sys</td>
    </tr>
    <tr>
      <td>module2</td>
      <td>Bio</td>
    </tr>
    <tr>
      <td>module2</td>
      <td>argparse</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>argparse</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>collections</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>os</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>glob</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>gzip</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>re</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>sys</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>bs4</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>pandas</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>sys</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>os</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>argparse</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>re</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>pandas</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>matplotlib</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>collections</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>pexpect</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>Bio</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>time</td>
    </tr>
    <tr>
      <td>module5</td>
      <td>sys</td>
    </tr>
    <tr>
      <td>module5</td>
      <td>re</td>
    </tr>
    <tr>
      <td>module7</td>
      <td>matplotlib</td>
    </tr>
    <tr>
      <td>module7</td>
      <td>pandas</td>
    </tr>
    <tr>
      <td>module7</td>
      <td>sys</td>
    </tr>
    <tr>
      <td>module7</td>
      <td>os</td>
    </tr>
    <tr>
      <td>module7</td>
      <td>copy</td>
    </tr>
    <tr>
      <td>module7</td>
      <td>itertools</td>
    </tr>
    <tr>
      <td>module7</td>
      <td>numpy</td>
    </tr>
    <tr>
      <td>module8</td>
      <td>sys</td>
    </tr>
    <tr>
      <td>module8</td>
      <td>re</td>
    </tr>
  </tbody>
</table>

</div>



## 2. Databases

本流程使用的R包如下：The R packages used in this workflow are listed below

<div align="center">

<table width="80%" align="center">
  <thead>
    <tr style="text-align: right;">
      <th>module</th>
      <th>R</th>
      <th>版本参数Version</th>
      <th>用途（统计/绘图）Purpose (Statistical Analysis/Visualization)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>module2</td>
      <td>ggplot2</td>
      <td>v3.4.2</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>dplyr</td>
      <td>v1.1.2</td>
      <td>统计Statistical</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>ggthemes</td>
      <td>v4.2.4</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>svglite</td>
      <td>v2.1.1</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>gridExtra</td>
      <td>v2.3</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>ggstance</td>
      <td>v0.3.6</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>UpSetR</td>
      <td>v1.4.0</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>RColorBrewer</td>
      <td>v1.1.3</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module3</td>
      <td>svglite</td>
      <td>v2.1.1</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>svglite</td>
      <td>v2.1.1</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>RColorBrewer</td>
      <td>v1.1.3</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>ggplot2</td>
      <td>v3.4.2</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>stringr</td>
      <td>v1.5.0</td>
      <td>统计Statistical</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>magrittr</td>
      <td>v2.0.3</td>
      <td>统计Statistical</td>
    </tr>
    <tr>
      <td>module4</td>
      <td>dplyr</td>
      <td>v1.1.2</td>
      <td>统计Statistical</td>
    </tr>
    <tr>
      <td>module5</td>
      <td>stringr</td>
      <td>v1.5.0</td>
      <td>统计Statistical</td>
    </tr>
    <tr>
      <td>module5</td>
      <td>circlize</td>
      <td>v0.4.15</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module5</td>
      <td>ComplexHeatmap</td>
      <td>v2.10.0</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module5</td>
      <td>grid</td>
      <td>v4.1.0</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module5</td>
      <td>svglite</td>
      <td>v2.1.1</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module6</td>
      <td>VennDiagram</td>
      <td>v1.7.3</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module6</td>
      <td>UpSetR</td>
      <td>v1.4.0</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module6</td>
      <td>RColorBrewer</td>
      <td>v1.1.3</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module6</td>
      <td>dplyr</td>
      <td>v1.1.2</td>
      <td>统计Statistical</td>
    </tr>
    <tr>
      <td>module6</td>
      <td>svglite</td>
      <td>v2.1.1</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module7</td>
      <td>ggplot2</td>
      <td>v3.4.2</td>
      <td>绘图Visualization</td>
    </tr>
    <tr>
      <td>module7</td>
      <td>reshape2</td>
      <td>v1.4.4</td>
      <td>统计Statistical</td>
    </tr>
  </tbody>
</table>

</div>


