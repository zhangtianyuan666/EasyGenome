# Easygenome易基因组

  # Author作者: Tianyuan Zhang(张天缘), Yong-xin Liu(刘永鑫), Yilin Li(李伊琳) Zhihao Zhu(朱志豪), Qingrun Xue(薛清润) et al.
  # Update更新时间: 2025-10-28
  # Version版本: 1.0 

## Step01.data_preparation_and_quanlity_control

  # 创建分析目录  建议创建在Easygenome目录下，以样本名命名 Create an analysis directory. It is recommended to create it in the Easygenome directory and name it after the sample.
  mkdir -p /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata  
  # 提取二代测序5000条reads Extract 5000 short-reads
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/sub_fq.py  /data6/zhangtianyuan/Pipeline/EasyGenome/data/SRR32313567_R1.fastq.gz 5000 /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/Read.5000.fq.gz
  # 使用kraken2对二代数据进行物种分类 - 判断是否有污染 Species classification of second-generation data using Kraken2 - Determine whether there is contamination.
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/kraken2_2.1.2-no-db.sif kraken2 --db /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/Kraken2/  --threads 48 --output /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567_kraken2.out.txt --report /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567_kraken2.report.txt --confidence 0.1 /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/Read.5000.fq.gz
  # 使用krona可视化kraken2分类结果 Visualize Kraken2 classification results using Krona.
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/kreport2krona.py -r /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567_kraken2.report.txt -o /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/krona
  cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata;singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/krona_2.7.1--e7615f7.sif ktImportText krona;cd ../../

  # 提取三代长读长测序数据，然后进行质控，使用kraken2物种分类 Extract third-generation long-read sequencing data, perform quality control, and conduct species classification using Kraken2.
  # 使用fastp过滤短读长数据并进行统计 Filter short-read data using fastp and perform statistical analysis.
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/fastp -i /data6/zhangtianyuan/Pipeline/EasyGenome/data/SRR32313567_R1.fastq.gz -I /data6/zhangtianyuan/Pipeline/EasyGenome/data/SRR32313567_R2.fastq.gz -w 16 -o /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567_R1.clean.fq.gz -O /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567_R2.clean.fq.gz -j /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.json -h /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.html
  # 进入分析目录后进行统计 Perform statistics after entering the analysis directory
  cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/fastp_stats.py
  cd ../../
  # 对三代长读长数据进行去接头处理 Perform adapter trimming on long-read data.
  singularity exec -B /data6/  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/porechop_v0.2.4.sif porechop -i /data6/zhangtianyuan/Pipeline/EasyGenome/data/SRR32313567.nanopore.fastq.gz  -o /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.Adapter.fastq -t 30
  # 三代长读长数据质控，保留1600bp及以上的reads、结果统计 Perform quality control on long-read data, retain reads of 1600bp or longer, and generate statistical results.
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/fastq_to_fastq.py --min_length 1600 /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.Adapter.fastq > /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered.min1600.fastq

  # 选取1G数据用于分析 Select 1G data for analysis  
  # 此处也可设置为500Mb，最佳组装深度为50-200X  This can also be set to 500Mb, with the optimal assembly depth being 50-200X.
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/select_1G_ont.py /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered.min1600.fastq   /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered.fastq 
  # 统计结果 Statistical results
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/seqkit stats /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered.fastq -a |awk -F' *' 'BEGIN{OFS="\t"}NR==1{print $1,$4,$5,$6,$7,$8,$13}NR==2{n=split($1,a,"[/]" );print a[n],$4,$5,$6,$7,$8,$13}'> /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/tmp_ONT_stats.xls
  # NanoPlot可视化 NanoPlot Visualization
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/nanoplot_latest.sif  NanoPlot --fastq /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered.fastq -o /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/fastq-plots --plots hex dot -t 48
  
  # 使用kraken2对三代数据进行物种分类 Kraken2 species classification long-read data (Extract 5000 reads)
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/sub_fq.py /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered.fastq 5000 /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered_5000.fastq
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/kraken2_2.1.2-no-db.sif kraken2 --db /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/Kraken2/  --threads 48 --output /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered_5000_kraken2.out.txt --report /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered_5000_kraken2.report.txt  --confidence 0.1 /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered_5000.fastq
  
  # krona可视化kraken2分类结果 krona visualization of kraken2 classification results
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/kreport2krona.py -r /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered_5000_kraken2.report.txt -o /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567_filtered_5000_krona
  cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata;singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/krona_2.7.1--e7615f7.sif ktImportText SRR32313567_filtered_5000_krona;cd ../../

## Step02.assembly_and_evaluation

  # 创建目录  Create directory
  mkdir -p /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/02.assembly;cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/02.assembly
  # 使用jellyfish进行kmer分析 Using Jellyfish for K-mer analysis
  # Generate configuration file
  echo "gzip -dc ../01.cleandata/SRR32313567_R1.clean.fq.gz" > generators
  echo "gzip -dc ../01.cleandata/SRR32313567_R2.clean.fq.gz" >> generators
  
  # 使用jellyfish计算count	Using Jellyfish to Calculate Count
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/jellyfish_v2.2.10.sif jellyfish count -t 30 -C -m 15 -s 1G -o kmer15 -Q 5 -g generators -G 4
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/jellyfish_v2.2.10.sif jellyfish histo -v -o kmer15.histo kmer15 -t 30 -h 100000
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/jellyfish_v2.2.10.sif jellyfish stats kmer15 -o kmer15.stat
  # 使用genomescope评估基因组大小 Assessing genome size using genomescope
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/ime_genomescope_v1.0.0.sif /opt/genomescope/genomescope.R kmer15.histo 15 150 ./genomescope 100000 verbose > genomescope_print.txt
  
  # 下面3个命令运行1个 The following three commands run 1
  # 若同时使用长短读长数据，则使用该命令拼接(For Hybrid Long read and Short read Data)
  # 使用unicycler 混合拼接Long-read 和 Short-read 数据 assembly genome using Unicycler by Hybrid Long read and Short read Data strategy
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/unicycler_v0.5.1.sif unicycler -1 ../01.cleandata/SRR32313567_R1.clean.fq.gz -2 ../01.cleandata/SRR32313567_R2.clean.fq.gz -l ../01.cleandata/SRR32313567.filtered.fastq -o unicycler_out -t 46 --keep 1 --mode conservative
  
  # 若仅使用短读长拼接，则使用该命令 (For Short-read Only data)
  # 使用unicycler对Short-read-Only数据进行拼接 Assembling Short-read Only data using unicycler
   singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/unicycler_v0.5.1.sif unicycler -1 ../01.cleandata/SRR32313567_R1.clean.fq.gz-2 ../01.cleandata/SRR32313567_R2.clean.fq.gz -o unicycler_short -t 46 --keep 1 --mode conservative
  
  # 若仅有Long-read数据，则使用该命令进行拼接 (For Long-read Only data)
  # 使用flye拼接基因组 Using Flye to assembly Genomes
   /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/flye_2.9.2.sif flye --nano-raw ../01.cleandata/SRR32313567.filtered.fastq --out-dir flye --threads 40 --iterations 5
  # 使用racon对拼接结果进行矫正  Correct the draft genome using racon
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/minimap2_2.24.sif minimap2 -d ref1.mmi flye/assembly.fasta
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/minimap2_2.24.sif minimap2 -ax map-ont -t 40 ref1.mmi ../01.cleandata/SRR32313567.filtered.fastq > map1.sam
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/canu-racon_2.0.sif racon -t 40 ../01.cleandata/SRR32313567.filtered.fastq map1.sam flye/assembly.fasta > flye/racon1.fa
  
  # 将目标结果进行链接，进行后续分析 Link the target assembly results for subsequent analysis
  # link assembly_result   下面2个命令运行1个 The following two commands run 1
  ln -s unicycler_out/assembly.fasta input.fa
  ln -s flye/racon1.fa input.fa
  
  # 使用pilon纠错  比完直接输出sorted.bam，纠错2轮 Output sorted.bam directly after comparison, and perform 2 rounds of error correction
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bwa-samtools_0.7.12_1.2.1.sif bwa index input.fa
  
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bwa-samtools_0.7.12_1.2.1.sif bwa mem -t 40 input.fa ../01.cleandata/SRR32313567_R1.clean.fq.gz ../01.cleandata/SRR32313567_R2.clean.fq.gz | singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools sort -@ 40 -T tmp.illumina -o illumina.sort.bam
  
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools index illumina.sort.bam
  java -jar /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/pilon-1.24.jar --fix bases --genome input.fa --frags illumina.sort.bam --diploid --outdir ./  --output genome 
  
  # 第二次使用pilon纠错 error correction again
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bwa-samtools_0.7.12_1.2.1.sif bwa index genome.fasta
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bwa-samtools_0.7.12_1.2.1.sif bwa mem -t 40 genome.fasta ../01.cleandata/SRR32313567_R1.clean.fq.gz ../01.cleandata/SRR32313567_R2.clean.fq.gz | singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools sort -@ 40 -T tmp.illumina -o illumina.sort.bam
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools index illumina.sort.bam
  java -jar /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/pilon-1.24.jar --fix bases --genome genome.fasta --frags illumina.sort.bam --diploid --outdir ./  --output assembly
  
  # 这里，我们得到了assembly.fasta（组装结果），用于后续分析。如果没有二代数据，则运行下一步 Here, we have obtained assembly.fasta (assembly results) for subsequent analysis. If there is no second-generation data, proceed to the next step
  # ln -s flye/racon1.fa assembly.fasta
  
  # 质粒鉴定 Plasmid identification
  mkdir PlasFlow;cd PlasFlow
  singularity exec  --env LC_ALL=C,LANG=C    -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plasflow_latest_zty.sif PlasFlow.py --input ../assembly.fasta --output SRR32313567
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/stat_plasflow.py SRR32313567
  cd ../
  
    
  # 确定组装版本  Determine assembly version
  mkdir -p assemble
  mv assembly.fasta assemble/assembly.fasta
  
  # 使用QUAST对组装结果进行评估与统计 Evaluate and analyze the assembly results using QUAST.
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/quast_5.2.0.sif quast.py -o quast -t 40 assemble/assembly.fasta
  # 使用CheckM对组装结果进行评估与统计 Evaluate and analyze the assembly results using CheckM. 
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/checkm.v1.1.3.sif checkm lineage_wf unicycler_out/ checkmout/ -x fasta -t 48  --pplacer_threads 8 --tab_table -f checkmout/checkm.txt
  # 使用BUSCO对组装结果进行评估与统计 Evaluate and analyze the assembly results using BUSCO.
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/staphb_busco_6.0.0-prok-bacteria_odb12_2024-11-14.sif  busco -o busco -i assemble/assembly.fasta -l /busco_downloads/lineages/bacteria_odb12  -m geno -c 40
  
  # QV准确度评估  QV accuracy evaluation.
  # 创建目录并进行分析
  mkdir merqury_qv ;cd merqury_qv
  singularity exec -B /data6/ --env MERQURY=/usr/local/share/merqury/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/merqury_v1.3.sif meryl count k=21 output read1.meryl ../../01.cleandata/SRR32313567_R1.clean.fq.gz
  singularity exec -B /data6/ --env MERQURY=/usr/local/share/merqury/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/merqury_v1.3.sif meryl count k=21 output read2.meryl ../../01.cleandata/SRR32313567_R2.clean.fq.gz
  singularity exec -B /data6/ --env MERQURY=/usr/local/share/merqury/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/merqury_v1.3.sif meryl count k=21 output read3.meryl ../../01.cleandata/SRR32313567.filtered.fastq
  singularity exec -B /data6/ --env MERQURY=/usr/local/share/merqury/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/merqury_v1.3.sif meryl  union-sum output illumina.meryl read1.meryl read2.meryl 
  singularity exec -B /data6/ --env MERQURY=/usr/local/share/merqury/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/merqury_v1.3.sif meryl  union-sum output all.meryl illumina.meryl read3.meryl
  singularity exec -B /data6/ --env MERQURY=/usr/local/share/merqury/   /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/merqury_v1.3.sif    bash /usr/local/share/merqury/merqury.sh all.meryl ../assemble/assembly.fasta  assembly
  cd ..
  
  # 深度统计 Depth statistics
  # 二代测序深度统计(如有二代数据) depth statistics for short-reads 
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools faidx assemble/assembly.fasta
  awk -F'\t' '{print $1"\t"$2}' assemble/assembly.fasta.fai  > ctg_length.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bwa-samtools_0.7.12_1.2.1.sif bwa index assemble/assembly.fasta 
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bwa-samtools_0.7.12_1.2.1.sif bwa mem -t 40 assemble/assembly.fasta ../01.cleandata/SRR32313567_R1.clean.fq.gz ../01.cleandata/SRR32313567_R2.clean.fq.gz > aln2.sam
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools sort aln2.sam > aln2sorted.sam
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bedtools_v2.30.0.sif bedtools makewindows -g ctg_length.xls  -w 2000 -i winnum > genome_windows_2k.bed
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bedtools_v2.30.0.sif bedtools coverage -mean -sorted -bed -g ctg_length.xls -a genome_windows_2k.bed -b aln2sorted.sam > aln2_2kbin.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools depth -a aln2sorted.sam > aln2_depth.xls
  
  # 三代测序深度统计（如有三代数据） Depth statistics of third-generation sequencing
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/minimap2_2.24.sif minimap2 -ax map-ont -t 40 assemble/assembly.fasta ../01.cleandata/SRR32313567.filtered.fastq -o aln3.sam
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools sort aln3.sam > aln3sorted.sam
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bedtools_v2.30.0.sif bedtools coverage -mean -sorted -bed -g ctg_length.xls -a genome_windows_2k.bed -b aln3sorted.sam > aln3_2kbin.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools depth -a aln3sorted.sam > aln3_depth.xls
  
  # GC含量及GC-skew计算 calculate GC content and GC-skew
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python3 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/GC.py -f assemble/assembly.fasta -w 2000 -s 2000 > GCstat.xls
  
  
  # GC_depth 图绘制   Plot GC_depth graph
  paste GCstat.xls aln2_2kbin.xls |awk '{print $1":"$2"-"$3"\t"$4"\t"$10}' >aln2_GC_depth.xls
  paste GCstat.xls aln3_2kbin.xls |awk '{print $1":"$2"-"$3"\t"$4"\t"$10}' >aln3_GC_depth.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/GC_depth.r aln2_GC_depth.xls aln2_GC_depth
  singularity exec -B  /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/GC_depth.r aln3_GC_depth.xls aln3_GC_depth
  cd ../


## Step03.Basic_Structual_and_Functional_Annotation

  # 创建目录 Create directory
  mkdir -p /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/03.anno;
  cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/03.anno/
  # 使用prokka进行基因组注释及统计  Using Prokka for genome annotation and statistics
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/prokka_1.14.5.sif prokka --force --outdir prokka_out --kingdom Bacteria --gcode 11 --cpus 40 --evalue 1e-6 --rnammer --rfam --genus Genus --species species --strain strain --compliant --locustag ctg --centre Bacteria --prefix SRR32313567 ../02.assembly/assemble/assembly.fasta
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python3 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/stat_prokka.py prokka_out/SRR32313567.gff prokka_out/sturcAnno.stat.xls
  
  # 使用Bakta注释 和prokka二选一即可，本示例用的是prokka的注释结果 You can choose between Bakta annotation and Prokka. This example uses the annotation result of Prokka
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bakta_v1.9.3.sif  bakta  --db /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/baketa/db-light  --verbose --threads 16 --output ./ ../02.assembly/assemble/assembly.fasta --output bakta_out 
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif  python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/stat_bakta.py  bakta_out/assembly.gff3 bakta_out/sturcAnno.stat.xls
  
  
  # 假基因注释 Pseudogene annotation
  mkdir pseudo_gene;cd pseudo_gene
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/pseudofinder-master.sif python /opt/pseudofinder-master/pseudofinder.py annotate -g ../prokka_out/SRR32313567.gbk --outprefix SRR32313567 -di -skpdb --threads SRR32313567 -db /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/nr/B_A_V_nr.fa.dmnd 
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/fasta_stat.py -all SRR32313567_pseudos.fasta --out SRR32313567.ppseq.stat
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/all_stat.py 
  cd ../
  # 重复序列注释 Repeat sequence annotation
  mkdir -p repeat;cd repeat
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/tetools_1.4.sif RepeatMasker -pa 10 -q -html -gff -dir ./ ../prokka_out/SRR32313567.fna
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/repeat_stat.py  SRR32313567.fna.out SRR32313567.fna.out.gff  SRR32313567  --genome ../prokka_out/SRR32313567.fna
  cd ../
  # eggnog注释 eggnog annotation for functional research
  mkdir eggnog;cd eggnog
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/eggnog-mapper_2.1.9.sif emapper.py -i ../prokka_out/SRR32313567.faa --cpu 40 --data_dir /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/eggnog/ --output SRR32313567
  cd ../
  # cog注释及作图 cog annotation and plot
  mkdir cog;cd cog
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/cogData.py ../eggnog/SRR32313567.emapper.annotations /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/cog/cog.txt /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/cog/CogClass.txt cog.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/cog.R cog.xls cog1
  cd ../
  # GO注释及作图 GO annotation and plot
  mkdir go;cd go
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/go.py /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/go/go.txt ../eggnog/SRR32313567.emapper.annotations go.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/go_plot.R go.xls
  cd ../
  # kegg注释及作图 kegg annotation and plot
  mkdir kegg;cd kegg
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/kegg.py /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/kegg/K_info.txt /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/kegg/PathwayHtext.txt ../eggnog/SRR32313567.emapper.annotations kegg.xls pathway.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/mk_id2K2ko.py kegg.xls /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/kegg/B_A_KoPathways.txt > id2K_ko.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/ko.R  pathway.xls B_A_V
  cd ../
  # Uniprot注释 Uniprot annotation
  mkdir uniprot;cd uniprot
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/diamond blastp --evalue 1e-05 --max-target-seqs 1 --db /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/uniprot/uniprot_B_A_V.fa.dmnd  --query ../prokka_out/SRR32313567.faa --threads 40 --outfmt 6 'qseqid' 'sseqid' 'pident' 'qcovhsp' 'length' 'mismatch' 'gapopen' 'qstart' 'qend' 'sstart' 'send' 'evalue' 'bitscore' 'stitle' --out uniprot_blast.xls
  sed -i '1i#qseqid\tsseqid\tpident\tqcovhsp\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle' uniprot_blast.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/uniprot.py uniprot_blast.xls uniprot.xls
  cd ../
  # Nr注释统计作图 Nr annotation, statistical analysis, and visualization
  mkdir nr;cd nr
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/diamond blastp --evalue 1e-05 --max-target-seqs 1 --db /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/nr/B_A_V_nr.fa.dmnd --query ../../03.anno/prokka_out/SRR32313567.faa --threads 40 --outfmt 6 'qseqid' 'sseqid' 'pident' 'qcovhsp' 'length' 'mismatch' 'gapopen' 'qstart' 'qend' 'sstart' 'send' 'evalue' 'bitscore' 'stitle' --out nr_blast.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/species_num.py nr_blast.xls
  head -n 11 Nr_Species_distribution.xls > head_10_Nr_Species_distribution.txt
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/nr_result.py nr_blast.xls nr.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/nr.R head_10_Nr_Species_distribution.txt
  cd ../
  
  # tigerfam注释及作图 tigerfam annotation and visualization
  mkdir tigerfam;cd tigerfam
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/eggnog-mapper_2.1.9.sif hmmscan --cpu 40 -E 0.01 --tblout tblout --domtblout domtblout /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/tigerfam/TIGRFAMs.hmm ../prokka_out/SRR32313567.faa 
  awk 'NR<4||$0 !~/^#/ ' domtblout > tigerfam_tmp_hmm.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/tigerfam-hmmscan.py tigerfam_tmp_hmm.xls 1e-5 0.35 >tigerfam_stringent_hmm.xls
  cat tigerfam_stringent_hmm.xls |sed 's/#Seq ID/Seq ID/'|awk -F'\t' '{print $1"\t"$3"\t"$4"\t"$12}' > tigerfam.xls
  sed '1d' tigerfam_tmp_hmm.xls | sed 's/^# //g' |sed 2d |sed 's/target name/target_name/g' |sed 's/query name/query_name/g' |sed 's/description of target/description_of_target/g' |sed 's/[[:space:]]\{1,100\}/\t/g' > tigerfam_hmm.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/tigerfam.R tigerfam.xls
  cd ../
  
  # 注释统计 Annotation statistics
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/stat_anno.py -seq prokka_out/SRR32313567.faa -kegg kegg/kegg.xls -ko kegg/pathway.xls -uniprot uniprot/uniprot.xls -go go/go.xls -cog cog/cog.xls -tigerfam tigerfam/tigerfam.xls -nr nr/nr.xls -o ./
  
  # upset图 upset graph
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/upset_annotation.R all_annotation.xls upset_annotation
  
  # antismash注释 antismas annotation
  mkdir antismash;cd antismash
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/antismash_standalone_8.0.2.sif antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go  --smcog-trees  --output-dir /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/03.anno/antismash/results_genome  ../prokka_out/SRR32313567.gbk
  cd ../
  
  # 使用BioMGCore统计antismash结果 Use BioMGCore to analyze and summarize antismash results.
  mkdir BioMGCore;cd BioMGCore
  python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/BioMGCore/antiSTAT.py -i /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/03.anno/antismash/  -o ./statistics.xlsx
  cd ../
  
  # 基因岛注释 gene island annotation
  mkdir island;cd island
  singularity exec /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/island.sif python /opt/Island/Dimob.py -i ../prokka_out/SRR32313567.gbk -o ./ -d /opt/Island/islandpath/Dimob.pl   
  cd ../
  
  # 前噬菌体注释 prophage annotaiton
  mkdir prophage;cd prophage
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/phispy_3.7.8.sif PhiSpy.py  ../prokka_out/SRR32313567.gbk -o ./
  cd ../
  
  # CRISPR注释 CRISPR annotation 
  mkdir CRISPR;cd CRISPR
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/minced_v0.4.2.sif minced   ../prokka_out/SRR32313567.fna  genome_CRISPR.txt genome_CRISPR.gff
  cd ../../


## Step04.Advanced_functional_annotation
  
  # 创建目录 Create directory
  mkdir -p /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/04.specific;
  cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/04.specific
  
  # cazy注释及作图 cazy annotation and visualization
  mkdir cazy;
  cd cazy
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/eggnog-mapper_2.1.9.sif hmmscan --domtblout cazy_tmp_hmm.xls --cpu 40 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/cazy/dbCAN.txt ../../03.anno/prokka_out/SRR32313567.faa
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif   python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/hmmscan-parser.py cazy_tmp_hmm.xls cazy_stringent_hmm.xls  1e-18  0.35
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/get_cazy.py /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/cazy/domin_type.txt cazy_stringent_hmm.xls cazy.xls cazy_class_stat.xls
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/CAZY_plot.R  cazy_class_stat.xls
  cd ../
  
  # tcdb注释及作图 tcdb annotation and visualization
  mkdir tcdb;cd tcdb
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/diamond blastp --evalue 1e-05 --max-target-seqs 1 --db /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/tcdb/tcdb.fa.dmnd --query ../../03.anno/prokka_out/SRR32313567.faa --threads 40 --outfmt 6 'qseqid' 'sseqid' 'pident' 'qcovhsp' 'length' 'mismatch' 'gapopen' 'qstart' 'qend' 'sstart' 'send' 'evalue' 'bitscore' 'stitle' --out tcdb_blast.xls
  sed -i '1i#qseqid\tsseqid\tpident\tqcovhsp\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle' tcdb_blast.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/get_tcdb.py tcdb_blast.xls tcdb.xls
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/TCDB_plot.R  tcdb.xls
  cd ../
  
  # phi注释 phi annotation
  mkdir phi;cd phi
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/ncbi-blast-2.17.0+/bin/blastp -evalue 1e-05 -max_target_seqs 1 -db /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/phi/phi-base.fa -query ../../03.anno/prokka_out/SRR32313567.faa -num_threads 40 -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -out phi_blast.xls
  sed -i '1i#qseqid\tsseqid\tpident\tqcovs\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle' phi_blast.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif  python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/phi.py -i phi_blast.xls -o phi.xls
  cd ../
  
  # card注释 card annotation
  mkdir card;cd card
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif  python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/ch_format_fasta.py ../../03.anno/prokka_out/*.faa tmp.SRR32313567.faa
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/rgi.sif rgi load --card_json /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/card/card.json --local 
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/rgi.sif rgi main -i tmp.SRR32313567.faa -o SRR32313567.rgi -t protein -a BLAST -n 16 --clean --local --include_loose
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/replace_empty.py SRR32313567.rgi.txt card_rgi.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/get_card.py card_rgi.xls card.xls
  
  # 画图 Plot
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/card.py
  cd ../
  
  # cyp450注释 cyp450 annotation
  # 创建分析目录
  mkdir cyped;cd cyped
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/diamond blastp --evalue 1e-05 --max-target-seqs 1 --db /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/cyped/cyped.fa.dmnd --query ../../03.anno/prokka_out/SRR32313567.faa --threads 40 --outfmt 6 'qseqid' 'sseqid' 'pident' 'qcovhsp' 'length' 'mismatch' 'gapopen' 'qstart' 'qend' 'sstart' 'send' 'evalue' 'bitscore' 'stitle' --out cyped_blast.xls
  sed -i '1i#qseqid\tsseqid\tpident\tqcovhsp\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle' cyped_blast.xls 
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/get_cyped.py cyped_blast.xls cyped.xls
  cd -
  
  # 毒力因子注释 VFDB annotation
  # 创建分析目录
  mkdir VFDB;cd VFDB
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/diamond blastp --evalue 1e-05 --max-target-seqs 1 --db /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/VFDB/VFDB_setB_pro.fas --query  ../../03.anno/prokka_out/SRR32313567.faa --threads 40 --outfmt 6 'qseqid' 'sseqid' 'pident' 'qcovhsp' 'length' 'mismatch' 'gapopen' 'qstart' 'qend' 'sstart' 'send' 'evalue' 'bitscore' 'stitle' --out  VFDB_blast.xls
   singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif  python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/VFDB.py 
  cd -
  
  # 信号肽预测 Signal peptide prediction
  mkdir SignalP;cd SignalP
  singularity exec /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/sravankrishna47_signalp-fast.sif signalp6 --fastafile ../../03.anno/prokka_out/SRR32313567.fna  --organism other --output_dir  ./ --format txt --mode fast
  cd ../
  
  # 跨膜蛋白预测 Transmembrane protein prediction
  mkdir  tmhmm;cd  tmhmm
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/tmhmm-2.0c/bin/tmhmm --short < ../../03.anno/prokka_out/SRR32313567.fna > tmhmm.out  2>tmhmm.err
  cd ..
  
  # Ⅲ 型分泌系统效应蛋白注释 Annotation of type III secretion system effector proteins
  mkdir EffectiveT3;cd EffectiveT3
  singularity exec /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/EffectiveT3.sif  java -jar /opt/TTSS_GUI-1.0.1.jar -m /opt/TTSS_STD-2.0.2.jar -f ../../03.anno/prokka_out/SRR32313567.faa -q  -t selective -o result.txt
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/filter_effective.py result.txt filter_result.txt
  cd ../
  
  # 启动子预测  基因组、窗口大小、软件路径 Promoter prediction genome, window size, software path
  mkdir PromPredict;cd PromPredict
  python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/PromPredict.py ../../03.anno/prokka_out/SRR32313567.fna 100 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software
  cd -
  
  # ICE预测 Integrative and Conjugative Elements
  mkdir ICEberg;cd ICEberg
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/diamond blastp --evalue 1e-05 --max-target-seqs 1 --db /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/icefinder/ICE_aa_all.fas.dmnd --query   ../../03.anno/prokka_out/SRR32313567.faa --threads 40 --outfmt 6 'qseqid' 'sseqid' 'pident' 'qcovhsp' 'length' 'mismatch' 'gapopen' 'qstart' 'qend' 'sstart' 'send' 'evalue' 'bitscore' 'stitle' --out  ICEberg_blast.xls

   
## Step05.circlize_visualization

  # 创建目录
  mkdir -p /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/05.circlize;cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/05.circlize
    
  # 基因组圈图 genome circle diagram
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/first.py ../02.assembly/GCstat.xls tmp_GCstat.xls '1_pilon_pilon'
  less ../03.anno/prokka_out/SRR32313567.gff|awk 'NF>2'|sed -E 's/gnl\|Bacteria\|ctg_([0-9]+)/\1_pilon_pilon/g' > SRR32313567.gff
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/Train_bac_genome/Train-Bacteria-genome/script/circle_plot.r SRR32313567 ./ 2000 SRR32313567.gff ../02.assembly/aln3_2kbin.xls ../02.assembly/aln3_depth.xls ../02.assembly/aln2_2kbin.xls ../02.assembly/aln2_depth.xls tmp_GCstat.xls '1_pilon_pilon'
  
  # 如果仅有三代数据（无二代数据)，则使用该流程绘制圈图 If there is only third-generation data (no second-generation data), use this process to draw a circle chart
  # 纯三代圈图绘制  three-generation-data-Only circle diagram drawing
  less ../03.anno/prokka_out/SRR32313567.gff|awk 'NF>2'|sed -E 's/gnl\|Bacteria\|ctg_([0-9]+)/\1_pilon_pilon/g' > SRR32313567_long-read.gff
   singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/Train_bac_genome/Train-Bacteria-genome/script/circle_plot_only_long.r SRR32313567 ./ 2000 SRR32313567_long-read.gff  ../02.assembly/aln3_2kbin.xls ../02.assembly/aln3_depth.xls tmp_GCstat.xls '1_pilon_pilon'
  
  # 绘制简易圈图 Draw a simple circle diagram
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/cgview_2.0.3.sif perl /usr/bin/cgview_xml_builder.pl -sequence ../03.anno/prokka_out/SRR32313567.gbk -gc_content T -gc_skew T -size large-v2 -tick_density 0.05 -draw_divider_rings T -custom showBorder=false title=Example map titleFontSize=200 -output map.xml 
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/cgview_2.0.3.sif  java -jar /usr/bin/cgview.jar -i map.xml -o map.png 
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/cgview_2.0.3.sif  java -jar /usr/bin/cgview.jar -i map.xml -o map.pdf
  cd ../


## Step06.Comparative genome

  # 创建并进入目录
  mkdir -p /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/06.comp;
  cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/06.comp
  
  #准备文件：按照分析类型分类，存放在/data6/zhangtianyuan/Pipeline/EasyGenome/input文件夹下面的不同的子目录中，基因组的fna文件存放在genome子文件夹下，运行ortho需要的蛋白文件存放在ortho文件夹的pep文件夹下、运行jcvi需要的gff文件存放在jcvi子文件夹下 Prepare files: Categorize them by analysis type and store them in different subdirectories under the /data6/zhangtianyuan/Pipeline/EasyGenome/input folder. The genome fna files are stored in the genome subfolder. Protein files required for running ortho are stored in the pep folder within the ortho folder. Gff files required for running jcvi are stored in the jcvi subfolder.
    
  #orthofinder 同源基因聚类分析 Homologous gene clustering analysis
  #创建目录
  mkdir ortho;cd ortho
  #注释分析OG的基因组得到的faa文件复制到此分析目录，用prokka或者Bateka注释得到的faa文件或者网上数据库中的faa文件都可 The faa files obtained by annotating the OG genome are copied to this analysis directory. The faa files annotated by prokka or Bateka or the faa files in the online database can be used.
  
  ln -s ../../../input/ortho/pep .   # 下载的链接时根据拉丁名直接改名成4字母，例如Escherichia coli ==>  Ecol When downloading a link, change the name to 4 letters according to the Latin name, for example, Escherichia coli ==> Ecol
    
  mkdir pep_dir  #创建目录 Create directory
  
  ls pep/*.faa|awk -F "[./]" '{print "/data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/orthomclSoftware-v2.0.9/bin/orthomclAdjustFasta "$2" "$0" 1;mv "$2".fasta pep_dir"}' > AdjustFasta.sh  #生成改名脚本 Generate rename script
  sh +x AdjustFasta.sh    #运行改名脚本 Run the rename script
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/orthofinder_v3.1.sif orthofinder  -og -f $PWD/pep_dir -t 20 -a 10 -o $PWD/ortho   # 运行orthofinder Run orthofinder
  cd ../  
  #venn 同源基因组基因韦恩图绘制 Drawing of Venn diagram of homologous genome genes
  mkdir venn;cd venn
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/venn.R ../ortho/ortho/*/Orthogroups/Orthogroups.GeneCount.tsv Pchl,Pcic,Pent,SRR32313567
  
  # jcvi共线性分析 JCVI collinearity analysis
  mkdir jcvi;cd jcvi
    
  # 提取gff中CDS的信息 Extract CDS information from gff
  mkdir input ;cd input
  # 将核酸序列中的cds提取出来，如果近缘下载的是cds可以忽略近缘的提取，直接将近缘的cds序列cp  到当前分析目录即可  Extract the cds from the nucleic acid sequence. If the closest relative is the cds downloaded, you can ignore the extraction of the closest relative and directly cp the cds sequence of the closest relative to the current analysis directory.
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/jcvi_w_last_latest.sif python -m jcvi.formats.gff bed --type=CDS  --key=ID ../../../../input/jcvi/Pchl.gff  -o Pchl.bed
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/jcvi_w_last_latest.sif python -m jcvi.formats.bed uniq Pchl.bed
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/jcvi_w_last_latest.sif python -m jcvi.formats.gff bed --type=CDS --key=ID ../../../../input/jcvi/Pcic.gff  -o Pcic.bed
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/jcvi_w_last_latest.sif python -m jcvi.formats.bed uniq Pcic.bed
  
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/jcvi_w_last_latest.sif python -m jcvi.formats.gff bed --type=CDS  --key=ID  ../../../03.anno/prokka_out/SRR32313567.gff  -o SRR32313567.bed
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/jcvi_w_last_latest.sif python -m jcvi.formats.bed uniq SRR32313567.bed 
  
  # 提取cds序列并去除cds的描述信息，只保留id，不然画图会报错，Extract the cds sequence and remove the description information of cds, leaving only the id, otherwise the drawing will give an error
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/seqkit grep -f <(cut -f 4 Pchl.uniq.bed) ../../../../input/jcvi/Pchl.cds |/data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/seqkit seq  -i  > Pchl.cds
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/seqkit grep -f <(cut -f 4 Pcic.uniq.bed) ../../../../input/jcvi/Pcic.cds -i |/data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/seqkit seq -i >Pcic.cds
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/seqkit grep -f <(cut -f 4 SRR32313567.uniq.bed) ../../../03.anno/prokka_out/SRR32313567.ffn -i|/data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/seqkit seq -i  >SRR32313567.cds
  cd ../
  ln -s  input/Pchl.uniq.bed Pchl.bed
  ln -s input/ Pcic.uniq.bed Pcic.bed
  ln -s input/SRR32313567.uniq.bed SRR32313567.bed
  ln -s input/*cds .
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/jcvi_w_last_latest.sif python -m jcvi.compara.catalog ortholog --dbtype nucl  Pchl SRR32313567  --no_strip_names
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/jcvi_w_last_latest.sif python -m jcvi.compara.catalog ortholog --dbtype nucl SRR32313567 Pcic --no_strip_names 
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/jcvi_w_last_latest.sif python -m jcvi.compara.synteny screen --simple Pchl.SRR32313567.anchors Pchl.SRR32313567.anchors.new
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/jcvi_w_last_latest.sif python -m jcvi.compara.synteny screen --simple SRR32313567.Pcic.anchors SRR32313567.Pcic.anchors.new
  echo "gnl|Bacteria|ctg_1" > seqids ## Pcic
  echo "gnl|Bacteria|ctg_1" >> seqids # SRR32313567
  echo "gnl|Bacteria|ctg_1" >> seqids # Pchl
  echo "#y,xstart,xend,rotation,color,label,va,bed" > layout
  echo "0.7,0.2,0.85,0,green,Pchl,top,Pchl.bed" >> layout
  echo "0.5,0.2,0.85,0,blue,SRR32313567,top,SRR32313567.bed" >> layout
  echo "0.3,0.2,0.85,0,red,Pcic, bottom,Pcic.bed" >> layout
  echo "# edges" >> layout
  echo "e,0,1,Pchl.SRR32313567.anchors.simple" >> layout
  echo "e,1,2,SRR32313567.Pcic.anchors.simple" >> layout
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/jcvi_w_last_latest.sif python -m jcvi.graphics.karyotype seqids layout --format=pdf  --figsize=15x10
  cd ../
  
  # gtdbtk系统发育学分类与注释 gtdbtk phylogenetic classification and annotation
  mkdir gtdbtk;cd gtdbtk
  cp ../../03.anno/prokka_out/SRR32313567.fna  ../../../input/genome
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/gtdbtk_2.3.0_r214.sif gtdbtk classify_wf --genome_dir ../../../input/genome  --out_dir output/classify_wf --extension fasta --prefix bac --cpu 40 --skip_ani_screen
  cd ../
 
  # pyani平均核苷酸相似度分析 Pyani average nucleotide similarity analysis
  mkdir pyani;cd pyani
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/pyani.sif average_nucleotide_identity.py -i ../../../input/genome   -o output -m ANIm -g
  cd ..
  
  # syri 基因组共线性分析 Genome collinearity analysis
  mkdir syri;cd syri
  
  # 提取基因组 Extract genome
  awk 'BEGIN{seq=0} /^>/{seq++; if(seq>1) exit; print ">ctg1"; next} {print}'   ../../03.anno/prokka_out/SRR32313567.fna >target.fasta
  awk 'BEGIN{seq=0} /^>/{seq++; if(seq>1) exit; print ">ctg1"; next} {print}' /data6/zhangtianyuan/Pipeline/EasyGenome/ref_genome/GCA_000412675.fna >ref.fasta
  
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/minimap2_2.24.sif minimap2   -ax asm5 --eqx ref.fasta target.fasta |    singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools view -bS  >out.bam
  singularity exec -B /data6 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/syri_v1.71.sif  syri -c out.sam -r ref.fasta -q target.fasta -k -F B 
  # 暂时不要运行
  #echo -e "#file\tname\ttags">genomes.txt  
  #echo -e "ref.fasta\tctg1\tlw:1.5">>genomes.txt
  #echo -e "target.fasta\tctg1\tlw:1.5">>genomes.txt
  #/data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/plotsr --sr syri.out --genomes genomes.txt     -o output_plot.png
  cd ..


## Step07.Pan-genome analysis
  # 泛基因组分析 Pan-genome analysis
  mkdir -p /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/07.roary;cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/07.roary
  cp ../03.anno/prokka_out/*gff ../../input/roary/
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/roary_v3.9.1_addlibbz2.sif  bash -c 'export LC_ALL=C && export LANG=C && roary   -e --mafft -p 40  ../../input/roary/*gff'
  # 结果统计图绘制 Result statistical graph drawing
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/roary_pie.py
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/stat_roary.py
    cd - 

   # 备选方案：如果没有GFF文件，需要进行注释，参考下列文档
      # mkdir -p /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/07.roary;cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/07.roary
      #复制注释的gff文件到此目录
      #cp ../03.anno/prokka_out/*gff .
      
      #对近缘基因组做prokka注释
      # 设置路径
      #FASTA_DIR="/data6/zhangtianyuan/Pipeline/EasyGenome/ref_genome/"
      #OUT_BASE="prokka_out"
      ## 创建输出目录
      #mkdir -p prokka_out
      #
      ## 循环处理所有fasta文件
      #for fasta in ${FASTA_DIR}/*.fasta; do
      #    # 提取样本名（不带路径与后缀）
      #    sample=$(basename "$fasta" .fasta)
      ##    sample=$(echo "$fasta" | cut -d '.' -f1)
      #
      #    echo "Processing $sample..."
      #
      #    singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/prokka_1.14.5.sif   prokka \
      #        --force \
      #        --outdir prokka_out \
      #        --kingdom Bacteria \
      #        --gcode 11 \
      #        --cpus 40 \
      #        --evalue 1e-6 \
      #        --rnammer \
      #        --rfam \
      #        --genus Genus \
      #        --species species \
      #        --strain strain \
      #        --compliant \
      #        --locustag ctg \
      #        --centre Bacteria \
      #        --prefix "$sample" \
      #        "$fasta"
      #done
    
      #把近缘基因组的prokka注释得到的gff文件都复制到/data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/03.anno/roary目录
      #singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/roary_v3.9.1.sif roary  --mafft -p 8  *.gff 

## Step08.MLST analysis
  # mlst分析 
  
  # 创建目录 Create a directory
  mkdir -p /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/08.MLST/mlst;cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/08.MLST/mlst
  # 运行MLST分析 Run MLST analysis
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/mlst_v2.9addlibz.sif  mlst /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/03.anno/prokka_out/SRR32313567.fna > tmp.txt 
  # 汇总MLST分析结果 Summarizing MLST analysis results
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/mlst2tsv.py tmp.txt MLST.xls
  # 返回主目录 Return to the main directory
  cd -
  
  # wgmlst分析
  
  # 如果有多个近缘基因组，下载近缘的cds，存放在input/cds目录  If there are multiple closely related genomes, download the CDS files of closely related genomes and store them in the input/cds directory.
  # 创建目录 Create a directory
  mkdir -p /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/08.MLST/wgmlst;cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/08.MLST/wgmlst
  # 拷贝测试数据或下载过的cds数据   Copy test data or downloaded CDS data
  cp -rf /data6/zhangtianyuan/Pipeline/EasyGenome/input/cds . 
  
  # 将目标基因组的cds拷贝到cds目录     Copy the target genome CDS to the CDS directory
  cp ../../03.anno/prokka_out/SRR32313567.ffn ./cds 
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/chewbbaca1_v3.4.sif chewBBACA.py CreateSchema -i  ./cds   -o schema_folder --cpu 16
  # 否则下载数据PubMLST 数据库现有的schema Otherwise download the existing schema of the data PubMLST database
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/chewbbaca1_v3.4.sif  chewBBACA.py AlleleCall -i ./cds -g schema_folder/schema_seed -o wgMLSTpath --cpu 20 --cds
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/chewbbaca1_v3.4.sif chewBBACA.py ExtractCgMLST -i wgMLSTpath/results_alleles.tsv -o cgMLSTPath  --threshold 0.95
  cd -




################################################  Thanks  #################################################################
#########################感谢使用，如有疑问，请联系邮箱:zhangtianyuan@caas.cn或微信号:allianzzhang############################
####Thank you for using. If you have any questions, please contact email: zhangtianyuan@caas.cn Or WeChat: allianzzhang###
#####################################################END##################################################################

