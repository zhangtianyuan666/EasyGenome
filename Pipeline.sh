# Easygenome易基因组

  # Author作者: Tianyuan Zhang(张天缘), Yong-xin Liu(刘永鑫), Yilin Li(李伊琳) Zhihao Zhu(朱志豪), Qingrun Xue(薛清润) et al.
  # Update更新时间: 2025-10-28
  # Version版本: 1.0 

## Step01.data_preparation_and_quanlity_control
  # 创建分析目录  建议创建在Easygenome目录下，以样本名命名 Create an analysis directory. It is recommended to create it in the Easygenome directory and name it after the sample.
  mkdir -p /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata  
  # 提取二代5000条reads Extract 5000 second-generation reads
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/sub_fq.py  /data6/zhangtianyuan/Pipeline/EasyGenome/data/SRR32313567_R1.fastq.gz 5000 /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/Read.5000.fq.gz
  # kraken2二代数据物种分类 - 判断是否有污染 Use Kraken2 for second-generation data species classification - determine whether there is pollution
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/kraken2_2.1.2-no-db.sif kraken2 --db /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/Kraken2/  --threads 48 --output /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567_kraken2.out.txt --report /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567_kraken2.report.txt --confidence 0.1 /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/Read.5000.fq.gz
  # krona可视化kraken2分类结果 krona visualization of kraken2 classification results
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/kreport2krona.py -r /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567_kraken2.report.txt -o /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/krona
  cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata;singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/krona_2.7.1--e7615f7.sif ktImportText krona;cd ../../

  # 提取三代数据，然后进行质控，使用kraken2物种分类 Extract three generations of data, perform quality control, and use kraken2 species classification
  # fastp过滤二代数据及结果统计 Fastp filtering of second-generation data and result statistics
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/fastp -i /data6/zhangtianyuan/Pipeline/EasyGenome/data/SRR32313567_R1.fastq.gz -I /data6/zhangtianyuan/Pipeline/EasyGenome/data/SRR32313567_R2.fastq.gz -w 16 -o /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567_R1.clean.fq.gz -O /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567_R2.clean.fq.gz -j /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.json -h /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.html
  # 进入分析目录后进行统计 Perform statistics after entering the analysis directory
  cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/fastp_stats.py
  cd ../../
  # 三代数据去接头 Three-generation data connector
  singularity exec -B /data6/  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/porechop_v0.2.4.sif porechop -i /data6/zhangtianyuan/Pipeline/EasyGenome/data/SRR32313567.nanopore.fastq.gz  -o /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.Adapter.fastq -t 30
  # 三代数据质控，保留1600bp及以上的reads、结果统计 Third-generation data quality control, retaining reads of 1600bp and above, and result statistics
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/fastq_to_fastq.py --min_length 1600 /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.Adapter.fastq > /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered.min1600.fastq

  # 选取1G数据用于分析 Select 1G data for analysis
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/select_1G_ont.py /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered.min1600.fastq   /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered.fastq 
  # 统计结果 Statistical results
  /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Software/seqkit stats /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered.fastq -a |awk -F' *' 'BEGIN{OFS="\t"}NR==1{print $1,$4,$5,$6,$7,$8,$13}NR==2{n=split($1,a,"[/]" );print a[n],$4,$5,$6,$7,$8,$13}'> /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/tmp_ONT_stats.xls
  # NanoPlot可视化 NanoPlotVisualization
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/nanoplot_latest.sif  NanoPlot --fastq /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered.fastq -o /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/fastq-plots --plots hex dot -t 48
  # kraken2物种分类 三代数据 Kraken2 species classification three generations of data
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/sub_fq.py /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered.fastq 5000 /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered_5000.fastq

  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/kraken2_2.1.2-no-db.sif kraken2 --db /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Database/Kraken2/  --threads 48 --output /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered_5000_kraken2.out.txt --report /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered_5000_kraken2.report.txt  --confidence 0.1 /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered_5000.fastq
  
  # krona可视化kraken2分类结果 krona visualization of kraken2 classification results
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/kreport2krona.py -r /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567.filtered_5000_kraken2.report.txt -o /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata/SRR32313567_filtered_5000_krona
  cd /data6/zhangtianyuan/Pipeline/EasyGenome/SRR32313567/01.cleandata;singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/krona_2.7.1--e7615f7.sif ktImportText SRR32313567_filtered_5000_krona;cd ../../

# Step02.assembly_and_evaluation

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
  
  # 根据质粒注释结果修改组装文件序列名 Modify the assembly file sequence name based on plasmid annotation results
  #singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/rename_fasta_by_plasflow.py ../assembly.fasta SRR32313567.txt rename_assembly.fasta
  
  
  # 确定组装版本  Determine assembly version
  mkdir -p assemble
  mv assembly.fasta assemble/assembly.fasta
  
  # Quast Quality Value
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/quast_5.2.0.sif quast.py -o quast -t 40 assemble/assembly.fasta
  # checkm  
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/checkm.v1.1.3.sif checkm lineage_wf unicycler_out/ checkmout/ -x fasta -t 48  --pplacer_threads 8 --tab_table -f checkmout/checkm.txt
  # busco Detecting genome integrity
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/staphb_busco_6.0.0-prok-bacteria_odb12_2024-11-14.sif  busco -o busco -i assemble/assembly.fasta -l /busco_downloads/lineages/bacteria_odb12  -m geno -c 40
  
  # qv base accuracy
  mkdir merqury_qv ;cd merqury_qv
  singularity exec -B /data6/ --env MERQURY=/usr/local/share/merqury/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/merqury_v1.3.sif meryl count k=21 output read1.meryl ../../01.cleandata/SRR32313567_R1.clean.fq.gz
  singularity exec -B /data6/ --env MERQURY=/usr/local/share/merqury/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/merqury_v1.3.sif meryl count k=21 output read2.meryl ../../01.cleandata/SRR32313567_R2.clean.fq.gz
  singularity exec -B /data6/ --env MERQURY=/usr/local/share/merqury/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/merqury_v1.3.sif meryl count k=21 output read3.meryl ../../01.cleandata/SRR32313567.filtered.fastq
  singularity exec -B /data6/ --env MERQURY=/usr/local/share/merqury/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/merqury_v1.3.sif meryl  union-sum output illumina.meryl read1.meryl read2.meryl 
  singularity exec -B /data6/ --env MERQURY=/usr/local/share/merqury/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/merqury_v1.3.sif meryl  union-sum output all.meryl illumina.meryl read3.meryl
  singularity exec -B /data6/ --env MERQURY=/usr/local/share/merqury/   /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/merqury_v1.3.sif    bash /usr/local/share/merqury/merqury.sh all.meryl ../assemble/assembly.fasta  assembly
  cd ..
  # 深度统计 Depth statistics
  # 二代测序深度统计 Next-generation sequencing depth statistics
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools faidx assemble/assembly.fasta
  awk -F'\t' '{print $1"\t"$2}' assemble/assembly.fasta.fai  > ctg_length.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bwa-samtools_0.7.12_1.2.1.sif bwa index assemble/assembly.fasta 
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bwa-samtools_0.7.12_1.2.1.sif bwa mem -t 40 assemble/assembly.fasta ../01.cleandata/SRR32313567_R1.clean.fq.gz ../01.cleandata/SRR32313567_R2.clean.fq.gz > aln2.sam
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools sort aln2.sam > aln2sorted.sam
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bedtools_v2.30.0.sif bedtools makewindows -g ctg_length.xls  -w 2000 -i winnum > genome_windows_2k.bed
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bedtools_v2.30.0.sif bedtools coverage -mean -sorted -bed -g ctg_length.xls -a genome_windows_2k.bed -b aln2sorted.sam > aln2_2kbin.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools depth -a aln2sorted.sam > aln2_depth.xls
  # 三代测序深度统计 Depth statistics of third-generation sequencing
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/minimap2_2.24.sif minimap2 -ax map-ont -t 40 assemble/assembly.fasta ../01.cleandata/SRR32313567.filtered.fastq -o aln3.sam
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools sort aln3.sam > aln3sorted.sam
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/bedtools_v2.30.0.sif bedtools coverage -mean -sorted -bed -g ctg_length.xls -a genome_windows_2k.bed -b aln3sorted.sam > aln3_2kbin.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/samtools_1.17.sif samtools depth -a aln3sorted.sam > aln3_depth.xls
  
  # GC含量及GC-skew GC content and GC-skew
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/python39pandas_pexpect_Bio_PromPredict_r.sif python3 /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/GC.py -f assemble/assembly.fasta -w 2000 -s 2000 > GCstat.xls
  
  
  # GC_depth 图 GC_depth graph
  paste GCstat.xls aln2_2kbin.xls |awk '{print $1":"$2"-"$3"\t"$4"\t"$10}' >aln2_GC_depth.xls
  paste GCstat.xls aln3_2kbin.xls |awk '{print $1":"$2"-"$3"\t"$4"\t"$10}' >aln3_GC_depth.xls
  singularity exec -B /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/GC_depth.r aln2_GC_depth.xls aln2_GC_depth
  singularity exec -B  /data6/ /data6/zhangtianyuan/Pipeline/EasyGenome/Public/Singularity/plot1.sif Rscript /data6/zhangtianyuan/Pipeline/EasyGenome/Public/script/GC_depth.r aln3_GC_depth.xls aln3_GC_depth
