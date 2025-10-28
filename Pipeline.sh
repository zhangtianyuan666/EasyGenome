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
