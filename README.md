# EasyGenome (易基因组）

a user-friendly and flexible pipeline for prokaryote genome analysis （一种用户友好且灵活的原核生物基因组分析管道）

Version：v1.00

Update：2025/10/27

# Pipeline manual and file description (流程使用和文件介绍)


# What can we do? (结果展示)
- EasyGenome is User-friendly pipeline for bacterial genome analysis, integrating eight essential steps.

-  Supports various sequencing platforms data (nanopore-only, PacBio-only, Hybrid long -read and short-read; short-read Only).

- Supports species classification, genome assembly, genome circular graph, comparative genome analysis, pan-genome visualizations, and MLST profiles.

- Fully open-source and customizable, enabling reproducible and high-quality research.

- Support 20+ analysis methods and publish-ready visualization;

- Chinese/English manual and video supported.

- The pipeline is maintained and updated regularly, and we encourage users to contribute appropriate code.

![annotation_1](https://github.com/user-attachments/assets/e446ca2a-0821-4c08-83af-3ca125a1c969)
**Figure 1. Summary of EasyGenome.**

# Install (安装)
#####################01.数据下载与解压  #####################
###############Data download and decompression##############
```bash
#下载，解压，删除压缩包 Download, decompress, and delete the compressed package
wget -c  -O /data6/zhangtianyuan/Pipeline/EasyGenome/data.tar.gz ftp://download.nmdc.cn/Easygenome/data.tar.gz   #下载测试数据到指定目录，注意修改/data6/zhangtianyuan/Pipeline/EasyGenome/为本服务器上的目录 #Download test data to the specified directory, please make sure to modify/data6/zhangtianyuan/Pipeline/EasyGenome/to the directory on this server
tar -zxvf data.tar.gz #解压目录 Unzip directory
rm data.tar.gz #删除压缩包 Delete compressed package

wget -c -O /data6/zhangtianyuan/Pipeline/EasyGenome/ ftp://download.nmdc.cn/Easygenome/input.tar.gz #下载准备数据、示例数据 Download preparation data and sample data 
tar -zxvf input.tar.gz #解压目录 Unzip directory
rm input.tar.gz #删除压缩包 Delete compressed package

wget -c -O /data6/zhangtianyuan/Pipeline/EasyGenome/ ftp://download.nmdc.cn/Easygenome/Public.tar.gz #下载软件、数据库 Download software and database 
tar -zxvf Public.tar.gz #解压目录 Unzip directory
rm Public.tar.gz  #删除压缩包 Delete compressed package

# 以上3个文件，内测时，若无法下载，可从/data6/zhangtianyuan/Pipeline/EasyGenome拷贝
```

#####################02.无root权限安装singularity#####################  
##############Installing singularity without root privileges#########
```bash
#下载GO并安装  Download GO and install it
#进到官网下载最新Linux版本安装，官网：https://go.dev/dl/在服务器上用wget工具下载或者手动进官网下载再传到服务器上都可
wget  https://go.dev/dl/go1.25.3.linux-amd64.tar.gz #下载软件压缩包Download the software compressed package
tar  -xzf   go1.25.3.linux-amd64.tar.gz  #解压 Unzip
rm go1.25.3.linux-amd64.tar.gz #删除压缩包 Delete compressed package
 
#将GO添加到环境变量中 Add to environment variables：
export PATH=/opt/soft/go/bin:$PATH  #在~/.bashrc文件中添加一行（/opt/soft/go/替换成你的go安装目录）  Add a line to the ~/.bashrc file (replace /opt/soft/go/ with your go installation directory)
source ~/.bashrc   # 添加完后运行这行，使重新登录用户会话或者运行也生效  After adding, run this line to make it effective to log in to the user session or run again

#下载singularity源代码包 Download the singularity source code package
#在官网https://github.com/sylabs/singularity/releases找想要安装的版手动下载或wget工具下载 

wget   -c https://github.com/sylabs/singularity/releases/download/v4.3.4/singularity-ce-4.3.4.tar.gz #下载 Download
tar  -xzf   singularity-ce-4.3.4.tar.gz  #解压 Unzip
cd  singularity-ce-4.3.4 #进到软件目录 Go to the software directory
./mconfig --prefix=./   #编译， --prefix=后面跟想要安装singularity软件的目录 Compile, --prefix= followed by the directory where you want to install the singularity software


#如果编译的时候和下面一样显示没有找到libseccomp包,可以直接禁用这个功能(不推荐，会导致对应的功能不可用) If the libseccomp package is not found during compilation as shown below, you can disable this feature directly(Not recommended, as it will cause the corresponding function to be unavailable)
#运行 You can disable this feature directly by running：
./mconfig --prefix=./   --without-seccomp
#或者让管理员安装 Or ask the administrator to install it：
apt-get install libseccomp-dev，
#安装完后再运行：
./mconfig --prefix=./

#mconfig运行完成后按照提示依次运行 After mconfig is finished, follow the prompts to run
cd /opt/soft/singularity-ce-4.3.4/builddir #进builddir目录  Enter the builddir directory 
make #编译 compile
make install  #安装 Install
```
 
#####################03.使用Root安装singularity##################
####################Install singularity using Root##############
```bash
#注意：如果上述02步骤无法安装，则使用本步骤
#NOTICE:If the above step 02 cannot be installed, use this step

#GO安装同上 GO installation is the same as above

#安装singularity参考官网：https://docs.sylabs.io/guides/latest/user-guide/quick_start.html#download
#安装依赖的基础软件 Install dependent software
sudo apt-get update  # Install debian packages for dependencies
sudo apt-get install -y \ 
		autoconf \
		automake \
	        cryptsetup \
		fuse2fs \
		git \
		fuse \
		libfuse-dev \
		libseccomp-dev \
		libtool \
		pkg-config \
		runc \
		squashfs-tools \
		squashfs-tools-ng \
		uidmap \
		wget \
		zlib1g-dev  #安装依赖 Install dependencies

#下载singularity源代码包  Download the singularity source code package

#在官网https://github.com/sylabs/singularity/releases找想要安装的版手动下载或wget工具下载 
#on the official website https://github.com/sylabs/singularity/releases Find the version you want to install and download it manually or with the wget tool

wget   -c https://github.com/sylabs/singularity/releases/download/v4.3.4/singularity-ce-4.3.4.tar.gz #download the package
tar  -xzf   singularity-ce-4.3.4.tar.gz  #解压缩目录	unzip the package
cd  singularity-ce-4.3.4 #进到软件目录  Enter the software directory
	./mconfig && make -C builddir   #生成配置文件	Generate configuration file
	    	sudo make -C builddir install	#安装软件	install the package

```

# Quick Start (快速运行)
Setup the work directory(wd), and  directory(db), then run each line by click run in top right corner

# Example dataset (示例数据)
Public目录包含了数据库，软件，镜像和脚本 
The Public directory contains databases, software, container images, and scripts.

input目录存放比较基因组和泛基因组分析的示例数据  
The input directory stores example data for comparative genomics and pan-genome analysis.

data目录存放测试数据，数据来源于Pseudomonas sp. KK18 (PRJNA1079050);其中, SRR32313567.nanopore.fastq.gz为纳米孔长读长测序数据;SRR32313567_R1.fastq.gz和SRR32313567_R2.fastq.gz为短读长测序数据    
The data directory contains test data sourced from Pseudomonas sp. KK18 (PRJNA1079050).  SRR32313567.nanopore.fastq.gz is Nanopore long-read data. SRR32313567_R1.fastq.gz and SRR32313567_R2.fastq.gz are short-read data.

# FAQ (常见问题)
Frequenty Asked Questions in pipeline.sh

Note: All the .sh script is writting in markdown format, using Notepad++ or VSCode for better reading experience.

# Citation (引文)
If used this script, please cited:

Copyright 2019-2025 Tianyuan Zhang zhangtianyuan@caas.cn, Yong-Xin Liu liuyongxin@caas.cn
