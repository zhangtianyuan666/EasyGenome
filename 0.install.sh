####
 # Author作者: Tianyuan Zhang(张天缘), Yong-xin Liu(刘永鑫), Yilin Li(李伊琳) Zhihao Zhu(朱志豪), Qingrun Xue(薛清润) et al.
    # Update更新时间: 2025-10-27
    # Version版本: 1.0

#####################01.数据下载与解压  #####################
###############Data download and decompression##############
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


#####################02.无root权限安装singularity#####################  
##############Installing singularity without root privileges#########
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
 
 
#####################03.使用Root安装singularity##################
####################Install singularity using Root##############
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



		
