####
 # Author作者: Tianyuan Zhang(张天缘), Yong-xin Liu(刘永鑫), Yilin Li(李伊琳) Zhihao Zhu(朱志豪), Qingrun Xue(薛清润) et al.
    # Update更新时间: 2025-10-28
    # Version版本: 1.0 

#####################################    此步骤根据服务器的存储情况可批量更改挂载目录和分析路径  #############################################
###############This step can batch change the mounting directory and analysis path based on the storage situation of the server##########


#!/bin/bash

# 批量路径替换脚本
# 用法：./00.prepare_batch_replace.sh

# 定义要替换的路径映射（旧路径 -> 新路径）
declare -A path_mapping=(
    ["/data6/zhangtianyuan/Pipeline/EasyGenome/"]="/home/EasyGenome/"
    ["-B /data6/"]="-B /mnt/"
    # 添加更多需要替换的路径...默认Database、script和Singularity目录在Public下面，如果另存目录请再这里按上面格式添加替换啊，注意路径最后的斜杠
)
# 要处理的脚本文件，支持通配符
script_files=("Pipeline.sh")

# 备份原文件并执行替换
for file in ${script_files[@]}; do
    if [[ -f "$file" ]]; then
        echo "处理文件: $file"
        cp "$file" "$file.bak"
        for old_path in "${!path_mapping[@]}"; do
            new_path="${path_mapping[$old_path]}"
            # 使用sed进行替换，注意转义路径中的斜杠
            escaped_old=$(printf '%s\n' "$old_path" | sed 's/[\/&]/\\&/g')
            escaped_new=$(printf '%s\n' "$new_path" | sed 's/[\/&]/\\&/g')
            sed -i "s|${escaped_old}|${escaped_new}|g" "$file"
        done
        
        echo " $file (备份保存为 $file.bak)"
    fi
done

