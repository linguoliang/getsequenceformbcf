# 从bcf文件获取序列
## 1. 软件说明
1. 本软件用于从bcf文件中提取序列，并提取出相应的位置信息。
2. 软件分为数据输入，获取序列，区间信息提取，数据输出四大模块。
3. 软件接受的参数 *-i/--input file.bcf -g/--gff3 file.gff3
   -o/--output-prefix outputprefix*
4. 

## 2. 数据输入模块
1. 输入文件bcf/vcf 前部分有以#开头的陈述信息用自定义的TrimHead() 函数去除
头部信息。
2. -i 参数可以接受管道操作符传递的数据。  
3.  


## 3. 数据输出模块
1. fasta 文件格式：  
> \>chromsome:POSstart-POSend|fulllength or intron
ACGTTACCGAGCGTGTGAGCTAGCAGTA  

   * 文件支持简并模式  
2. 统计文件