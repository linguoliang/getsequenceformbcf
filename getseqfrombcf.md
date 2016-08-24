# 从bcf文件获取序列

## 1. 软件说明

1. 本软件用于从bcf文件中提取序列，并提取出相应的位置信息。
2. 软件分为数据输入，获取序列，区间信息提取，数据输出四大模块。
3. 软件接受的参数 *-i/--input file.bcf -g/--gff3 file.gff3
4. -o/--output-prefix outputprefix*
5. 4.

## 2. 数据输入模块

1. 输入文件bcf/vcf 前部分有以#开头的陈述信息用自定义的TrimHead()
2. 函数去除头部信息。
3. -i 参数可以接受管道操作符传递的数据。  

## 3.数据处理模块

### 1.碱基处理模块(base_merge.py)

 碱基处理模块是用于处理碱基的简并性的一个模块。碱基的简并性

### 2.Gtf文件处理模块

## 4. 数据输出模块

- fasta 文件格式：  

  > >chromsome:POSstart-POSend|gene-0/1

  > ACGTTACCGAGCGTGTGAGCTAGCAGTA  

  - 文件支持简并模式  
  - gene 可以是多个，出现多个的时候为gene1_gene2...
  - 0/1/2/3 代表是否为全长1表示是全长，0/1也可出现多个

1. 统计文件
