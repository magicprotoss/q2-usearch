---
editor_options: 
  markdown: 
    wrap: 72
---

# q2-usearch

这是一个把 [USEARCH](https://github.com/rcedgar/usearch12/) 装进
[QIIME2](https://qiime2.org) 的插件

## 关于本插件

[USEARCH](https://drive5.com/usearch/)
作为扩增子数据上游分析的主流软件（[QIIME](https://qiime2.org)、[USEARCH](https://drive5.com/usearch/)、[MOTHUR](https://mothur.org/)）之一,
在过去十多年间广受欢迎 (也是我们的首选🙈).
我们最初写这个插件是为了组内使用方便, 现在 [USEARCH
开源](https://mp.weixin.qq.com/s/i0zzOP5IRNdY9PfqHbpDEQ)之后我们决定共享这个插件.
以下是本插件打包到 [QIIME2](https://qiime2.org/) 中的流程/命令:

-   使用unoise3算法对有效数据去噪生成zOTU表及zOTUs
-   使用uparse算法对有效数据聚类生成OTU表及97%OTUs
-   使用unoise3算法对有效数据去噪生成zOTUs，随后使用uclust算法将zOTUs聚类到OTUs，最后生成OTU表
-   合并双端序列

\*插件还在早起开发过程中，后续更新接口可能会有调整

## 插件安装

第一步: 激活我们想要安装到的 QIIME2 conda 环境

``` bash
# conda activate qiime2-amplicon-2024.5
conda activate <目标q2环境名称>
```

如果找不到环境名的话可以这样看一下

``` bash
conda env list | grep qiime2
```

这样当前conda下所有环境名中带有qiime2字符的虚拟环境都会被检索出来

``` bash
# qiime2-amplicon-2024.5     /home/navi/miniconda3/envs/qiime2-amplicon-2024.5
```

第二步: 下载插件

``` bash
git clone https://github.com/magicprotoss/q2-usearch .
```

第三步: 切换路径到插件文件夹下并安装插件

``` bash
cd q2-usearch # 切换路径
python ./setup.py install # 安装插件
qiime dev refresh-cache # 刷新QIIME2命令行界面缓存
```

第四步:
使用mamba或者conda安装[seqkit2](https://mp.weixin.qq.com/s/ucCA4-eZINkppdyQFItXHw)和[usearch12](https://mp.weixin.qq.com/s/i0zzOP5IRNdY9PfqHbpDEQ)

``` bash
# 建议使用mamba
# 国内网络不稳conda圈圈转到一半容易报断连错误
mamba install -c bioconda seqkit">=2.0.0" usearch
# conda install -c bioconda seqkit">=2.0.0" usearch
```

测试一下usearch是否可以被正常调用

``` bash
usearch
# usearch v12.0 [b1d935b], 132Gb RAM, 24 cores
# (C) Copyright 2013-24 Robert C. Edgar.
# https://drive5.com/usearch
```

第五步（可选）: 清理残留文件

``` bash
cd .. && rm -rf q2-usearch
```

## 方法

### 1. denoise-no-primer-pooled

-   封装流程概览

    ![](docs/denoise_no_primer_pooled.png)

-   [API文档](docs/denoise_no_primer_pooled.md)

### 2. cluster-no-primer-pooled

-   封装流程概览

    ![](docs/cluster_no_primer_pooled.png)

-   [API文档](docs/cluster_no_primer_pooled.md)

### 3. denoise-then-cluster-no-primer-pooled

-   封装流程概览

    ![](docs/denoise_then_cluster_no_primer_pooled.png)

-   [API文档](docs/denoise_then_cluster_no_primer_pooled.md)

## 教程：使用插件去噪并生成100% zOTU

### 处理测序公司返回的有效数据

1.  准备序列文件

    现在在测序公司返回的原始数据一般会直接打包成一个压缩包，我们来看一下里面的内容:

    ![](docs/folder_structure_zh_cn.png)

    一般情况下，使用 "Valid data" 文件夹里的 "有效数据" 会比较省时省力.
    接下来我们把序列文件文件名中不是样本ID的部分删掉.
    （公司那个压缩包别删, 万一重命名的时候手滑了...）
    如果您用的是windows的电脑, 微软出的
    [PowerToys](https://learn.microsoft.com/en-us/windows/powertoys/)
    中的那个 PowerRename 批量重命名文件超级好用:

    ![](docs/power_rename.png)

2.  将序列文件导入到QIIME2中

    首先激活 QIIME2 的 conda 环境

    ``` bash
    conda activate qiime2-amplicon-2024.5
    ```

    我们在上一步已经把文件名都换成样本ID了,
    现在只需要做一个manifest文件（和QIIME里面那个mapping文件一个作用,
    用来对应样本ID和文件名）就能将序列导入到QIIME2中了.
    我们可以使用[纯bash脚本的方案](https://library.qiime2.org/plugins/qiime2-manifest-metadata-generator/23/),
    或者也可以用我们写的python脚本（还需要装一个python包,
    同样使用mamba安装）：

    ``` bash
    # conda install xlsxwriter
    # 建议使用mamba
    # 国内网络不稳conda圈圈转到一半容易报断连错误
    mamba install xlsxwriter
    ```

    装好后（只需要装一次），运行脚本生成manifest文件和metadata表格

    ``` bash
    python generate_metadata.py --input_path <path-to-your-valid-data> --from_filename
    ```

    脚本输出的 'manifest.tsv'
    中记录了所有序列文件的样本ID和绝对路径，这个待会导入序列文件的时候用.

    ![](docs/manifest_example.png)

    1.  脚本输出的 'metadata.xlsx' 是已经按照 QIIME2
        元数据文件格式格式化好的表格，样本ID一列已经自动填充好了.

    ![](docs/metadata_example.png)

    我们建议大家将这个文件下载下来填好,
    一是因为QIIME2自带可交互的可视化模块,
    方便按分组初步探索数据；二是这个去除第二行注释其实就是分组文件，后续R中调包分析也少不了.
    关于QIIME2元数据的具体填法可以看[这里](https://docs.qiime2.org/2024.5/tutorials/metadata/#metadata-formatting-requirements)
    .
    接下来我们可以运行以下命令直接将该Excel表格转换为QIIME2兼容的tsv格式，以便后续使用.

    ``` bash
    python generate_metadata.py --to_tsv
    ```

    上述准备工作完成之后，运行下方命令将序列文件导入到QIIME2中.
    注意我们使用的是经过公司拼接预处理之后的 '有效数据',
    导入命令中的输入文件格式和类型都应该使用单端序列所对应的参数,
    无论我们使用的是单端还是双端测序策略.

    ``` bash
    qiime tools import \
        --input-path manifest.tsv \
        --input-format 'SingleEndFastqManifestPhred33V2' \
        --type 'SampleData[SequencesWithQuality]' \
        --output-path fastq-seqs.qza
    ```

3.  去噪生成zOTU表与zOTUs

    以上准备工作完成后，我们就可以运行插件中的denoise-no-primer-pooled流程进行去噪生成zOTU表与zOTUs了.
    在参数设置方面，这里选择将 `--p-min-size` 设置成 4 来提高灵敏度.

    ``` bash
    qiime usearch denoise-no-primer-pooled \
        --i-demultiplexed-sequences fastq-seqs.qza \
        --p-min-size 4 \
        --o-representative-sequences rep-seqs-unoise3.qza \
        --o-table table-unoise3.qza \
        --o-denoising-stats stats-unoise3.qza
    ```

    如果您这边的下机数据还没有释放, 我们准备了[(Dong, Guo et al.
    2021)](https://doi.org/10.1016/j.aquaculture.2020.736199)
    发表在Aquaculture杂志上文章的数据用于测试,
    下载链接[点我](https://zenodo.org/records/12655543?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjQyMzZlOTlkLWFkZmEtNDM3NS1hMzUzLTFlNWMxZjJlM2Q1YSIsImRhdGEiOnt9LCJyYW5kb20iOiI2NzJmNDlkYTRmNTUwZDdlYzA4ZjllMmUyYmU2MTYzNCJ9.U9S4WNLhCgTYMz37j2QUatQ2F-HYzl9suSWlFr11nMdKslJNTraD5IuxvDmBgsAzUc9s8wTsm_9UJhTs6VZ9wA)（百度云[点我](https://pan.baidu.com/s/1C_TwzBLh85k16JxEGEkDjQ?pwd=r6y8)）.

    ``` bash
    qiime usearch denoise-no-primer-pooled \
        --i-demultiplexed-sequences ddbj_dl.qza \
        --p-min-size 4 \
        --o-representative-sequences rep-seqs-unoise3.qza \
        --o-table table-unoise3.qza \
        --o-denoising-stats stats-unoise3.qza
    ```

### General Use Cases

1.  单端测序数据处理

    这里使用 QIIME2官方教程
    ["人体各部位微生物组"](https://mp.weixin.qq.com/s/p2Snx0v8Fh_BOY-z2vVUCg)
    , 请耐心观看教学视频到 '**a文件准备和样本拆分**' 部分.
    在6分:50秒（运行DADA2之前），除了选项 **1: DADA2** 和选项 **2:
    Deblur** 之外, 我们现在又多了第三种选择.

    #### **选项 3: UNOISE3**

    [unoise3](https://drive5.com/usearch/manual/cmd_unoise3.html)
    命令首先使用
    [UNOISE](https://drive5.com/usearch/manual/unoise_algo.html)算法对序列文件进行去噪得到扩增子.
    随后使用经过进一步优化的[uchime2
    algorithm](https://drive5.com/usearch/manual/uchime2_algo.html)算法去除嵌合体得到zOTUs.
    `q2-usearch` 插件将 [(Yan, Lin et al.
    2024)](https://link.springer.com/article/10.1186/s13717-023-00480-7)
    发表在Ecological Processes杂志上文章中所用的方法封装进了
    `denoise-no-primer-pooled` 方法中.

    对单端测序下机序列来说, 进行 [global
    trimming](https://drive5.com/usearch/manual/global_trimming.html)
    很有必要.
    去噪生成zOTU过程中两条序列只要出现单碱基的不同就会被识别为不同zOTU,
    而单端测序每条序列3'端终止的位置有可能不同，就会导致来自同一生物模版的序列出现一长一短的现象.
    这批教程中数据处理过程中通过设置参数 `--p-trunc-len 120`
    来实现全局剪切对齐序列. 设置该参数后将输入文件中所有序列在 `120` bp
    的位置截断，并丢弃长度 \< `120` bp的序列.

    ``` bash
    qiime usearch denoise-no-primer-pooled \
        --i-demultiplexed-sequences demux.qza \
        --p-min-size 4 \
        --p-trunc-len 120 \
        --o-representative-sequences rep-seqs-unoise3.qza \
        --o-table table-unoise3.qza \
        --o-denoising-stats stats-unoise3.qza
    ```

    ``` bash
    qiime metadata tabulate \
        --m-input-file stats-unoise3.qza \
        --o-visualization stats-unoise3.qzv
    ```

    我们这里重命名一下文件以便使用UNOISE3结果继续进行下游分析

    ``` bash
    mv rep-seqs-unoise3.qza rep-seqs.qza
    mv table-unoise3.qza table.qza
    ```

2.  双端测序数据处理

    对于双端序列文件，在进行去噪之前我们需要先将他们合并为单端序列.
    这里使用 QIIME2官方教程
    ["沙漠土壤分析"](https://mp.weixin.qq.com/s/xU-3RGzdgH8rIimkDrHIIA).
    请耐心观看教学视频到8分:30秒（运行DADA2之前）我们需要合并双端序列文件.

    `merge-pairs` 方法是从官方插件 `q2-vsearch` 里⌘+c ⌘+v过来的,
    我们把它的默认参数调成了usearch11里的默认值.
    关于双端序列文件合并的参数怎样设置请移步这里 [usearch
    manual](https://rcedgar.github.io/usearch12_documentation/cmd_fastq_mergepairs.html).

    ``` bash
    qiime usearch merge-pairs \
        --i-demultiplexed-seqs demux.qza \
        --o-merged-sequences merged.qza \
        --o-unmerged-sequences unmerged.qza \
        --verbose
    ```

    接下来就可以去噪生成zOTU表与zOTUs啦

    ``` bash
    qiime usearch denoise-no-primer-pooled \
        --i-demultiplexed-sequences merged.qza \
        --p-min-size 4 \
        --o-representative-sequences rep-seqs-unoise3.qza \
        --o-table table-unoise3.qza \
        --o-denoising-stats stats-unoise3.qza
    ```

## 这是一张大饼

这座💩⛰️如果我们还更的动的话（划掉）:

方法:

-   使用SINTAX算法注释特征序列
-   对DADA2输出的ASV进行聚类生成OTUs，然后使用Pacbio CCS 数据生成OTU表.

流程:

-   对去除barcodes的原始序列进行：双端序列合并(仅双端) ➡️ 去除引物 ➡️
    去噪/OTU聚类（一步到位）

-   使用 global search
    将特征序列与参考数据库中的序列进行100%匹配比对，随后使用SINTAX算法对未匹配上的序列进行物种注释

如果使用中遇到任何问题，可以在GitHub里面直接用中文发issue😉

祝大家终端里跑代码再也不报错（不存在的） 🎉🎉🎉
