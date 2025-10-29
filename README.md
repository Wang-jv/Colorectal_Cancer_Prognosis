# 结直肠癌TLS微环境与预后分析项目

## 项目概述

本项目是一个基于多组学数据整合的结直肠癌预后分析研究。项目利用MOVICS (Multi-Omics Integration and Visualization in Cancer Subtyping)分析平台，结合转录组、甲基化、体细胞突变等多维度数据，对结直肠癌患者进行分子分型并分析其与预后的关系。特别关注了肿瘤相关淋巴样结构(Tertiary Lymphoid Structures, TLS)在结直肠癌预后中的作用。

## 数据来源

- TCGA-COAD队列数据
- GEO验证队列(GSE39582)数据
- MSK-IMPACT结直肠癌队列数据

## 分析内容

### 1. 多组学整合分析
- 基因表达谱(mRNA)分析
- 甲基化数据分析
- 体细胞突变分析
- lncRNA表达分析

### 2. 分子分型
- 使用10种算法进行分子分型
- 整合分型结果获得consensus分型
- 分型结果可视化与验证

### 3. 生存分析
- 不同分子亚型间生存差异分析
- 构建6基因预后模型
  - PDZD4
  - PPP1R1A
  - PCOLCE2
  - ACSL6
  - CALB2
  - PTH1R

### 4. 功能分析
- 差异表达基因分析
- KEGG通路富集分析
- 免疫浸润分析
- 药物敏感性分析

### 5. 临床相关性分析
- 临床病理特征关联分析
- 多因素Cox回归分析
- 列线图预测模型构建

## 关键结果

1. 确定了基于TLS特征的两个分子亚型
2. 构建了6基因预后预测模型
3. 发现了潜在的治疗靶点和生物标志物
4. 验证了模型在独立队列中的预测效能

## 代码结构

- `Colorectal_MOVICS_TLS.R`: 主要分析流程代码
- `supplementary_code.R`: 补充分析代码
- `data/`: 原始数据和处理后的数据文件

## 使用说明

### 环境要求
```R
- R >= 4.0
- MOVICS
- survival
- survminer
- ggplot2
- dplyr
- 其他依赖包
```

### 运行步骤
1. 数据预处理和标准化
2. 多组学整合分析
3. 分子分型
4. 生存分析
5. 功能注释
6. 结果可视化

## 结果文件
- 分子分型结果
- 生存分析图
- 免疫组化验证结果
- 富集分析结果
- 预后模型性能评估结果

## 参考文献
项目分析流程参考MOVICS教程：https://xlucpu.github.io/MOVICS/MOVICS-VIGNETTE.html

## 联系方式
[18135079495@163.com]
