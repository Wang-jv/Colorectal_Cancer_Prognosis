library(MOVICS)
library(dplyr)


#原教程链接：https://xlucpu.github.io/MOVICS/MOVICS-VIGNETTE.html

#一、整理数据为标准格式————"2.COAD.tcga.already_match_samples.RData" ———

##这一部分是完全匹配好样本的原始数据，还未进行TLS相关基因的筛选
#————————mRNA.expr————————————————————————
# 选择您要筛选的行名
rows_to_select <- c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19",
                    "CCL21","CXCL9","CXCL10","CXCL11","CXCL13","CXCL13",
                    "CD200","FBLN7","ICOS","SGPP2","SH2D1A","TIGIT","PDCD1",
                    "CD4","CCR5","CXCR3","CSF2","IGSF6","IL2RA","CD38","CD40",
                    "CD5","MS4A1","SDC1","GFI1","IL1R1","IL1R2","IL10","CCL20",
                    "IRF4","TRAF6","STAT5A","TNFRSF17")


# 筛选特定的28行
fpkm2 <-fpkm
mRNA.expr <- fpkm2[rownames(fpkm2) %in% rows_to_select, ]

mRNA.expr <- coad.tcga$mRNA.expr
# 移除批次后的PCA
pca.after <- PCAplot(mRNA.expr, 
                     coad.tcga$clin.info$CANCER_TYPE_ACRONYM, 
                     palette = "npg",
                     legend.title = "TCGA cohort",
                     geom2 = "point") +
  ggtitle("TCGA Remove batch effects PCA") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 14, hjust = 0.5),
    legend.text = element_text(size = 10)
  )
ggsave(
  filename = paste0("Result/", "Batch rectification after pca.png"),
  plot = pca.after,
  width = 6,
  height = 5,
  device = "png",
  dpi = 300
)
ggsave(
  filename = paste0("Result/", "Batch rectification after pca.pdf"),
  plot = pca.after,
  width = 6,
  height = 5,
  device = "pdf"
)

#提取高变异信息lnc ===============

########用01.cor.py做两个表格的相关性，筛选出｜cor｜最大的前40个lncRNA，再导入
write.table(lncRNA.expr,"lncRNA.expr.txt",sep="\t")
write.table(mRNA.expr,"mRNA.expr.txt",sep="\t")
name <- read.table("lncRNA.expr_cor_result_top284.txt",sep="\t",header=T,row.names = 1)
# Get matches 
lncRNA.expr2 <- lncRNA.expr[rownames(lncRNA.expr) %in% 
                              rownames(name),]
####重命名
lncRNA.expr <-lncRNA.expr2


#—————————————————提取高变异信息———————meth.beta——————————————————————————————————————————————————————————————————
meth.beta2 <-meth.beta
# Get matches 
meth.beta3 <- meth.beta2[rownames(meth.beta2) %in% 
                           rownames(mRNA.expr),]
meth.beta <-meth.beta3
#—————————————————提取高变异信息———————mut.status——————————————————————————————————————————————————————————————————
mut.status2 <- mut.status

# Get matches 
mut.status3 <- mut.status2[rownames(mut.status2) %in% 
                             rownames(mRNA.expr),]
mut.status <-mut.status3
#——————————————————————————————————————————————————————————————————————————————————————
coad.tcga <- list(mRNA.expr = mRNA.expr,
                  lncRNA.expr =lncRNA.expr ,
                  meth.beta=meth.beta,
                  mut.status=mut.status ,
                  count=count ,
                  fpkm=fpkm ,
                  maf=maf  ,
                  segment=segment,
                  clin.info=clin.info)

save(coad.tcga, file = "3.COAD.tcga.final.RData")  ##这部分是已经筛选了TLS相关基因作为分群数据的数据

# GEO数据准备==================
load("COAD.GSE39582.mtx.meta.RData")
xdata <- geo_xdata.Scale.raw%>%t()%>%as.data.frame()
ydata <- geo_ydata.raw[,2:5]%>%as.data.frame()
colnames(ydata) <- c("futime","fustat","rf_stime","rfs_status")
ydata$futime<- 30*ydata$futime

coad.yau <- list(mRNA.expr = xdata,
                 clin.info =ydata)
save(coad.yau, file = "3.coad.yau.final.RData")  ##这部分是已经筛选了TLS相关基因作为分群数据的数据


#二、GET Module(获取模块)———————————————————————————————————————————————————————————————————————————————
load("3.COAD.tcga.final.RData")
load("3.coad.yau.final.RData")


#1）获取数据
# print name of example data
names(coad.tcga)
#> [1] "mRNA.expr"   "lncRNA.expr" "meth.beta"   "mut.status"  "count"      
#> [6] "fpkm"        "maf"         "segment"     "clin.info"
names(coad.yau)
#> [1] "mRNA.expr" "clin.info"

# extract multi-omics data
mo.data   <- coad.tcga[1:4]

# extract raw count data for downstream analyses
count    <- coad.tcga$count

# extract fpkm data for downstream analyses
fpkm     <- coad.tcga$fpkm

# extract maf for downstream analysis
maf      <- coad.tcga$maf

# extract segmented copy number for downstream analyses
segment   <- coad.tcga$segment

# extract survival information
surv.info <- coad.tcga$clin.info


#3）获得最佳聚类数
#在任何聚类研究中，最重要的参数是估计数据的最佳聚类数k，其中k需要足够小以减少噪声，但又足够大以保留重要信息。
#此处，MOVICS指的是使用getClustNum()函数通过CPI5和Gaps-statistics6来估计聚类数。

## FIG1A. ==========
# identify optimal clustering number (may take a while)
optk.coad <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,T), #是否是二进制，mo.data共有四个数据，其中第四个是二进制的
                         #（注：第4个数据是体细胞突变，它是一个二进制矩阵） note: the 4th data is somatic mutation which is a binary matrix
                         try.N.clust = 2:8, # try cluster number from 2 to 8
                         fig.name    = "3.CLUSTER NUMBER OF TCGA-coad")

#以上对聚类数的估计给出了3的任意k。然而，用于癌症的流行PAM50分类器有5个分类，如果仔细观察描述性图，
#可以注意到CPI和Gaps-统计数据在k为5时都没有下降太多。因此，考虑到这些知识，选择5的k作为进一步分析的最优聚类数。
#4）从单一算法得到结果
#在这一部分中，我将首先展示如何使用 MOVICS 通过指定一种具有详细参数的算法来执行多组学整合聚类。例如，让我们尝试如下所示的 iClusterBayes 方法：
# perform iClusterBayes (may take a while)
iClusterBayes.res <- 
  getiClusterBayes(data    = mo.data,
                   N.clust     = 5,
                   type        = c("gaussian","gaussian","gaussian","binomial"), #Data type corresponding to the list of matrics, which can be gaussian, binomial or possion.
                   #与矩阵列表相对应的数据类型，可以是高斯、二项式或可能性。
                   n.burnin    = 1800,  #一个整数值，用于指示MCMC burn-in次数。
                   #MCMC burn-in是一个概念，最早来自于统计学和机器学习中的马尔可夫链蒙特卡洛（MCMC）采样。
                   #在开始采样之前，通常会运行一段时间的马尔可夫链以达到平稳分布。这段时间被称为“burn-in”期。
                   #在RNN中，“burn-in”过程是指在开始学习之前先通过网络运行一段时间的序列数据，使得RNN的隐藏状态有机会收敛到某种有意义的状态。
                   #这样，当我们开始学习时，RNN的隐藏状态就已经包含了一些有用的历史信息，而不仅仅是初始的零状态。
                   #这种方法尤其对处理长序列的任务有帮助，因为这些任务可能需要网络记住距离当前时刻较远的历史信息。
                   n.draw      = 1200,    #一个整数值，用于指示MCMC绘制的数量。
                   # MCMC draw的数量是指从后验分布中抽取的样本数量。在MCMC采样中，我们使用马尔可夫链来生成后验分布的样本。
                   #MCMC采样的目的是生成足够多的样本，以便我们可以对后验分布进行准确的估计。
                   #在MCMC采样中，我们通常会运行多个独立的马尔可夫链，以确保我们获得的样本是从后验分布中独立抽取的。
                   #因此，MCMC draw的数量是指从每个独立的马尔可夫链中抽取的样本数量.
                   prior.gamma = c(0.5, 0.5, 0.5, 0.5),  #A numerical vector to indicate the prior probability for the indicator variable gamma of each subdataset.
                   #在统计学中，gamma分布是一种连续概率分布，通常用于建模正值的随机变量。
                   #在贝叶斯统计中，gamma分布通常用作参数的先验分布。在这种情况下，
                   #gamma分布的参数可以被视为超参数，因为它们控制了参数的分布。
                   #在某些情况下，我们可能会使用gamma分布作为指示变量的先验分布。
                   #指示变量是一个二元变量，它的值为0或1，通常用于表示某个事件是否发生。
                   #在这种情况下，gamma分布的参数被解释为指示变量为1的先验概率。
                   sdev        = 0.05,         #A numerical value to indicate the standard deviation of random walk proposal for the latent variable.
                   thin        = 3)            #为了减少自相关，使MCMC链变细的数值。


#如果同时在getMOIC（）中为methodslist参数指定一个算法列表，它将自动逐个执行具有默认参数的每个算法，
#并最终返回从指定算法派生的结果列表。既然iClusterBayes已经完成，让我们同时尝试其他9种默认参数算法。
#这需要一些时间，所以休息一下喝杯咖啡。

# perform multi-omics integrative clustering with the rest of 9 algorithms
pdf("5.moic.res.list.pdf")
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", 
                                            "LRAcluster", "ConsensusClustering",
                                            "IntNMF", "CIMLR", "MoCluster"),
                         N.clust     = 5,
                         type        = c("gaussian", "gaussian", "gaussian", "binomial"))
dev.off()

# attach iClusterBayes.res as a list using append() to moic.res.list with 9 results already
#使用append（）将iClusterBayes.res作为list附加到moic.res.list，已经有9个结果
moic.res.list <- append(moic.res.list, list("iClusterBayes" = iClusterBayes.res))
# save moic.res.list to local path
save(iClusterBayes.res,moic.res.list, file = "4-5.moic.res.list.rda")
#6) get consensus from different algorithms
## FIG1B. ==========
load("4-5.moic.res.list.rda")

cmoic.coad <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "6.CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")

###保存数据
save(iClusterBayes.res,moic.res.list,cmoic.coad, file = "4-6.cmoic.coad.rda")

#※导出10种聚类方法得到的综合聚类的结果
clust <- cmoic.coad$clust.res
write.csv(clust,"clust.csv")
write.csv(fpkm,"fpkm_final.csv")
write.csv(count,"count_final.csv")
write.csv(surv.info,"surv.info.csv")

#7.1) get quantification of similarity using silhoutte
## FIG1C. ==========
load("4-6.cmoic.coad.rda")
getSilhouette(sil      = cmoic.coad$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "7.1.SILHOUETTE",
              height   = 5.5,
              width    = 5)
#7.2) get multi-omics heatmap based on clustering result
## FIG1D. ========
load("4-6.cmoic.coad.rda")



# convert beta value to M value for stronger signal
indata <- mo.data
indata$meth.beta <- log2(indata$meth.beta / (1 - indata$meth.beta))

# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,F), # F:no center for mutation
                     scaleFlag  = c(T,T,T,F)) #F: no scale for mutation


feat   <- iClusterBayes.res$feat.res
feat1  <- feat[which(feat$dataset == "mRNA.expr"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA.expr"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "meth.beta"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "mut.status"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4)

# set color for each omics data
# if no color list specified all subheatmaps will be unified to green and red color pattern
mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, meth.col, mut.col)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = cmoic.coad$clust.res, # cluster results，如果此处改为：iClusterBayes.res$clust.res，则意为用单个iClusterBayes聚类的结果来作为分群数据
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "7.2COMPREHENSIVE HEATMAP OF CONSENSUSMOIC_无表型")

#现在回到cmoic.coad的共识结果，它集成了10个算法，这次还提供了样本的注释来生成热图。
#由于getMoHeatmap（）的核心函数基于ComplexHeatmap R包，因此在创建注释时，
#应始终使用circize:：colorRamp2（）函数为连续变量（例如，本例中的年龄）生成颜色映射函数。

# extract PAM50, pathologic stage and age for sample annotation
annCol    <- surv.info[,c("M_STAGE", "N_STAGE","pstage", "age","SEX"), drop = FALSE]

# generate corresponding colors for sample annotation
annColors <- list(age    = circlize::colorRamp2(breaks = c(min(annCol$age),
                                                           median(annCol$age),
                                                           max(annCol$age)), 
                                                colors = c("#0000AA", "#555555", "#AAAA00")),
                  pstage = c("T1"    = "green",
                             "T2"    = "blue",
                             "T3"    = "red",
                             "T4"    = "yellow", 
                             "T4A"    = "yellow", 
                             "T4B"    = "yellow", 
                             
                             "TX"    = "black",
                             "TIS"    = "black"),
                  
                  M_STAGE = c("M0"    = "green",
                              "M1"    = "blue",
                              "M1A"    = "blue",
                              "M1B"    = "blue",
                              "MX"    = "black"),
                  N_STAGE = c("N0"    = "green",
                              "N1"    = "blue",
                              "N1A"    = "blue",
                              "N1B"    = "blue",
                              "N1C"    = "blue",
                              "N2"    = "red",
                              "N2A"    = "red",
                              "N2B"    = "red",
                              "NX"    = "black"),
                  SEX = c("Male"    = "green",
                          "Female"    = "blue")
)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = cmoic.coad$clust.res, # consensusMOIC results
             clust.dend    = NULL, # show no dendrogram for samples
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(F,F,F,F), # show no dendrogram for features
             color         = col.list,
             annCol        = annCol, # annotation for samples
             annRow        = annRow, # mark selected features
             annColors     = annColors, # annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "7.2COMPREHENSIVE HEATMAP OF CONSENSUSMOIC")





#===========三、 COMP Module=======
#1) compare survival outcome

## FIG1E. ============
load("4-6.cmoic.coad.rda")

surv.coad <- compSurv(moic.res         = cmoic.coad,
                      surv.info        = surv.info,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(5,10), # estimate 5 and 10-year survival
                      fig.name         = "1.KAPLAN-MEIER CURVE OF CONSENSUSMOIC")
#> --a total of 643 samples are identified.
#> --removed missing values.
#> --leaving 642 observations.
#> --cut survival curve up to 10 years.

mut.coad <- compMut(moic.res     = cmoic.coad,
                    mut.matrix   = coad.tcga$mut.status, # binary somatic mutation matrix
                    doWord       = TRUE, # generate table in .docx format
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.01, # keep those genes that mutated in at least 5% of samples
                    p.cutoff = 1, 
                    p.adj.cutoff = 0.5, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                    innerclust   = TRUE, # perform clustering within each subtype
                    annCol       = annCol, # same annotation for heatmap
                    annColors    = annColors, # same annotation color for heatmap
                    width        = 8, 
                    height       = 5,
                    fig.name     = "./result/ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")

# compare TMB
tmb.coad <- compTMB(moic.res     = cmoic.coad,
                    maf          = maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "./result/DISTRIBUTION OF TMB AND TITV")

# change column names of segment data
colnames(segment) <- c("sample","chrom","start","end","value")
fga.coad <- compFGA(moic.res     = cmoic.coad,
                    segment      = segment,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.2, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "./result/BARPLOT OF FGA")


#————整理数据——————————————————————————————————————————————————————————————————

xdata.NoneScale.raw_t <- read.csv("fpkm_final.csv",header = T,row.names = 1)
colnames(xdata.NoneScale.raw_t) <- gsub("\\.", "-", colnames(xdata.NoneScale.raw))
xdata.NoneScale.raw <- xdata.NoneScale.raw_t %>% t()%>% as.data.frame()
ydata.raw <- read.csv("surv.info.csv",header = T,row.names = 1)

##均一化操作
xdata.Scale.raw <- xdata.NoneScale.raw %>% 
  { (apply(., 2, sd) != 0) } %>% 
  { xdata.NoneScale.raw[, .] } %>% 
  scale %>% as.data.frame()


#完全整理好数据之后运行这行代码，保存数据为.RData
save(xdata.NoneScale.raw, xdata.Scale.raw, ydata.raw, file = "COAD.TCGA.mtx_TLS.RData")
load("E:/pulis/202510/2023.9.19_胃肠道癌生信分析_MOVICS_TLS/2023.9.19_胃肠道癌生信分析_MOVICS_TLS/2.建模/COAD.TCGA.mtx_TLS.RData")

# ============= step1:用diff.py来筛选差异 ########
##先在R上整理好格式：
load("4-6.cmoic.coad.rda")
clust <- cmoic.coad$clust.res
load("3.COAD.tcga.final.RData")
fpkm<-coad.tcga$fpkm
fpkm_t <- fpkm %>% t() %>% as.data.frame()
fpkm_t$group <- clust$clust
fpkm_t$group <- fpkm_t$group %>% as.character()
fpkm_t2 <- fpkm_t[, c(32359, 1:32358)]

####删除每组中全为0的列
fpkm_t3 <- fpkm_t2 %>% 
  group_by(group) %>% 
  dplyr::select(where(~!all(. == 0)))
rownames(fpkm_t3)<-rownames(fpkm_t2)
#####删除0占比>20%的列
cols_to_remove <- which(colMeans(fpkm_t3 == 0) > 0.2)
fpkm_t4 <- fpkm_t3[, -cols_to_remove]
rownames(fpkm_t4)<-rownames(fpkm_t3)

fpkm4 <- fpkm_t4 %>%
  t() %>% 
  as.data.frame()

#######写出制表符分隔的表格时，尽量用write.csv，而不是write.table，然后指定分隔符为"\t"
write.table(fpkm4,"fpkm_cluster.csv",sep="\t")

# step2:在本地把各组间差异显著的gene和normal_tumor差异显著的基因取交集，然后用交集基因建模########
###overlap_韦恩图
# 清除当前环境中的变量
rm(list=ls())
# 读取测试数据
data<- read.table("venn.txt",header = TRUE,sep="\t")
head(data)

venn.plot <- venn.diagram(
  x = list(
    deg_our=data$cluster_top, normal_tumor=data$normal_tumor
  ),
  filename = NULL,
  col = "black",#边框颜色
  lwd = 2, # 边框线的宽度
  fill = c("#f5cac3","#588b8b"),
  alpha = 0.70,
  label.col = "black",
  cex = 2.0,
  cat.col = c("#f5cac3","#588b8b"),
  cat.cex = 1.5,#类名字体大小
  margin = 0.04, #边际距离
  scaled = F
)
pdf(file="venn.pdf",width = 5,height =5)
grid.draw(venn.plot)
dev.off()





#======step3:=cox单因素回归==================
##筛选出相关基因的列
xdata.NoneScale.raw.select <- xdata.NoneScale.raw%>% 
  as.data.frame() %>% 
  dplyr::select(
    COMP,
    NPTXR,
    CILP,
    PRELP,
    HAND2,
    ADGRD1,
    NRXN2,
    PTGIS,
    THBS4,
    AOX1,
    L1CAM,
    RNF150,
    OGN,
    MAP6,
    CSDC2,
    FBXO17,
    CHRDL1,
    ADAMTSL3,
    WASF3,
    PDZRN4,
    TNXB,
    ANK2,
    `HAND2-AS1`,
    CASQ2,
    NGFR,
    SYPL2,
    CLDN11,
    ANGPTL1,
    PLA2G5,
    PDZD4,
    GAP43,
    BHMT2,
    AGTR1,
    PLIN4,
    SCN7A,
    SFRP1,
    CADM3,
    ATP1A2,
    CRTAC1,
    NECAB1,
    LDB3,
    CPEB1,
    CCBE1,
    ABCA6,
    PI16,
    GPM6A,
    FABP4,
    BCHE,
    NOVA1,
    MASP1,
    PPP1R1A,
    ALDH1L1,
    TRIM9,
    NPTX1,
    PLP1,
    FMN2,
    TMEM59L,
    LONRF2,
    SLC8A2,
    PRIMA1,
    DPP6,
    SORCS1,
    RBFOX3,
    MAPK4,
    PTPRZ1,
    ABCA8,
    ABCA9,
    ABI3BP,
    `ADAMTS9-AS1`,
    ADH1B,
    BNC2,
    C1QTNF7,
    C2orf74,
    CCL18,
    CD163,
    CD36,
    CHL1,
    CLIP4,
    COL10A1,
    CP,
    CXCL10,
    CXCL11,
    CXCL9,
    DPY19L2,
    F13A1,
    FAM110B,
    FGF7,
    FGFBP2,
    FGL2,
    GALNT15,
    GAS1,
    GFPT2,
    GLI3,
    `IGLV3-21`,
    KCNMA1,
    LIFR,
    LILRB5,
    LRCH2,
    LYVE1,
    `MAGI2-AS3`,
    MAL,
    MAMDC2,
    MFAP5,
    MMRN1,
    PCOLCE2,
    PDE1B,
    PRDM8,
    RIC3,
    SIGLEC1,
    SLIT2,
    STEAP4,
    SYNE1,
    TMOD1,
    TMTC1,
    TPO,
    TWIST2,
    VSIG4,
    ZBTB16,
    ACSL6,
    `ADAMTS9-AS2`,
    BEND5,
    BOC,
    BVES,
    C14orf132,
    C16orf89,
    C7,
    C8orf88,
    CALB2,
    CELP,
    CHRDL2,
    CNR1,
    DNER,
    FAM107A,
    FGF2,
    FILIP1,
    GALNT16,
    GRIK5,
    GRM8,
    GSTM5,
    HSPA7,
    JAM2,
    JPH2,
    KCNMB1,
    LMO3,
    LY6G6D,
    LYNX1,
    MEIS2,
    MITF,
    MPDZ,
    NEXN,
    NTSR1,
    OSR1,
    PDE2A,
    PEG10,
    PHLDB2,
    PHYHIP,
    PLN,
    PRPH,
    PTGER3,
    PTH1R,
    PYGM,
    RBP2,
    SCHIP1,
    SLIT3,
    SNCA,
    SPEG,
    STON1,
    TAGLN3,
    TMEM100,
    TTLL7,
    UST) %>%
  as.data.frame() 


##导入time\status
xdata.NoneScale.raw.select$time <- ydata.raw$time
xdata.NoneScale.raw.select$status <- ydata.raw$status

#调整顺序
cols <- colnames(xdata.NoneScale.raw.select)
new_cols <- c(cols[length(cols)],cols[length(cols)-1], cols[1:(length(cols) - 2)])

# 然后将 dataframe 按照新的列名顺序排列
xdata.NoneScale.raw.select2 <- xdata.NoneScale.raw.select[, new_cols]

##计算
pFilter=1
outResult=data.frame()
sigGenes=c("time","status")
for(i in colnames(xdata.NoneScale.raw.select2[,3:ncol(xdata.NoneScale.raw.select2)])){
  mycox <- coxph(Surv(time, status) ~ xdata.NoneScale.raw.select2[,i], data = xdata.NoneScale.raw.select2)
  mycoxSummary = summary(mycox)
  pvalue=mycoxSummary$coefficients[,"Pr(>|z|)"]
  if(pvalue<pFilter){
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=mycoxSummary$conf.int[,"exp(coef)"],
                          L95CI=mycoxSummary$conf.int[,"lower .95"],
                          H95CI=mycoxSummary$conf.int[,"upper .95"],
                          pvalue=mycoxSummary$coefficients[,"Pr(>|z|)"])
    )
    
  }
}

##写出结果文件
write.csv(outResult, file='outResult_建模.csv')

# 单因素森林图 ==================
### FIG3B. ==========
outResult<-read.csv("outResult.csv", header = T, sep = ",")
outResult %>%
  dplyr::mutate(p.adjust = signif(p.adjust(pvalue, method = "BH"),digits = 3),
                lower = signif(lower, digits = 3),
                upper = signif(upper, digits = 3),
                HR = signif(HR, digits = 3),
                pvalue = signif(pvalue, digits = 3),
                `HR(95% CI)` = paste0(
                  HR, " (",
                  lower, "-", upper, ")"
                )) %>% 
  write.csv("single_cox_outResult.csv",row.names = F)
#森林图绘制：
photo <- outResult %>%
  dplyr::mutate(p.adjust = signif(p.adjust(pvalue, method = "BH"),digits = 3),
                lower = signif(lower, digits = 3),
                upper = signif(upper, digits = 3),
                HR = signif(HR, digits = 3),
                pvalue = signif(pvalue, digits = 3),
                `HR(95% CI)` = paste0(
                  HR, " (",
                  lower, "-", upper, ")"
                )) %>%
  forestplot(labeltext = c(id,`HR(95% CI)`,p.adjust), 
             clip = c(0.5, 3),#森林图中置信区间上下限设置
             xlog = F,
             zero=2)%>%
  fp_set_style(box = "#588b8b", #权重方块颜色
               line = "#588b8b") %>%#是否对X轴取对数
  fp_set_zebra_style("#ffebe8") %>% #更改为斑马样式，并修改汇总菱形的填充颜色
  fp_add_header(id = c("","gene"),
                `HR(95% CI)` =c("","HR(95% CI)"),
                p.adjust = c("","p.adjust"))
pdf(file="single_cox_forestPlot.pdf",width = 7,height =8)
print(photo)
dev.off()

#############单因素森林图模型里6个基因
## FIG3C. ==========
outResult<-read.table("out_result_6gene.txt", header = T, sep = "\t")
#森林图绘制：
photo <- outResult %>%
  dplyr::mutate(p.adjust = signif(p.adjust(pvalue, method = "BH"),digits = 3),
                lower = signif(lower, digits = 3),
                upper = signif(upper, digits = 3),
                HR = signif(HR, digits = 3),
                pvalue = signif(pvalue, digits = 3),
                `HR(95% CI)` = paste0(
                  HR, " (",
                  lower, "-", upper, ")"
                )) %>%
  dplyr::filter(id %in% c("PDZD4",
                          "PPP1R1A",
                          "PCOLCE2",
                          "ACSL6",
                          "CALB2",
                          "PTH1R")) %>%
  forestplot(labeltext = c(id,`HR(95% CI)`,p.adjust), 
             clip = c(0.5, 3),#森林图中置信区间上下限设置
             xlog = F,
             zero=2)%>%
  fp_set_style(box = "#588b8b", #权重方块颜色
               line = "#588b8b") %>%#是否对X轴取对数
  fp_set_zebra_style("#ffebe8") %>% #更改为斑马样式，并修改汇总菱形的填充颜色
  fp_add_header(id = c("","gene"),
                `HR(95% CI)` =c("","HR(95% CI)"),
                p.adjust = c("","p.adjust"))
pdf(file="6genes_cox_forestPlot.pdf",width = 7,height =3)
print(photo)
dev.off()           

#============step4:lasso:建模用n个cox显著基因======================

small.subset.ll <- c(
  "PCOLCE2",
  "PPP1R1A",
  "PTH1R",
  "PDE1B",
  "PDZD4",
  "CALB2",
  "CLDN11",
  "PLIN4",
  "FABP4",
  "TMEM59L",
  "PRELP",
  "NPTX1",
  "OSR1",
  "CHRDL1",
  "ADGRD1",
  "CILP",
  "PLA2G5",
  "COMP",
  "SYPL2",
  "RIC3",
  "PTGIS",
  "CD36",
  "CSDC2",
  "MEIS2",
  "HSPA7",
  "NGFR",
  "GRIK5",
  "ABCA9",
  "THBS4",
  "NOVA1",
  "CADM3",
  "CLIP4",
  "SIGLEC1",
  "PHYHIP",
  "ABCA6",
  "MAP6",
  "KCNMA1",
  "GPM6A",
  "SPEG",
  "SLC8A2",
  "LONRF2",
  "DNER",
  "NEXN",
  "PTPRZ1",
  "TMOD1",
  "CRTAC1",
  "GAP43",
  "BCHE",
  "MAMDC2",
  "BOC",
  "BHMT2",
  "JPH2",
  "ABI3BP",
  "GAS1",
  "LDB3",
  "PRPH",
  "ACSL6",
  "PLN",
  "PEG10",
  "CCBE1",
  "ANK2",
  "PDZRN4",
  "BVES",
  "KCNMB1",
  "CASQ2",
  "MAL",
  "C7",
  "RNF150",
  "ATP1A2",
  "SYNE1",
  "DPY19L2",
  "AGTR1",
  "C2orf74",
  "GALNT16",
  "ABCA8",
  "NPTXR",
  "HAND2",
  "VSIG4",
  "`HAND2-AS1`",
  "PTGER3",
  "UST")




##整理数据
xdata.ll <- xdata.Scale.raw[, small.subset.ll[small.subset.ll %in% colnames(xdata.Scale.raw)]]
ydata.ll <- ydata.raw %>% dplyr::select(time, status)
ydata.ll <- ydata.ll[rownames(xdata.ll), ]


xdata.ll <-xdata.ll %>%as.matrix()
ydata.ll <-ydata.ll %>%as.data.frame()

## 建模 ================================
set.seed(123456)

fitted.ll <- cv.glmHub(xdata.ll, Surv(ydata.ll$time, ydata.ll$status),
                       family  = 'cox',
                       lambda = buildLambda(1),
                       network = 'correlation', 
                       network.options = networkOptions(cutoff = .001, #这里可以修改参数
                                                        min.degree = .5))
coefs.ll.v <- coef(fitted.ll, s = 'lambda.min')[,1] %>% { .[. != 0]}
summary(fitted.ll)
lassoLambdat(model = fitted.ll, file.out = "./")
laassDeviance(model = fitted.ll, file.out = "./")

coefs.ll.v %>% { 
  data.frame(gene.name   = names(.),
             coefficient = .,
             stringsAsFactors = FALSE)
} %>%
  arrange(gene.name) %>%
  write.csv("lasso_cox_coefficients.csv", row.names = FALSE)

#####找到这个比较好的模型之后，要作生存分析的图
##作生存分析图
pdf("polyII-related-survival-in-COAD.pdf")
separate2GroupsCox(as.vector(coefs.ll.v), 
                   xdata.ll[, names(coefs.ll.v)], 
                   ydata.ll, 
                   plot.title = 'Full dataset', 
                   legend.outside = FALSE)###结果输出的$pvalue值越小越好。
dev.off()

##保存以上模型信息
save(xdata.ll, ydata.ll, fitted.ll, coefs.ll.v, file = "5geneSignature-20240111COAD_1.537021e-05.RData")


load("E:/pulis/202510/2023.9.19_胃肠道癌生信分析_MOVICS_TLS/2023.9.19_胃肠道癌生信分析_MOVICS_TLS/2.建模/5geneSignature-20240111COAD_1.537021e-05.RData")
#========================== 基因多因素cox ========================
########确定好模型后，分别在Scale和Nonescale数据表中加入risk分群列，然后保存分群后的数据信息
# CNOT1 CPM EPAS1 LRPPRC NCKAP1 NDUFB6 OXSM POU4F1 RPN1 SLC7A11 ZCCHC14 ZNF484
####Scale数据
##转换为数据框
xdata.tmp <- xdata.Scale.raw %>% as.data.frame()
##用得到的系数计算指数
index <- -xdata.tmp$ACSL6 *0.0015881+
  xdata.tmp$CALB2*0.0368085 +
  xdata.tmp$PCOLCE2  *0.1360931 +
  xdata.tmp$PDZD4 *0.0444171 +
  xdata.tmp$PPP1R1A*0.0500747 +
  xdata.tmp$PTH1R *0.0488492
##把生成的指数这一列加到原来的表格中
xdata.tmp$index <- index
##根据index排序,并分为2群
xdata.Scale.risk<-xdata.tmp[order(xdata.tmp$index,decreasing = T),]
group <- c(rep("High-risk",241),rep("Low-risk",242))
cluster <- c(rep("1",241),rep("2",242))

xdata.Scale.risk$risk <- group
xdata.Scale.risk$cluster <- cluster

####NoneScale数据
##转换为数据框
xdata.NoneScale.tmp <- xdata.NoneScale.raw %>% as.data.frame()
##把生成的指数这一列加到原来的表格中
xdata.NoneScale.tmp $index <- index
##根据index排序,并分为2群
xdata.NoneScale.risk<-xdata.NoneScale.tmp[order(xdata.NoneScale.tmp$index,decreasing = T),]
group <- c(rep("High-risk",241),rep("Low-risk",242))
cluster <- c(rep("1",241),rep("2",242))

xdata.NoneScale.risk$risk <- group
xdata.NoneScale.risk$cluster <- cluster

####ydata.risk数据
##重新对xdata和ydata的3个表格排序，使其sample 编号顺序一致，方便后续操作
xdata.NoneScale.risk <- xdata.NoneScale.risk[rownames(xdata.NoneScale.risk) %in% 
                                               rownames(ydata.raw),]
xdata.Scale.risk <- xdata.Scale.risk[rownames(xdata.Scale.risk) %in% 
                                       rownames(ydata.raw),]
ydata.raw_order  <- ydata.raw[rownames(xdata.Scale.risk), ]

##把分群的这一列信息加到ydata表格中，方便后续操作。
ydata.risk<-ydata.raw_order%>%
  as.data.frame()
ydata.risk$risk <- xdata.Scale.risk$risk
ydata.risk$index <- xdata.Scale.risk$index
ydata.risk$cluster <-xdata.Scale.risk$cluster

#读入包含add表型的表格
clin.add <-read.table("clin_add.txt",sep="\t",header=T,row.names = 1)
clin.add2     <- clin.add[rownames(ydata.risk), ]
#####合并ydata
ydata.risk_add <- merge(ydata.risk,clin.add2, by = "row.names")
rownames(ydata.risk_add)<-ydata.risk_add$Row.names
ydata.risk_add<-ydata.risk_add[,-1]
####调整列顺序
cols <- colnames(ydata.risk_add)
new_cols <- c(cols[2],cols[1],cols[3:35],cols[38],cols[39:46],cols[36:37])

# 然后将 dataframe 按照新的列名顺序排列
ydata.risk_add2 <-ydata.risk_add[, new_cols]


######修改名字
ydata.risk <- ydata.risk_add2

######再保存一份ydata里离散指标转换为数字的表格。
#导出数据
#write.csv(ydata.risk, file='ydata.risk.number.csv')
##这一步需要整理好一个把离散参数变为数字的ydata表格：ydata.raw.number_forma.cs，再导入

ydata.raw.number_forma <- read.csv("ydata.raw.number_forma.csv",header=T,row.names = 1)
#####加入index和risk
ydata.raw.number_forma <- ydata.raw.number_forma[rownames(ydata.risk), ]
ydata.raw.number_forma$risk<-ydata.risk$risk
ydata.raw.number_forma$index<-ydata.risk$index
ydata.raw.number_forma$age <-ydata.risk$age
ydata.risk.number_forma<-ydata.raw.number_forma

clin.add2     <- clin.add[rownames(ydata.risk.number_forma), ]
#####合并ydata
ydata.risk.number_forma_add <- merge(ydata.risk.number_forma,clin.add2, by = "row.names")
rownames(ydata.risk.number_forma_add )<-ydata.risk.number_forma_add $Row.names
ydata.risk.number_forma_add <- ydata.risk.number_forma_add[,-1]

####调整列顺序
cols <- colnames(ydata.risk.number_forma_add)
new_cols <- c(cols[2],cols[1],cols[3:35],cols[38],cols[39:46],cols[36:37])

# 然后将 dataframe 按照新的列名顺序排列
ydata.risk.number_forma_add2 <-ydata.risk.number_forma_add[, new_cols]


######修改名字
ydata.risk.number_forma <- ydata.risk.number_forma_add2

#完全整理好数据之后运行这行代码，保存包含分群信息的数据为.RData
save(xdata.NoneScale.raw, xdata.Scale.raw, ydata.raw, xdata.NoneScale.risk, 
     xdata.Scale.risk, ydata.risk,ydata.risk.number_forma, 
     file = "COAD.TCGA.mtx.meta_already_cluster_TLS.RData")

# 列线图====================
## FIG4D. ==========
## 添加变量标签
surv.info$pstage <- factor(surv.info$pstage, levels = c("T1","T2","T3","T4","T4A","T4B"))
surv.info$SEX <- factor(surv.info$SEX, levels = c("Male","Female"))
surv.info$M_STAGE <- factor(surv.info$M_STAGE, levels = c("MX","M0","M1","M1X","M1B"))
surv.info$N_STAGE <- factor(surv.info$N_STAGE, levels = c("NX","N0","N1","N1A","N1B","N1C","N2","N2A","N2B"))
surv.info2<-surv.info
# Get matches 
surv.info2     <- surv.info2 [rownames(ydata.risk), ]

surv.info2$index<-ydata.risk$index
#data(surv.info2)
head(surv.info2)
## 根据nomogram要求处理数据
dd=datadist(surv.info2)
options(datadist="dd")
f2 <- psm(Surv(futime,fustat) ~ index+age+pstage+M_STAGE+N_STAGE, data =  surv.info2, dist='lognormal') 
med <- Quantile(f2) # 计算中位生存时间
surv <- Survival(f2) # 构建生存概率函数

## 绘制COX回归中位生存时间的Nomogram图
nom <- nomogram(f2, fun=function(x) med(lp=x),
                funlabel="Median Survival Time")
plot(nom)

## 绘制COX回归生存概率的Nomogram图
## 注意surv.info2数据的time是以’天‘为单位
nom <- nomogram(f2, fun=list(function(x) surv(365, x),
                             function(x) surv(1825, x),
                             function(x) surv(3650, x)),
                funlabel=c("1-year Survival Probability",
                           "5-year Survival Probability",
                           "10-year Survival Probability"))
plot(nom, xfrac=.6)


## 评价COX回归的预测效果
## 计算c-index
rcorrcens(Surv(futime,fustat) ~ predict(f2), data =  surv.info2)

# Somers' Rank Correlation for Censored Data    Response variable:Surv(futime, fustat)
# 
#                C   Dxy  aDxy    SD    Z      P   n
# predict(f2) 0.59 0.179 0.179 0.087 2.06 0.0391 287
## 重新调整模型函数f2，也即添加x=T, y=T
f2 <- psm(
  Surv(futime,fustat) ~age+SEX+pstage+M_STAGE+N_STAGE+index,
  data =  surv.info2, x=T, y=T, dist='lognormal') 

## 构建校正曲线
cal1 <- calibrate(f2, 
                  cmethod='KM', 
                  method="boot", 
                  u=365,# 需要与之前模型中定义好的time.inc一致，即365或730；
                  m=80, #每次抽样的样本量，
                  B=1000,
                  smoother="x") 
plot(cal1)

cal2 <- calibrate(f2, 
                  cmethod=c('hare', 'KM'),
                  method="boot", 
                  u=365, 
                  m=150,  ##m要根据样本量来确定，由于标准曲线一般将所有样本分为3组（在图中显示3个点）
                  B=40, 
                  bw=FALSE, 
                  rule="aic", 
                  type="residual", sls=0.05, aics=0, force=NULL,
                  estimates=TRUE,
                  pr=FALSE, what="observed-predicted", tol=1e-12, maxdim=5)

## FIG4E. 绘制校正曲线 =====
plot(cal1,lwd=2,lty=1,
     conf.int=T,# 是否显示置信区间
     errbar.col="blue",#直线曲线bar颜色
     col="red", # 曲线颜色
     xlim=c(0.6,1),ylim=c(0.6,1),
     xlab="Nomogram-Predicted Probability of 1-Year DFS",
     ylab="Actual 1-Year DFS (proportion)",
     subtitles = F)#不显示副标题

#——方法2————————————————————————————————————————————————————————————————————————————————————————————

library(riskRegression)
library(survival)       
library(prodlim)
library(pec)
# 导入mgus2数据集
dt<- surv.info[,1:7]           
# 数据转换
dt <- na.omit(dt) # 删除缺失值          
dt$etime <-dt$futime       
dt$event <- dt$fustat     
dt$event <- factor(dt$event, 0:2)
#2、构建竞争风险模型
f1 <- CSC(Hist(etime,event)~age+pstage+SEX,
          data = dt)          
f2 <- CSC(Hist(etime,event)~age+mspike,
          data = dt)

#3、绘制ROC曲线
x_roc <- Score(list(model1=f1, model2=f2),          
               Hist(etime,event)~1,          
               data=dt,          
               cause= 1,     # 结局事件
               times=360,          
               plots="roc",          
               metrics = "auc")
plotROC(x_roc)
#4、绘制校准曲线

#4.1 多个竞争风险模型校准曲线

x_cal <- Score(list(model1=f1,model2=f2),          
               Hist(etime,event)~1,          
               data=dt,          
               cause= 1,          
               times=360,          
               plots="cal",          
               metrics = "auc")         
plotCalibration(x_cal)


#4.2 不同时间点的校准曲线

calPlot(f1,time = 120,          
        xlim = c(0,0.8),ylim = c(0,0.8),          
        col = "darkcyan",legend = F)          
calPlot(f1,time = 240,          
        xlim = c(0,0.8),ylim = c(0,0.8),          
        add = T,col = "tomato")          
calPlot(f1,time = 360,          
        xlim = c(0,0.8),ylim = c(0,0.8),          
        add = T,col = "purple")          
legend("bottomright",          
       legend=c("10年","20年","30年"),          
       col=c("darkcyan","tomato","purple"),          
       lwd=2,cex=0.7,bty='n')









#生存状态填充柱状图-原模型图
## FIGS1C. =========

ydata.risk2<-ydata.risk %>%
  ggplot(aes(x=risk,fill=factor(status)))+
  geom_bar(stat="count",position='fill',width = .6) +
  geom_text(stat='count',aes(label=after_stat(count)), color="white", size=7,position=position_fill(0.5))+
  scale_fill_manual(values=c("#f5cac3","#588b8b"))+
  coord_flip()+
  theme(panel.background=element_rect(fill='transparent'),
        axis.line = element_line(arrow = arrow(length = unit(0.15, 'cm')),
                                 colour = "black"),panel.grid =element_blank(),
        text=element_text(size= 50),axis.text = element_text(colour="black"))
ggsave("生存状态填充柱状图.PDF",ydata.risk2,width=14,height=5.5)
write.csv(ydata.risk, file='生存状态填充柱状图.csv')


#—————————————————————TLS热图—————————————————————————————
# 导入数据并调整样本顺序
df2   <- read.csv("TLS_ssgsea.csv", header = T, sep = ",")
df <- t(df2)
colnames(df) <- df[1,]
df <- df[-1,]
rownames(df) <- gsub("\\.", "-", rownames(df))
rownames(df) <- paste0(rownames(df), "A")
df<-df%>%as.data.frame()
# Get matches 
df2 <- df [rownames(df) %in% 
             rownames(xdata.NoneScale.risk),]
df2     <- df2 [rownames(xdata.NoneScale.risk), ]
df2 <- df2 %>% as.data.frame()
########排序
df2$risk <- xdata.NoneScale.risk$risk
df_order <-df2[order(df2$risk,decreasing = F),]
sample_info<- df_order $risk %>%as.data.frame()
sample_info$. <- df_order $risk
rownames(sample_info) <- rownames(df_order)
df_order <- df_order[,-3]

df_t <- t(df_order)
df_t <- df_t%>%as.data.frame()

#可视化-热图
#####热图数据变为数值型
df_t2 = as.data.frame(lapply(df_t,as.numeric))
colnames(df_t2)=colnames(df_t)

##作图，加了annotation_col 参数，设置 annotation_col =sample_info
#————————————————————————列不聚类——————————————————————————————
tu<-
  pheatmap(df_t2,scale="row",border=NA,
           color=colorRampPalette(c("#004444","white","#c10017"))(50),
           cluster_row=T,cluster_cols=F,
           gaps_col = c(251),
           angle_col=45, show_rownames=T,show_colnames =T ,
           treeheight_col =10,treeheight_row =15,legend = T,
           annotation_col =sample_info )

pdf(file="pheatmap_TLS_nocluster.pdf",width = 7,height =2)
print(tu)
dev.off()

#————————————————————————列聚类——————————————————————————————
tu<-
  pheatmap(df_t2,scale="row",border=NA,
           color=colorRampPalette(c("#004444","white","#c10017"))(50),
           cluster_row=T,cluster_cols=T,
           # gaps_col = c(338,452),
           angle_col=45, show_rownames=T,show_colnames =T ,
           treeheight_col =10,treeheight_row =15,legend = T,
           annotation_col =sample_info )

pdf(file="pheatmap_TLS.pdf",width = 7,height =2)
print(tu)
dev.off()

#TLS分数差异box图
# 导入数据并调整样本顺序
df2   <- read.csv("TLS_ssgsea.csv", header = T, sep = ",")
df <- t(df2)
colnames(df) <- df[1,]
df <- df[-1,]
rownames(df) <- gsub("\\.", "-", rownames(df))
rownames(df) <- paste0(rownames(df), "A")
df<-df%>%as.data.frame()
# Get matches 
df2 <- df [rownames(df) %in% 
             rownames(xdata.NoneScale.risk),]
df2     <- df2 [rownames(xdata.NoneScale.risk), ]
df2 <- df2 %>% as.data.frame()
df2$risk <- xdata.NoneScale.risk$risk
df2$TLS<- as.numeric(df2$TLS)
pdf("TLS_box.pdf",width=4,height= 4)
ggplot(df2, aes(x =factor(risk), y=log2(TLS) ,fill =factor(risk))) +
  geom_boxplot(size=0.5,varwidth = F, alpha=0.8) + 
  geom_jitter(position=position_jitter(width=0.15,height=0.2), size=0.8,alpha=0.3,colour="#bfdbf7") + 
  # facet_wrap(cell~.,scales="free",ncol = 10)+
  scale_fill_manual(values=c("#588b8b","#f5cac3")) +
  theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
        axis.text.x = element_text(angle = 30,vjust = 0.9,hjust = 1), panel.grid =element_blank(),axis.text = element_text(colour="black"),
        strip.background = element_rect(colour="black", fill="white"))+
  geom_signif(
    comparisons = list( c("High-risk","Low-risk")),
    map_signif_level = F, textsize = 6,step_increase = 0.1,test=t.test
  ) 
dev.off()

#——————确认完之后，把分组信息加到原来表格————————————————————————————————————————————————————————————————————————————————————————————————————
clust <- cmoic.coad$clust.res
ydata.risk$cluster <- as.numeric(ydata.risk$cluster)
# Get matches
ydata.risk.tmp <- ydata.risk[rownames(ydata.risk) %in% 
                               rownames(clust),]
ydata.risk.tmp     <- ydata.risk.tmp [rownames(clust), ]
#####把分组信息加到cmoic.coad并重命名为cmoic.coad_TLS
clust$cluster <- ydata.risk.tmp$cluster
cmoic.coad_TLS<-cmoic.coad
cmoic.coad_TLS$clust.res$clust<-ydata.risk.tmp$cluster

###保存数据
save(iClusterBayes.res,moic.res.list,cmoic.coad,cmoic.coad_TLS, file = "4-6.cmoic.coad_TLS.rda")


#7.2) get multi-omics heatmap based on clustering result


load("4-6.cmoic.coad_TLS.rda")

# convert beta value to M value for stronger signal
indata <- mo.data
indata$meth.beta <- log2(indata$meth.beta / (1 - indata$meth.beta))

# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,F), # F:no center for mutation
                     scaleFlag  = c(T,T,T,F)) #F: no scale for mutation


feat   <- iClusterBayes.res$feat.res
feat1  <- feat[which(feat$dataset == "mRNA.expr"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA.expr"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "meth.beta"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "mut.status"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4)

# set color for each omics data
# if no color list specified all subheatmaps will be unified to green and red color pattern
mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, meth.col, mut.col)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = cmoic.coad_TLS$clust.res, # cluster results，如果此处改为：iClusterBayes.res$clust.res，则意为用单个iClusterBayes聚类的结果来作为分群数据
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "7.2COMPREHENSIVE HEATMAP OF CONSENSUSMOIC_无表型")

# extract PAM50, pathologic stage and age for sample annotation
annCol    <- surv.info[,c("M_STAGE", "N_STAGE","pstage", "age","SEX"), drop = FALSE]

# generate corresponding colors for sample annotation
annColors <- list(age    = circlize::colorRamp2(breaks = c(min(annCol$age),
                                                           median(annCol$age),
                                                           max(annCol$age)), 
                                                colors = c("#0000AA", "#555555", "#AAAA00")),
                  pstage = c("T1"    = "green",
                             "T2"    = "blue",
                             "T3"    = "red",
                             "T4"    = "yellow", 
                             "T4A"    = "yellow", 
                             "T4B"    = "yellow", 
                             
                             "TX"    = "black",
                             "TIS"    = "black"),
                  
                  M_STAGE = c("M0"    = "green",
                              "M1"    = "blue",
                              "M1A"    = "blue",
                              "M1B"    = "blue",
                              "MX"    = "black"),
                  N_STAGE = c("N0"    = "green",
                              "N1"    = "blue",
                              "N1A"    = "blue",
                              "N1B"    = "blue",
                              "N1C"    = "blue",
                              "N2"    = "red",
                              "N2A"    = "red",
                              "N2B"    = "red",
                              "NX"    = "black"),
                  SEX = c("Male"    = "green",
                          "Female"    = "blue")
)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = cmoic.coad_TLS$clust.res, # consensusMOIC results
             clust.dend    = NULL, # show no dendrogram for samples
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(F,F,F,F), # show no dendrogram for features
             color         = col.list,
             annCol        = annCol, # annotation for samples
             annRow        = annRow, # mark selected features
             annColors     = annColors, # annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "7.2COMPREHENSIVE HEATMAP OF CONSENSUSMOIC")


# 1) survival comparison
load("4-6.cmoic.coad_TLS.rda")
surv.coad <- compSurv(moic.res         = cmoic.coad_TLS,
                      surv.info        = surv.info,
                      convt.time       = "m", # convert day unit to month
                      clust.col = rev(c("#2EC4B6", "#E71D36")),
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(5,10), # estimate 5 and 10-year survival
                      fig.name         = "1.KAPLAN-MEIER CURVE OF CONSENSUSMOIC")
#> --a total of 643 samples are identified.
#> --removed missing values.

print(surv.coad)
dev.copy2pdf(file = "1.KAPLAN-MEIER CURVE OF CONSENSUSMOIC.PDF",width=5,height=6)


#2)比较临床特征
#然后比较不同亚型之间的临床特征。MOVICS 提供的功能compClinvar()可以总结连续变量和
#分类变量并执行适当的统计测试。该函数可以给出.docx格式的表格，易于在医学研究论文中使用。
surv.info2<-surv.info[,1:7]
clin.coad <- compClinvar(moic.res      = cmoic.coad_TLS,
                         var2comp      = surv.info2, # data.frame needs to summarize (must has row names of samples)
                         strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         factorVars    = c("pstage","fustat","SEX","M_STAGE","N_STAGE"), # features that are considered categorical variables
                         nonnormalVars = c("futime","age"), # feature(s) that are considered using nonparametric test
                         exactVars     = c("pstage","M_STAGE","N_STAGE"), # feature(s) that are considered using exact test
                         doWord        = TRUE, # generate .docx file in local path
                         tab.name      = "2.SUMMARIZATION OF CLINICAL FEATURES")
#> --all samples matched.
print(clin.coad$compTab)

#4）比较总突变负荷
#不用说，免疫疗法正在成为现代癌症治疗的支柱。最近的分析将肿瘤基因组景观与抗肿瘤免疫联系起来。
#特别是，一项新的研究表明，肿瘤特异性基因组病变与免疫检查点激活以及患者对免疫疗法的反应程度和持续时间有关。
#这些病变包含高突变负荷8和非整倍体9。为了量化这些可能影响免疫治疗的基因组改变，
#MOVICS提供了两种功能来计算总突变负荷（TMB）和基因组改变部分（FGA）
#。具体来说，TMB是指在肿瘤基因组中发现的突变数量，而FGA是受拷贝数增加或减少影响的基因组百分比。
#这两个属性对基因研究人员都很有用，因为它们为他们提供了关于肿瘤基因组构成的更深入的信息。
#让我从compTMB（）开始，向您展示如何使用这两个函数。首先，此函数的输入maf数据必须至少具有以下10列：
head(maf)
#然后运行compTMB（）。默认情况下，此函数仅在计算体细胞突变频率时考虑非同义变体，
#包括Frame_Shift_Del、Frame_Shift_Ins、Splice_Site、Translation_Start_Site、Nonsense_mutation、
#Nonstop_mutation、In_Frame_Del、In_FFrame_Ins和Missense_mutation。
#除了计算TMB，该函数还将单核苷酸变体分为转变和转变（TiTv），并描述TMB和TiTv的分布。
# compare TMB

tmb.coad <- compTMB(moic.res     = cmoic.coad_TLS,
                    maf          = maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "4.DISTRIBUTION OF TMB AND TITV")

head(tmb.coad$TMB.dat)
TMB.dat <- tmb.coad$TMB.dat
write.csv(TMB.dat,"4.Demo of comparison of TMB among 5 identified subtype of breast cancer in TCGA-COAD cohort.csv")

#5) compare fraction genome altered
#接下来，compFGA（）不仅计算FGA，还计算每个亚型中每个样本的特定增益（FGG）或损失（FGL）。
#为了实现此功能，应使用完全相同的列名来准备符合条件的分段副本编号输入，如下所示：      
# change column names of segment data
colnames(segment) <- c("sample","chrom","start","end","value")
#Let’s see how the input looks like:

head(segment)
#   Then this function can be run without any difficulties. Notably, if your CNA calling procedure did not provide a segmented copy number as value column but the original copy number, argument of iscopynumber must be switched to TRUE instead.
#那么这个功能就可以毫无困难地运行了。值得注意的是，如果您的CNA调用过程没有提供分段副本编号作为值列，
#而是提供原始副本编号，则iscopynumber的参数必须切换为TRUE。   
# compare FGA, FGG, and FGL
fga.coad <- compFGA(moic.res     = cmoic.coad_TLS,
                    segment      = segment,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.2, # threshold to determine CNA gain or loss
                    test.method  = "parametric", # statistical testing method
                    fig.name     = "5.BARPLOT OF FGA")
#> --2 samples mismatched from current subtypes.
#> 5% 10% 15% 21% 26% 31% 36% 41% 46% 51% 57% 62% 67% 72% 77% 82% 88% 93% 98%
head(fga.coad$summary)
summary <- fga.coad$summary
write.csv(summary,"5.Demo of comparison of fraction genome altered among 5 identified subtype of breast cancer in TCGA-COAD cohort.csv")
#6) compare drug sensitivity
drug.coad <- compDrugsen(moic.res    = cmoic.coad_TLS,
                         norm.expr   = fpkm[,cmoic.coad_TLS$clust.res$samID], # double guarantee sample order
                         drugs       = c(  "Erlotinib",        #####引号里面的药物不要有空格
                                           # "Rapamycin",
                                           # "Sunitinib",
                                           # "PHA-665752",
                                           #"MG-132",
                                           # "Paclitaxel",
                                           # "Cyclopamine",
                                           # "AZ628",
                                           # "Sorafenib",
                                           # "VX-680",
                                           "Imatinib",
                                           "Cisplatin",
                                           # "TAE684",
                                           "Crizotinib",
                                           # "Saracatinib",
                                           "S-Trityl-L-cysteine",
                                           #                 "Z-LLNle-CHO",
                                           #                 "Dasatinib",
                                           #                 "GNF-2",
                                           #                 "CGP-60474",
                                           #                 "CGP-082996",
                                           #                 "A-770041",
                                           #                 "WH-4-023",
                                           #                 "WZ-1-84",
                                           #                 "BI-2536",
                                           #                 "BMS-536924",
                                           #                 "BMS-509744",
                                           "CMK",
                                           #                 "Pyrimethamine",
                                           #                 "JW-7-52-1",
                                           #                 "A-443654",
                                           #                 "GW843682X",
                                           #                 "MS-275",
                                           #                 "Parthenolide",
                                           #                 "KIN001-135",
                                           #                 "TGX221",
                                           #                 "Bortezomib",
                                           #                 "XMD8-85",
                                           #                 "Roscovitine",
                                           #                 "Salubrinal",
                                           #                 "Lapatinib",
                                           #                 "GSK269962A",
                                           "Doxorubicin",
                                           #                 "Etoposide",
                                           "Gemcitabine",#
                                           # "MitomycinC",
                                           # "Vinorelbine",
                                           # "NSC-87877",
                                           # "Bicalutamide",
                                           # "QS11",
                                           # "CP466722",
                                           "Midostaurin",
                                           # "CHIR-99021",
                                           # "AP-24534",
                                           # "AZD6482",
                                           # "JNK-9L",
                                           # "PF-562271",
                                           # "HG-6-64-1",
                                           # "JQ1",
                                           # "JQ12",
                                           "DMOG",
                                           # "FTI-277",
                                           # "OSU-03012",
                                           "Shikonin",
                                           "Embelin",
                                           # "FH535",
                                           # "PAC-1",
                                           # "IPA-3",
                                           # "GSK-650394",
                                           # "5-Fluorouracil",
                                           "Thapsigargin",
                                           # "BMS-754807",
                                           # "Lisitinib",
                                           "Bexarotene"
                                           # "Bleomycin",
                                           # "LFM-A13",
                                           # "GW-2580",
                                           # "AUY922",
                                           # "Phenformin",
                                           # "Pazopanib",
                                           # "LAQ824",
                                           # "GSK1904529A",
                                           # "BMS345541",
                                           # "Tipifarnib",
                                           # "BMS-708163",
                                           # "Ruxolitinib",
                                           # "AS601245",
                                           # "TL-2-105",
                                           # "AT-7519",
                                           # "TAK-715",
                                           # "BX-912",
                                           # "ZSTK474",
                                           # "AS605240",
                                           # "GSK1070916",
                                           # "KIN001-102",
                                           # "LY317615",
                                           # "GSK429286A",
                                           # "FMK",
                                           # "QL-XII-47",
                                           # "CAL-101",
                                           # "UNC0638",
                                           # "XL-184",
                                           # "WZ3105",
                                           # "XMD14-99",
                                           # "AC220",
                                           # "CP724714",
                                           # "JW-7-24-1",
                                           # "NPK76-II-72-1",
                                           # "STF-62247",
                                           # "NG-25",
                                           # "TL-1-85",
                                           # "VX-11e",
                                           #  "FR-180204",
                                           #  "Zibotentan",
                                           # "YM155",
                                           # "NSC-207895",
                                           # "AR-42",
                                           # "CUDC-101",
                                           # "Belinostat",
                                           # "I-BET-762",
                                           # "CAY10603",
                                           # "BIX02189",
                                           # "CH5424802",
                                           # "EKB-569",
                                           # "GSK2126458",
                                           # "KIN001-236",
                                           #  "KIN001-244",
                                           # "KIN001-055",
                                           # "KIN001-260",
                                           # "KIN001-266",
                                           # "Masitinib",
                                           # "MP470",
                                           # "MPS-1-IN-1",
                                           # "BHG712",
                                           # "OSI-930",
                                           # "OSI-027",
                                           # "CX-5461",
                                           # "PHA-793887",
                                           #  "PI-103",
                                           # "PIK-93",
                                           # "SB52334",
                                           # "TPCA-1",
                                           # "TG101348",
                                           # "Foretinib",
                                           # "Y-39983",
                                           #  "YM201636",
                                           #  "Tivozanib",
                                           # "GSK690693",
                                           # "SNX-2112",
                                           #  "QL-XI-92",
                                           # "XMD13-2",
                                           # "QL-X-138",
                                           # "XMD15-27"
                         ), # a vector of names of drug in GDSC
                         clust.col = rev(c("#2EC4B6", "#E71D36")),
                         tissueType  = "digestive_system", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "6.BOXVIOLIN OF ESTIMATED IC50") 

drug.coad <- compDrugsen(moic.res    = cmoic.coad_TLS,
                         norm.expr   = fpkm[,cmoic.coad_TLS$clust.res$samID], # double guarantee sample order
                         drugs       = c(  "5-Fluorouracil",
                                           "AC220",
                                           "BMS-536924",
                                           "BMS-708163",
                                           #"BMS-754806",
                                           "CH5424802",
                                           "Cyclopamine",
                                           "EKB-569",
                                           "GSK690693",
                                           "GW843682X",
                                           "LY317615",
                                           "Lapatinib",
                                           "Lisitinib",
                                           "MG-132",
                                           "MPS-1-IN-1",
                                           "Masitinib",
                                           "Phenformin",
                                           "Roscovitine",
                                           "Sorafenib",
                                           "TGX221",
                                           "VX-11e",
                                           "XMD8-85"
                         ), # a vector of names of drug in GDSC
                         clust.col = rev(c("#2EC4B6", "#E71D36")),
                         tissueType  = "digestive_system", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "6.BOXVIOLIN OF ESTIMATED IC50_high") 

## FIG8A. =============
drug.df <- do.call(rbind, lapply(names(drug.coad), function(nm) {
  cbind(Drug = nm, drug.coad[[nm]])
}))
drug.df %>% 
  dplyr::filter(drug.df$Drug %in% c("Imatinib", "CMK", "Crizotinib","Embelin","Cisplatin","Doxorubicin")) %>% 
  grouped_ggbetweenstats(., 
                         x = "Subtype", 
                         y = "Est.IC50",
                         grouping.var    = Drug,
                         palette = "Set1", #rev(c("#2EC4B6", "#E71D36"))
                         plot.type = "box",
                         type = "nonparametric", ## type of test
                         bf.message = FALSE,
                         pairwise.comparisons = TRUE,
                         pairwise.display = "all",
                         outlier.tagging      = F,
                         p.adjust.method = "BH",
                         effsize.type = "d",
                         var.equal = TRUE,
                         tr = 0.2,
                         k = 3, # 小数点位数
                         centrality.plotting = T,
                         centrality.label.args = list(size = 3),
                         plotgrid.args = list(ncol = 3),
                         ggsignif.args = list(test.args = c(alternative = "two.sided")),
                         ylab = "",
                         xlab = "",
                         ggtheme = theme_bw(base_size = 18) +
                           theme(axis.text = element_text(face = "plain", size = 16),
                                 axis.text.y = element_text(face = "plain", size = 16),
                                 axis.title = element_text(face = "plain", size = 18),
                                 plot.title = element_text(size = 14),
                                 plot.subtitle = element_text(size = 6.5),
                                 plot.caption = element_text(size = 12)))

ggsave(
  filename = "药物_IC50_COAD-BoxPlot.png",
  # plot = p2,
  path = "./",
  width = 12,
  height = 10,
  device = "png",
  dpi = 300
)
ggsave(
  filename = "药物_IC50-BoxPlot.pdf",
  # plot = p2,
  path = "./",
  width = 15,
  height = 10,
  device = "pdf"
)

########将药物敏感性的结果变为dataframe,其中行名中有药物信息
lala <- do.call(rbind, lapply(drug.coad, as.data.frame))
save(drug.coad,file = "drug.coad.Rdata")
write.csv(drugData2016,"drugData2016.csv")

######将行名中的药物提取出来
lala2 <- lala %>% 
  rownames_to_column(var = "rowname") %>% 
  separate(rowname, into = c("drug", "sample"), sep = "\\.")
#####先变为宽格式
lala_wider <- tidyr::pivot_wider(lala2, 
                                 names_from=drug,
                                 values_from =Est.IC50)
#####加入index信息
rownames(lala_wider)<-lala_wider$sample
lala_wider    <- lala_wider[rownames(ydata.risk), ]
lala_wider$index<-ydata.risk$index
#####变为长格式
colnames(lala_wider)

lala_wider_long <- 
  tidyr::pivot_longer(lala_wider,
                      cols = c(
                        `5-Fluorouracil`,
                        AC220,
                        `BMS-536924`,
                        `BMS-708163`,
                        #`BMS-754806`,
                        CH5424802,
                        Cyclopamine,
                        `EKB-569`,
                        GSK690693,
                        GW843682X,
                        LY317615,
                        Lapatinib,
                        Lisitinib,
                        `MG-132`,
                        `MPS-1-IN-1`,
                        Masitinib,
                        Phenformin,
                        Roscovitine,
                        Sorafenib,
                        TGX221,
                        `VX-11e`,
                        `XMD8-85` ), names_to = "drug", values_to = "value")
########做散点图
cell <- "drug"
pdf_file <- paste(cell ,"index_药物敏感性_COAD-pointPlot_high.pdf", sep="_")
pdf(pdf_file,width=5,height= 10)

ggplot(lala_wider_long , aes(x=log2(index),y=value))+
  geom_point(size=0.8,colour ="#588b8b" )+
  geom_smooth(method = "lm",se = T,colour="#588b8b",fill ="#588b8b")+
  # stat_poly_eq(aes(label = paste(after_stat(eq.label),
  #                                after_stat(p.value.label),
  #                                sep = "~~~")),
  #              formula = y~x,
  #              parse = TRUE, 
  #              size=2,colour ="black")+
  stat_correlation(mapping = use_label("R"), 
                   method = "pearson",
                   label.x="left",
                   # npcy = 0.8,
                   # size=2.5,
                   colour ="#588b8b")+
  labs(x=" ",y=NULL)+
  facet_wrap(.~drug,scales="free",ncol=3)+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
        panel.grid =element_blank(),
        axis.text = element_text(colour="black"),legend.key = element_blank(),
        legend.title=element_blank())+
  labs(fill=' ')
dev.off()



#==================== 表型cox=====================

head(ydata.risk.number_forma,1)
#####转换为长数据
ydata.risk.number_forma_long <- 
  tidyr::pivot_longer(ydata.risk.number_forma,
                      cols = c(
                        age,
                        SEX ,
                        pstage,
                        N_stage,
                        M_stage,
                        # tumor_size,
                        # Aneuploidy_Score, 
                        # MSI_MANTIS_Score,
                        MSIsensor_Score,
                        # Mutation_Count, 
                        # TMB_.nonsynonymous.,
                        Buffa_Hypoxia_Score, 
                        Ragnum_Hypoxia_Score, 
                        # Winter_Hypoxia_Score
                      ), names_to = "bx", values_to = "value")

##作图-如果有显著离群值，可删除后重新导入数据作图
#导入数据
formula.trans <- y ~ I(x^2)
pdf_file <- paste("E:/pulis/202510/2023.9.19_胃肠道癌生信分析_MOVICS_TLS/" ,"模型分群后_index_表型_COAD-pointPlot.pdf", sep="_")
pdf(pdf_file,width=10,height= 5)

ggplot(ydata.risk.number_forma_long, aes(x=index,y=value))+
  geom_point(size=1,colour ="#219ebc" )+
  geom_smooth(method = "lm",se = T,colour="#219ebc",fill = "#219ebc")+
  # stat_poly_eq(aes(label = paste(after_stat(eq.label),
  #                                after_stat(p.value.label),
  #                                sep = "~~~")),
  #              formula = y~x,
  #              parse = TRUE,
  #              size=2,
  #              colour ="#219ebc")+
  # stat_poly_eq(aes(label =  paste(after_stat(rr.label),
  #                                 after_stat(n.label), sep = "*\", \"*")),
  #              formula = formula.trans) +
  stat_correlation(mapping = use_label("R", "P", "n"),
                   method = "pearson",
                   label.x="left",
                   npcy = 0.8,
                   size=3,
                   parse = TRUE,  # 禁用表达式解析
                   colour ="black")+
  # stat_correlation(aes(label = paste(after_stat(r.label),
  #                                    after_stat(p.value.label),
  #                                    after_stat(n.label),
  #                                    sep = "*\", \"*"))) +
  labs(x=" ",y=NULL)+
  facet_wrap(bx~.,scales="free",ncol = 4)+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
        panel.grid =element_blank(),
        axis.text = element_text(colour="black"),legend.key = element_blank(),
        legend.title=element_blank())+
  labs(fill=' ')
dev.off()


# 表型的cox单因素及与index的多因素回归================
##调整顺序
ydata.risk.number_forma <- ydata.risk.number_forma[rownames(ydata.risk), ]
##导入index列
ydata.risk.number_forma$index <- ydata.risk$index
#导出数据，将后1列移到最前面，并删除无关列（如全为0的列、行号列、risk、pvalue等列）再导入

#调整顺序
cols <- colnames(ydata.risk.number_forma)
new_cols <- c(cols[length(cols)-1],cols[length(cols)] ,cols[1:2],cols[36:44],cols[3:35])

# 然后将 dataframe 按照新的列名顺序排列
ydata.risk.number_forma2 <- ydata.risk.number_forma[, new_cols]


#write.csv(ydata.risk.number_forma, file='ydata.risk.number_forma.csv')
#ydata.risk_forma <-read.csv('ydata.risk.number_forma.csv', header = T, row.names = 1, sep=",")

#============= 表型cox单因素回归 ====================
##计算
pFilter=1
outResult=data.frame()
sigGenes=c("time","status")
for(i in colnames(ydata.risk.number_forma2[,c(2,5:18)])){
  mycox <- coxph(Surv(time, status) ~ ydata.risk.number_forma2[,i], data = ydata.risk.number_forma2)
  mycoxSummary = summary(mycox)
  pvalue=mycoxSummary$coefficients[,"Pr(>|z|)"]
  if(pvalue<pFilter){
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=mycoxSummary$conf.int[,"exp(coef)"],
                          L95CI=mycoxSummary$conf.int[,"lower .95"],
                          H95CI=mycoxSummary$conf.int[,"upper .95"],
                          pvalue=mycoxSummary$coefficients[,"Pr(>|z|)"])
    )
    
  }
}


##写出结果文件
write.csv(outResult, file='单因素COX_outResult.csv')

#############单因素森林图
## FIG4B. 单因素森林图=================
outResult<-read.csv("单因素_outResult.csv", header = T, sep = ",")
outResult %>%
  dplyr::mutate(p.adjust = signif(p.adjust(pvalue, method = "BH"),digits = 3),
                L95CI = signif(as.numeric(lower), digits = 3),
                H95CI = signif(as.numeric(upper), digits = 3),
                HR = signif(as.numeric(HR), digits = 3),
                pvalue = signif(as.numeric(pvalue), digits = 3),
                `HR(95% CI)` = paste0(
                  HR, " (",
                  L95CI, "-", H95CI, ")"
                )) %>% 
  write.csv("单因素_outResult_forForest.csv")
#森林图绘制：
photo <- outResult %>%
  dplyr::mutate(p.adjust = signif(p.adjust(pvalue, method = "BH"),digits = 3),
                L95CI = signif(as.numeric(lower), digits = 3),
                H95CI = signif(as.numeric(upper), digits = 3),
                HR = signif(as.numeric(HR), digits = 3),
                pvalue = signif(as.numeric(pvalue), digits = 3),
                `HR(95% CI)` = paste0(
                  HR, " (",
                  L95CI, "-", H95CI, ")"
                )) %>% 
  forestplot(labeltext = c(id,`HR(95% CI)`,p.adjust), 
             clip = c(0.5, 3),#森林图中置信区间上下限设置
             xlog = F,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1.2),     # 坐标刻度字体大小
               xlab  = gpar(cex = 1.3),     # 坐标轴标题字体
               label = gpar(cex = 1.0)      # 表格文字大小
             ),
             zero=2) %>%
  fp_set_style(box = "#588b8b", #权重方块颜色
               line = "#588b8b") %>%#是否对X轴取对数
  fp_set_zebra_style("#ffebe8") %>% #更改为斑马样式，并修改汇总菱形的填充颜色
  fp_add_header(id = c("","gene"),
                `HR(95% CI)` =c("","HR(95% CI)"),
                p.adjust = c("","p.adjust"))

pdf(file="表型_单因素_forestPlot.pdf",width = 6,height =2.7)
print(photo)
dev.off()

##################模型分群后做表型和index的cox多因素回归。
## FIG4C. 多因素森林图=================
head(ydata.risk.number_forma2,1)
###计算五因素回归-CNR1_index+P_stage+radiation_therapy+age——————————因为都不显著，就不往下做了。
#age pstage N_stage M_stage SEX tumor_size 
res.cox <- coxph(Surv(time, status) ~ index+age+pstage+ N_stage+ M_stage ,data = ydata.risk.number_forma2)
summary(res.cox)

#森林图绘制：
outResult <- read.csv("多因素_outResult.csv", header = T, sep = ",")
outResult %>%
  dplyr::mutate(p.adjust = signif(p.adjust(pvalue, method = "BH"),digits = 3),
                L95CI = signif(as.numeric(lower), digits = 3),
                H95CI = signif(as.numeric(upper), digits = 3),
                HR = signif(as.numeric(HR), digits = 3),
                pvalue = signif(as.numeric(pvalue), digits = 3),
                `HR(95% CI)` = paste0(
                  HR, " (",
                  L95CI, "-", H95CI, ")"
                )) %>% 
  write.csv("多因素_outResult_forForest.csv")

photo<-outResult %>%
  dplyr::mutate(p.adjust = signif(p.adjust(pvalue, method = "BH"),digits = 3),
                L95CI = signif(as.numeric(lower), digits = 3),
                H95CI = signif(as.numeric(upper), digits = 3),
                HR = signif(as.numeric(HR), digits = 3),
                pvalue = signif(as.numeric(pvalue), digits = 3),
                `HR(95% CI)` = paste0(
                  HR, " (",
                  L95CI, "-", H95CI, ")"
                )) %>% 
  forestplot(labeltext = c(id,`HR(95% CI)`,p.adjust), 
             clip = c(0.5, 3),#森林图中置信区间上下限设置
             xlog = F,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1.2),     # 坐标刻度字体大小
               xlab  = gpar(cex = 1.3),     # 坐标轴标题字体
               label = gpar(cex = 1.0)      # 表格文字大小
             ),
             zero=2) %>%
  fp_set_style(box = "#588b8b", #权重方块颜色
               line = "#588b8b") %>%#是否对X轴取对数
  fp_set_zebra_style("#ffebe8") %>% #更改为斑马样式，并修改汇总菱形的填充颜色
  fp_add_header(id = c("","gene"),
                `HR(95% CI)` =c("","HR(95% CI)"),
                p.adjust = c("","p.adjust"))
pdf(file="多因素_forestPlot.pdf",width = 6,height =2.7)
print(photo)
dev.off()




#7）与其他亚型比较一致性
# customize the factor level for pstage
surv.info$pstage <- factor(surv.info$pstage, levels = c("T1","T2","T3","T4","T4A","T4B"))
surv.info$SEX <- factor(surv.info$SEX, levels = c("Male","Female"))
surv.info$M_STAGE <- factor(surv.info$M_STAGE, levels = c("MX","M0","M1","M1X","M1B"))
surv.info$N_STAGE <- factor(surv.info$N_STAGE, levels = c("NX","N0","N1","N1A","N1B","N1C","N2","N2A","N2B"))

# agreement comparison (support up to 6 classifications include current subtype)
agree.coad <- compAgree(moic.res  = cmoic.coad_TLS,
                        subt2comp = surv.info[,c("SEX","pstage","M_STAGE","N_STAGE")],
                        doPlot    = TRUE,
                        box.width = 0.2,
                        fig.name  = "7.AGREEMENT OF CONSENSUSMOIC WITH PAM50 AND PSTAGE")
#> --all samples matched.

print(agree.coad)
write.csv(agree.coad,"7.Agreement of 5 identified subtypes with PAM50 classification and pathological stage in TCGA-BRCA cohort.csv")

# ============= 四、 RUN Module===============
#1） 运行差分表达式分析

# run DEA with edgeR
runDEA(dea.method = "edger",      ##三选一：edger和deseq2，以及limma。
       expr       = count, # raw count data
       moic.res   = cmoic.coad_TLS,


# run DEA with DESeq2
runDEA(dea.method = "deseq2",
       expr       = count,
       moic.res   = cmoic.coad_TLS,
       prefix     = "TCGA-COAD")


# run DEA with limma
runDEA(dea.method = "limma",
       expr       = fpkm, # normalized expression data
       moic.res   =cmoic.coad_TLS,
       prefix     = "TCGA-COAD")
#> --all samples matched.
#> --you choose limma and please make sure a microarray profile or a normalized expression data [FPKM or TPM without log2 transformation is recommended] was provided.
#> --log2 transformation done for expression data.
#> limma of CS1_vs_Others done...
#> limma of CS2_vs_Others done...

#每个已识别的癌症亚型将与其他（其他）亚型进行比较，并根据res.path的参数存储相应的.txt文件。
#默认情况下，这些文件将保存在当前工作目录下。

#2) run biomarker identification procedure

# choose edgeR result to identify subtype-specific up-regulated biomarkers
marker.up <- runMarker(moic.res      = cmoic.coad,
                       dea.method    = "limma",   ##三选一：edger和deseq2，以及limma。
                       prefix        = "TCGA-COAD", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = fpkm, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "2.UPREGULATED BIOMARKER HEATMAP_limma")
#> --all samples matched.
#> --log2 transformation done for expression data.
# check the upregulated biomarkers
head(marker.up$templates)
templates <- marker.up$templates
write.csv(templates,"2.Demo of subtype-specific upregulated biomarkers for 5 identified subtypes of breast cancer in TCGA-COAD cohort.csv")


#然后尝试limma的结果来鉴定亚型特异性下调的生物标志物。请随时尝试DESeq2。
# choose limma result to identify subtype-specific down-regulated biomarkers
marker.dn <- runMarker(moic.res      = cmoic.coad,
                       dea.method    = "limma",        ##三选一：edger和deseq2，以及limma。
                       prefix        = "TCGA-COAD",
                       dirct         = "down",
                       n.marker      = 100, # switch to 50
                       doplot        = TRUE,
                       norm.expr     = fpkm,
                       annCol        = annCol,
                       annColors     = annColors,
                       fig.name      = "2.DOWNREGULATED BIOMARKER HEATMAP_limma")
#> --all samples matched.
#> --log2 transformation done for expression data.
head(marker.dn$templates)
templates <- marker.dn$templates

write.csv(templates,"2.Demo of subtype-specific downregulated biomarkers for 5 identified subtypes of breast cancer in TCGA-COAD cohort.csv")

#3) run gene set enrichment analysis
#类似地，基于其相应的DEA结果对每个亚型运行GSEA，以识别亚型特异性功能途径。
#为此，我准备了一个基因集背景，该背景包括来自分子特征数据库（MSigDB，https://www.gsea-msigdb.org/gsea/msigdb/index.jsp).你可以下载其他感兴趣的背景资料供自己学习。

# MUST locate ABSOLUTE path of msigdb file
MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)
#"/Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library/MOVICS/extdata/c5.bp.v7.1.symbols.xls"

#同样，这些确定的特定途径应通过显著性阈值（例如，标称P值<0.05和调整后P值<0.025），
#并且不得与其他亚型的任何途径重叠。在具有亚型特异性途径后，检索途径内的基因，
#通过使用GSVA R包计算单个样本富集分数。随后，亚型特异性富集分数将由亚型内的平均值或中值
#（默认为平均值）表示，并将通过对角线热图进一步可视化。

# run GSEA to identify up-regulated GO pathways using results from edgeR
gsea.up <- runGSEA(moic.res     = cmoic.coad_TLS,
                   dea.method   = "limma",   ##三选一：edger和deseq2，以及limma。
                   prefix       = "TCGA-COAD", # MUST be the same of argument in runDEA()
                   dat.path     = "2023.9.19_胃肠道癌生信分析_MOVICS_TLS/4.其他分析/1.差分表达式分析/", # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = "E:/pulis/202510/2023.9.19_胃肠道癌生信分析_MOVICS_TLS/2023.9.19_胃肠道癌生信分析_MOVICS_TLS/4.其他分析/4.ssGSEA分析/c2.cp.kegg.v2022.1.Hs.symbols.gmt", # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = fpkm, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "3.UPREGULATED PATHWAY HEATMAP_limma")
#> --all samples matched.
#> GSEA done...
#> --log2 transformation done for expression data.
#> Estimating GSVA scores for 50 gene sets.
#> Estimating ECDFs with Gaussian kernels


#   |===================================================================== |  98%
# |                                                                            
#   |======================================================================| 100%
#> gsva done...
#> heatmap done...

#Check some columns of GSEA results for the first subtype (CS1).#
print(gsea.up$gsea.list$CS1[1:6,3:6])

#Also check results of subtype-specific enrichment scores.

head(round(gsea.up$grouped.es,3))

#然后尝试来源于DESeq2的结果来鉴定亚型特异性下调途径。请随意尝试limma。

# run GSEA to identify down-regulated GO pathways using results from DESeq2
gsea.dn <- runGSEA(moic.res     = cmoic.coad_TLS,
                   dea.method   = "limma",      ##三选一：edger和deseq2，以及limma。
                   prefix       = "TCGA-COAD",
                   msigdb.path  = "E:/pulis/202510/2023.9.19_胃肠道癌生信分析_MOVICS_TLS/2023.9.19_胃肠道癌生信分析_MOVICS_TLS/4.其他分析/4.ssGSEA分析/c2.cp.kegg.v2022.1.Hs.symbols.gmt",
                   norm.expr    = fpkm,
                   dirct        = "down",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method  = "gsva", # switch to ssgsea/gsva
                   norm.method  = "mean", # switch to median/mean
                   fig.name     = "3.DOWNREGULATED PATHWAY HEATMAP_limma") 

# =========  4)run gene set variation analysis ================
#对于所有新定义的分子亚型，描述其通过基因集的不同特征验证的特征是至关重要的。
#MOVICS提供了一个简单的函数，该函数使用基因集变异分析来基于给定的感兴趣基因集列表计算
#每个亚型中每个样本的富集分数。首先，我们必须准备一份感兴趣的基因列表，保存为GMT格式。

# MUST locate ABSOLUTE path of gene set file
GSET.FILE <- 
  system.file("extdata", "gene sets of interest.gmt", package = "MOVICS", mustWork = TRUE)
GSET.FILE <- 
  system.file("extdata", "AUTOPHAGY_geneset.gmt", package = "MOVICS", mustWork = TRUE)
GSET.FILE <- 
  system.file("extdata", "AUTOPHAGY_geneset_select.gmt", package = "MOVICS", mustWork = TRUE)




GSET.FILE <- 
  system.file("extdata", "c2.cp.kegg.v2022.1.Hs.symbols.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "h.all.v2022.1.Hs.symbols.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "h.all.v2022.1.Hs.symbols_select.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "REACTOME_metabolism.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "28_immune_cell_GeneSignatures.gmt", package = "MOVICS", mustWork = TRUE)
GSET.FILE <- 
  system.file("extdata", "CMS1.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "CMS3.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "CMS4.gmt", package = "MOVICS", mustWork = TRUE)


GSET.FILE <- 
  system.file("extdata", "TF.gmt", package = "MOVICS", mustWork = TRUE)
GSET.FILE <- 
  system.file("extdata", "TF_wx.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "immune_gene.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "chromatin.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "signature_metabolism.gmt", package = "MOVICS", mustWork = TRUE)

GSET.FILE <- 
  system.file("extdata", "signature_metabolism_select.gmt", package = "MOVICS", mustWork = TRUE)

#然后我们可以根据指定的方法和给定的感兴趣的基因集来计算单个样本的富集分数。
# run GSVA to estimate single sample enrichment score based on given gene set of interest
gsva.res <- 
  runGSVA(moic.res      = cmoic.coad_TLS,
          norm.expr     = fpkm,
          gset.gmt.path = "E:/pulis/202510/2023.9.19_胃肠道癌生信分析_MOVICS_TLS/2023.9.19_胃肠道癌生信分析_MOVICS_TLS/4.其他分析/4.ssGSEA分析/c2.cp.kegg.v2022.1.Hs.symbols.gmt", # ABSOLUTE path of gene set file
          gsva.method   = "ssgsea", # method to calculate single sample enrichment score
          annCol        = annCol,
          annColors     = annColors,
          fig.path      = getwd(),
          fig.name      = "4.GENE SETS OF INTEREST HEATMAP_hallmarker_TLS",
          height        = 7,
          width         = 8)

# check raw enrichment score
print(gsva.res$raw.es[1:3,1:3])

signature_metabolism114 <- gsva.res$raw.es
write.csv(signature_metabolism114 ,"4.GENE SETS OF INTEREST HEATMAP_hallmarker_TLS.csv")


#———免疫基因box图———————————————————————————————————————————————————————————————————————————————————————————————————  
##————————————————————作高低风险——免疫gene的差异box图————————————————————————————————————

###########facet图
#PDCD1, CD274,CTLA4, BCL2, HLA-E, GZMA, GZMB, GZMK, GZMM, GZMH, CD83, CD8A, TGFB1, TGFB2, TGFB3, TNF, IL6, IL27, IL27RA
#免疫检查点：PDCD1, CD274,CTLA4, BCL2,IL6, IL27, IL27RA,
#其他：GZMA, GZMB, GZMK, GZMM, GZMH, CD83, CD8A, TGFB1, TGFB2, TGFB3, TNF
#HLA:HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB, HLA-DPA1, HLA-DPB1, HLA-DPB2, HLA-DQA1, HLA-DQB1, HLA-DQB2, HLA-DRA, HLA-DRB1, HLA-DRB5, HLA-DRB6, HLA-E
# #`HLA-DMA`,
# `HLA-DMB`,
# `HLA-DOA`,
# `HLA-DOB`,
# `HLA-DPA1`,
# `HLA-DPB1`,
# `HLA-DPB2`,
# `HLA-DQA1`,
# `HLA-DQB1`,
# `HLA-DQB2`,
# `HLA-DRA`,
# `HLA-DRB1`,
# `HLA-DRB5`,
# `HLA-DRB6`,
# `HLA-E`,

# 免疫基因 =============================================
## FIG5C. 免疫基因-box图==================
##筛选出risk及免疫相关基因的列
xdata.NoneScale.risk.select <- xdata.NoneScale.risk %>% 
  as.data.frame() %>% 
  dplyr::select(CD8A, GZMK, GZMM, GZMH, 
                TGFB1, TGFB2, TGFB3, TNF,  #CD83,GZMA, GZMB,
                risk) %>%
  as.data.frame()

ydata.risk<-ydata.risk%>%
  as.data.frame()

#加入index列
xdata.NoneScale.risk.select$index <-ydata.risk$index
# 转换数据框为长格式
xdata.NoneScale.risk.select_long2 <- tidyr::pivot_longer(xdata.NoneScale.risk.select, 
                                                         cols = c(CD8A, GZMK, GZMM, GZMH, 
                                                                  TGFB1, TGFB2, TGFB3, TNF,  #CD83,GZMA, GZMB,
                                                         ), names_to = "gene", values_to = "value")

p1 <- grouped_ggbetweenstats(xdata.NoneScale.risk.select_long2, 
                             x = "risk", 
                             y = "value",
                             grouping.var    = gene,
                             palette = "Set1", #rev(c("#2EC4B6", "#E71D36"))
                             plot.type = "box",
                             type = "parametric", ## type of test
                             bf.message = FALSE,
                             pairwise.comparisons = TRUE,
                             pairwise.display = "all",
                             outlier.tagging      = F,
                             p.adjust.method = "BH",
                             effsize.type = "d",
                             var.equal = TRUE,
                             tr = 0.2,
                             k = 3, # 小数点位数
                             centrality.plotting = T,
                             centrality.label.args = list(size = 3),
                             plotgrid.args = list(ncol = 4),
                             ggsignif.args = list(test.args = c(alternative = "two.sided")),
                             ylab = "",
                             xlab = "",
                             ggtheme = theme_bw(base_size = 18) +
                               theme(axis.text = element_text(face = "plain", size = 16),
                                     axis.text.y = element_text(face = "plain", size = 16),
                                     axis.title = element_text(face = "plain", size = 18),
                                     plot.title = element_text(size = 14),
                                     plot.subtitle = element_text(size = 6.5),
                                     plot.caption = element_text(size = 12)))
ggsave(
  filename = "免疫基因_COAD-BoxPlot2.png",
  plot = p1,
  path = "./",
  width = 15,
  height = 10,
  device = "png",
  dpi = 300
)
ggsave(
  filename = "免疫基因_COAD-BoxPlot2.pdf",
  plot = p1,
  path = "./",
  width = 15,
  height = 10,
  device = "pdf"
)

# 免疫基因-box图
pdf_file <- paste("免疫基因_COAD-BoxPlot.pdf")
dev.off()
pdf(pdf_file,width=9,height= 4)
ggplot(xdata.NoneScale.risk.select_long2, aes(x =factor(risk), y=log2(value) ,fill =factor(risk))) +
  geom_boxplot(size=0.5,varwidth = F, alpha=0.8) + 
  geom_jitter(position=position_jitter(width=0.15,height=0.2), size=0.8,alpha=0.3,colour="#bfdbf7") + 
  facet_wrap(gene~.,scales="free",ncol = 6)+
  scale_fill_manual(values=rev(c("#2EC4B6", "#E71D36"))) +
  theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
        axis.text.x = element_text(angle = 30,vjust = 0.9,hjust = 1), panel.grid =element_blank(),axis.text = element_text(colour="black"),
        strip.background = element_rect(colour="black", fill="white"))+
  geom_signif(
    comparisons = list( c("Low-risk", "High-risk")),
    map_signif_level = TRUE, textsize = 6,step_increase = 0.1
  ) 
dev.off()  

#免疫检查点基因==================
##  FIG5D. 免疫检查点基因-box图==================
##筛选出risk及免疫相关基因的列
xdata.NoneScale.risk.select <- xdata.NoneScale.risk%>% 
  as.data.frame() %>% 
  dplyr::select(PDCD1, CD274,CTLA4, IL27,# IL6,  IL27RA,BCL2,
                risk) %>%
  as.data.frame() 
ydata.risk<-ydata.risk%>%
  as.data.frame() 
#加入index列
xdata.NoneScale.risk.select$index <-ydata.risk$index
# 转换数据框为长格式
xdata.NoneScale.risk.select_long2 <- tidyr::pivot_longer(xdata.NoneScale.risk.select, 
                                                         cols = c(PDCD1, CD274,CTLA4, IL27,
                                                                  # IL6,  IL27RA,BCL2,
                                                         ), names_to = "gene", values_to = "value")
p2 <- grouped_ggbetweenstats(xdata.NoneScale.risk.select_long2, 
                             x = "risk", 
                             y = "value",
                             grouping.var    = gene,
                             palette = "Set1", #rev(c("#2EC4B6", "#E71D36"))
                             plot.type = "box",
                             type = "parametric", ## type of test
                             bf.message = FALSE,
                             pairwise.comparisons = TRUE,
                             pairwise.display = "all",
                             outlier.tagging      = F,
                             p.adjust.method = "BH",
                             effsize.type = "d",
                             var.equal = TRUE,
                             tr = 0.2,
                             k = 3, # 小数点位数
                             centrality.plotting = T,
                             centrality.label.args = list(size = 3),
                             plotgrid.args = list(ncol = 4),
                             ggsignif.args = list(test.args = c(alternative = "two.sided")),
                             ylab = "",
                             xlab = "",
                             ggtheme = theme_bw(base_size = 18) +
                               theme(axis.text = element_text(face = "plain", size = 16),
                                     axis.text.y = element_text(face = "plain", size = 16),
                                     axis.title = element_text(face = "plain", size = 18),
                                     plot.title = element_text(size = 14),
                                     plot.subtitle = element_text(size = 6.5),
                                     plot.caption = element_text(size = 12)))
ggsave(
  filename = "免疫检查点基因_COAD-BoxPlot.png",
  plot = p2,
  path = "./",
  width = 15,
  height = 5,
  device = "png",
  dpi = 300
)
ggsave(
  filename = "免疫检查点基因_COAD-BoxPlot.pdf",
  plot = p2,
  path = "./",
  width = 15,
  height = 5,
  device = "pdf"
)

###########免疫基因-box图 ==========
pdf_file <- paste("免疫检查点基因_COAD-BoxPlot.pdf")

pdf(pdf_file,width=9,height= 4)
ggplot(xdata.NoneScale.risk.select_long2, aes(x =factor(risk), y=log2(value) ,fill =factor(risk))) +
  geom_boxplot(size=0.5,varwidth = F, alpha=0.8) + 
  geom_jitter(position=position_jitter(width=0.15,height=0.2), size=0.8,alpha=0.3,colour="#bfdbf7") + 
  facet_wrap(gene~.,scales="free",ncol = 6)+
  scale_fill_manual(values=c("#588b8b","#f5cac3")) +
  theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
        axis.text.x = element_text(angle = 30,vjust = 0.9,hjust = 1), panel.grid =element_blank(),axis.text = element_text(colour="black"),
        strip.background = element_rect(colour="black", fill="white"))+
  geom_signif(
    comparisons = list( c("Low-risk", "High-risk")),
    map_signif_level = TRUE, textsize = 6,step_increase = 0.1
  ) 
dev.off()  


# HLA基因============
## FIG5E. HLA基因-box图==================
##筛选出risk及免疫相关基因的列
xdata.NoneScale.risk.select <- xdata.NoneScale.risk %>% 
  as.data.frame() %>% 
  dplyr::select(`HLA-DMA`,
                `HLA-DMB`,
                `HLA-DOA`,
                `HLA-DOB`,
                # `HLA-DPA1`,
                # `HLA-DPB1`,
                # `HLA-DPB2`,
                # `HLA-DQA1`,
                # `HLA-DQB1`,
                # `HLA-DQB2`,
                # `HLA-DRA`,
                # `HLA-DRB1`,
                # `HLA-DRB5`,
                # `HLA-E`,
                risk) %>%
  as.data.frame() 
ydata.risk<-ydata.risk%>%
  as.data.frame() 
#加入index列
xdata.NoneScale.risk.select$index <-ydata.risk$index
# 转换数据框为长格式
xdata.NoneScale.risk.select_long2 <- tidyr::pivot_longer(xdata.NoneScale.risk.select, 
                                                         cols = c(`HLA-DMA`,
                                                                  `HLA-DMB`,
                                                                  `HLA-DOA`,
                                                                  `HLA-DOB`,
                                                                  # `HLA-DPA1`,
                                                                  # `HLA-DPB1`,
                                                                  # `HLA-DPB2`,
                                                                  # `HLA-DQA1`,
                                                                  # `HLA-DQB1`,
                                                                  # `HLA-DQB2`,
                                                                  # `HLA-DRA`,
                                                                  # `HLA-DRB1`,
                                                                  # `HLA-DRB5`,
                                                                  # `HLA-E`
                                                         ), names_to = "gene", values_to = "value")

p3 <- grouped_ggbetweenstats(xdata.NoneScale.risk.select_long2, 
                             x = "risk", 
                             y = "value",
                             grouping.var    = gene,
                             palette = "Set1", #rev(c("#2EC4B6", "#E71D36"))
                             plot.type = "box",
                             type = "parametric", ## type of test
                             bf.message = FALSE,
                             pairwise.comparisons = TRUE,
                             pairwise.display = "all",
                             outlier.tagging      = F,
                             p.adjust.method = "BH",
                             effsize.type = "d",
                             var.equal = TRUE,
                             tr = 0.2,
                             k = 3, # 小数点位数
                             centrality.plotting = T,
                             centrality.label.args = list(size = 3),
                             plotgrid.args = list(ncol = 4),
                             ggsignif.args = list(test.args = c(alternative = "two.sided")),
                             ylab = "",
                             xlab = "",
                             ggtheme = theme_bw(base_size = 18) +
                               theme(axis.text = element_text(face = "plain", size = 16),
                                     axis.text.y = element_text(face = "plain", size = 16),
                                     axis.title = element_text(face = "plain", size = 18),
                                     plot.title = element_text(size = 14),
                                     plot.subtitle = element_text(size = 6.5),
                                     plot.caption = element_text(size = 12)))
ggsave(
  filename = "HLA基因_COAD-BoxPlot.png",
  plot = p3,
  path = "./",
  width = 15,
  height = 5,
  device = "png",
  dpi = 300
)
ggsave(
  filename = "HLA基因_COAD-BoxPlot.pdf",
  plot = p3,
  path = "./",
  width = 15,
  height = 5,
  device = "pdf"
)

###########免疫基因-box图
pdf_file <- paste("HLA基因_COAD-BoxPlot.pdf")

pdf(pdf_file,width=9,height= 5.9)
ggplot(xdata.NoneScale.risk.select_long2, aes(x =factor(risk), y=log2(value) ,fill =factor(risk))) +
  geom_boxplot(size=0.5,varwidth = F, alpha=0.8) + 
  geom_jitter(position=position_jitter(width=0.15,height=0.2), size=0.8,alpha=0.3,colour="#bfdbf7") + 
  facet_wrap(gene~.,scales="free",ncol = 6)+
  scale_fill_manual(values=c("#588b8b","#f5cac3")) +
  theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
        axis.text.x = element_text(angle = 30,vjust = 0.9,hjust = 1), panel.grid =element_blank(),axis.text = element_text(colour="black"),
        strip.background = element_rect(colour="black", fill="white"))+
  geom_signif(
    comparisons = list( c("Low-risk", "High-risk")),
    map_signif_level = TRUE, textsize = 6,step_increase = 0.1
  ) 
dev.off()  
#==========================免疫浸润堆积柱状图 ==================
##  FIG5A  ================
# 导入数据并调整样本顺序
df   <- read.table("immune_cell.txt", header = T, sep = "\t")
row.names(df) <- df$X
df<-df[-1,]
colnames(df) <- gsub("\\.", "-", colnames(df))

df<-df%>%t()%>%as.data.frame()

# Get matches 
df <- df [rownames(df) %in% 
            rownames(ydata.risk),]
df     <- df [rownames(ydata.risk), ]

df2 <- df
df2 <- df2[,-1]
df3 <- df2
df3 = as.data.frame(lapply(df3,as.numeric))
rownames(df3)=rownames(df2)
df2 = as.data.frame(lapply(df2,as.numeric))
rownames(df2)=rownames(df3)

for (i in colnames(df2)) {
  #i <- colnames(ssgsea)[1]
  df3[,i] <- df2[,i] /apply(df2,1,sum)
}

df3$risk <- ydata.risk$risk
df3$index<-ydata.risk$index
df3$X <- rownames(df3)
colnames(df3)  

# 转换数据框为长格式
df_long <- tidyr::pivot_longer(df3, cols = c(Central_memory_CD8_T_cell       ,Effector_memeory_CD8_T_cell,    
                                             Activated_CD4_T_cell, Central_memory_CD4_T_cell,      
                                             Effector_memeory_CD4_T_cell,     T_follicular_helper_cell,       
                                             Gamma_delta_T_cell,   Type_1_T_helper_cell,           
                                             Type_17_T_helper_cell,Type_2_T_helper_cell,           
                                             Regulatory_T_cell,    Activated_B_cell,               
                                             Immature__B_cell,     Memory_B_cell,                  
                                             Natural_killer_cell,  CD56bright_natural_killer_cell, 
                                             CD56dim_natural_killer_cell,     Myeloid_derived_suppressor_cell,
                                             Natural_killer_T_cell,Activated_dendritic_cell,       
                                             Plasmacytoid_dendritic_cell,     Immature_dendritic_cell,        
                                             Macrophage,           Eosinophil,                     
                                             Mast_cell,            Monocyte,                       
                                             Neutrophil), names_to = "variable", values_to = "value")


df_long.low <- df_long %>% filter(risk == "Low-risk" )   #数据分成两部分分别展示
df_long.high <- df_long %>% filter(risk == "High-risk" )

p0 <- grouped_ggbetweenstats(df_long, 
                             x = "risk", 
                             y = "value",
                             grouping.var    = variable,
                             palette = "Set1", #rev(c("#2EC4B6", "#E71D36"))
                             plot.type = "box",
                             type = "nonparametric", ## type of test
                             bf.message = FALSE,
                             pairwise.comparisons = TRUE,
                             pairwise.display = "all",
                             outlier.tagging      = F,
                             p.adjust.method = "BH",
                             effsize.type = "d",
                             var.equal = TRUE,
                             tr = 0.2,
                             k = 3, # 小数点位数
                             centrality.plotting = T,
                             centrality.label.args = list(size = 3),
                             plotgrid.args = list(ncol = 4),
                             ggsignif.args = list(test.args = c(alternative = "two.sided")),
                             ylab = "",
                             xlab = "",
                             ggtheme = theme_bw(base_size = 16) +
                               theme(axis.text = element_text(face = "plain", size = 14),
                                     axis.text.y = element_text(face = "plain", size = 14),
                                     axis.title = element_text(face = "plain", size = 18),
                                     plot.title = element_text(size = 14),
                                     plot.subtitle = element_text(size = 6.5),
                                     plot.caption = element_text(size = 12)))
ggsave(
  filename = "免疫细胞丰度-BoxPlot.png",
  plot = p0,
  path = "./",
  width = 17,
  height = 25,
  device = "png",
  dpi = 300
)
ggsave(
  filename = "免疫细胞丰度_COAD-BoxPlot.pdf",
  plot = p0,
  path = "./",
  width = 17,
  height = 25,
  device = "pdf"
)

# 绘制堆积柱状图
pdf("low_免疫浸润堆积图.pdf",width=30,height=5)
ggplot(df_long.low, aes(x = X, y = value, fill = variable)) +
  geom_bar(stat="identity",position="stack", width=0.7,linewidth=0.25)+
  geom_col(position = "stack", color = "white") +
  labs(title = "Stacked Bar Chart Example", x = "", y = "Value") +
  scale_y_continuous(expand = c(0, 0),breaks=seq(0,100,25),labels=c('0','25%','50%','75%','100%'))+
  scale_fill_manual(values = c(brewer.pal(6,"Blues"),brewer.pal(8,"Set2"),brewer.pal(7,"Set3"),"#31572c","#226f54","#4f772d","#90a955","#a40404","#5c0000","#43291f"))+
  theme(axis.text.x = element_blank()) +
  xlab("")
dev.off()

pdf("high_免疫浸润堆积图.pdf",width=30,height=5)
ggplot(df_long.high, aes(x = X, y = value, fill = variable)) +
  geom_bar(stat="identity",position="stack", width=0.7,linewidth=0.25)+
  geom_col(position = "stack", color = "white") +
  labs(title = "Stacked Bar Chart Example", x = "", y = "Value") +
  scale_y_continuous(expand = c(0, 0),breaks=seq(0,100,25),labels=c('0','25%','50%','75%','100%'))+
  scale_fill_manual(values = c(brewer.pal(6,"Blues"),brewer.pal(8,"Set2"),brewer.pal(7,"Set3"),"#31572c","#226f54","#4f772d","#90a955","#a40404","#5c0000","#43291f"))+
  theme(axis.text.x = element_blank()) +
  xlab("")
dev.off()


#免疫细胞与index相关性散点图————————————————————————————————————————————————————————————————————————————————————————  
cell <- "celllist"
pdf_file <- paste(cell ,"index_免疫浸润_COAD-pointPlot.pdf", sep="_")
pdf(pdf_file,width=8,height= 11)

ggplot(df_long , aes(x=index,y=value))+
  geom_point(size=0.8,colour ="#588b8b" )+
  geom_smooth(method = "lm",se = T,colour="#588b8b",fill ="#588b8b")+
  # stat_poly_eq(aes(label = paste(after_stat(eq.label),after_stat(p.value.label),sep = "~~~")),formula = y~x,parse = TRUE, size=2,colour ="black")+
  # stat_correlation(mapping = use_label("R", "P", "n"),
  #                  method = "pearson",
  #                  label.x="left",
  #                  npcy = 0.8,
  #                  size=3,
  #                  parse = T,  # 禁用表达式解析
  #                  colour ="black")+
  # stat_correlation(aes(label = paste(after_stat(r.label),
  #                                    after_stat(p.value.label),
  #                                    after_stat(n.label),
  #                                    sep = " *\",\"* ")),
  #                  parse = TRUE) +
  labs(x=" ",y=NULL)+
  facet_wrap(.~variable,scales="free",ncol=4)+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
        panel.grid =element_blank(),
        axis.text = element_text(colour="black"),legend.key = element_blank(),
        legend.title=element_blank())+
  labs(fill=' ')
dev.off()


##筛选出risk及免疫相关基因的列
xdata.NoneScale.risk.select <- xdata.NoneScale.risk%>% 
  as.data.frame() %>% 
  dplyr::select(CNR1,GZMA, GZMB, GZMK, GZMM, GZMH, CD83, CD8A, TGFB1, TGFB2, TGFB3, TNF, 
                risk) %>%
  as.data.frame() 
ydata.risk<-ydata.risk%>%
  as.data.frame() 
#加入index列
xdata.NoneScale.risk.select$index <-ydata.risk$index
# 转换数据框为长格式
xdata.NoneScale.risk.select_long2 <- tidyr::pivot_longer(xdata.NoneScale.risk.select, cols = c(GZMA, GZMB, GZMK, GZMM, GZMH, CD83, CD8A, TGFB1, TGFB2, TGFB3, TNF), names_to = "gene", values_to = "value")

###########免疫基因-box图
pdf_file <- paste("免疫基因_COAD-BoxPlot.pdf")

pdf(pdf_file,width=9,height= 4)
ggplot(xdata.NoneScale.risk.select_long2, aes(x =factor(risk), y=log2(value) ,fill =factor(risk))) +
  geom_boxplot(size=0.5,varwidth = F, alpha=0.8) + 
  geom_jitter(position=position_jitter(width=0.15,height=0.2), size=0.8,alpha=0.3,colour="#bfdbf7") + 
  facet_wrap(gene~.,scales="free",ncol = 6)+
  scale_fill_manual(values=c("#588b8b","#f5cac3")) +
  theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
        axis.text.x = element_text(angle = 30,vjust = 0.9,hjust = 1), panel.grid =element_blank(),axis.text = element_text(colour="black"),
        strip.background = element_rect(colour="black", fill="white"))+
  geom_signif(
    comparisons = list( c("Low-risk", "High-risk")),
    map_signif_level = TRUE, textsize = 6,step_increase = 0.1
  ) 
dev.off()
#————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————


gene <- "genelist"
pdf_file <- paste(gene ,"index_免疫浸润_COAD-pointPlot.pdf", sep="_")
pdf(pdf_file,width=8,height= 5)

ggplot( xdata.NoneScale.risk.select_long2 , aes(x=index,y=value))+
  geom_point(size=0.8,colour ="#219ebc" )+
  geom_smooth(method = "lm",se = T,colour="#219ebc",fill = "#219ebc")+
  stat_poly_eq(aes(label = paste(after_stat(eq.label),after_stat(p.value.label),sep = "~~~")),formula = y~x,parse = TRUE, size=2,colour ="#219ebc")+
  stat_correlation(mapping = use_label("R", "P", "n"),
                   method = "pearson",
                   label.x="right",
                   npcy = 0.8,
                   size=2.5,
                   colour ="#219ebc")+
  labs(x=" ",y=NULL)+
  facet_wrap(.~gene,scales="free",ncol=4)+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
        panel.grid =element_blank(),
        axis.text = element_text(colour="black"),legend.key = element_blank(),
        legend.title=element_blank())+
  labs(fill=' ')
dev.off()  





#=========原癌基因抑癌基因===================

##筛选出risk及原癌基因抑癌基因相关基因的列
###正图
xdata.NoneScale.risk.select <- xdata.NoneScale.risk%>% 
  as.data.frame() %>% 
  dplyr::select(
    FLT3,ATM,NRAS,BRCA2,PALB2,PDGFRA,ESR1,RET,BRAF,NTRK3,CHEK2,index,
    risk) %>%
  as.data.frame() 


#####附图
xdata.NoneScale.risk.select <- xdata.NoneScale.risk%>% 
  as.data.frame() %>% 
  dplyr::select(  ABL1,ALK,BARD1,BRCA1,BRIP1,CDK12,CHEK1,EGFR,ERBB2,
                  EZH2,FANCL,FGFR1,FGFR2,FGFR3,IDH1,IDH2,KIT,KRAS,MET,NF1,
                  NTRK1,NTRK2,PDGFB,PDGFRB,PIK3CA,RAD51C,
                  RAD54L,ROS1,SMARCB1,TSC1,TSC2,index,
                  risk) %>%
  as.data.frame() 


ydata.risk<-ydata.risk%>%
  as.data.frame() 

# 转换数据框为长格式
xdata.NoneScale.risk.select_long2 <- tidyr::pivot_longer(
  xdata.NoneScale.risk.select, 
  cols = c( ABL1, ALK,BARD1,BRCA1,BRIP1,CDK12,CHEK1,EGFR,ERBB2,
            EZH2,FANCL,FGFR1,FGFR2,FGFR3,IDH1,IDH2,KIT,KRAS,MET,NF1,
            NTRK1,NTRK2,PDGFB,PDGFRB,PIK3CA,RAD51C,
            RAD54L,ROS1,SMARCB1,TSC1,TSC2
  ), names_to = "gene", values_to = "value")

###########基因-box图
pdf_file <- paste("原癌基因抑癌基因_COAD-BoxPlot.pdf")

pdf(pdf_file,width=9,height= 4.3)
ggplot(xdata.NoneScale.risk.select_long2, aes(x =factor(risk), y=log2(value) ,fill =factor(risk))) +
  geom_boxplot(size=0.5,varwidth = F, alpha=0.8) + 
  geom_jitter(position=position_jitter(width=0.15,height=0.2), size=0.8,alpha=0.3,colour="#bfdbf7") + 
  facet_wrap(gene~.,scales="free",ncol = 6)+
  scale_fill_manual(values=c("#588b8b","#f5cac3")) +
  theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
        axis.text.x = element_text(angle = 30,vjust = 0.9,hjust = 1), panel.grid =element_blank(),axis.text = element_text(colour="black"),
        strip.background = element_rect(colour="black", fill="white"))+
  geom_signif(
    comparisons = list( c("Low-risk", "High-risk")),
    map_signif_level = TRUE, textsize = 6,step_increase = 0.1
  ) 
dev.off()

#————————————————————————————————————————————————————————————————————————————————————
##筛选出risk及原癌基因抑癌基因相关基因的列
xdata.NoneScale.risk.select <- xdata.NoneScale.risk%>% 
  as.data.frame() %>% 
  dplyr::select(ABL1,ALK,ATM,BARD1,BRAF,BRCA1,BRCA2,BRIP1,CDK12,CHEK1,CHEK2,EGFR,ERBB2,
                ESR1,EZH2,FANCL,FGFR1,FGFR2,FGFR3,FLT3,IDH1,IDH2,KIT,KRAS,MET,NF1,NRAS,
                NTRK1,NTRK2,NTRK3,PALB2,PDGFB,PDGFRA,PDGFRB,PIK3CA,RAD51C,
                RAD54L,RET,ROS1,SMARCB1,TSC1,TSC2) %>%
  as.data.frame() 
ydata.risk<-ydata.risk%>%
  as.data.frame() 
xdata.NoneScale.risk.select_t <- xdata.NoneScale.risk.select %>% t()%>%as.data.frame()
xdata.NoneScale.risk.select_t = as.data.frame(lapply(xdata.NoneScale.risk.select_t,as.numeric))
rownames(xdata.NoneScale.risk.select_t)=colnames(xdata.NoneScale.risk.select)
xdata.NoneScale.risk.select_t <- xdata.NoneScale.risk.select_t %>%as.matrix()
colnames(xdata.NoneScale.risk.select_t ) <- gsub("\\.", "-", colnames(xdata.NoneScale.risk.select_t ))

#可视化-热图
library(pheatmap)
library(RColorBrewer)
#####先加好risk分组信息：
##把risk列提出来变为sample_info
sample_info<-xdata.NoneScale.risk$risk%>%as.data.frame()
##将ssgsea.1的列名加到sample_info的行名上
rownames(sample_info) = colnames(xdata.NoneScale.risk.select_t)
##作图，加了annotation_col 参数，设置 annotation_col =sample_info
pdf(file="pheatmap_hallmark.pdf",width = 10,height =10)
pheatmap(xdata.NoneScale.risk.select_t,scale="row",border=NA,
         color=colorRampPalette(c("#023e8a","white","#d1495b"))(500),
         cluster_row=T,cluster_cols=F,gaps_col =241,
         show_rownames=T,show_colnames = F,treeheight_col =10,treeheight_row =15,legend = T,
         annotation_col =sample_info )
dev.off()




gene <- "genelist"
pdf_file <- paste(gene ,"index_原癌基因抑癌基因_COAD-pointPlot_附图.pdf", sep="_")
pdf(pdf_file,width=10,height= 9.5)

ggplot( xdata.NoneScale.risk.select_long2 , aes(x=index,y=log(value)))+
  geom_point(size=0.8,colour ="#219ebc" )+
  geom_smooth(method = "lm",se = T,colour="#219ebc",fill = "#219ebc")+
  stat_poly_eq(aes(label = paste(after_stat(eq.label),after_stat(p.value.label),sep = "~~~")),formula = y~x,parse = TRUE, size=2,colour ="#219ebc")+
  stat_correlation(mapping = use_label("R", "P", "n"),
                   method = "pearson",
                   label.x="right",
                   npcy = 0.8,
                   size=2.5,
                   colour ="#219ebc")+
  labs(x=" ",y=NULL)+
  facet_wrap(.~gene,scales="free",ncol=6)+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
        panel.grid =element_blank(),
        axis.text = element_text(colour="black"),legend.key = element_blank(),
        legend.title=element_blank())+
  labs(fill=' ')
dev.off()  

#=========代谢：hallmarker——ssgsea热图===================
library(ComplexHeatmap) 

m = matrix(rnorm(10000), nrow = 100)
rownames(m) = 1:100

ha = rowAnnotation(foo = anno_mark(at = c(1:4, 20, 60, 97:100), labels = month.name[1:10]))

pdf("heatmap.pdf",width = 6,height = 7)

Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha,
        
        row_names_side = "left", row_names_gp = gpar(fontsize = 4))

dev.off() 

#————————————————————————————————————————————————————————
A <- xdata.NoneScale.risk[,-32359]
A <-A[,-32359]
A <-A[,-32359]
A<-A%>%t()
for (i in 1:nrow(A)) A[i, ] <- scale(log(unlist(A[i, ] + 1), 2)) #对数据进行标准化处理
samples <- rep(c('High-risk', 'Low-risk'), c(241, 242)) #定义样本分组信息  
genes <- c("A2M","AAAS")
genes <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "HALLMARK_KRAS_SIGNALING_UP",
           "HALLMARK_MYOGENESIS",
           "HALLMARK_ANGIOGENESIS",
           "HALLMARK_APICAL_JUNCTION",
           "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
           "HALLMARK_APOPTOSIS",
           "HALLMARK_ALLOGRAFT_REJECTION",
           "HALLMARK_INFLAMMATORY_RESPONSE",
           "HALLMARK_PEROXISOME",
           "HALLMARK_COMPLEMENT",
           "HALLMARK_IL2_STAT5_SIGNALING",
           "HALLMARK_INTERFERON_GAMMA_RESPONSE",
           "HALLMARK_UV_RESPONSE_DN",
           "HALLMARK_IL6_JAK_STAT3_SIGNALING",
           "HALLMARK_HYPOXIA",
           "HALLMARK_HEDGEHOG_SIGNALING")
genes <- as.data.frame(genes)

Heatmap(A,#表达矩阵
        col = colorRampPalette(c("#0000AA", "#555555", "#AAAA00"))(600),#颜色定义
        show_row_names = F,#不展示行名，
        use_raster=TRUE,
        right_annotation = rowAnnotation(link = anno_mark(at = which(rownames(A) %in% genes$genes),
                                                          labels = genes$genes, labels_gp = gpar(fontsize = 10))),
        top_annotation = HeatmapAnnotation(Group = samples, 
                                           simple_anno_size = unit(2, 'mm'), 
                                           col = list(Group = c('High-risk' = "#00DAE0",  'Low-risk' = "#FF9289")),
                                           show_annotation_name = FALSE))#分组注释

#print(B)



pdf("hallmarker_heatmap.pdf",width=8,height= 8)
#  B + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% genes$genes), 
#                                                                         labels = genes$genes))


dev.off()







# =============== kegg—ssgsea===============
ssgsea<- read.csv("4.GENE SETS OF INTEREST HEATMAP_c2.cp.kegg_long.csv",header=T,row.names = 1)
ssgsea<-ssgsea[-1,]
colnames(ssgsea ) <- gsub("\\.", "-", colnames(ssgsea ))

ssgsea.t<-ssgsea%>% t()%>%as.data.frame()
ssgsea.t    <- ssgsea.t [rownames(xdata.NoneScale.risk), ]

ssgsea.t$risk<-xdata.NoneScale.risk$risk   
#########筛选显著的通路
ssgsea.t.select <- ssgsea.t%>% 
  as.data.frame() %>% 
  dplyr::select(risk,
                KEGG_BUTANOATE_METABOLISM,
                KEGG_PEROXISOME,
                KEGG_GLYCOSYLPHOSPHATIDYLINOSITOL_GPI_ANCHOR_BIOSYNTHESIS,
                KEGG_PYRUVATE_METABOLISM,
                KEGG_STEROID_BIOSYNTHESIS,
                KEGG_FATTY_ACID_METABOLISM,
                KEGG_LYSINE_DEGRADATION,
                KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION) %>%
  as.data.frame() 

# 转换数据框为长格式
df_long <- tidyr::pivot_longer(ssgsea.t.select, cols = c( KEGG_BUTANOATE_METABOLISM,
                                                          KEGG_PEROXISOME,
                                                          KEGG_GLYCOSYLPHOSPHATIDYLINOSITOL_GPI_ANCHOR_BIOSYNTHESIS,
                                                          KEGG_PYRUVATE_METABOLISM,
                                                          KEGG_STEROID_BIOSYNTHESIS,
                                                          KEGG_FATTY_ACID_METABOLISM,
                                                          KEGG_LYSINE_DEGRADATION,
                                                          KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION),
                               names_to = "variable", values_to = "value")

##做box图——facet
#value列变为数值型
df_long$value <- as.numeric(df_long$value)

kegg<- "kegglist"
pdf_file <- paste(kegg,"kegg_COAD-BoxPlot.pdf", sep="_")
pdf(pdf_file,width=6,height=4)
ggplot(df_long, aes(x =factor(risk), y=log2(value) ,fill =factor(risk))) +
  geom_boxplot(size=0.5,varwidth = F, alpha=0.8) + 
  geom_jitter(position=position_jitter(width=0.15,height=0.2), size=0.8,alpha=0.3,colour="#bfdbf7") + 
  facet_wrap(variable~.,scales="free",ncol = 4)+
  scale_fill_manual(values=c("#588b8b","#f5cac3")) +
  theme(legend.position="bottom",panel.background=element_rect(fill='transparent',color='black'),
        axis.text.x = element_text(angle = 30,vjust = 0.9,hjust = 1), panel.grid =element_blank(),axis.text = element_text(colour="black"),
        strip.background = element_rect(colour="black", fill="white"))+
  theme(strip.text.x = element_text(size =4)) +
  geom_signif(
    comparisons = list( c("Low-risk", "High-risk")),
    map_signif_level = TRUE, textsize = 4,step_increase = 0.1
  )
dev.off()  










# ==================== keggtreemap-2 ===================
# devtools::install_github("m-jahn/WeightedTreemaps")
library(WeightedTreemaps)
library(RColorBrewer)

GO2<-read.table("High_risk_up_treemap.txt", header = TRUE, sep = "\t")
oct_coord <- list(
  x = sin(seq(0, 2, 2/8)*pi) * 1000 + 1000,
  y = cos(seq(0, 2, 2/8)*pi) * 1000 + 1000
)
GO2<- voronoiTreemap(
  data = GO2,
  levels = c("Input.number", "Term"),
  cell_size = "Input.number",
  shape = oct_coord,
  seed = 123
)
##绘制树状图
#创建一个pdf文件，准备写入
pdf("kegg_high_risk_up.PDF",width = 5,height = 5)
#画图，写入到pdf
##如果画不出图，就新建画板

drawTreemap(GO2, label_size = 2.5, label_color = "black",
            color_palette = c(brewer.pal(5,"Set3"),brewer.pal(5,"Pastel1")) ,color_level = 2 )
##关闭画板(关闭pdf，解除占用）
dev.off()

#—————————————————————————————————————— 
GO2<-read.table("Low_risk_up_treemap.txt", header = TRUE, sep = "\t")
oct_coord <- list(
  x = sin(seq(0, 2, 2/8)*pi) * 1000 + 1000,
  y = cos(seq(0, 2, 2/8)*pi) * 1000 + 1000
)
GO2<- voronoiTreemap(
  data = GO2,
  levels = c("Input.number", "Term"),
  cell_size = "Input.number",
  shape = oct_coord,
  seed = 123
)
##绘制树状图
#创建一个pdf文件，准备写入
pdf("kegg_low_risk_up.PDF",width = 5,height = 5)
#画图，写入到pdf
##如果画不出图，就新建画板

drawTreemap(GO2, label_size = 2.5, label_color = "black",
            color_palette = c(brewer.pal(5,"Set3"),brewer.pal(5,"Pastel1")) ,color_level = 2 )
##关闭画板(关闭pdf，解除占用）
dev.off()


#  ============= 5) run nearest template prediction in external cohort =============
#哦，等等，我们忘了什么吗？是的，有一个数据集还没有使用，所以让我们看看我们是
#否可以使用这些亚型特异性生物标志物来验证外部Yau队列中当前的乳腺癌症亚型。
#在这一部分中，我们的核心目的是预测外部数据集中每个样本的可能亚型。
#在大多数情况下，这是一个多分类问题，并且所识别的生物标志物可能很难在外部队列
#中进行整体匹配，因此使用基于模型的预测算法可能不可靠。因此，MOVICS为验证队列
#中的亚型预测提供了两种无模型的方法。首先，MOVICS切换到最近模板预测（NTP），
#它可以灵活地应用于跨平台、跨物种和多类预测，而无需对分析参数进行任何优化。
#只需要做一件事就是生成一个模板文件，幸运的是，该文件已经准备好了。

# run NTP in Yau cohort by using up-regulated biomarkers
yau.ntp.pred <- runNTP(expr       = coad.yau$mRNA.expr,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "5.NTP HEATMAP FOR YAU")
#> --original template has 500 biomarkers and 272 are matched in external expression profile.
#> cosine correlation distance
#> 682 samples; 5 classes; 39-66 features/class
#> serial processing; 1000 permutation(s)...
#> predicted samples/class (FDR<0.05)

#>  CS1  CS2  CS3  CS4  CS5 <NA> 
#>  104   85   90  120  140  143
head(yau.ntp.pred$ntp.res)

ntp.res  <-yau.ntp.pred$ntp.res
write.csv(ntp.res,"5.Demo of predicted subtypes in Yau cohort by NTP using subtype-specific upregulated biomarkers identified from TCGA-BRCA cohort.csv")

#上面准备了一个对象yau.ntp.pred，该对象在结构上与getMOIC（）返回的对象相似，
#但只存储cluster.res，如果有额外的数据可用，这些cluster.res可以传递给COMP模块内的函数。
#例如，我在此首先比较Yau队列中预测的5种癌症亚型的生存结果。

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = coad.yau$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "h", # switch to both
                     fig.name         = "5.KAPLAN-MEIER CURVE OF NTP FOR YAU",
                     #p.adjust.method='BY'   #c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none' 
) 
#> --a total of 682 samples are identified.
#> --cut survival curve up to 10 years.
print(surv.yau)
#> $fitd
#> Call:
#> survdiff(formula = Surv(futime, fustat) ~ Subtype, data = mosurv.res, 
#>     na.action = na.exclude)
#> 
#>               N Observed Expected (O-E)^2/E (O-E)^2/V
#> Subtype=CS1 136       51     43.8     1.169     1.452
#> Subtype=CS2 117       51     38.5     4.037     4.876
#> Subtype=CS3 119       33     42.3     2.045     2.517
#> Subtype=CS4 159       43     57.3     3.590     4.820
#> Subtype=CS5 151       50     46.0     0.351     0.441
#> 
#>  Chisq= 11.2  on 4 degrees of freedom, p= 0.02 
#> $fit
#> Call: survfit(formula = Surv(futime, fustat) ~ Subtype, data = mosurv.res, 
#>     na.action = na.exclude, error = "greenwood", type = "kaplan-meier", 
#>     conf.type = "plain")
#> 
#>       n events median 0.95LCL 0.95UCL
#> CS1 136     51     NA      NA      NA
#> CS2 117     51    205     109      NA
#> CS3 119     33    236      NA      NA
#> CS4 159     43    222     191      NA
#> CS5 151     50     NA      NA      NA
#> 
#> $xyrs.est
#> [1] "[Not Available]: argument of xyrs.est was not specified."
#> 
#> $overall.p
#> [1] 0.02410534
#> 
#> $pairwise.p
#> 
#>  Pairwise comparisons using Log-Rank test 
#> 
#> data:  mosurv.res and Subtype 
#> 
#>     CS1   CS2   CS3   CS4  
#> CS2 0.699 -     -     -    
#> CS3 0.158 0.059 -     -    
#> CS4 0.100 0.039 0.901 -    
#> CS5 0.804 0.504 0.244 0.158
#> 
#> P value adjustment method: BH

# 我进一步检查了预测的亚型与PAM50分类之间的一致性。

# compare agreement in Yau cohort
agree.yau <- compAgree(moic.res  = yau.ntp.pred,
                       subt2comp = coad.yau$clin.info[, "PAM50", drop = FALSE],
                       doPlot    = TRUE,
                       fig.name  = "5.YAU PREDICTEDMOIC WITH PAM50")
#> --all samples matched.

print(agree.yau)

#很明显，与TCGA-BRCA相比，Yau队列中预测的癌症亚型可以很好地区分预后，并且在一定程度上与PAM50分类显示出相似的一致性模式。

# 6） 围绕类药物分类器运行分区
# 6) run partition around medoids classifier
# 除了NTP，MOVICS还提供了另一种无模型的方法来预测亚型。具体而言，runPAM（）首先在发现
#（训练）队列（即TCGA-BRCA）中围绕类药物（PAM）分类器训练一个分区，以预测外部验证
#（测试）队列（如BRCA-Yau）中患者的亚型，验证队列中的每个样本都被分配到一个亚型标签，
#其质心与样本具有最高的Pearson相关性17。最后，将进行组内比例（IGP）统计，
#以评估发现和验证队列之间获得的亚型的相似性和再现性18。
yau.pam.pred <- runPAM(train.expr  = fpkm,
                       moic.res    = cmoic.coad,
                       test.expr   = coad.yau$mRNA.expr)
#> --all samples matched.
#> --a total of 7303 genes shared and used.
#> --log2 transformation done for training expression data.
#> --testing expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.

#  yau.pam.pred对象还存储了cluster.res，这些cluster.res可以传递给其他函数，用户可以检查IGP信息，如下所示：
print(yau.pam.pred$IGP)
#>       CS1       CS2       CS3       CS4       CS5 
#> 0.4545455 0.6416667 0.5616438 0.7814570 0.9225806


#    7) run consistency evaluation using Kappa statistics  
# 7） 使用Kappa统计进行一致性评估

# 想知道NTP或PAM在使用发现队列时的准确性吗？想知道不同的预测结果有多一致吗？然后使用runKappa（）：

# predict subtype in discovery cohort using NTP
tcga.ntp.pred <- runNTP(expr      = fpkm,
                        templates = marker.up$templates,
                        doPlot    = FALSE) 
#> --original template has 500 biomarkers and 500 are matched in external expression profile.
#> cosine correlation distance
#> 643 samples; 5 classes; 100-100 features/class
#> serial processing; 1000 permutation(s)...
#> predicted samples/class (FDR<0.05)
#> 
#>  CS1  CS2  CS3  CS4  CS5 <NA> 
#>   99  105  138  155  107   39

# predict subtype in discovery cohort using PAM
tcga.pam.pred <- runPAM(train.expr  = fpkm,
                        moic.res    = cmoic.coad,
                        test.expr   = fpkm)
#> --all samples matched.
#> --a total of 13771 genes shared and used.
#> --log2 transformation done for training expression data.
#> --log2 transformation done for testing expression data.

# check consistency between current and NTP-predicted subtype in discovery TCGA-BRCA
runKappa(subt1     = cmoic.coad$clust.res$clust,
         subt2     = tcga.ntp.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "NTP",
         fig.name  = "7.CONSISTENCY HEATMAP FOR TCGA between CMOIC and NTP")
# check consistency between current and PAM-predicted subtype in discovery TCGA-BRCA
runKappa(subt1     = cmoic.coad_TLS$clust.res$clust,
         subt2     = tcga.pam.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "PAM",
         fig.name  = "7.CONSISTENCY HEATMAP FOR TCGA between CMOIC and PAM")

# check consistency between NTP and PAM-predicted subtype in validation Yau-BRCA
runKappa(subt1     = yau.ntp.pred$clust.res$clust,
         subt2     = yau.pam.pred$clust.res$clust,
         subt1.lab = "NTP",
         subt2.lab = "PAM",
         fig.name  = "7.CONSISTENCY HEATMAP FOR YAU")

#5. Little Trick

# include original clinical information as `clust.res` and a string value for `mo.method` to a list
pseudo.moic.res                 <- list("clust.res" = surv.info,
                                        "mo.method" = "PAM50")

# make pseudo samID
pseudo.moic.res$clust.res$samID <- rownames(pseudo.moic.res$clust.res)

# make pseudo clust using a mapping relationship
pseudo.moic.res$clust.res$clust <- sapply(pseudo.moic.res$clust.res$PAM50,
                                          switch,
                                          "Basal"   = 1, # relabel Basal as 1
                                          "Her2"    = 2, # relabel Her2 as 2
                                          "LumA"    = 3, # relabel LumA as 3
                                          "LumB"    = 4, # relabel LumnB as 4
                                          "Normal"  = 5) # relabel Normal as 5
# Let’s check how the pseudo moic.res object looks like.

head(pseudo.moic.res$clust.res)

#  好吧，一切都准备好了，只需记住您感兴趣的子类型是如何映射的，并可以使用这样的伪对象执行下游分析，如compSurv（），如下所示：
# survival comparison
pam50.brca <- compSurv(moic.res         = pseudo.moic.res,
                       surv.info        = surv.info,
                       convt.time       = "y", # convert day unit to year
                       surv.median.line = "h", # draw horizontal line at median survival
                       fig.name         = "KAPLAN-MEIER CURVE OF PAM50 BY PSEUDO")
#> --a total of 643 samples are identified.
#> --removed missing values.
#> --leaving 642 observations.
#> --cut survival curve up to 10 years.
















#================== cluster_筛选出要进行TF和目标基因相关性计算的基因========================
TF  <- read.table("TF_uniq.txt",header=T,sep="\t") 

fpkm_t.select <-fpkm_t%>% 
  as.data.frame() %>% 
  dplyr::select(
    AEBP2,
    AFF4,
    AHR,
    AR,
    ARID3A,
    ARNT,
    #    ARNTL,
    ARRB1,
    ASCL1,
    ASCL2,
    ATF1,
    ATF2,
    ATF3,
    ATF5,
    ATF7IP,
    ATOH1,
    ATRX,
    BACH1,
    BACH2,
    BARHL1,
    BARX2,
    BATF,
    BCL11A,
    BCL3,
    BCL6,
    BCOR,
    BDP1,
    BHLHE40,
    BMI1,
    BRCA1,
    BRD2,
    BRD3,
    BRD4,
    BRD7,
    BRF1,
    BRF2,
    BTAF1,
    # C10orf12,
    #C17orf49,
    #C17orf96,
    CAMTA2,
    CARM1,
    CASP8AP2,
    CASZ1,
    CBFB,
    CBX2,
    CBX3,
    CBX4,
    CBX8,
    # CCDC101,
    CCNT2,
    CDK2,
    CDK7,
    CDK8,
    CDK9,
    CDX2,
    CEBPA,
    CEBPB,
    CEBPD,
    CEBPZ,
    CENPA,
    CHD1,
    CHD3,
    CHD4,
    CHD7,
    CIITA,
    CLOCK,
    CREB1,
    CREB3,
    CREBBP,
    CREM,
    CTBP2,
    CTCF,
    CTCFL,
    CTNNB1,
    CUX1,
    DBP,
    DDX5,
    DLX2,
    DLX4,
    DMRT1,
    DUX4,
    E2F1,
    E2F2,
    E2F3,
    E2F4,
    E2F6,
    E2F7,
    E2F8,
    EBF1,
    EED,
    EGLN2,
    EGR1,
    EGR2,
    EHF,
    EHMT2,
    ELF1,
    ELF2,
    ELF3,
    ELF5,
    ELK1,
    ELK3,
    ELK4,
    ELL2,
    EMX1,
    EOMES,
    EP300,
    EPAS1,
    ERF,
    ERG,
    ERMAP,
    ESR1,
    ESR2,
    ESRRA,
    ETS1,
    ETS2,
    ETV1,
    ETV4,
    ETV5,
    ETV6,
    ETV7,
    EWSR1,
    EZH1,
    EZH2,
    FANCD2,
    FLI1,
    FOS,
    FOSL1,
    FOSL2,
    FOXA1,
    FOXA2,
    FOXD2,
    FOXG1,
    FOXH1,
    FOXK1,
    FOXM1,
    FOXO3,
    FOXO4,
    FOXP1,
    FOXP2,
    FOXP3,
    GABPA,
    GATA1,
    GATA2,
    GATA3,
    GATA4,
    GATA6,
    GATAD1,
    GFI1B,
    GLI1,
    GLI2,
    GLIS1,
    GLIS2,
    GLIS3,
    GLYR1,
    GMEB1,
    GMEB2,
    GRHL1,
    GRHL2,
    GRHL3,
    GTF2B,
    GTF2I,
    GTF3C1,
    GTF3C2,
    HAND1,
    HAND2,
    HBP1,
    HCFC1,
    HDAC1,
    HDAC2,
    HDAC3,
    HDAC6,
    HDAC8,
    HES1,
    HES4,
    HEY1,
    HIF1A,
    HINFP,
    HIRA,
    HIVEP1,
    HMGN1,
    HMGN3,
    HNF1B,
    HNF4A,
    HNF4G,
    HNRNPL,
    HOXA1,
    HOXA13,
    HOXA4,
    HOXA5,
    HOXA6,
    HOXA7,
    HOXA9,
    HOXB13,
    HOXB4,
    HOXB7,
    HOXC11,
    HOXC6,
    HOXC8,
    HSF1,
    HSF2,
    ICE1,
    ICE2,
    ID3,
    IKZF1,
    IRF1,
    IRF2,
    IRF3,
    IRF4,
    IRF5,
    ISL1,
    JAZF1,
    JMJD6,
    JUN,
    JUNB,
    JUND,
    KAT2A,
    KAT2B,
    KAT5,
    KAT8,
    KDM1A,
    KDM4A,
    KDM4C,
    KDM5A,
    KDM5B,
    KDM5C,
    KLF1,
    KLF10,
    KLF15,
    KLF3,
    KLF4,
    KLF5,
    KLF9,
    KMT2A,
    KMT2B,
    L3MBTL2,
    LEF1,
    LHX2,
    LIN9,
    LMNB1,
    LMO2,
    LYL1,
    MAF,
    MAFB,
    MAFF,
    MAFG,
    MAFK,
    MAX,
    MAZ,
    MBD2,
    MBD3,
    MBD4,
    MECOM,
    MECP2,
    MED1,
    MED12,
    MEF2A,
    MEF2C,
    MEIS1,
    MEIS2,
    MITF,
    MNT,
    # MRE11A,
    MSX1,
    MTA3,
    MXI1,
    MYB,
    MYBL1,
    MYBL2,
    MYC,
    MYH11,
    MYOD1,
    MZF1,
    # NANOG,
    NCOR1,
    NCOR2,
    NFAT5,
    NFATC1,
    NFE2,
    NFE2L2,
    NFIC,
    NFKB2,
    NFYA,
    NFYB,
    NFYC,
    NIPBL,
    `NKX2-1`,
    `NKX2-5`,
    `NKX3-1`,
    NOTCH1,
    NPAT,
    NR1H3,
    NR2C2,
    NR2F1,
    NR2F2,
    NR3C1,
    NR4A1,
    NR5A2,
    NRF1,
    NRIP1,
    OGT,
    ONECUT1,
    ORC1,
    OTX2,
    PALB2,
    PAX3,
    PAX5,
    PAX6,
    PBX1,
    PBX2,
    PBX3,
    PCGF2,
    PCGF6,
    PDX1,
    # PGBD3,
    PGR,
    PHF20,
    PHF8,
    PIAS1,
    PIAS4,
    PITX1,
    PML,
    POLR2A,
    POLR3A,
    POLR3G,
    POU2F1,
    POU2F2,
    POU5F1,
    PPARD,
    PPARG,
    PPARGC1A,
    PRAME,
    PRDM1,
    PRDM14,
    PRKDC,
    RAC3,
    RAD21,
    RARA,
    RARG,
    RB1,
    RBBP5,
    RBCK1,
    RBL1,
    RBL2,
    RBPJ,
    RCOR1,
    RELA,
    REPIN1,
    REST,
    RFX2,
    RFX3,
    RFX5,
    RNF2,
    RORA,
    RUNX1,
    RUNX1T1,
    RUNX2,
    RUNX3,
    RXRA,
    RXRG,
    RYBP,
    SALL4,
    SAP30,
    SETDB1,
    SFMBT1,
    SIN3A,
    SIRT6,
    SIX5,
    SMAD1,
    SMAD2,
    SMAD2,
    SMAD3,
    SMAD4,
    SMARCA4,
    SMARCB1,
    SMARCC1,
    SMARCC2,
    SMC1A,
    SMC3,
    # SNAI2,
    # SNAPC1,
    #SNAPC2,
    #SNAPC4,
    #SNAPC5,
    SOX10,
    SOX17,
    SOX2,
    SOX3,
    SOX4,
    SOX9,
    SP1,
    SP2,
    SP3,
    SP4,
    SPDEF,
    SPI1,
    SPIB,
    SRC,
    SREBF1,
    SREBF2,
    SRF,
    STAG1,
    STAT1,
    STAT2,
    STAT3,
    STAT4,
    STAT5A,
    STAT5B,
    STAT6,
    SUMO1,
    SUMO2,
    SUMO2,
    SUPT20H,
    SUZ12,
    TAF1,
    TAF3,
    TAF7,
    TAL1,
    TBL1X,
    TBL1XR1,
    TBP,
    TBX21,
    TBX3,
    TCF12,
    TCF21,
    TCF3,
    TCF4,
    TCF7L2,
    TEAD4,
    TERF1,
    TET2,
    TF,
    TFAP2A,
    TFAP2C,
    TFAP4,
    TFDP1,
    TFE3,
    THAP1,
    THAP11,
    TOP1,
    TP53,
    TP63,
    TP73,
    TRIM24,
    TRIM28,
    TRRAP,
    TTF2,
    UBE2I,
    UBN1,
    UBP1,
    UBTF,
    USF1,
    USF2,
    VDR,
    VEZF1,
    WDR5,
    #  WHSC1,
    WRNIP1,
    XBP1,
    XRCC4,
    XRN2,
    YY1,
    ZBED4,
    ZBTB10,
    ZBTB17,
    ZBTB2,
    ZBTB33,
    ZBTB7A,
    ZBTB7B,
    ZC3H8,
    ZEB1,
    ZFAT,
    ZFHX3,
    ZFP42,
    ZFX,
    ZKSCAN1,
    ZMIZ1,
    ZNF12,
    ZNF143,
    ZNF217,
    ZNF236,
    ZNF250,
    ZNF263,
    ZNF266,
    ZNF280D,
    ZNF281,
    ZNF384,
    ZNF395,
    ZNF706,
    ZNF711,
    ZNF76,
    ZNF83,
    ZNF84,
    ZNF92,
    ZZZ3) %>%
  as.data.frame()     

TLS.gene <-fpkm_t%>% 
  as.data.frame() %>% 
  dplyr::select(
    CCL2,CCL3,CCL4,CCL5,CCL8,CCL18,CCL19,CCL21,CXCL9,CXCL10,CXCL11,CXCL13,CXCL13,
    CD200,FBLN7,ICOS,SGPP2,SH2D1A,TIGIT,PDCD1,CD4,CCR5,CXCR3,CSF2,IGSF6,IL2RA,CD38,
    CD40,CD5,MS4A1,SDC1,GFI1,IL1R1,IL1R2,IL10,CCL20,IRF4,TRAF6,STAT5A,TNFRSF17
  ) %>%
  as.data.frame()   

fpkm.select<-fpkm_t.select%>%t()%>%as.data.frame()
TLS.gene<-TLS.gene%>%t()%>%as.data.frame()

write.table(fpkm.select,"fpkm.select.txt",sep = "\t")
write.table(TLS.gene,"TLS.gene.txt",sep = "\t")

ydata.tmp<-ydata.risk[,36:38]
write.csv(ydata.tmp,"ydata.tmp.csv")


# ICGC ============
sample_info <- read.delim("E:/database/ICGC/sp_donor_family.all_projects_transfer_specimen.gz",header=T,row.names = NULL)
table(sample_info$project_code)
