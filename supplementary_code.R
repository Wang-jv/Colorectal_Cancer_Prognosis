library(clusterProfiler)
library(fgsea)
library(org.Hs.eg.db)
# supplementary Table14 Pathway enrichment using shuffled gene labels ==========================
deg_res <- read.delim("2023.9.19_胃肠道癌生信分析_MOVICS_TLS/4.其他分析/1.差分表达式分析/consensusMOIC_TCGA-COAD_deseq2_test_result.CS1_vs_Others.txt")

up_degs <- deg_res %>% 
  dplyr::filter(log2fc > 0.35, padj < 0.05)
down_degs <- deg_res %>% 
  dplyr::filter(log2fc < -0.35, padj < 0.05)

upIDs <- bitr(up_degs$id,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = "org.Hs.eg.db"
)
downIDs <- bitr(down_degs$id,
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = "org.Hs.eg.db")

bg_genes <- keys(org.Hs.eg.db, keytype = "ENTREZID")

# 模拟置换（例如1000次）
set.seed(1234)
id_list <- list(
  Upregulated = sample(bg_genes, length(up_degs$id)/20),
  Downregulated = sample(bg_genes, length(down_degs$id))
)
kegg_res <- compareCluster(
  id_list,
  fun = "enrichKEGG",
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  pAdjustMethod = "fdr"
)
kegg_res
write.csv(summary(kegg_res), 
          file = "KEGG随机富集结果.csv",row.names = F)

enrichplot::dotplot(enrich_rand,
                    showCategory = 20,
                    color = "pvalue",  
                    label_format = 60)


kegg <- enrichKEGG(gene = entrezIDs$ENTREZID, 
                   organism = "hsa", 
                   pvalueCutoff = 0.2,
                   pAdjustMethod = "BH")

## Fig.9C Immunohistochemical ===========
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(tidyr)
sheets <- getSheetNames("./免疫组化.xlsx")

df <- map_dfc(set_names(sheets), function(sheet) {
  read.xlsx("./免疫组化.xlsx", sheet = sheet,rowNames = F)$`area%`
})
df2 <- df %>% 
  dplyr::mutate(RiskScore  =  -ACSL6 *0.0015881+
                  CALB2*0.0368085 +
                  PCOLCE2  *0.1360931 +
                  PDZD4 *0.0444171 +
                  PPP1R1A*0.0500747 +
                  PTH1R *0.0488492) %>% 
  dplyr::mutate(risk = ifelse(RiskScore >= median(RiskScore), "High-Risk", "Low-Risk"))
write.csv(df2, file = "免疫组化丰度及风险评分.csv", row.names = T)

df_long <- df2 %>% 
  tibble::rownames_to_column(var = "SampleID") %>% 
  pivot_longer(cols = -c(SampleID, risk), 
               names_to = "variable", 
               values_to = "value")

p0 <- grouped_ggbetweenstats(df_long, 
                             x = "risk", 
                             y = "value",
                             grouping.var    = variable,
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
                             ylab = "Positive area(%)",
                             xlab = "",
                             ggtheme = theme_bw(base_size = 16) +
                               theme(axis.text = element_text(face = "plain", size = 14),
                                     axis.text.y = element_text(face = "plain", size = 14),
                                     axis.title = element_text(face = "plain", size = 18),
                                     plot.title = element_text(size = 22),
                                     plot.subtitle = element_text(size = 9),
                                     plot.caption = element_text(size = 12)))
ggsave(
  filename = "免疫组化丰度-BoxPlot.png",
  plot = p0,
  path = "./",
  width = 20,
  height = 5,
  device = "png",
  dpi = 300
)
ggsave(
  filename = "免疫组化丰度_COAD-BoxPlot.pdf",
  plot = p0,
  path = "./",
  width = 20,
  height = 5,
  device = "pdf"
)

## Fig3.D rectal msk Survival========
library(survival)
library(survminer)

clinical_data <- read.delim("rectal_msk_2022/rectal_msk_2022_clinical_data.tsv")
expr_data <- read.delim("rectal_msk_2022/data_mrna_seq_expression.txt", row.names = NULL, check.names = F)
xdata.tmp <- expr_data %>% 
  dplyr::filter(Hugo_Symbol %in% c("ACSL6","CALB2","PCOLCE2","PDZD4","PPP1R1A","PTH1R")) %>%
  tibble::column_to_rownames(var = "Hugo_Symbol") %>%
  sjmisc::rotate_df()

index <- -xdata.tmp$ACSL6 *0.0015881+
  xdata.tmp$CALB2*0.0368085 +
  xdata.tmp$PCOLCE2  *0.1360931 +
  xdata.tmp$PDZD4 *0.0444171 +
  xdata.tmp$PPP1R1A*0.0500747 +
  xdata.tmp$PTH1R *0.0488492
xdata.tmp$index <- index



surv_data <- clinical_data %>%
  dplyr::select(1:3, starts_with("Overall"), starts_with("Race")) %>% 
  dplyr::mutate(Time = Overall.Survival.Neo,
                Event = ifelse(Overall.Survival.Status == "1:DECEASED", 1, 0))

train.data.plot <- xdata.tmp %>%
  dplyr::select(index) %>%
  # dplyr::mutate(risk = ifelse(index >= median(index), "High TRCMLS", "Low TRCMLS")) %>%
  tibble::rownames_to_column(var = "Sample.ID") %>%
  dplyr::left_join(surv_data, by = "Sample.ID") %>%
  na.omit()

res.cut <- surv_cutpoint(train.data.plot,
                         time = "Time",
                         event = "Event",
                         "index",
                         minprop = .1,
                         progressbar = TRUE)
res.cat <- surv_categorize(res.cut)
my.surv <- Surv(res.cat$Time, res.cat$Event)
fit <- survfit(my.surv ~ index, data = res.cat)

summary(fit)
ggsurvplot(fit, 
           data = res.cat,
           conf.int=F,
           pval=TRUE,
           risk.table.col = "strata",
           tables.height = 0.5,
           #中位生存是否标识，无："none", 横线加垂线："hv", 横线："h", 垂线："v"
           surv.median.line = 'h',
           xscale = 1,
           palette = rev(c("#2EC4B6", "#E71D36")),
           conf.int.style = "step",
           break.x.by = 20,
           title="Survival probability (%)", 
           xlab = paste0("Time (", "Months", ")"),
           risk.table=F,
           risk.table.y.text = FALSE,
           risk.table.height=.2)
ggsave("result/Survival_curve_TRCMLS.pdf", width = 5, height = 5)

## FIG.s1D TCGA-CRC randomly 6 gene set Survival========
set.seed(1)
random_genes <- sample(colnames(xdata.Scale.raw), 6)
train.data <- xdata.Scale.raw %>%
  as.data.frame() %>% 
  dplyr::select(all_of(random_genes)) %>%
  merge(ydata.ll, by = 0) %>%
  na.omit()

head(train.data)
multi_formulas <- as.formula(paste0('Surv(time,status)~ ', paste0(random_genes, collapse = '+')))
multi_models <- coxph(multi_formulas, data = train.data)
summary(multi_models)
riskScore <- predict(multi_models, type = "risk", newdata = train.data) # 选择作为生存分析的模型
train.data.plot <- cbind(train.data, riskScore)
res.cut <- surv_cutpoint(train.data.plot,
                         time = "time",
                         event = "status",
                         "riskScore",
                         minprop = .5,
                         progressbar = TRUE)
res.cat <- surv_categorize(res.cut)
my.surv <- Surv(res.cat$time, res.cat$status)
fit <- survfit(my.surv ~ riskScore, data = res.cat)
summary(fit)
data.survdiff <- survdiff(my.surv~ riskScore,data = res.cat)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
p.val
ggsurvplot(fit, 
           data = res.cat,
           conf.int=F,
           pval=TRUE,
           risk.table.col = "strata",
           tables.height = 0.5,
           #中位生存是否标识，无："none", 横线加垂线："hv", 横线："h", 垂线："v"
           surv.median.line = 'h',
           xscale = 30.4,
           palette = rev(c("#2EC4B6", "#E71D36")),
           conf.int.style = "step",
           break.x.by = 20*30.4,
           title="Survival probability (%)", 
           xlab = paste0("Time (", "Months", ")"),
           risk.table=F,
           risk.table.y.text = FALSE,
           risk.table.height=.2)
ggsave("result/randomly_6Genes_Survival_curve_TRCMLS.pdf", width = 5, height = 5)
