## 1. DESeq2
```
library(DESeq2)
# 设置分组信息
group1 = c(rep('Old',5),rep('Young',5))
colData1 <- data.frame(row.names=colnames(data1),group_list=group1)
dds <- DESeqDataSetFromMatrix(countData = data1,colData = colData1,design = ~group_list)
dds <- dds[rowSums(counts(dds)) >= 100,]
dds <- DESeq(dds)
res <- results(dds,contrast = c("group_list","Old","Young"))
resOrdered <- res[order(res$padj),]
DEG <- as.data.frame(resOrdered)
DEG <- na.omit(DEG)
# 筛选差异基因
diff <- DEG[DEG$pvalue<0.05&abs(DEG$log2FoldChange)>1&DEG$baseMean>30,]
diff$change <- "up"
diff[diff$log2FoldChange<0,]$change <- "down"
table(diff$change)
```
## 2. Pathway enrichment analysis
```
library(clusterProfiler)
# GO富集
go_Enrich <- function(Gene_list,Group){
  tmp <- enrichGO(gene = Gene_list,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = 'BP',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.1)
  tmp <- data.frame(tmp)[,c(2,6,8,9)]
  tmp$Group <- Group
  return(tmp)
}
t <- go_Enrich(row.names(up), 'up_gene')
t$LogP <- -log10(t$p.adjust)
t$Description <- factor(t$Description, levels = rev(t$Description))
# kegg富集
kegg_Enrich <- function(genes,Group){
  id <- bitr(genes, fromType = 'SYMBOL',
             toType = c('ENTREZID'), OrgDb = 'org.Hs.eg.db')
  kk <- enrichKEGG(gene = id$ENTREZID, organism = "hsa", keyType = "kegg",
                   pAdjustMethod = "BH", pvalueCutoff = 0.05)
  kk <- data.frame(kk)[,c(2,5,6)]
  kk$Group <- Group
  return(kk)
}
tt <- kegg_Enrich(unique(df$V5), 'NED_gene')
tt$LogP <- -log10(tt$p.adjust)
tt$Description <- factor(tt$Description, levels = rev(tt$Description))
# 可视化
p <- ggplot(t, aes(x =LogP, y = Description, fill = Group)) + geom_col()+theme_bw()
mytheme <- theme(
  legend.position = 'none',
  axis.text.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = 'black',size = 0.75),
  axis.line.y = element_line(color = 'black',size = 0.75),
  axis.text = element_text(size = 12)
)
p1 <- p + mytheme+scale_x_continuous(expand = c(0,0))
p2 <- p1 +
  geom_text(data = t,
            aes(x = 0.05, y = Description, label = Description),
            size = 4.5,
            hjust = 0)+
  labs(x = '-log10(p-adj)', y = ' ') + 
  scale_fill_manual(values = c('#deb941','#92BD4B'))
p2
```
