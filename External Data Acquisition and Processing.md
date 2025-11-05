# External Data Acquisition and Processing

## Single-cell data procession
The rawdata was downloaded from GEO database(GSE252646), five samples(SRR27442072 SRR27442073 SRR27442074 SRR27442075 SRR27442076) rawdata were downlaoded.
### rawdata filtered
cat rawdata/sample.txt |while read id;
  do (fastp  -i rawdata/${id}_1.fastq.gz -I rawdata/${id}_2.fastq.gz -o cleandata/${id}_R1.fastq.gz -O cleandata/${id}_R2.fastq.gz &);
done
### cellranger files generated
#!/bin/bash
bin=../pipeline/cellranger-9.0.1/bin/cellranger
db=../pipeline/refdata-gex-GRCm39-2024-A 
ls $bin; ls $db 

fq_dir=/home/data/project/10x/raw/PBMC1 # Change to your own path
$bin count --id= PBMC1_outs \
--localcores= 8 \
--transcriptome= $db \
--fastqs= $fq_dir \
--sample= PBMC1   \
--localmem=64
### matrix filtered, data intergrated and cell annotation
The codes in this step could be referred to the code in file soupx.md, data integration, cluster analysis, cell annotation.md in this project.
### cell to cell communcation
This codes could be found in file cellchat.md in this peoject.

## The transcriptome data of rat model
### gene expression file 
The expression file was downloaded from GEO database(GSE133750).
### diff-expression: limma
df <- read.table("expr_data.txt", header = T, sep = "\t", row.names = 1, check.names = F)
head(df)
list <- c(rep("Control", 3), rep("Early neuritis",3)) %>% factor(., levels = c("Control", "Early neuritis"), ordered = F)

list <- model.matrix(~factor(list)+0)  #把group设置成一个model matrix
colnames(list) <- c("CK", "Treat")
df.fit <- lmFit(df, list)  ## 数据与list进行匹配
df.matrix <- makeContrasts(Early neuritis - Control, levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
head(tempOutput)
nrDEG = na.omit(tempOutput) 
diffsig <- nrDEG  
write.csv(diffsig, "all.limmaOut.csv") 
sigdiff genes: |logFC|>1&adj.P.Val<0.05 # cutoff

### STEM(Science, Technology, Engineering, Mathematics)
This was analysis on the OECloud tools at https://cloud.oebiotech.com. This online analyses software need to register before use, and it is easy to register and use.
The results could be found at link: https://data.mendeley.com/datasets/yxd283n6xf/2

### KEGG analysis
data <- read_excel("path/to/your/7_Profile_42.xlsx", sheet=1) # the file was from the STEM, 7_Profile_42.xls, this file could be found in 'stem_report.zip' file, link: https://data.mendeley.com/datasets/yxd283n6xf/2
gene.df <- bitr(gene$symbol,fromType="ENTREZID",toType="SYMBOL", OrgDb = org.Mm.eg.db)
gene <- gene.df$ENTREZID
head(gene)
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 1,
                 qvalueCutoff = 1)
library(enrichplot)
kkx <- pairwise_termsim(kk)
treeplot(kkx, hclust_method = "average")
ggsave("KEGG.tree.pdf", width = 12, height = 6)

### heatmap
library(pheatmap)
test <- read.csv("target_genes.csv" , row.names = 1) # target_genes.csv could be found in 'stem_report.zip' file, link: https://data.mendeley.com/datasets/yxd283n6xf/2
annotation_col = read.csv("group.csv" , row.names = 1)
ann_colors = list(
  Group = c(SLE = "#D05146",  Health = "#6091AB" )
                 )
pheatmap(test, scale = "row", annotation_col = annotation_col,
         annotation_colors = ann_colors, show_rownames=T)








