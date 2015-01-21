# Project Oyster

Analysis of heat shocked and non heat shocked oysters using R

Start by adding an information library

> library(DESeq2)

Then get the data file

> data <- read.table("/Users/alumnomatlab/desktop/project/data/Cgigas-HS-count.txt", header = T, sep = "\t")
> rownames(data) <- data$Feature
> data <- data[,-1]

Define what each column of data represents

> deseq2.colData <- data.frame(condition=factor(c(rep("Control", 3), rep("Treated", 3))), 
+                              type=factor(rep("single-read", 6)))
> rownames(deseq2.colData) <- colnames(data)
> deseq2.dds <- DESeqDataSetFromMatrix(countData = data,
+                                      colData = deseq2.colData, 
+                                      design = ~ condition)

Use the library to compare with the data

> deseq2.dds <- DESeq(deseq2.dds)

Readout

estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing

Create the results file and give it a look

> deseq2.res <- results(deseq2.dds)
> deseq2.res <- deseq2.res[order(rownames(deseq2.res)), ]
> head(deseq2.res)

Readout

log2 fold change (MAP): condition Treated vs Control 
Wald test p-value: condition Treated vs Control 
DataFrame with 6 rows and 6 columns
               baseMean log2FoldChange     lfcSE       stat       pvalue         padj
              <numeric>      <numeric> <numeric>  <numeric>    <numeric>    <numeric>
CGI_10000001  31.577608      0.6109184 0.5772895  1.0582530 2.899401e-01 5.004287e-01
CGI_10000002  31.411455      0.2376639 0.5525684  0.4301077 6.671173e-01 8.108845e-01
CGI_10000003  27.843547      1.0103915 0.7971988  1.2674273 2.050025e-01 4.046120e-01
CGI_10000004   1.634132      0.5839041 0.7983546  0.7313844 4.645444e-01           NA
CGI_10000005   0.000000             NA        NA         NA           NA           NA
CGI_10000009 736.886210      4.3733499 0.3622502 12.0727342 1.471621e-33 2.084212e-29

Dimensions of the results file

> dim(deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ])

Readout

[1] 3155    6

Assign the results to a temporary file and graph first layer

> tmp <- deseq2.res
> plot(tmp$baseMean, tmp$log2FoldChange, pch=20, cex=0.45, ylim=c(-3, 3), log="x", col="black",
+      main="DEG Virus Exposure  (pval <= 0.05)",
+      xlab="mean of normalized counts",
+      ylab="Log2 Fold Change")

Readout

Mensajes de aviso perdidos
In xy.coords(x, y, xlabel, ylabel, log) :
  1357 x values <= 0 omitted from logarithmic plot
  
Graph second layer

> tmp.sig <- deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ]
> points(tmp.sig$baseMean, tmp.sig$log2FoldChange, pch=20, cex=0.45, col="green")
> abline(h=c(-1,1), col="red")
> write.table(tmp.sig, "/Users/alumnomatlab/desktop/project/output/", row.names = T)




Generate the results table

> write.table(tmp.sig, "/Users/alumnomatlab/desktop/project/output/oyster.tab", row.names = T)



