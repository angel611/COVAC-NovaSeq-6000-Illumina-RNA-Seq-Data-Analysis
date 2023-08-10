# COVAC-NovaSeq-6000-Illumina-RNA-Seq-Data-Analysis

## Installation

We install most tools from conda and some from bioconductor.

### Required tools

+ **Fastqc** for checking the quality of the data
+ Trimmomatic: 
+ star
+ sva
+ edgeR
+ DESeq2
+ multiqc
+ fgsea
+ ComplexHeatmap
+ volcanoPlot
+ ggplot2
+ RColorBrewer

Adding a package to apath
```
source ~/.bash_profile
```

Install trimmomatic using `conda` via `biconda`

```
conda install trimmomatic -c bioconda
```

checking the quality of fastq

```
fastqc demodata/sample5_1.*fq -o qc
```
open the html file
```
open qc/sample5_1000_1_fastqc.html
```
Run multiqc to aggragate all the qc results for all samples into a single file

```
multiqc qc/* -o qc
```

Trimming with Trimmomatic

```
trimmomatic PE demodata/sample5_1000_1.fq demodata/sample5_1000_2.fq paired_trimmed_1.fq unpaired_trimmed_1.fq paired_trimmed_2.fq unpaired_trimmed_2.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:15 TRAILING:20 SLIDINGWINDOW:4:10 MINLEN:100 AVGQUAL:20 -phred33

```
Trimming many files at once. Create a loop.Put the file in the same directory as the fastq.files and run the command below in terminal.
```
bash trimmomatic-batch.sh
```

  Generate genome index

```
star --runMode genomeGenerate --genomeDir . --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile  gencode.v29.annotation.gtf --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 
```

Map reads onto the reference genome and generate the read GeneCounts-applies to one file
```
star --readFilesIn paired_trimmed_1.fq paired_trimmed_2.fq --genomeDir .  --outSAMtype BAM SortedByCoordinate  --outFileNamePrefix SRR629561 --quantMode GeneCounts
```
Referencing to the human reference genome for many files

Mapping many files at once. Create a loop.Put the file in the same directory as the fastq.files and run the command below in terminal.

```
bash star-batch.sh
```

Count the number of reads

```
for i in `ls *.fastq`; do echo $(cat ${i} | wc -l)/4|bc; 
```
Look at the bam file

```/usr/local/Cellar/samtools/1.17/bin/samtools view -h sample5Aligned.sortedByCoord.out.bam | less```

Basic alignment/mapping summary

```/usr/local/Cellar/samtools/1.17/bin/samtools flagstat sample5Aligned.sortedByCoord.out.bam```

Subasample FASTQ files to pick a Million reads incase sub-sampling is needed

```seqtk sample -s100 sample5_1.fq 1000000 > sample5_sub_1.fq```
```seqtk sample -s100 sample5_2.fq 1000000 > sample5_sub_2fq```

Counting paired reads after trimming

```for i in `ls paired_trimmed*.fq`; do echo $(cat ${i} | wc -l)/4|bc; done```

Counting unpaired reads after trimming

```for i in `ls unpaired_trimmed*.fq`; do echo $(cat ${i} | wc -l)/4|bc; done```

Run HTSEQ - for counting gene abundance

```htseq-count --stranded=yes --idattr=ID --type=gene --mode=union --format=bam SRR629561Aligned.sortedByCoord.out.bam gencode.v29.annotation.gff3 > sample5_pergene_counts.txt```

Counting the gene abundance from different files. Create a loop.Put the file in the same directory as the bam files and run the command below in terminal.

```
bash htseq-batch.sh
```
Proceed to R and import the gene counts in R

#Load the packages

```
pacman::p_load(sva, edgeR, tidyverse,ggplot2,limma,ComplexHeatmap,circlize, volcanoPlot,tidyverse,dplyr,RColorBrewer,ggrepel,fgsea)
```

#Load the count table data and combine the HTSeq outputs

```
files <- list.files(path=".", pattern="*.txt")
# using perl to manpulate file names by trimming file extension
labs <- paste("", gsub("\\.txt", "", files, perl=TRUE), sep="")
cov <- list()
for (i in labs) {
  filepath <- file.path(".",paste(i,".txt",sep=""))
  cov[[i]] <- read.table(filepath,sep = "\t", header=F, stringsAsFactors=FALSE)
  colnames(cov[[i]]) <- c("ENSEMBL_GeneID", i)
}
## construct one data frame from list of data.frames using reduce function
df <-Reduce(function(x,y) merge(x = x, y = y, by ="ENSEMBL_GeneID"), cov)
```

#Removing batch effects from the data using Combat function from sva R-Bioconductor package

```
# 'batch_variable' is the name of the column representing batch information

# Convert data to a matrix if it's a data frame
if (is.data.frame(df)) {
    data_matrix <- as.matrix(df)
}

# Run Combat
Batcheffectsremoved <- ComBat(df, batch = batch_variable)
```


#Normalising data using edgeR

```
d <- calcNormFactors(df)
# Adjust raw counts to normalized counts using normalization factors
normalized_counts <- cpm(df, normalized.lib.sizes = d$norm.factors)
#The cpm() function calculates CPM-normalized counts using the provided normalization factors.
```

#Filtering the data #First get rid of genes which do not occur frequently enough. We can choose this cutoff by saying we must have at least 100 counts per million (calculated with cpm() in R on any particular gene that we want to keep. #Here, we're only keeping a gene if it has a cpm of 100 or greater for at least two samples.
#Filtering reduces the dataset so there is very little power to detect differential expression, so little information is lost by filtering. After filtering, it is a good idea to reset the library sizes

```
dim(d)
head(cpm(d))
apply(d$counts, 2, sum) # total gene counts per sample
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)
```

#Determine the differentially expressed genes between visits and between groups using edgeR

#Create a DGE list\*\* #edgeR works on a table of integer read counts, with rows corresponding to genes and columns to independent libraries. edgeR stores data in a simple list-based data object called a DGEList. This type of object is easy to use because it can be manipulated like any list in R. You can make this in R by specifying the counts and the groups in the function DGEList().


```
# Assuming you have sample information stored in a data frame 'sample_info'

# Create a DGEList object from normalized counts
dge <- DGEList(counts = normalized_counts, group = sample_info$group_column)
# Perform differential expression analysis
dge <- calcNormFactors(dge)
design <- model.matrix(~group_column, data = sample_info)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
results <- glmQLFTest(fit, contrast = "group_column")
# Access differential expression results
top_diff_genes <- topTags(results, n = 10)
```

#Second method of differential expression

#Differential Expression #The function exactTest() conducts tagwise tests using the exact negative binomial test. The test results for the n most significant tags are conveniently displayed by the topTags() function. By default, Benjamini and Hochberg's algorithm is used to control the false discovery rate (FDR).

```
# Perform exact test for differential expression
exact_test_results <- exactTest(dge)
# Display top significant genes
top_diff_genes <- topTags(exact_test_results, n = 10)
```

#Identify differentially expressed genes (DEGs) based on the results of a differential expression analysis. 
#"BH" refers to the Benjamini-Hochberg method, which controls the false discovery rate (FDR) when dealing with multiple hypothesis tests

```
# Assuming 'exact_test_results' contains the results of a differential expression analysis
DEG <- decideTestsDGE(exact_test_results, adjust.method = "BH", p.value = 0.05)
# Summary of the results
summary(DEG)
```
#Visualizing differentially expressed genes (DEGs)

```
de1tags12 <- rownames(df)[as.logical(DEG)] 
plotSmear(exact_test_results, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")
```

#Adding more information to gene expression matrix

#Heat map

```
Heatmap(rawcounts, name = "Expression", cluster_columns = TRUE, show_column_dend = FALSE, cluster_column_slices = TRUE, column_title_gp = gpar(fontsize = 8), column_gap = unit(0.5, "mm"), cluster_rows = TRUE, show_row_dend = FALSE,row_names_gp = gpar(fontsize = 4), column_title_rot = 90, top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))), show_column_names = FALSE, use_raster = TRUE, raster_quality = 4)
```

#Correcting p-values(FDR-corrected) for multiple comparisons using Benjamin and Hochberg method
#This is used to identify genes that are both upregulated and statistically significant after FDR correction.

```
library(dplyr)
alpha <- .05
diffexp <- df %>%
   mutate(padj = p.adjust(p_value, method="BH")) %>%
   filter(log2FC > 0, padj < alpha)
```
```
# Create a basic volcano plot
ggplot(data = diff_data, aes(x = log2fc, y = -log10(pval))) +
  geom_point() +
  labs(x = "Log2 Fold Change", y = "-log10(p-value)",
       title = "Volcano Plot of Differentially Expressed Genes")
```
