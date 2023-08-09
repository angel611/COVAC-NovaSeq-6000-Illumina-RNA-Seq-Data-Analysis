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
```

#Determine the differentially expressed genes between visits and between groups using edgeR

#Create a DGE list\*\* #edgeR works on a table of integer read counts, with rows corresponding to genes and columns to independent libraries. edgeR stores data in a simple list-based data object called a DGEList. This type of object is easy to use because it can be manipulated like any list in R. You can make this in R by specifying the counts and the groups in the function DGEList().

```
d <- DGEList(counts=as.matrix(df))
d$samples
```
#Filtering the data #First get rid of genes which do not occur frequently enough. We can choose this cutoff by saying we must have at least 100 counts per million (calculated with cpm() in R on any particular gene that we want to keep. #Here, we're only keeping a gene if it has a cpm of 100 or greater for at least two samples.

```
dim(d)
head(cpm(d))
apply(d$counts, 2, sum) # total gene counts per sample
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)
```

#Filtering reduces the dataset so there is very little power to detect differential expression, so little information is lost by filtering. After filtering, it is a good idea to reset the library sizes

```{r}
d$samples$lib.size <- colSums(d$counts)
d$samples
```


#Differential Expression #The function exactTest() conducts tagwise tests using the exact negative binomial test. The test results for the n most significant tags are conveniently displayed by the topTags() function. By default, Benjamini and Hochberg's algorithm is used to control the false discovery rate (FDR).

```
design <- ~ Batch
d$samples$group <- rep(c("group1", "group2"))
DG <- estimateDisp(d,robust=TRUE)
Group <- exactTest(DG, pair=c("group1","group2")) 
```

#The total number of differentially expressed genes at FDR\< 0:05 is:

```
DEG <- decideTestsDGE(Group, adjust.method="BH", p.value=0.05)
summary(DEG)
```

# differentially expressed tags from the naive method in d1

```
de1tags12 <- rownames(d)[as.logical(DEG)] 
plotSmear(Group, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")
```

#Plot MA

```
plotMA(d, ylim=c(-2,2))
```
#Adding more information to gene expression matrix

#Heat map

```{r}
Heatmap(rawcounts, name = "Expression", cluster_columns = TRUE, show_column_dend = FALSE, cluster_column_slices = TRUE, column_title_gp = gpar(fontsize = 8), column_gap = unit(0.5, "mm"), cluster_rows = TRUE, show_row_dend = FALSE,row_names_gp = gpar(fontsize = 4), column_title_rot = 90, top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))), show_column_names = FALSE, use_raster = TRUE, raster_quality = 4)
```

# Correcting p-values for multiple comparisons using Benjamin and Hochberg method

```{r}
library(dplyr)
alpha <- .05

diffexp <- df %>%
   mutate(padj = p.adjust(0.05, method="BH")) %>%
   filter(log2FC > 0, padj < alpha)
```

```{r}
# Create a basic volcano plot
ggplot(data = severevshealthy_degresults_1, aes(x = log2fc, y = -log10(pval))) +
         geom_point()
```

#Pathway enrivhment analysis

```{r}
GO_file = "/Users/admin/Desktop/Gene\ counts/h.all.v2023.1.Hs.symbols.gmt.txt "
GSEA = function(df, GO_file, pval = 0.05) {
  set.seed(54321)
  library(dplyr)
  library(fgsea)

  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)

  fgRes <- fgsea::fgsea(pathways = myGO,
                           stats = gene_list,
                           minSize=15, ## minimum gene set size
                           maxSize=400, ## maximum gene set size
                           nperm=10000) %>% 
                  as.data.frame() %>% 
                  dplyr::filter(padj < !!pval) %>% 
                  arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))

  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGO,
                                      stats = gene_list)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))

  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))

  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  header = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")

  colos = setNames(c("firebrick2", "dodgerblue2"),
                 c("Up-regulated", "Down-regulated"))

g1= ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
  geom_point( aes(fill = Enrichment, size = size), shape=21) +
  scale_fill_manual(values = colos ) +
  scale_size_continuous(range = c(2,10)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=header) + 
        th

  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}
```


```{r}
library(dplyr)

S4table = read.csv("/Users/admin/Desktop/Gene\ counts/Copy\ of\ pone.0145322.s006.csv", header=TRUE, skip =1) %>%
  filter(Gene.Symbol != "")
gene_list = S4table$DESeq2.Log2.Fold.Change
names(gene_list) = S4table$Gene.Symbol
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
head(gene_list)
GO_file = "/Users/admin/Desktop/Gene/counts/h.all.v2023.1.Hs.symbols.gmt.txt "
res = GSEA(gene_list, GO_file, pval = 0.05)
dim(res$Results)
## [1] 159   9
res$Plot
plot_geneset_clusters = function( gs_results, GO_file, min.sz = 4, main="GSEA clusters"){
  library(ggplot2)
  library(ggrepel)
  library(stringr)

  myGO = fgsea::gmtPathways(GO_file)
  df = matrix(nrow=nrow(gs_results), ncol = nrow(gs_results), data = 0)
  rownames(df) = colnames(df) = gs_results$pathway

  for ( i in 1:nrow(gs_results)) {
    genesI =  unlist(myGO[names(myGO) == gs_results$pathway[i] ])
    for (j in 1:nrow(gs_results)) {
      genesJ = unlist(myGO[names(myGO) == gs_results$pathway[j] ])
      ## Jaccards distance  1 - (intersection / union )
      overlap = sum(!is.na(match(genesI, genesJ )))
      jaccards = overlap / length(unique(c(genesI, genesJ) ))
      df[i,j] = 1-jaccards
    }
  }

  ## Cluster nodes using dynamic tree cut, for colors
  distMat = as.dist(df)
  dendro = hclust(distMat, method = "average" )
  clust = dynamicTreeCut::cutreeDynamicTree( dendro, minModuleSize = min.sz )
  ## Note: in dynamicTreeCut, cluster 0, is a garbage cluster for things that dont cluster, so we remove it

  gs_results$Cluster = clust
  gs_results = gs_results[gs_results$Cluster != 0, ]

  ## select gene sets to label for each clusters
  bests = gs_results %>%  
    group_by( Cluster ) %>% 
    top_n(wt = abs(size), n = 1) %>% 
    .$pathway
  ## determine cluster order for plotting
  clust_ords = gs_results %>% 
    group_by( Cluster ) %>% 
    summarise("Average" = NES ) %>% 
    arrange(desc(Average)) %>% 
    .$Cluster %>% 
    unique

  gs_results$Cluster = factor(gs_results$Cluster, levels = clust_ords)

  gs_results$Label = ""
  gs_results$Label[gs_results$pathway %in% bests ] = gs_results$pathway[gs_results$pathway %in% bests ]
  gs_results$Label = str_remove(gs_results$Label, "GO_")
  gs_results$Label = tolower(gs_results$Label)

  g1 = ggplot(gs_results, aes(x = Cluster, y = NES, label = Label )) +
    geom_jitter( aes(color = Cluster,  size = size), alpha = 0.8, height = 0, width = 0.2 ) +
    scale_size_continuous(range = c(0.5,5)) +
    geom_text_repel( force = 2, max.overlaps = Inf) +
    ggtitle(main) +
    th

return(g1)
}
plot_geneset_clusters( gs_results = res$Results[res$Results$NES > 0, ], 
                       main = "Up-regulated GSEA clusters",
                                  GO_file = GO_file,
                                  min.sz = 4 )

plot_geneset_clusters( gs_results = res$Results[res$Results$NES < 0, ], 
                       main = "Down-regulated GSEA clusters",
                                  GO_file = GO_file,
                                  min.sz = 4 )
```
#Data exploration

```{r}
plotMDS(d, method="bcv",labels = NULL, pch = NULL, cex = 1,
     dim.plot = c(1,2), gene.selection = "pairwise",xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE)
```

