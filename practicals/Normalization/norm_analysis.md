### Load packages

First of all we will load *scran* and *scater* as well as some other R
packages which we are going to use in this lab:

``` r
suppressMessages(require(scran))
suppressMessages(require(scater))
suppressMessages(require(edgeR))
suppressMessages(require(matrixStats))
suppressMessages(require(mixOmics))
suppressMessages(require(RColorBrewer))
suppressMessages(require(Rtsne))
suppressMessages(require(limSolve))
suppressMessages(require(DESeq2))
suppressMessages(require(tweeDEseq))
```

### Load Expression Values and Metadata

Let us start with loading the matrix of raw gene expression counts and
filtering away lowly expressed genes. Lowly expressed genes should be
excluded from the downstream analysis as they might lead to spurious
results. To correctly compare expression of different genes one needs to
normalize the raw counts by gene length, i.e. the RPKM normalization
does. Therefore we will use RPKM normalized data set in order to exclude
lowly expressed genes which we define as genes with mean RPKM count
across all cells below 1.

``` r
unzip("data/ILC/ensembl_countvalues_ILC.zip", exdir = "data/ILC/")
unzip("data/ILC/ensembl_rpkmvalues_ILC.zip", exdir = "data/ILC/")
D <- read.table("data/ILC/ensembl_countvalues_ILC.csv",sep=",",header=T,row.names=1)
rpkm <- read.table("data/ILC/ensembl_rpkmvalues_ILC.csv",sep=",",header=T,row.names=1)
D<-D[rowMeans(rpkm)>=1,]
rpkm<-rpkm[rowMeans(rpkm)>=1,]
D[1:4,1:4]
```

    ##                 T74_P1_A9_ILC1 T74_P1_B4_NK T74_P1_B7_ILC2 T74_P1_B9_NK
    ## ENSG00000000003              0            0              0            0
    ## ENSG00000000419              0            0              0            0
    ## ENSG00000000457              0            0              0            0
    ## ENSG00000000460              1            0              1            0

``` r
dim(D)
```

    ## [1] 15454   648

The rows of the matrix represent Ensembl gene IDs (you can convert them
to gene symbols using biomaRt package) from 15454 genes, and the columns
are IDs from 648 cells from different individuals sequenced at different
plates. To see how many individuals and plates we have let us load the
meta-information and have a look:

``` r
M <- read.table("data/ILC/Metadata_ILC.csv",sep=",",header=T,row.names=1)
M$Plate<-matrix(unlist(strsplit(as.character(M$Plate),"_")),byrow=TRUE,ncol=2)[,2]
head(M)
```

    ##                 Plate Donor Celltype
    ## T74_P1_A9_ILC1     P1   T74     ILC1
    ## T74_P1_B4_NK       P1   T74       NK
    ## T74_P1_B7_ILC2     P1   T74     ILC2
    ## T74_P1_B9_NK       P1   T74       NK
    ## T74_P1_D10_ILC2    P1   T74     ILC2
    ## T74_P1_E1_ILC3     P1   T74     ILC3

Thus we have cells from 3 individuals with IDs:

``` r
levels(factor(M$Donor))
```

    ## [1] "T74" "T75" "T86"

that were pooled together and sequenced at 4 plates with IDs:

``` r
levels(factor(M$Plate))
```

    ## [1] "P1" "P2" "P3" "P4"

and finally we have 4 cell-types with the following IDs:

``` r
levels(factor(M$Celltype))
```

    ## [1] "ILC1" "ILC2" "ILC3" "NK"

We will color cells on the PCA, MDS and tSNE plots by their cell type,
therfore let us add a column of individual colors for each cell type:

``` r
Celltype_Colors<-data.frame(Celltype=levels(factor(M$Celltype)),Color=brewer.pal(length(levels(factor(M$Celltype))),'Dark2'))
M$Celltype_Color<-as.character(Celltype_Colors$Color[match(M$Celltype,Celltype_Colors$Celltype)])
```

Finally, before we start comparing different normalization strategies,
let us check how many spikes and how many endogenous genes we have after
we have filtered the unnormalized data set:

``` r
print(paste0("There are ",sum(grepl("ERCC_",rownames(D)))," spikes in the scRNAseq data set after filtering"))
```

    ## [1] "There are 71 spikes in the scRNAseq data set after filtering"

``` r
D_ERCC<-D[grepl("ERCC_",rownames(D)),]
D_ERCC[1:5,1:5]
```

    ##                                         T74_P1_A9_ILC1 T74_P1_B4_NK
    ## ERCC_0.11444092:mix1_0.11444092:mix2                 0            1
    ## ERCC_0.22888184:mix1_0.22888184:mix2..2              4            2
    ## ERCC_0.22888184:mix1_0.34332275:mix2..1              0            0
    ## ERCC_0.45776367:mix1_0.11444092:mix2..2              0            0
    ## ERCC_0.45776367:mix1_0.45776367:mix2                 0            0
    ##                                         T74_P1_B7_ILC2 T74_P1_B9_NK
    ## ERCC_0.11444092:mix1_0.11444092:mix2                 1            1
    ## ERCC_0.22888184:mix1_0.22888184:mix2..2              4            4
    ## ERCC_0.22888184:mix1_0.34332275:mix2..1              0            0
    ## ERCC_0.45776367:mix1_0.11444092:mix2..2              0            0
    ## ERCC_0.45776367:mix1_0.45776367:mix2                 0            0
    ##                                         T74_P1_D10_ILC2
    ## ERCC_0.11444092:mix1_0.11444092:mix2                  0
    ## ERCC_0.22888184:mix1_0.22888184:mix2..2               1
    ## ERCC_0.22888184:mix1_0.34332275:mix2..1               0
    ## ERCC_0.45776367:mix1_0.11444092:mix2..2               0
    ## ERCC_0.45776367:mix1_0.45776367:mix2                  0

``` r
print(paste0("There are ",sum(grepl("ERCC_",rownames(D))==FALSE)," endogenous genes in the scRNAseq data set after filtering"))
```

    ## [1] "There are 15383 endogenous genes in the scRNAseq data set after filtering"

This unnormalized and filtered data set we will use later for applying
different normalization methods. Note that RPKM normalized data set
contains exactly the same genes since we used mean RPKM count above or
equalt to 1 as a criterion for filtering.

### Defining Visualization Functions

In this section we will define a few custom functions which will be
useful later for comparison of normalization techniques. Since we are
going to repeat the same visualization commands such as plotting CV^2
vs. mean expression, PCA, tSNE for data sets normalized by different
methods it makes sense to put those commands in separate functions. The
functions we define here take a data set (normalized or unnormalized) as
*df* and a short custom description of the normalization method as
*description*, and produce a bunch of informatice plots. There is no
need to go through the code of the functions unless you really want it,
the code is not difficult but bulky and boring. Below, *df* is a matrix
with samples as columns and features as rows, and *description* is
supposed to say something about the normalization used, i.e. “RAW”,
“RPKM”, “DECONVOLUTION” or “SPIKE LIB SIZE” normalization; *layout* and
*print* define how many plots we want to combine and if the
visualization should be displayed, respectively. First function plots
the squared coefficient of variation against mean gene expression for
endogenous genes and spikes, fit the spike points and determines
endogeneous genes with variation above spikes level:

``` r
cv2_vs_mean_function<-function(df,description,layout,print)
{
if(layout==2){par(mfrow=c(layout,layout))}
df_ERCC<-df[grepl("ERCC_",rownames(df)),]
#PLOT CV^2 VS. MEAN PLOT FOR ENDOGENOUS GENES
mean_expr<-as.numeric(rowMeans(df,na.rm=TRUE))
sd_expr<-rowSds(as.matrix(df),na.rm=TRUE)
cv_squared_expr<-(sd_expr/mean_expr)^2
if(print==TRUE)
{
plot(log10(cv_squared_expr)~log10(mean_expr),pch=20,cex=0.5,xlab=paste0("log10 ( MEAN ",description," COUNT )"),ylab="log10 ( CV² )",main=paste0(description," COUNTS")) 
}
#PLOT CV^2 VS. MEAN PLOT FOR SPIKES
mean_expr_ercc<-as.numeric(rowMeans(df_ERCC,na.rm=TRUE))
sd_expr_ercc<-rowSds(as.matrix(df_ERCC),na.rm=TRUE)
cv_squared_expr_ercc<-(sd_expr_ercc/mean_expr_ercc)^2
mean_expr_ercc_filtered<-mean_expr_ercc
cv_squared_expr_ercc_filtered<-cv_squared_expr_ercc
if(print==TRUE)
{
points(log10(cv_squared_expr_ercc_filtered)~log10(mean_expr_ercc_filtered),col="red",pch=20,cex=1.5)
}
#FIT SPIKES WITH A CURVE
fit_expr_ercc<-loess(log10(cv_squared_expr_ercc_filtered)[is.finite(log10(mean_expr_ercc_filtered))]~log10(mean_expr_ercc_filtered)[is.finite(log10(mean_expr_ercc_filtered))],span=1)
j<-order(log10(mean_expr_ercc_filtered)[is.finite(log10(mean_expr_ercc_filtered))])
if(print==TRUE)
{
lines(fit_expr_ercc$fitted[j]~log10(mean_expr_ercc_filtered)[is.finite(log10(mean_expr_ercc_filtered))][j],col="red",lwd=3)
}
pred_expr<-predict(fit_expr_ercc,log10(mean_expr))
#DETERMINE VARIABLE GENES THAT ARE ABOVE THE SPIKEINS CURVE
filtered_expr<-df[log10(cv_squared_expr)>=pred_expr,]
filtered_expr<-filtered_expr[grepl("NA",rownames(filtered_expr))==FALSE,]
if(print==TRUE)
{
text(3,2.5,paste0(dim(filtered_expr)[1]," variable genes"))
}
if(print==FALSE)
{
return(filtered_expr)
}
}
```

Second function accepts the filtered data set (only endogenous genes
above spike level) and makes a PCA plot using *prcomp* R function:

``` r
pca_function <- function(df, description, layout)
{
  if(layout == 2) { par(mfrow = c(layout, layout)) }
  
  filtered_expr <- cv2_vs_mean_function(df, description, layout, FALSE)
  
  # PCA PLOT WITH PRCOMP
  PC <- prcomp(log10(t(filtered_expr + 1)), center = TRUE, scale = FALSE)
  
  plot(PC$x[,1:2], col = M$Celltype_Color, main = paste0("PCA PLOT: ", description, " COUNTS"))
  
  # Get unique cell types and their corresponding colors
  cell_types <- unique(M$Celltype)
  cell_colors <- unique(M$Celltype_Color)
  
  legend("topleft", 
         legend = cell_types, 
         cex = 1, 
         fill = cell_colors, 
         inset = 0.02)
}
```

Third function accepts the filtered data set (only endogenous genes
above spikes level) and makes a tSNE plot using *Rtsne* R function:

``` r
tsne_function <- function(df, description, layout)
{
  if(layout == 2) { par(mfrow = c(layout, layout)) }
  
  filtered_expr <- cv2_vs_mean_function(df, description, layout, FALSE)
  
  # TSNE PLOT
  library("Rtsne")
  set.seed(1)
  tsne.out_expr <- Rtsne(t(log10(filtered_expr + 1)), initial_dims = 20, verbose = TRUE, perplexity = 30)
  
  plot(tsne.out_expr$Y, main = paste0("tSNE: ", description, " COUNTS"), 
       col = M$Celltype_Color, xlab = "tSNE1", ylab = "tSNE2")
  
  # Get unique cell types and colors (same fix as PCA)
  cell_types <- unique(M$Celltype)
  cell_colors <- unique(M$Celltype_Color)
  
  legend("topright", 
         legend = cell_types, 
         cex = 1, 
         fill = cell_colors, 
         inset = 0.02)
}
```

### Deconvolution Normalization

Until now we have been dealing with two data sets: 1) raw counts, 2)
RPKM counts. Despite RPKM is satisfactory for scRNAseq data in some
situations, it over-corrects the sequencing depth bias by equalizing
library sizes between the cells.

Here we will apply a pooling across cell, or deconvolution,
normalization method developed in the lab of John Marioni (EBI). For
this purpose we need to create a *scran* object and mark spikes inside
the object:

``` r
library(SingleCellExperiment)
sce_deconv <- SingleCellExperiment(assays = list(counts = as.matrix(D)))

library(scater)
qc_stats <- perCellQCMetrics(sce_deconv, subsets = list(MySpikes = rep(c(FALSE, TRUE), c(dim(D)[1] - dim(D_ERCC)[1], dim(D_ERCC)[1]))))
colData(sce_deconv) <- cbind(colData(sce_deconv), qc_stats)

is_spike <- rep(c(FALSE, TRUE), c(dim(D)[1] - dim(D_ERCC)[1], dim(D_ERCC)[1]))
spike_sce <- sce_deconv[is_spike, ]
sce_deconv <- sce_deconv[!is_spike, ]
altExp(sce_deconv, "MySpikes") <- spike_sce

# Extract the spike-in SingleCellExperiment object
spike_sce <- altExp(sce_deconv, "MySpikes")
# Get the count matrix (first 5 rows and columns)
counts(spike_sce)[1:5, 1:5]
```

    ##                                         T74_P1_A9_ILC1 T74_P1_B4_NK
    ## ERCC_0.11444092:mix1_0.11444092:mix2                 0            1
    ## ERCC_0.22888184:mix1_0.22888184:mix2..2              4            2
    ## ERCC_0.22888184:mix1_0.34332275:mix2..1              0            0
    ## ERCC_0.45776367:mix1_0.11444092:mix2..2              0            0
    ## ERCC_0.45776367:mix1_0.45776367:mix2                 0            0
    ##                                         T74_P1_B7_ILC2 T74_P1_B9_NK
    ## ERCC_0.11444092:mix1_0.11444092:mix2                 1            1
    ## ERCC_0.22888184:mix1_0.22888184:mix2..2              4            4
    ## ERCC_0.22888184:mix1_0.34332275:mix2..1              0            0
    ## ERCC_0.45776367:mix1_0.11444092:mix2..2              0            0
    ## ERCC_0.45776367:mix1_0.45776367:mix2                 0            0
    ##                                         T74_P1_D10_ILC2
    ## ERCC_0.11444092:mix1_0.11444092:mix2                  0
    ## ERCC_0.22888184:mix1_0.22888184:mix2..2               1
    ## ERCC_0.22888184:mix1_0.34332275:mix2..1               0
    ## ERCC_0.45776367:mix1_0.11444092:mix2..2               0
    ## ERCC_0.45776367:mix1_0.45776367:mix2                  0

Further, the most optimal way to apply the deconvolution method is when
the pooling across cells is done for each cell sub-population
separately. Therefore we will pre-cluster all the cells with
hierarchical clustering and dynamic tree cut. After this we will proceed
with normalizing endogenous genes with deconvolution method. Spikes will
be normalized separately by traditional library size normalization since
spike counts are independent on the total mRNA content. The library size
normalization is appropriate for spikes since we do no expect any
differentially expressed spikes between the cells.

``` r
sce_deconv <- logNormCounts(sce_deconv)
clusters <- quickCluster(sce_deconv, min.size = 199, d = NA)


sce_deconv <- computeSumFactors(sce_deconv,cluster=clusters)
summary(sizeFactors(sce_deconv))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.1548  0.5034  0.9275  1.0000  1.3812  3.8531

``` r
spike_factors <- librarySizeFactors(altExp(sce_deconv, "MySpikes"))
colData(sce_deconv)$sizeFactor_Spikes <- spike_factors
colnames(colData(sce_deconv))
```

    ## [1] "sum"                       "detected"                 
    ## [3] "subsets_MySpikes_sum"      "subsets_MySpikes_detected"
    ## [5] "subsets_MySpikes_percent"  "total"                    
    ## [7] "sizeFactor"                "sizeFactor_Spikes"

``` r
head(colData(sce_deconv)$sizeFactor)         # Endogenous size factors
```

    ## [1] 0.2610978 2.5647336 1.1807359 1.0575056 0.4481234 0.8592482

``` r
head(colData(sce_deconv)$sizeFactor_Spikes)  # Spike-in size factors
```

    ##  T74_P1_A9_ILC1    T74_P1_B4_NK  T74_P1_B7_ILC2    T74_P1_B9_NK T74_P1_D10_ILC2 
    ##       0.6519482       2.7998218       2.5414232       2.3246439       0.9562505 
    ##  T74_P1_E1_ILC3 
    ##       0.7434053

``` r
par(mfrow=c(1,2))

# 1. Plot Endogenous size factors
hist(colData(sce_deconv)$sizeFactor, 
     xlab = "Endogenous genes deconvolution size factors", 
     main = "")

# 2. Plot Spike-in size factors
hist(colData(sce_deconv)$sizeFactor_Spikes, 
     xlab = "Spikes library size factors", 
     main = "")
```

![](norm_analysis_files/figure-markdown_github/deconvolution-1.png)

``` r
sce_deconv <- logNormCounts(sce_deconv)

# 1. Endogenous genes normalized by deconvolution factors
# Note: Using colData(sce_deconv)$sizeFactor instead of @phenoData@data$size_factor
expr_deconv <- as.data.frame(t(apply(D[grepl("ERCC_", rownames(D)) == FALSE, ], 1, "/", colData(sce_deconv)$sizeFactor)))

# 2. ERCC genes normalized by spike-in factors
# Note: Using colData(sce_deconv)$sizeFactor_Spikes instead of @phenoData@data$size_factor_MySpikes
expr_deconv_ercc <- as.data.frame(t(apply(D_ERCC, 1, "/", colData(sce_deconv)$sizeFactor_Spikes)))

# 3. Merge them
expr_deconv_merge <- rbind(expr_deconv, expr_deconv_ercc)
expr_deconv_merge[1:5,1:5]
```

    ##                 T74_P1_A9_ILC1 T74_P1_B4_NK T74_P1_B7_ILC2 T74_P1_B9_NK
    ## ENSG00000000003       0.000000        0.000      0.0000000    0.0000000
    ## ENSG00000000419       0.000000        0.000      0.0000000    0.0000000
    ## ENSG00000000457       0.000000        0.000      0.0000000    0.0000000
    ## ENSG00000000460       3.829982        0.000      0.8469295    0.0000000
    ## ENSG00000000938       0.000000     2755.452      0.8469295    0.9456215
    ##                 T74_P1_D10_ILC2
    ## ENSG00000000003        0.000000
    ## ENSG00000000419        0.000000
    ## ENSG00000000457      147.280860
    ## ENSG00000000460       89.261128
    ## ENSG00000000938        2.231528

Here *expr_deconv_merge* contains endogenous genes normalized by
deconvolution normalization and spikes normalized by library size
normalization. We will proceed with detection of highly variable genes
by decomposing biological and technical variation, the latter is
determined as the variation of spikes. We will also correct for donor
batch on this step.

``` r
# Detecting highly variable genes
# 1. Define your batch design (same as before)
batch <- M$Donor
design <- model.matrix(~batch)

# 2. Calculate variance and fit the trend (accounts for donor batch)
fit <- modelGeneVar(sce_deconv, design = design)

# 3. Extract the top highly variable genes
top.hvgs <- order(fit$bio, decreasing = TRUE)
head(fit[top.hvgs, ])
```

    ## DataFrame with 6 rows and 6 columns
    ##                      mean     total      tech       bio      p.value
    ##                 <numeric> <numeric> <numeric> <numeric>    <numeric>
    ## ENSG00000188404   6.91576   20.2677 10.115828  10.15190  1.81946e-03
    ## ENSG00000090104   7.11316   18.8717  9.683338   9.18838  2.98485e-03
    ## ENSG00000170345  11.83839   10.0436  0.931944   9.11166 7.50561e-177
    ## ENSG00000118503   8.31849   15.2407  6.929209   8.31151  2.54879e-04
    ## ENSG00000160593   5.64920   20.6866 13.029557   7.65700  4.43043e-02
    ## ENSG00000110848   9.90708   10.1568  3.105617   7.05117  2.37204e-11
    ##                          FDR
    ##                    <numeric>
    ## ENSG00000188404  8.80214e-01
    ## ENSG00000090104  8.80214e-01
    ## ENSG00000170345 1.15459e-172
    ## ENSG00000118503  3.56436e-01
    ## ENSG00000160593  8.80214e-01
    ## ENSG00000110848  6.08152e-08

``` r
# 4. Plot the top 7 highly variable genes
par(mfrow = c(1, 1))
barplot(as.numeric(fit$bio[top.hvgs])[1:7], 
        names.arg = rownames(fit)[top.hvgs][1:7], 
        las = 1, 
        ylim = c(0, max(fit$bio[top.hvgs][1:7]) * 1.2), 
        main = "Highly Variable Genes", 
        cex.names = 0.5, 
        col = "brown", 
        ylab = "BIOLOGICAL VARIATION",las=2)
```

![](norm_analysis_files/figure-markdown_github/highly%20variable%20genes-1.png)

``` r
# 5. Plot the Mean-Variance trend
plot(fit$mean, fit$total, 
     xlab = "Mean log-expression", 
     ylab = "Variance")

# Add the fitted technical trend line (stored in fit$tech)
o <- order(fit$mean)
lines(fit$mean[o], fit$tech[o], col = "red", lwd = 2)
points(fit$mean, fit$tech, col = "red", pch = 16)
```

![](norm_analysis_files/figure-markdown_github/highly%20variable%20genes-2.png)

Further we will rank gene pairs by their par-wise Spearman correlation.
To test significance of the correlation the observed Spearman rho is
compared with permuted Spearman rho (null hupothesis). Again
batch-effects are explicitly taken into account:

``` r
#Detecting correlated genes
#BiocManager::install("ResidualMatrix")

cor.pairs <- correlatePairs(sce_deconv[top.hvgs[1:200], ], design = design)

# View the top correlated pairs
head(cor.pairs)
```

    ## DataFrame with 6 rows and 5 columns
    ##             gene1           gene2       rho      p.value          FDR
    ##       <character>     <character> <numeric>    <numeric>    <numeric>
    ## 1 ENSG00000167286 ENSG00000160654  0.827789 2.67927e-164 5.33175e-160
    ## 2 ENSG00000170345 ENSG00000120129  0.726744 1.86861e-107 1.85926e-103
    ## 3 ENSG00000110002 ENSG00000162594  0.680174  3.46581e-89  2.29899e-85
    ## 4 ENSG00000162594 ENSG00000104951  0.670494  8.28114e-86  4.11987e-82
    ## 5 ENSG00000110002 ENSG00000157404  0.659814  3.16318e-82  1.25895e-78
    ## 6 ENSG00000110002 ENSG00000186074  0.651316  1.76501e-79  5.85396e-76

``` r
# Create the barplot for the top 4 pairs
barplot(abs(as.numeric(cor.pairs$rho))[1:4],
        names.arg = paste(cor.pairs$gene1[1:4], cor.pairs$gene2[1:4], sep = "_"),
        las = 1,
        main = "Correlated Gene Pairs",
        cex.names = 0.6,
        col = "darkgreen",
        ylab = "SPEARMAN RHO")
```

![](norm_analysis_files/figure-markdown_github/correlated%20genes-1.png)

Finally we will perform cell cycle phase assignment using a pre-trained
classifier. The classifier was trained by looking at the pairs of genes
for which the difference in their expression changes its sign from cell
phase to cell phase.

``` r
#ASSIGN CELL CYCLE STATUS
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran")) #launch classifier for human data
assigned <- cyclone(sce_deconv, pairs=hs.pairs) #classification
rownames(assigned$scores)<-colnames(sce_deconv)
head(assigned$scores)
```

    ##                    G1     S   G2M
    ## T74_P1_A9_ILC1  0.182 0.670 0.050
    ## T74_P1_B4_NK    0.237 0.928 0.000
    ## T74_P1_B7_ILC2  0.352 0.330 0.436
    ## T74_P1_B9_NK    0.041 0.353 0.822
    ## T74_P1_D10_ILC2 0.372 0.312 0.233
    ## T74_P1_E1_ILC3  0.870 0.218 0.049

``` r
phase <- rep("S", ncol(sce_deconv))
phase[assigned$scores$G1 > 0.5] <- "G1"
phase[assigned$scores$G2M > 0.5] <- "G2M"
phase[assigned$scores$G1 > 0.5 & assigned$scores$G2M > 0.5] <- "unknown"
table(phase)
```

    ## phase
    ##      G1     G2M       S unknown 
    ##     287      22     334       5

``` r
assigned$scores$phase<-phase
color<-ifelse(phase=="G1","blue","red")
color[phase=="S"]<-"green"
color[phase=="unknown"]<-"black"
assigned$scores$color<-color
head(assigned$scores)
```

    ##                    G1     S   G2M phase color
    ## T74_P1_A9_ILC1  0.182 0.670 0.050     S green
    ## T74_P1_B4_NK    0.237 0.928 0.000     S green
    ## T74_P1_B7_ILC2  0.352 0.330 0.436     S green
    ## T74_P1_B9_NK    0.041 0.353 0.822   G2M   red
    ## T74_P1_D10_ILC2 0.372 0.312 0.233     S green
    ## T74_P1_E1_ILC3  0.870 0.218 0.049    G1  blue

``` r
#PCA PLOT
PC<-prcomp(log10(t(expr_deconv_merge+1)),center=TRUE,scale=FALSE)
plot(PC$x[,1:2],col=color,main=paste0("PCA PLOT: DECONVOLUTION COUNTS"),pch=19)
legend("topleft",c("G1","G2M","S","unknown"),cex=1,fill=c("blue","red","green","black"),inset=0.02)
```

![](norm_analysis_files/figure-markdown_github/cell%20cycle-1.png)

``` r
#tSNE PLOT
set.seed(1)
tsne.out_expr<-Rtsne(t(log10(expr_deconv_merge+1)),initial_dims=20,verbose=TRUE,perplexity=30)
```

    ## Performing PCA
    ## Read the 648 x 20 data matrix successfully!
    ## OpenMP is working. 1 threads.
    ## Using no_dims = 2, perplexity = 30.000000, and theta = 0.500000
    ## Computing input similarities...
    ## Building tree...
    ## Done in 0.02 seconds (sparsity = 0.179393)!
    ## Learning embedding...
    ## Iteration 50: error is 59.335185 (50 iterations in 0.04 seconds)
    ## Iteration 100: error is 57.383215 (50 iterations in 0.03 seconds)
    ## Iteration 150: error is 57.213989 (50 iterations in 0.03 seconds)
    ## Iteration 200: error is 57.137256 (50 iterations in 0.03 seconds)
    ## Iteration 250: error is 57.104999 (50 iterations in 0.03 seconds)
    ## Iteration 300: error is 0.849042 (50 iterations in 0.03 seconds)
    ## Iteration 350: error is 0.778409 (50 iterations in 0.02 seconds)
    ## Iteration 400: error is 0.758987 (50 iterations in 0.02 seconds)
    ## Iteration 450: error is 0.751841 (50 iterations in 0.02 seconds)
    ## Iteration 500: error is 0.744384 (50 iterations in 0.02 seconds)
    ## Iteration 550: error is 0.740466 (50 iterations in 0.02 seconds)
    ## Iteration 600: error is 0.738461 (50 iterations in 0.02 seconds)
    ## Iteration 650: error is 0.736865 (50 iterations in 0.02 seconds)
    ## Iteration 700: error is 0.735535 (50 iterations in 0.02 seconds)
    ## Iteration 750: error is 0.733632 (50 iterations in 0.02 seconds)
    ## Iteration 800: error is 0.731711 (50 iterations in 0.02 seconds)
    ## Iteration 850: error is 0.730194 (50 iterations in 0.02 seconds)
    ## Iteration 900: error is 0.728687 (50 iterations in 0.02 seconds)
    ## Iteration 950: error is 0.727766 (50 iterations in 0.02 seconds)
    ## Iteration 1000: error is 0.727104 (50 iterations in 0.02 seconds)
    ## Fitting performed in 0.52 seconds.

``` r
plot(tsne.out_expr$Y,main="tSNE: DECONVOLUTION COUNTS",col=color,xlab="tSNE1",ylab="tSNE2",pch=19)
legend("topright",c("G1","G2M","S","unknown"),cex=1,fill=c("blue","red","green","black"),inset=0.02)
```

![](norm_analysis_files/figure-markdown_github/cell%20cycle-2.png)

We conclude that there are a lot of cycling cells in the data set.
However they are spread more or less homogeneously over all clusters,
and importantly do not form a separate cluster.

### Spike Library Size Normalization

As we saw in the previous section we can also normalize both endogenous
genes and spikes by spikes library size normalization, i.e. size factors
are calculated assuming that there is no differential expression between
spikes. This method will preserve library size differences between cells
for endogenous genes. Below, spikes are normalized by their own library
sizes, and those same factors are applied to the endogenous genes to
preserve cell-to-cell library size differences.

``` r
### Spike Library Size Normalization

# 1. Create the object and calculate QC metrics
sce_spike <- SingleCellExperiment(assays = list(counts = as.matrix(D)))
qc_spike <- perCellQCMetrics(sce_spike, subsets = list(MySpikes = rep(c(FALSE, TRUE), c(dim(D)[1] - dim(D_ERCC)[1], dim(D_ERCC)[1]))))
colData(sce_spike) <- cbind(colData(sce_spike), qc_spike)

# 2. Split the spikes into an altExp
is_spike <- rep(c(FALSE, TRUE), c(dim(D)[1] - dim(D_ERCC)[1], dim(D_ERCC)[1]))
spike_sce <- sce_spike[is_spike, ]
sce_spike <- sce_spike[!is_spike, ]
altExp(sce_spike, "MySpikes") <- spike_sce

# 3. Verify the spike matrix extraction (Extracts first 5 rows and columns of spike counts)
assay(altExp(sce_spike, "MySpikes"))[1:5, 1:5] 
```

    ##                                         T74_P1_A9_ILC1 T74_P1_B4_NK
    ## ERCC_0.11444092:mix1_0.11444092:mix2                 0            1
    ## ERCC_0.22888184:mix1_0.22888184:mix2..2              4            2
    ## ERCC_0.22888184:mix1_0.34332275:mix2..1              0            0
    ## ERCC_0.45776367:mix1_0.11444092:mix2..2              0            0
    ## ERCC_0.45776367:mix1_0.45776367:mix2                 0            0
    ##                                         T74_P1_B7_ILC2 T74_P1_B9_NK
    ## ERCC_0.11444092:mix1_0.11444092:mix2                 1            1
    ## ERCC_0.22888184:mix1_0.22888184:mix2..2              4            4
    ## ERCC_0.22888184:mix1_0.34332275:mix2..1              0            0
    ## ERCC_0.45776367:mix1_0.11444092:mix2..2              0            0
    ## ERCC_0.45776367:mix1_0.45776367:mix2                 0            0
    ##                                         T74_P1_D10_ILC2
    ## ERCC_0.11444092:mix1_0.11444092:mix2                  0
    ## ERCC_0.22888184:mix1_0.22888184:mix2..2               1
    ## ERCC_0.22888184:mix1_0.34332275:mix2..1               0
    ## ERCC_0.45776367:mix1_0.11444092:mix2..2               0
    ## ERCC_0.45776367:mix1_0.45776367:mix2                  0

``` r
# 4. Calculate spike-in size factors and APPLY them globally (general.use = TRUE)
# We calculate library sizes from the altExp and store them in the main object
spike_factors <- librarySizeFactors(altExp(sce_spike, "MySpikes"))
colData(sce_spike)$sizeFactor <- spike_factors

# 5. Check the size factors
summary(sizeFactors(sce_spike))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.1466  0.5450  0.9735  1.0000  1.3861  3.1264

``` r
head(colData(sce_spike)$sizeFactor) # These are now the spike factors
```

    ##  T74_P1_A9_ILC1    T74_P1_B4_NK  T74_P1_B7_ILC2    T74_P1_B9_NK T74_P1_D10_ILC2 
    ##       0.6519482       2.7998218       2.5414232       2.3246439       0.9562505 
    ##  T74_P1_E1_ILC3 
    ##       0.7434053

``` r
# 6. Normalize the endogenous genes using these spike factors
# Use colData() to extract the actual numbers instead of the string name
sce_spike <- logNormCounts(sce_spike, size.factors = colData(sce_spike)$sizeFactor)

# 7. Extract the normalized (linear) expression matrix
altExp(sce_spike, "MySpikes") <- logNormCounts(altExp(sce_spike, "MySpikes"))
expr_spike_endo <- as.data.frame(2^logcounts(sce_spike))
expr_spike_ercc <- as.data.frame(2^logcounts(altExp(sce_spike, "MySpikes")))
expr_spike <- rbind(expr_spike_endo, expr_spike_ercc)
expr_spike[1:5,1:5]
```

    ##                 T74_P1_A9_ILC1 T74_P1_B4_NK T74_P1_B7_ILC2 T74_P1_B9_NK
    ## ENSG00000000003       1.000000        1.000        1.00000     1.000000
    ## ENSG00000000419       1.000000        1.000        1.00000     1.000000
    ## ENSG00000000457       1.000000        1.000        1.00000     1.000000
    ## ENSG00000000460       2.533864        1.000        1.39348     1.000000
    ## ENSG00000000938       1.000000     2525.089        1.39348     1.430173
    ##                 T74_P1_D10_ILC2
    ## ENSG00000000003        1.000000
    ## ENSG00000000419        1.000000
    ## ENSG00000000457       70.019573
    ## ENSG00000000460       42.830044
    ## ENSG00000000938        2.045751

### Comparing Normalization Methods

Finally, let us put all main plots for different normalization
stratageies together:

``` r
cv2_vs_mean_function(D,"RAW",2,TRUE)
cv2_vs_mean_function(rpkm,"RPKM",1,TRUE)
cv2_vs_mean_function(expr_deconv_merge,"DECONVOLUTION",1,TRUE)
cv2_vs_mean_function(expr_spike,"SPIKE LIB SIZE",1,TRUE)
```

![](norm_analysis_files/figure-markdown_github/compare%20normalizations-1.png)

``` r
pca_function(D,"RAW",2)
pca_function(rpkm,"RPKM",1)
pca_function(expr_deconv_merge,"DECONVOLUTION",1)
pca_function(expr_spike,"SPIKE LIB SIZE",1)
```

![](norm_analysis_files/figure-markdown_github/compare%20normalizations-2.png)

``` r
tsne_function(D,"RAW",2)
```

    ## Performing PCA
    ## Read the 648 x 20 data matrix successfully!
    ## OpenMP is working. 1 threads.
    ## Using no_dims = 2, perplexity = 30.000000, and theta = 0.500000
    ## Computing input similarities...
    ## Building tree...
    ## Done in 0.02 seconds (sparsity = 0.193725)!
    ## Learning embedding...
    ## Iteration 50: error is 59.632217 (50 iterations in 0.04 seconds)
    ## Iteration 100: error is 56.690742 (50 iterations in 0.03 seconds)
    ## Iteration 150: error is 56.437020 (50 iterations in 0.03 seconds)
    ## Iteration 200: error is 56.317623 (50 iterations in 0.02 seconds)
    ## Iteration 250: error is 56.295228 (50 iterations in 0.02 seconds)
    ## Iteration 300: error is 0.790524 (50 iterations in 0.03 seconds)
    ## Iteration 350: error is 0.722773 (50 iterations in 0.02 seconds)
    ## Iteration 400: error is 0.708921 (50 iterations in 0.02 seconds)
    ## Iteration 450: error is 0.699436 (50 iterations in 0.02 seconds)
    ## Iteration 500: error is 0.693221 (50 iterations in 0.02 seconds)
    ## Iteration 550: error is 0.688829 (50 iterations in 0.02 seconds)
    ## Iteration 600: error is 0.685454 (50 iterations in 0.02 seconds)
    ## Iteration 650: error is 0.682069 (50 iterations in 0.02 seconds)
    ## Iteration 700: error is 0.673707 (50 iterations in 0.02 seconds)
    ## Iteration 750: error is 0.670312 (50 iterations in 0.02 seconds)
    ## Iteration 800: error is 0.669339 (50 iterations in 0.02 seconds)
    ## Iteration 850: error is 0.668377 (50 iterations in 0.02 seconds)
    ## Iteration 900: error is 0.667583 (50 iterations in 0.02 seconds)
    ## Iteration 950: error is 0.667410 (50 iterations in 0.02 seconds)
    ## Iteration 1000: error is 0.665854 (50 iterations in 0.02 seconds)
    ## Fitting performed in 0.50 seconds.

``` r
tsne_function(rpkm,"RPKM",1)
```

    ## Performing PCA
    ## Read the 648 x 20 data matrix successfully!
    ## OpenMP is working. 1 threads.
    ## Using no_dims = 2, perplexity = 30.000000, and theta = 0.500000
    ## Computing input similarities...
    ## Building tree...
    ## Done in 0.02 seconds (sparsity = 0.176826)!
    ## Learning embedding...
    ## Iteration 50: error is 63.147194 (50 iterations in 0.04 seconds)
    ## Iteration 100: error is 59.539731 (50 iterations in 0.03 seconds)
    ## Iteration 150: error is 59.266316 (50 iterations in 0.03 seconds)
    ## Iteration 200: error is 59.253743 (50 iterations in 0.03 seconds)
    ## Iteration 250: error is 59.250104 (50 iterations in 0.03 seconds)
    ## Iteration 300: error is 1.011642 (50 iterations in 0.02 seconds)
    ## Iteration 350: error is 0.955185 (50 iterations in 0.02 seconds)
    ## Iteration 400: error is 0.940126 (50 iterations in 0.02 seconds)
    ## Iteration 450: error is 0.932540 (50 iterations in 0.02 seconds)
    ## Iteration 500: error is 0.924758 (50 iterations in 0.02 seconds)
    ## Iteration 550: error is 0.920378 (50 iterations in 0.02 seconds)
    ## Iteration 600: error is 0.916017 (50 iterations in 0.02 seconds)
    ## Iteration 650: error is 0.912934 (50 iterations in 0.02 seconds)
    ## Iteration 700: error is 0.909739 (50 iterations in 0.02 seconds)
    ## Iteration 750: error is 0.908415 (50 iterations in 0.02 seconds)
    ## Iteration 800: error is 0.907452 (50 iterations in 0.02 seconds)
    ## Iteration 850: error is 0.906041 (50 iterations in 0.02 seconds)
    ## Iteration 900: error is 0.905559 (50 iterations in 0.02 seconds)
    ## Iteration 950: error is 0.903892 (50 iterations in 0.02 seconds)
    ## Iteration 1000: error is 0.902754 (50 iterations in 0.02 seconds)
    ## Fitting performed in 0.48 seconds.

``` r
tsne_function(expr_deconv_merge,"DECONVOLUTION",1)
```

    ## Performing PCA
    ## Read the 648 x 20 data matrix successfully!
    ## OpenMP is working. 1 threads.
    ## Using no_dims = 2, perplexity = 30.000000, and theta = 0.500000
    ## Computing input similarities...
    ## Building tree...
    ## Done in 0.02 seconds (sparsity = 0.183523)!
    ## Learning embedding...
    ## Iteration 50: error is 60.608479 (50 iterations in 0.04 seconds)
    ## Iteration 100: error is 57.528961 (50 iterations in 0.03 seconds)
    ## Iteration 150: error is 57.380074 (50 iterations in 0.03 seconds)
    ## Iteration 200: error is 57.346575 (50 iterations in 0.03 seconds)
    ## Iteration 250: error is 57.282000 (50 iterations in 0.03 seconds)
    ## Iteration 300: error is 0.898894 (50 iterations in 0.03 seconds)
    ## Iteration 350: error is 0.833777 (50 iterations in 0.02 seconds)
    ## Iteration 400: error is 0.817838 (50 iterations in 0.02 seconds)
    ## Iteration 450: error is 0.809622 (50 iterations in 0.02 seconds)
    ## Iteration 500: error is 0.803634 (50 iterations in 0.02 seconds)
    ## Iteration 550: error is 0.801737 (50 iterations in 0.02 seconds)
    ## Iteration 600: error is 0.799093 (50 iterations in 0.02 seconds)
    ## Iteration 650: error is 0.796596 (50 iterations in 0.02 seconds)
    ## Iteration 700: error is 0.794356 (50 iterations in 0.02 seconds)
    ## Iteration 750: error is 0.792353 (50 iterations in 0.02 seconds)
    ## Iteration 800: error is 0.790805 (50 iterations in 0.02 seconds)
    ## Iteration 850: error is 0.788944 (50 iterations in 0.02 seconds)
    ## Iteration 900: error is 0.788625 (50 iterations in 0.02 seconds)
    ## Iteration 950: error is 0.787893 (50 iterations in 0.02 seconds)
    ## Iteration 1000: error is 0.787305 (50 iterations in 0.02 seconds)
    ## Fitting performed in 0.49 seconds.

``` r
tsne_function(expr_spike,"SPIKE LIB SIZE",1)
```

    ## Performing PCA
    ## Read the 648 x 20 data matrix successfully!
    ## OpenMP is working. 1 threads.
    ## Using no_dims = 2, perplexity = 30.000000, and theta = 0.500000
    ## Computing input similarities...
    ## Building tree...
    ## Done in 0.02 seconds (sparsity = 0.188910)!
    ## Learning embedding...
    ## Iteration 50: error is 60.881411 (50 iterations in 0.04 seconds)
    ## Iteration 100: error is 56.785143 (50 iterations in 0.03 seconds)
    ## Iteration 150: error is 56.314284 (50 iterations in 0.03 seconds)
    ## Iteration 200: error is 56.267560 (50 iterations in 0.03 seconds)
    ## Iteration 250: error is 56.288395 (50 iterations in 0.03 seconds)
    ## Iteration 300: error is 0.823930 (50 iterations in 0.03 seconds)
    ## Iteration 350: error is 0.751203 (50 iterations in 0.02 seconds)
    ## Iteration 400: error is 0.734555 (50 iterations in 0.02 seconds)
    ## Iteration 450: error is 0.726368 (50 iterations in 0.02 seconds)
    ## Iteration 500: error is 0.720548 (50 iterations in 0.02 seconds)
    ## Iteration 550: error is 0.717255 (50 iterations in 0.02 seconds)
    ## Iteration 600: error is 0.714131 (50 iterations in 0.02 seconds)
    ## Iteration 650: error is 0.711349 (50 iterations in 0.02 seconds)
    ## Iteration 700: error is 0.709299 (50 iterations in 0.02 seconds)
    ## Iteration 750: error is 0.707931 (50 iterations in 0.02 seconds)
    ## Iteration 800: error is 0.707534 (50 iterations in 0.02 seconds)
    ## Iteration 850: error is 0.706738 (50 iterations in 0.02 seconds)
    ## Iteration 900: error is 0.705620 (50 iterations in 0.02 seconds)
    ## Iteration 950: error is 0.706150 (50 iterations in 0.02 seconds)
    ## Iteration 1000: error is 0.704617 (50 iterations in 0.02 seconds)
    ## Fitting performed in 0.51 seconds.

![](norm_analysis_files/figure-markdown_github/compare%20normalizations-3.png)

### Compare Size Factors for Deconvolution with TMM, DESEQ and RPKM

Last thing we will do, we will compare scRNAseq specific normalization
such as deconvolution with main bulk RNAseq normalization startegies.
For this purpose we will plot deconvolution size factors against TMM,
DEseq and RPKM size factors:

``` r
# DECONVOLUTION NORMALIZATION
# Assuming 'sce_deconv' is your SingleCellExperiment object and 'clusters' are defined
sce_deconv <- computeSumFactors(sce_deconv, clusters = clusters)
out.facs <- sizeFactors(sce_deconv)

# DESEQ2 NORMALIZATION
library(DESeq2)
# Create DESeq2 object (condition is a dummy variable here)
dds <- DESeqDataSetFromMatrix(countData = D, colData = data.frame(condition = rep("1", ncol(D))), design = ~ 1)
dds <- estimateSizeFactors(dds)
head(sizeFactors(dds))
```

    ##  T74_P1_A9_ILC1    T74_P1_B4_NK  T74_P1_B7_ILC2    T74_P1_B9_NK T74_P1_D10_ILC2 
    ##       0.5879694       3.0117471       2.7712479       1.7430712       0.6802755 
    ##  T74_P1_E1_ILC3 
    ##       1.0773749

``` r
# TMM NORMALIZATION
library(edgeR)
cds_TMM <- DGEList(counts = D)
cds_TMM <- calcNormFactors(cds_TMM)

# CENTER SIZE FACTORS AROUND MEDIAN AND PLOT
# Extracting the effective library sizes for TMM (norm.factors * lib.size)
tmm_lib_sizes <- cds_TMM$samples$norm.factors * cds_TMM$samples$lib.size

# Centering all factors
med_out.facs <- as.numeric(out.facs) / median(as.numeric(out.facs))
med_sizeFactors_DESeq <- sizeFactors(dds) / median(sizeFactors(dds))
med_sizeFactors_TMM <- tmm_lib_sizes / median(tmm_lib_sizes)
med_sizeFactors_RPKM <- cds_TMM$samples$lib.size / median(cds_TMM$samples$lib.size)

# PLOTTING
par(mfrow = c(2, 2))

plot(log10(med_sizeFactors_TMM) ~ log10(med_sizeFactors_DESeq), 
     xlim = c(-1, 1), ylim = c(-1, 1), 
     xlab = "DESEQ", ylab = "TMM", main = "TMM vs DESeq")
abline(a = 0, b = 1, col = "red")

plot(log10(med_out.facs) ~ log10(med_sizeFactors_DESeq), 
     xlim = c(-1, 1), ylim = c(-1, 1), 
     xlab = "DESEQ", ylab = "DECONVOLUTION", main = "Deconvolution vs DESeq")
abline(a = 0, b = 1, col = "red")

plot(log10(med_out.facs) ~ log10(med_sizeFactors_TMM), 
     xlim = c(-1, 1), ylim = c(-1, 1), 
     xlab = "TMM", ylab = "DECONVOLUTION", main = "Deconvolution vs TMM")
abline(a = 0, b = 1, col = "red")

plot(log10(med_out.facs) ~ log10(med_sizeFactors_RPKM), 
     xlim = c(-1, 1), ylim = c(-1, 1), 
     xlab = "RPKM", ylab = "DECONVOLUTION", main = "Deconvolution vs RPKM")
abline(a = 0, b = 1, col = "red")
```

![](norm_analysis_files/figure-markdown_github/compare%20size%20factors-1.png)

Finally here is the details on the system on which this document was
compiled:

``` r
sessionInfo()
```

    ## R version 4.3.3 (2024-02-29)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Riga
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] tweeDEseq_1.48.0            DESeq2_1.42.1              
    ##  [3] limSolve_2.0.3              Rtsne_0.17                 
    ##  [5] RColorBrewer_1.1-3          mixOmics_6.26.0            
    ##  [7] lattice_0.22-5              MASS_7.3-60.0.1            
    ##  [9] edgeR_4.0.16                limma_3.58.1               
    ## [11] scater_1.30.1               ggplot2_4.0.3              
    ## [13] scran_1.30.2                scuttle_1.12.0             
    ## [15] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
    ## [17] Biobase_2.62.0              GenomicRanges_1.54.1       
    ## [19] GenomeInfoDb_1.38.8         IRanges_2.36.0             
    ## [21] S4Vectors_0.40.2            BiocGenerics_0.48.1        
    ## [23] MatrixGenerics_1.14.0       matrixStats_1.5.0          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.1-0              gridExtra_2.3.1          
    ##  [3] rlang_1.3.0               magrittr_2.0.5           
    ##  [5] otel_0.2.0                compiler_4.3.3           
    ##  [7] DelayedMatrixStats_1.24.0 vctrs_0.7.3              
    ##  [9] reshape2_1.4.5            quadprog_1.5-8           
    ## [11] stringr_1.6.0             pkgconfig_2.0.3          
    ## [13] crayon_1.5.3              fastmap_1.2.0            
    ## [15] XVector_0.42.0            rmarkdown_2.32           
    ## [17] ggbeeswarm_0.7.3          purrr_1.2.2              
    ## [19] xfun_0.62                 bluster_1.12.0           
    ## [21] zlibbioc_1.48.2           beachmat_2.18.1          
    ## [23] DelayedArray_0.28.0       BiocParallel_1.36.0      
    ## [25] irlba_2.4.1               parallel_4.3.3           
    ## [27] cluster_2.1.6             R6_2.6.1                 
    ## [29] stringi_1.8.9             Rcpp_1.1.2               
    ## [31] knitr_1.52                splines_4.3.3            
    ## [33] Matrix_1.6-5              igraph_2.3.4             
    ## [35] tidyselect_1.2.1          rstudioapi_0.19.0        
    ## [37] abind_1.4-8               yaml_2.3.12              
    ## [39] viridis_0.6.5             codetools_0.2-19         
    ## [41] tibble_3.3.1              plyr_1.8.9               
    ## [43] withr_3.0.3               rARPACK_0.11-0           
    ## [45] S7_0.2.2                  evaluate_1.0.5           
    ## [47] lpSolve_5.6.23            pillar_1.11.1            
    ## [49] ellipse_0.5.0             generics_0.1.4           
    ## [51] RCurl_1.98-1.20           cqn_1.48.0               
    ## [53] sparseMatrixStats_1.14.0  scales_1.4.0             
    ## [55] glue_1.8.1                metapod_1.10.1           
    ## [57] tools_4.3.3               BiocNeighbors_1.20.2     
    ## [59] ScaledMatrix_1.10.0       RSpectra_0.16-2          
    ## [61] locfit_1.5-9.12           grid_4.3.3               
    ## [63] tidyr_1.3.2               GenomeInfoDbData_1.2.11  
    ## [65] nor1mix_1.3-3             beeswarm_0.4.0           
    ## [67] BiocSingular_1.18.0       vipor_0.4.7              
    ## [69] cli_3.6.6                 rsvd_1.0.6               
    ## [71] S4Arrays_1.2.1            viridisLite_0.4.3        
    ## [73] dplyr_1.2.1               ResidualMatrix_1.12.0    
    ## [75] corpcor_1.6.10            gtable_0.3.6             
    ## [77] digest_0.6.39             SparseArray_1.2.4        
    ## [79] ggrepel_0.9.8             dqrng_0.4.1              
    ## [81] farver_2.1.2              htmltools_0.5.9          
    ## [83] lifecycle_1.0.5           statmod_1.5.2
