# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Sat Oct 19 09:42:40 EDT 2019

################################################################
#   Differential expression analysis with limma
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install('biomaRt')

library(Biobase)
library(GEOquery)
library(limma)
library(biomaRt)

# load series and platform data from GEO

gset <- getGEO("GSE7621", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "0000000001111111111111111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf) #Benjamini Hochberg multiple testing correction

topTableCtrlVsPD <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

table <- write.table(topTableCtrlVsPD, file=stdout(), row.names=F, sep="\t")
########################################AKT1
# gene list from initial exploratory analysis DISEASES
geneSymbolsShared <- c("AKT1","APOE","IL1B","MAPK8","SLC30A10","SREBF1","TNF","VPS13C")
#gene list from diffusionAPP on shared PPI STRING network
geneSymbolsDiffusion <-c("IL1B",
  "ADAMTS9",
  "KIF11",
  "CXCR4",
  "CACNG2",
  "KLHL42",
  "NTN1",
  "PRKD1",
  "SPTA1",
  "ANK1",
  "GABRG3",
  "HNF4A",
  "RHOU",
  "THSD4",
  "TMEM163",
  "IDE",
  "BNC2",
  "IGF1",
  "LAMA1",
  "LPL",
  "ADIPOQ",
  "ZBTB16",
  "CADPS",
  "MAP3K1",
  "AP3B1",
  "APOB",
  "TGFBR2",
  "CDK1",
  "VEGFA",
  "SCARB2",
  "SYK",
  "VAV2",
  "COL4A1",
  "SYT1")
geneSymbolsTopTable <- topTableCtrlVsPD$Gene.symbol
matches <- is.element(geneSymbolsDiffusion,geneSymbolsTopTable)
diffusedAndDEG <- topTableCtrlVsPD[topTableCtrlVsPD$Gene.symbol %in% geneSymbolsDiffusion,]
#unable to get this code running - difficulty defining "mart = "
# add more stable and Cytoscape ready IDs to table
#useMart(biomart = "ENSEMBL_MART_ENSEMBL")
#diffusedAndDEG <- getBM(filters = "ensembl_peptide_id", 
#                attributes = c("ensembl_peptide_id", "entrezgene", "description"),
#                values = genes, mart = mart)
#export data to file
write.table(diffusedAndDEG, "Gene List T2DM diffused with DE values", sep="\t")
################################################################
#   Boxplot for selected GEO samples
#################################################
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE7621", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- "0000000001111111111111111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  #set group names

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("control","PD")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE7621", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")

