# You will have to analyze the RNA-seq data presented in: Henn, A. D. et al. High-resolution temporal response patterns to influenza vaccine reveal a distinct human plasma cell gene signature. Scientific Reports 3, 2327 (2013).
# Use voom and limma to find genes that are differentially expressed at each time point compared to baseline (day 0). Use an FDR cutoff of 0.01. Display your results using pheatmap showing the log fold-change of the differentially expressed genes grouped by time point.
# Perform a GSEA analysis using camera and the MSigDB Reactome pathway gene signatures. Display your results using pheatmap, again group by timepoint. This is similar to what we've done in class.


#' Creates the R markdown files.
# library("knitr")
# opts_knit$set(progress = FALSE, verbose = FALSE, message=FALSE)
# spin(hair = "CobbMareaHW3.R", format = "Rmd")
# file.rename("CobbMareaHW3.md", "CobbMareaHW3.Rmd")


#Loads bioconductor packages
source("http://bioconductor.org/biocLite.R")
# biocLite()
library(limma)
library(BiocInstaller)
biocLite("GSEABase")
library(GSEABase)

#Loads GEO libraries and specific dataset
library(GEOmetadb)
library(GEOquery)

# Download database only if it's not done already
if (!file.exists("GEOmetadb.sqlite")) {
  getSQLiteFile()
}

#http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1002366
gse <- getGEO("GSE45735", GSEMatrix=TRUE, destdir = "Data/GEO")
pd <- pData(gse[[1]])
getGEOSuppFiles("GSE45735", makeDirectory=FALSE, baseDir = "Data/GEO/")


# Note the regular expression to grep file names
files <- list.files(path = "Data/GEO/", pattern = "GSE45735_T.*.gz", full.names = TRUE)

# Read in gzip-compressed, tab-delimited files
file_list <- lapply(files, read.table, sep='\t', header=TRUE)

# Subset to only those rows where Gene contains only non-space characters
# This addresses problems with T14 file containing 28 invalid rows at end of file
file_list <- lapply(file_list, function(file_list)subset(file_list, grepl('^[^[:space:]]+$', Gene)))

# Remove duplicated rows
file_list_unique <- lapply(file_list, function(x){x<-x[!duplicated(x$Gene),]; 
                                                  x <- x[order(x$Gene),]; 
                                                  rownames(x) <- x$Gene;
                                                  x[,-1]})
# Take the intersection of all genes
gene_list <- Reduce(intersect, lapply(file_list_unique, rownames))
file_list_unique <- lapply(file_list_unique, "[", gene_list,)
matrix <- as.matrix(do.call(cbind, file_list_unique))

# Clean up the pData
pd_small <- pd[!grepl("T13_Day8",pd$title),]
pd_small$Day <- sapply(strsplit(gsub(" \\[PBMC\\]", "", pd_small$title),"_"),"[",2)
pd_small$subject <- sapply(strsplit(gsub(" \\[PBMC\\]", "", pd_small$title),"_"),"[",1)
colnames(matrix) <- rownames(pd_small)


# Note that I add one to the count
new_set <- ExpressionSet(assayData = matrix+1)
pData(new_set) <- pd_small

design <- model.matrix(~subject+Day, new_set)
new_set_voom <- voom(new_set,design = design)

lm <- lmFit(new_set_voom, design)
eb <- eBayes(lm)

cont_matrix <- makeContrasts(timeD7-timeD3,levels=mm_TIV_08)
fit2 <- contrasts.fit(fit_TIV_08, cont_matrix)
fit2 <- eBayes(fit2)
topTable(fit2, adjust = "fdr")

# Look at the other time-points
# topTable(eb, coef = "DayDay1", number = 5)

# ordinary_t <- lm$coef / lm$stdev.unscaled / lm$sigma
# head(ordinary_t)
# ordinary_t <- ordinary_t[,"DayDay1"]
# 
# head(ordinary_t)
# 
# ordinary_p = p.adjust(2*pnorm(abs(ordinary_t), lower.tail=FALSE), method="BH")



c2_set <- getGmt("GSEA-sets/c2.all.v4.0.symbols.gmt")
gene_ids <- geneIds(c2_set)

# Camera requires gene-indices. 
# Which function to use will depend on which version of limma you have.
#     http://bioconductor.org/packages/release/bioc/news/limma/NEWS
#     "symbols2indices() renamed to ids2indices()."
if (exists("ids2indices")) { 
  sets_indices <- ids2indices(gene_ids, rownames(new_set))
}
if (exists("symbols2indices")) {
  sets_indices <- symbols2indices(gene_ids, rownames(new_set))    
}


res <- vector("list",length = 10)
for(i in 1:10)
{
  contrast <- paste0("DayDay",i)
  cont_matrix <- makeContrasts(contrast, levels=design)
  res[[i]] <- camera(new_set_voom, sets_indices, design=design, contrast=cont_matrix, sort=FALSE)
}


library(pheatmap)
PValue <- sapply(res, function(x){ifelse(x$Direction=="Up", -10*log10(x$PValue), 10*log10(x$PValue))})
rownames(PValue) <- rownames(res[[1]])
PValue_max <- rowMax(abs(PValue))
PValue_small <- PValue[PValue_max>30, ]
anno <- data.frame(Time=paste0("Day",1:10))
rownames(anno) <- colnames(PValue_small)  <- paste0("Day",1:10)

pheatmap(PValue_small, cluster_cols=FALSE,fontsize_row = 5)

