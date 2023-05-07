#!/usr/bin/Rscript
suppressPackageStartupMessages({
  library("data.table")
  library(DESeq2)
  library(ggplot2)
  library(stringr)
  library(dplyr)
  library(ggpubr)
  library(ggthemes)
  library(tximport)
  library(ggrepel)
  library("tibble")
  library("optparse")
  library("BiocParallel")
  
})

option_list <- list(
  make_option(c("-i", "--inputfile"), dest = "Rinput", default = "", 
              help="[opt] inputdir"),
  make_option(c("-n", "--namefle"), dest = "namefile", default = "",
              help="[opt] namefile"),
  make_option(c("-c","--t2g"),dest = "t2gfile",default = "",
              help = "[opt] t2g config file"),
  make_option(c("-t","--contrast"),dest = "contrastfile",default = "",
              help = "[opt] contrast config file"),
  make_option(c("-p", "--padj"), dest = "p", default = "",
              help="[opt] padj cutoff "),
  make_option(c("-f", "--foldchange"), dest = "f", default = "0",
              help="[opt] log2foldchange defalut 0"),
  make_option(c("-x", "--threads"), dest = "x", default = "4",
              help="[opt] threads defalut 4"),
  make_option(c("-o", "--outprefix"), dest = "outfile", default = "",
              help="[opt] DEseq2 out file")
) 

parser <- OptionParser(usage = "%prog [options] file",
     option_list=option_list, description = " 2019-3-24 \
Description: \
  Wrapper for DEseq2.\
Example: \
  DEseqForKallisto.R -i 01.kallisto   -n design.conf -c t2g -t contrast.conf -x 1  -o out "
)
#



arguments <- parse_args(parser)
opt <- arguments$options
test <- arguments$Rinput
files = file.path(test ,list.files(test),"abundance.h5")

names(files) = list.files(test)

register(MulticoreParam(arguments$x))

if(test == "") {
  print_help(parser)
  stop(sprintf("input file is not specified"))
  test = file("stdin")  
}

designfile <- arguments$namefile

designfile =  read.table(designfile,header=T,stringsAsFactors = F)

t2gfile <- arguments$t2gfile

tx2gene <- read.table(t2gfile,header=T)



contrastfile <- arguments$contrastfile

contrastfile = read.csv(contrastfile,sep="\t",header = FALSE,stringsAsFactors = F)



P <- as.numeric(arguments$p)
F <- as.numeric(arguments$f)
o <- arguments$outfile



txi = tximport(files, type = "kallisto", tx2gene = tx2gene, 
               txIn = TRUE, txOut = FALSE, countsFromAbundance = "no")


if(  unique(colnames(txi$abundance) == designfile$sample) ) {
  print("all  samples matched")
} else {
  print( "rename file from :" )
  colnames(txi$abundance)
  print("to :")
  designfile$sample
  }

dds <- DESeqDataSetFromTximport(txi,colData = designfile,design= ~ condition)

dds2 <- DESeq(dds,parallel = TRUE)

mergetable <- data.frame(row.names =row.names(dds2))

for (i in 1:nrow(contrastfile)) {
    contrast_name <- paste0(as.character(unlist(contrastfile[i,]))[2:3],collapse = "vs")
    contrast_vector <- as.vector(unlist(contrastfile[i,]))
    results <- results(dds2, contrast = contrast_vector,parallel = T)
    results <- lfcShrink(dds2, contrast = contrast_vector, res=results,parallel = T)
    if( arguments$p  == "") {
        diff <- subset(results,(log2FoldChange > F|log2FoldChange < -F)) %>% as.data.frame()
    }else {
    diff <- subset(results, padj<P & (log2FoldChange > F|log2FoldChange < -F)) %>% as.data.frame()
    }
    diff2 <- subset(results) %>% as.data.frame()
    diff2$logQ <- -log10(diff2$padj)
    diff2$Group <- "Not-changed"
    diff2$Label <- ""
    diff2 <- diff2[order(diff2$padj),] 
    diff2$Group[which( (diff2$padj<0.05) & (diff2$log2FoldChange > F))] = "Up-regulated"
    diff2$Group[which( (diff2$padj<0.05) & (diff2$log2FoldChange < -F) )] = "Down-regulated"
    upgenes <- head(rownames(diff2[which(diff2$Group =="Up-regulated"),]),10)
    downgenes <- head(rownames(diff2[which(diff2$Group =="Down-regulated"),]),10)
    top10genes <- c(as.character(upgenes),as.character(downgenes))
    diff2$Label[match(top10genes,rownames(diff2)) ] <- top10genes
    Vplot <- ggscatter(diff2,x = "log2FoldChange",y = "logQ", color = "Group" , palette=c("#2f5688","#BBBBBB","#CC0000"),size=1,label =diff2$Label,font.label = 8,repel = T,xlab = "Log2FC",ylab="-Log10(padj)") + theme_base() +
      geom_hline(yintercept = -log10(P), linetype="dashed")+
      geom_vline(xintercept = -c(-F, F), linetype="dashed")
    ggsave(paste0("./",contrast_name,"vp.pdf"),plot = Vplot)
    colnames(diff) <-  paste(contrast_name,colnames(diff),sep = "_")
    mergetable <- merge(mergetable,diff[,c(2,6)],by="row.names",sort=FALSE,all=TRUE)
    rownames(mergetable) <- mergetable$Row.names
    mergetable$Row.names <- NULL
    final <- merge(diff,as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
    outfile <- paste(o,contrast_name,sep = ".")
    write.table(final,file = outfile,row.names = FALSE,sep = "\t",quote = FALSE)
}


mergetable <- merge(mergetable,counts(dds2,normalize=TRUE),by="row.names",sort=FALSE,all=TRUE)



write.table(mergetable,file = paste0(o,"all"),row.names = FALSE,sep = "\t",quote = FALSE)

