#!/bin/Rscript
packages <- c("stringr","dplyr","ggplot2","cowplot","pheatmap","reshape","tibble","ggpubr","plotly","RColorBrewer","optparse")
invisible(lapply(packages, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE))))
#lapply(packages, function(x)suppressPackageStartupMessages(x))

option_list <- list(
  make_option(c("-m","--martix"), dest = "martix", default = "",
              help="   "),
  make_option(c("-t","--title"), dest = "title", default = "",
              help="   "),
  make_option(c("-a","--ann"), dest = "ann", default = "",
              help="    "),
  make_option(c("-s","--subsample"), dest = "subsample", default = "",
              help="   "),
  make_option(c("-b","--breaks"), dest = "breaks", default = "",type = "double",
              help=" "),
  make_option(c("-o", "--out"), dest = "out", default = ".",
              help="[opt] out folder")
  
)
#
parser <- OptionParser(usage = "%prog [options] file",
                       option_list=option_list, description = " 2022-01-28 \
Description: \
  plot gene reads depth \
Example: \
   /usr/bin/Rscript plotcoverage.R  -m ./TraesCS..._2000_coverage.martix  -a  331.ann -s subID -b 2000 -o PIF4  \ "
)
arguments <- parse_args(parser)
opt <- arguments$options

martix <- arguments$martix  
ann <- arguments$ann  
subsample  <- arguments$subsample 
breaks  <- arguments$breaks   
out <- arguments$out 
title <- arguments$title 


infile= read.table(martix,check.names = F,header = T)
infile  <- t(infile)
colnames(infile) <- seq(1,ncol(infile)) 
infile <- infile[-1,]

infile <- as.data.frame(infile)

anncol =  read.table(ann,check.names = F,header = T)
anncol <- arrange(anncol,DEL)
anncol=column_to_rownames(anncol,var = "gvcfid")

ann_colors = list(
  lumpy = c(L = "orange", S = "red2",N = "slateblue3",unknow="white"),
  DEL = c(L = "orange", S = "red2",N = "slateblue3"),
  indel = c(Ref = "slateblue3", Alt = "red2",hetero = "grey60",unknow="white"),
  indel_beagle = c(Ref = "slateblue3", Alt = "red2",hetero = "grey60"),
  hap = c(Ref = "black", Alt = "grey"),
  GROUP =  c(ASC = "#E83A38", ASL = "#F19595",EUM = "#6AB82D",EXM ="#304EA0")
)

dim(infile)

pheatmap(log2(infile+1),cluster_cols=FALSE,
         cluster_rows=T,display_numbers=F,legend = T,
         annotation_row=anncol,gaps_col = c(breaks,dim(infile)[2]-breaks), 
         annotation_colors = ann_colors,show_colnames = F,main = title,
         annotation_names_col=FALSE,angle_col ="45",filename = paste0(out,"all_sample.pdf") )

if (subsample != "" ) {
  subsample <- read.table(subsample)
  as.vector(subsample$V1)
  infile3 <- dplyr::filter(infile,row.names(infile) %in% subsample$V1 )
  
  pheatmap(log2(infile3+1),cluster_cols=FALSE,
           cluster_rows=T,display_numbers=F,legend = T,
           annotation_row=anncol,gaps_col = c(breaks,dim(infile)[2]-breaks), 
           annotation_colors = ann_colors, show_colnames = F,main = title,
           annotation_names_col=FALSE,angle_col ="45",filename = paste0(out,"sub_sample.pdf") )
  
}


