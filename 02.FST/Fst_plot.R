#!/usr/bin/Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(dplyr))
load("~/R/theme_qz.RData")

# Arguments
option_list <- list(
  make_option(c("-f", "--fst"), dest = "fst", default = "",
              help="input fst file"),
  make_option(c("-p", "--piinfo"), dest = "piinfo", default = "",
              help="input piinfo file"),
  make_option(c("-t", "--tajimasD"), dest = "tajimasD", default = "",
              help="input tajimasD file"),
  make_option(c("-x", "--xpclr"), dest = "xplcr", default = "",
              help="input xplcr file"),
  make_option(c("-a", "--ann"), dest = "ann", default = "",
              help="input ann file"),
  make_option(c("-w", "--width"), dest = "width", type = "double",
              help="plot width"),
  make_option(c("-l", "--left"), dest = "left",type = "double",default = 0,
              help="left grey"),
  make_option(c("-r", "--right"), dest = "right", type = "double",default = 0,
              help="right grey"),
  make_option(c("-v", "--height"), dest = "height", type = "double",
              help="plot height"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "out",
              help = " output plot name prefix. ")
)
parser <- OptionParser(usage = "Fst_plot.R -f chr5B_1000k.fst.windowed.weir.fst -p C_L.1000k.windowed.pi -t chr5B_1000k.TD -x ~/qinz/caojie/FST/oldVCF/chinasample/XPCLR/1000k/chr5B.pyxpclr -a FstAnnForRplot.txt -w 8 -v 6 -l 546 -r 576 -o Fst_plot",
                       description = 
                         "qinz001@qq.com",
                       option_list = option_list)

arguments <- parse_args(parser)


left_x <- arguments$left
right_x <- arguments$right
outfile <- arguments$outfile

width_p <- arguments$width
height_p <- arguments$height

FSTinfo <- arguments$fst
piinfo <- arguments$piinfo
TDinfo <- arguments$tajimasD
XPCLR <- arguments$xplcr
Ann <- arguments$ann

# 1Mb bin 
#FSTinfo <- read.table("~/qinz/caojie/FST/filterVCF/1000k/chr5B_1000k.fst.windowed.weir.fst",header =T)
#FSTinfo <- read.table("~/qinz/ProfessorXIN/GenomeScan/fst/chr1A_100k.fst.windowed.weir.fst",header =T)

FSTinfo = read.table(FSTinfo,header =T)
FSTinfo[FSTinfo$WEIGHTED_FST < 0 ,5] = 0 

piinfo = read.table(piinfo,header =F)
#piinfo <- read.table("~/qinz/caojie/FST/filterVCF/1000k/C_L.1000k.windowed.pi",header =F)
#piinfo <- read.table("~/qinz/ProfessorXIN/GenomeScan/fst/chr1A_20svs80s.100k.windowed.pi",header =F)

colnames(piinfo) = c("c","p","p2","PC","PL","C/L","L/C")

#TDinfo <- read.table("~/qinz/caojie/FST/filterVCF/1000k/chr5B_1000k.TD",header =T)
#TDinfo <- read.table("~/qinz/ProfessorXIN/GenomeScan/fst/chr1A_20svs80s_100k.TD",header =T)

TDinfo = read.table(TDinfo,header =T)

colnames(TDinfo) = c("c","p","s","TC","TL")
TDinfo$p  = TDinfo$p +1 


# Ann <- read.table("~/qinz/caojie/FST/filterVCF/FstAnnForRplot.txt",header =T); VRN-B1 TraesCS5B02G396600 ;PHYC TraesCS5B02G396200 
# Ann <- read.table("~/qinz/ProfessorXIN/GenomeScan/fst/chr1A_FstAnnForRplot.txt",header =T) 

Ann <- read.table(Ann,header =T)


if (XPCLR!="") {
  #~/qinz/caojie/FST/oldVCF/chinasample/XPCLR/1000k/chr5B.pyxpclr
  XPCLRinfo <- read.table(XPCLR,sep="\t",header =T)
  XPCLRinfo <- XPCLRinfo[,c(2,3,4,12)] 
  colnames(XPCLRinfo) = c("c","p","s","XPCLR")
  
  #XPCLRinfo$p  = XPCLRinfo$p +1 
  allinfo = base::merge(FSTinfo,piinfo,by.x=2,by.y=2) %>% 
    dplyr::select(BIN_START,WEIGHTED_FST,PC,PL) %>% 
    base::merge(.,TDinfo,by.x=1,by.y=2) %>% 
    base::merge(.,XPCLRinfo,by.x=1,by.y=2) %>% 
    dplyr::select(BIN_START,WEIGHTED_FST,PC,PL,TC,TL,XPCLR) %>%
    #dplyr::select(BIN_START,WEIGHTED_FST,PC,PL,TC,TL) %>%
    na.omit() %>% 
    reshape2::melt(id.vars="BIN_START")
  
} else  {
  allinfo = base::merge(FSTinfo,piinfo,by.x=2,by.y=2) %>% 
    dplyr::select(BIN_START,WEIGHTED_FST,PC,PL) %>% 
    base::merge(.,TDinfo,by.x=1,by.y=2) %>% 
    #base::merge(.,XPCLRinfo,by.x=1,by.y=2) %>% 
    #dplyr::select(BIN_START,WEIGHTED_FST,PC,PL,TC,TL,XPCLR) %>%
    dplyr::select(BIN_START,WEIGHTED_FST,PC,PL,TC,TL) %>%
    na.omit() %>% 
    reshape2::melt(id.vars="BIN_START")
  
}

allinfo2 <- dplyr::mutate(allinfo, group = case_when( (variable == "WEIGHTED_FST") ~  "FST",
                                               (variable == "X0") ~  "pi-ratio(L/C)",
                                               (variable == "TL") ~  "tajima's D",
                                               (variable == "TC") ~  "tajima's D",
                                               (variable == "PC") ~  "pi",
                                               (variable == "PL") ~  "pi",
                                               (variable == "XPCLR") ~  "XPCLR")) 

Fst_cutoff <-  dplyr::arrange(FSTinfo,desc(WEIGHTED_FST))[round(nrow(FSTinfo) * 0.05,0),5]
print(Fst_cutoff)
p1 <- ggplot(dplyr::filter(allinfo2),
             aes(x=BIN_START/1000000,
                 y=value,
                 group = variable)) + 
  geom_line(aes(group = variable,color = variable ),size =0.5) +
  scale_color_manual(values=c("black","#ED1C28","#3C499F","#ED1C28","#3C499F","black"))+
  facet_wrap(~ group,ncol = 1,scales = "free") +
  xlab(unique(TDinfo$c)) +
  theme_tufte() +
  theme(axis.line=element_line(),legend.position = "bottom")

## Ann FST and grey area
Fst_cutoff_Df <- data.frame(group = c("FST"), Z = Fst_cutoff)

p2 <- p1 + geom_hline(data = Fst_cutoff_Df, 
                aes(yintercept = Z),
                linetype="dashed",
                color="#8B1D22") 


if (left_x != 0 && right_x != 0) {
  p2 <- p1 + geom_rect(aes(xmin=left_x,xmax=right_x,ymin=-Inf,ymax=Inf),fill="grey",alpha=0.01)
}

## Ann line and geneID 
for (i in seq(1,nrow(Ann)) ) {
  
  p2 <- p2 + geom_vline(xintercept = Ann[i,1],
                        size=0.7,
                        linetype="dashed",
                        color=Ann[i,3]) +
    annotate("text",label = Ann[i,2],
             x = Ann[i,1],y=Inf,vjust = 1,hjust = 1, angle = 90,
             size =3, colour = "red")
  
}

p2 <- p2+
  scale_x_continuous(expand = c(0,0,0,0)) 



ggsave(p2,filename = paste0(outfile,".pdf"),width = width_p,height = height_p )

# TODO :  Zoom a input region  when bin size is smaller enough
if (FALSE) {
  

p3 <-  ggplot(dplyr::filter(allinfo2,group=="tajima's D"|group=="pi"|group=="XPCLR"),
              aes(x=BIN_START/1000000,
                  y=value,
                  group = variable)) + 
  geom_line(aes(group = variable,color = variable ),size =0.5) +
  scale_color_manual(values=c("#ED1C28","#3C499F","#ED1C28","#3C499F","black"))+
  geom_vline(xintercept = 558.1,linetype="dashed",size=0.7,color="seagreen4") +
  geom_vline(xintercept = 573.816,linetype="dashed",size=0.7,color="seagreen4") +
  
  #geom_hline(yintercept = 0.348143,linetype="dashed",size=0.7,color="seagreen4") +
  facet_wrap(~ group,ncol = 1,scales = "free") +
  theme_tufte() +
  theme(axis.line=element_line(),legend.position = "bottom") +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(546, 580)) +
  scale_y_continuous(expand = c(0, 0))

}
