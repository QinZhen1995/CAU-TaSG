#!/bin/bash
 
SM=${1} #  bam path 
GeneID=${2} #  02g
Flank=${3} #  2000
WdSize=${4} # 1 

## define POS 
#if false;then
grep  -w ${GeneID}    ~/genome/index/IWGSC_v1_wzh/IWGSC_v1.1_HC_LC.genebed2  \
    | awk  -F"\t" -vOFS="\t" -vFlank=${Flank} '{print $2,$3-Flank,$4+Flank}' > ${GeneID}_${Flank}.bed 

CHR=`cat  ${GeneID}_${Flank}.bed  | awk '{print $1}'`
POS=`cat  ${GeneID}_${Flank}.bed  | awk '{print $1":"$2"-"$3}'`
Strand=`grep  -w ${GeneID}    ~/genome/index/IWGSC_v1_wzh/IWGSC_v1.1_HC_LC.genebed2  \
   | awk '{print $5}'  `


## make bams link 
mkdir -p ${GeneID}_${Flank}_bams


for i in ` cat  ${SM}    `   ;do 
    ln -s  ~/project/reseq_data/workflowfile/${i}/05.mergeAsplit/${CHR}.dedup.bam ./${GeneID}_${Flank}_bams/${i}.bam
    ln -s  ~/project/reseq_data/workflowfile/${i}/05.mergeAsplit/${CHR}.dedup.bam.bai ./${GeneID}_${Flank}_bams/${i}.bam.bai
done 
wait 

echo "can't find below bams, skip them"
find ./${GeneID}_${Flank}_bams -xtype l   
find ./${GeneID}_${Flank}_bams -xtype l  -exec rm {} \;

## make window  
bedtools makewindows -b ${GeneID}_${Flank}.bed  -w ${WdSize}  > ${GeneID}_${Flank}_wd${WdSize}.bed
wait 



## bedtools coverage    
mkdir -p ${GeneID}_${Flank}_coverage

for i in ./${GeneID}_${Flank}_bams/*.bam ;do
    name=`basename ${i}`
    samtools view -b -h ${i} ${POS} | bedtools coverage -split  -a ${GeneID}_${Flank}_wd${WdSize}.bed  -b - | \
    awk -F"\t" -vOFS="\t" -v name=${name%%.bam} '{print $1,$2,$3,$4,name}' > ./${GeneID}_${Flank}_coverage/${name%%.bam}.bed
done
wait
header=$(find ./${GeneID}_${Flank}_coverage/ -type "f" -name "*.bed" | xargs cat | datamash  -s -g 1,2,3  collapse 5 | awk '{print $NF}' | uniq | sed 's/,/\t/g' )

find ./${GeneID}_${Flank}_coverage/ -type "f" -name "*.bed" | xargs cat  | datamash  -s -g 1,2,3  collapse 4 | sed 's/,/\t/g' |  \
    awk -F"\t" -vOFS="\t"  '{print $0}' > ./${GeneID}_${Flank}_coverage/${GeneID}_${Flank}_coverage.bed 

sed  "1iCHR\tSTART\tEND\t\\$header"   ./${GeneID}_${Flank}_coverage/${GeneID}_${Flank}_coverage.bed \
    | csvtk cut -t -f -1,-3 > ${GeneID}_${Flank}_wd${WdSize}_coverage.martix

wait

bks=$(expr ${Flank} / ${WdSize} )


/usr/bin/Rscript plotcoverage.R  \
    -m ${GeneID}_${Flank}_wd${WdSize}_coverage.martix  \
    -a ~/qinz/caojie/10DFYC_bam2fasta/331cailiao.gvcfid5 \
    -s /data/user/qinz/qinz/caojie/Bin2_promoter/DFYC.gvcfid \
    -t "${GeneID} ${POS}(${Strand}) " \
    -b ${bks} -o ${GeneID}_${Flank}_wd${WdSize}
