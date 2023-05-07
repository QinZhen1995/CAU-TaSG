#!/bin/bash


F=${3}
#chr1A.ann.bcf.gz
CHR=`basename ${F} | cut -d"." -f1 `

#S1=./CNC_group
S1=${4}

#S2=./CNL_group
S2=${5}

size=${1}

step=${2}
# fst 
#if false ;then
vcftools --bcf  ${F}   --weir-fst-pop  ${S1}  --weir-fst-pop  ${S2}  --fst-window-size ${size}000 --fst-window-step ${step}000  --out ${CHR}_${size}k.fst &

vcftools --bcf $F --keep ${S1} --TajimaD  ${size}000 --out ${CHR}_${S1}_${size}k.tajimaD & 

vcftools --bcf $F --keep ${S2} --TajimaD  ${size}000 --out ${CHR}_${S2}_${size}k.tajimaD &

wait
#awk -vOFS="\t" -F"\t" '{split($1,a,".");if(a[2] == "1" || a[2] == ""){$1=a[1] ; print} else {$1=a[1];$2 = $2+400000000;print} }' ${CHR}_${size}k.fst.windowed.weir.fst | sponge chr5B_${size}k.fst.windowed.weir.fst

paste ${CHR}_${S1}_${size}k.tajimaD.Tajima.D ${CHR}_${S2}_${size}k.tajimaD.Tajima.D | cut -f1,2,3,4,8  > ${CHR}_${S1}vs${S2}_${size}k.TD 

#awk -vOFS="\t" -F"\t" '{split($1,a,".");if(a[2] == "1" || a[2] == ""){$1=a[1] ; print} else {$1=a[1];$2 = $2+400000000;print} }'  ${CHR}_${size}k.TD  | sponge  ${CHR}_${size}k.TD 


# pi1 pi2 
vcftools --bcf $F --keep ${S1} --window-pi ${size}000 --window-pi-step ${step}000 --out ${CHR}_${S1} &

vcftools --bcf $F --keep ${S2} --window-pi ${size}000 --window-pi-step ${step}000 --out ${CHR}_${S2} &

wait

# pi1 in pi2 
gawk -vOFS="\t" 'ARGIND==1{A[$1":"$2]=$5;} ARGIND==2&&(FNR>1){if($1":"$2 in A){print $1, $2, $3, A[$1":"$2], $5, A[$1":"$2]/$5, $5/A[$1":"$2];}}' \
    ${CHR}_${S1}.windowed.pi ${CHR}_${S2}.windowed.pi  > ${CHR}_${S1}vs${S2}.${size}k.windowed.pi

wait 

# pi ratio 
#awk -vOFS="\t" -F"\t" '{split($1,a,".");if(a[2] == "1" || a[2] == ""){$1=a[1] ; print} else {$1=a[1];$2 = $2+400000000;print} }' ${CHR}_res.${size}k.windowed.pi | sponge ${CHR}_res.${size}k.windowed.pi 

gawk -F"\t" -vOFS="\t" '{print $1,$2,$3,$7}' ${CHR}_${S1}vs${S2}.${size}k.windowed.pi  > ${CHR}_${S1}vs${S2}.${size}k.windowed.piratio

sed -i '1i CHROM\tSTART\tEND\tPiDO/PiDE' ${CHR}_${S1}vs${S2}.${size}k.windowed.piratio

# log
gawk -vOFS="\t" '{$4=log($4+1)/log(10);print}' ${CHR}_${S1}vs${S2}.${size}k.windowed.piratio > ${CHR}_${S1}vs${S2}.${size}k.windowed.logpiratio


wait 
# Plot
wait
Fst_plot.R -f ${CHR}_${size}k.fst.windowed.weir.fst \
    -p ${CHR}_${S1}vs${S2}.${size}k.windowed.pi \
    -t ${CHR}_${S1}vs${S2}_${size}k.TD \
    -a ${CHR}_FstAnnForRplot.txt \
    -w 8 -v 6 \
    -o ${CHR}_Fst_plot_${size}k




wait 
