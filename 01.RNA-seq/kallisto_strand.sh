#!/bin/bash

pwd; hostname; mydate

mkdir -p 00.clean
mkdir -p 01.kallisto

for i in `echo ${1}`
    do 
        (R1=`echo ${i} | awk -F "|" '{print$1}'`
        R2=`echo ${i} | awk -F "|" '{print$2}'`
        name=`echo ${i} | awk -F "|" '{print$1}' | cut -d '_' -f1`
        echo "R1 fastq file:" ${R1}
        echo "R2 fastq file:" ${R2}
        echo "   out prefix:" ${name}
        echo "------fastp START-----"
        fastp -i ./${R1} -I ./${R2} \
            -o ./00.clean/${name}_R1_clean.fq.gz \
            -O ./00.clean/${name}_R2_clean.fq.gz \
            -3 -5 -W 6 -l 30 -c -h ./00.clean/${name}.html -z 6 -w 10
        wait
        echo `mydate` "------fastp END-----"
        echo "------quant START-----"
        kallisto quant \
            -i ~/genome/index/IWGSC_v1_part/IWGSC_part_v1.1_HC_LC_kallisto  \
            -o ./01.kallisto/${name}  \
      		--bias -t 1 \
      		-b 2 \
      		--single-overhang \
            --rf-stranded \
            --genomebam \
            -g  ~/genome/index/IWGSC_v1_part/IWGSC_part_v1.1_HC_LC.gtf \
            -c  ~/genome/index/IWGSC_v1_part/chrlength.tab \
            ./00.clean/${name}_R1_clean.fq.gz ./00.clean/${name}_R2_clean.fq.gz ) &
done
wait
mydate
echo '------ALL FINISH-----'
