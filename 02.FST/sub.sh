#!/bin/bash

for i in ../VCFs/chr1A.ann.bcf.gz ../VCFs/chr1B.ann.bcf.gz ../VCFs/chr1D.ann.bcf.gz ../VCFs/chr2A.ann.bcf.gz ../VCFs/chr2B.ann.bcf.gz ../VCFs/chr2D.ann.bcf.gz ../VCFs/chr3A.ann.bcf.gz ../VCFs/chr3B.ann.bcf.gz ../VCFs/chr3D.ann.bcf.gz ../VCFs/chr4A.ann.bcf.gz ../VCFs/chr4B.ann.bcf.gz ../VCFs/chr4D.ann.bcf.gz ../VCFs/chr5A.ann.bcf.gz ../VCFs/chr5B.ann.bcf.gz ../VCFs/chr5D.ann.bcf.gz ../VCFs/chr6A.ann.bcf.gz ../VCFs/chr6B.ann.bcf.gz ../VCFs/chr6D.ann.bcf.gz ../VCFs/chr7A.ann.bcf.gz ../VCFs/chr7B.ann.bcf.gz ../VCFs/chr7D.ann.bcf.gz ;do
    ./FST.sh 100 100 ${i} CNC CNL  & 
done
wait
