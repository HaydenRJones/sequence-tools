#! /bin/bash

bcftools query -f '%CHROM %DP [%AD]' $1 > ${1%.vcf}_stats.tsv

#echo 'CHROM DP AD' | cat - tmp > ${1%.vcf}_stats.vcf
#rm -rf tmp

python /mnt/g/ubuntu/plotVCF.py ${1%.vcf}_stats.tsv

