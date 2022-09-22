#!/bin/bash

#PBS -N cc_seqkit_paired
#PBS -l nodes=1:ppn=16
#PBS -q cv2
#PBS -j oe
#PBS -l walltime=72:00:00
#PBS -V


#cat all the lanes BJ333D30 library
#cd /p200/home/chengchen_gibh/rawReads/V350092063/
#mkdir -p BJ333D30
#cat rawBJ333D30/V350092063_L01_7_1.fq.gz rawBJ333D30/V350092063_L02_7_1.fq.gz rawBJ333D30/V350092063_L03_7_1.fq.gz rawBJ333D30/V350092063_L04_7_1.fq.gz > ./BJ333D30/BJ333D30_1_fq.gz

#cat rawBJ333D30/V350092063_L01_7_2.fq.gz rawBJ333D30/V350092063_L02_7_2.fq.gz rawBJ333D30/V350092063_L03_7_2.fq.gz rawBJ333D30/V350092063_L04_7_2.fq.gz > ./BJ333D30/BJ333D30_2_fq.gz

#/public/home/chengchen_gibh/software/DNBelab_C_Series_HT_scRNA-analysis-software/software/STAR --runThreadN 16 --runMode alignReads --genomeLoad LoadAndKeep --limitBAMsortRAM 64000000000  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --genomeDir /public/home/chengchen_gibh/refdata-gex-GRCh38-2020-A/customHg38BGI --outFileNamePrefix ./BJ333D30/BJ333D30 --readFilesIn ./BJ333D30/BJ333D30_1_fq.gz ./BJ333D30/BJ333D30_2_fq.gz

#extract unmapped read1 and read2 from bam files
#samtools view -b -f 13 ./BJ333D30/BJ333D30Aligned.sortedByCoord.out.bam > ./BJ333D30/BJ333D30_oligo.bam
#echo "extracted paired unmapped reads"
#bedtools bamtofastq -i ./BJ333D30/BJ333D30_oligo.bam -fq ./BJ333D30/BJ333D30_oligo_1.fq -fq2 ./BJ333D30/BJ333D30_oligo_2.fq
#echo "etracted oligo read1 and read2,needs further quality check"


#convert oligo files
cd ./BJ333D30
gzip BJ333D30_oligo_1.fq
gzip BJ333D30_oligo_2.fq

python3 ../huada_fqsplit.py -1 BJ333D30_oligo_1.fq.gz -2 BJ333D30_oligo_2.fq.gz -c BJ333D30_oligo_1_format.fq.gz -d BJ333D30_oligo_2_format.fq.gz


#extract cDNAs
samtools view -bf 69 BJ333D30Aligned.sortedByCoord.out.bam > BJ333D30_cDNA_1.bam
bedtools bamtofastq -i BJ333D30_cDNA_1.bam -fq BJ333D30_cDNA_1.fq

samtools view -bf 137 BJ333D30Aligned.sortedByCoord.out.bam > BJ333D30_cDNA_2.bam
bedtools bamtofastq -i BJ333D30_cDNA_2.bam -fq BJ333D30_cDNA_2.fq

samtools cat BJ333D30_cDNA_1.bam BJ333D30_cDNA_2.bam -o BJ333D30_cDNA.bam
samtools sort -n BJ333D30_cDNA.bam -o BJ333D30_cDNA_sorted.bam
bedtools bamtofastq -i BJ333D30_cDNA_sorted.bam -fq BJ333D30_cDNA_1.fq -fq2 BJ333D30_cDNA_2.fq
#success with bedtools

#gzip BJ333D30_cDNA_1.fq
#gzip BJ333D30_cDNA_2.fq

#quality check by fastp
fastp -i BJ333D30_cDNA_1.fq -I BJ333D30_cDNA_2.fq -o BJ333D30_cDNA_1_trimmed.fq.gz -O  BJ333D30_cDNA_2_trimmed.fq.gz -h BJ333D30_cDNA.html -A

fastp -i BJ333D30_oligo_1.fq.gz -I BJ333D30_oligo_2.fq.gz -o BJ333D30_oligo_1_trimmed.fq.gz -O BJ333D30_oligo_2_trimmed.fq.gz -A -h BJ333D30_oligo.html

