#!/bin/bash
for i in /WES_fastq/Sample*/*_R1.fastq.gz
do
j=${i##*/WES_fastq/Sample*/}
java -jar ~/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 16 $i ${i%R1.fastq.gz}R2.fastq.gz -baseout /cleandata/${j%_R1.fastq.gz}.fastq.gz ILLUMINACLIP:~/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 &
done
wait
exit 0
