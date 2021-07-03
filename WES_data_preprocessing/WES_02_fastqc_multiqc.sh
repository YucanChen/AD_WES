#!/bin/bash
fastqc -o /clean_fastqc/ -t 16 /cleandata/*.fastq.gz
multiqc /clean_fastqc/ -o /clean_fastqc/ -i WES_batch2
exit 0
