#!/bin/bash
cat Dementia_gvcf_list_batch1.txt| while read id; do cp /gvcf/$id* /gvcf/; done
cat Dementia_gvcf_list_batch2.txt| while read id; do cp /gvcf/R*$id* /gvcf/; done

#recall file for varifybamid
ls /recall/|while read id;do mv /recall/$id /rm_recall/;done
ls /recall/|while read id;do mv /recall/$id /rm_recall/;done
cat Dementia_gvcf_list_batch1.txt| while read id; do mv /rm_recall/$id* /recall/; done
cat Dementia_gvcf_list_batch2.txt| while read id; do mv /rm_recall/R*$id* /recall/; done

exit 0
