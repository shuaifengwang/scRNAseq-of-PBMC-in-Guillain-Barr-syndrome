#!/bin/bash
bin=../pipeline/cellranger-7.2.0/bin/cellranger
db=../pipeline/refdata-gex-GRCh38-2020-A  # 
ls $bin; ls $db 

fq_dir=/home/data/project/10x/raw/PBMC1
$bin count --id= PBMC1_outs \
--localcores= 8 \
--transcriptome= $db \
--fastqs= $fq_dir \
--sample= PBMC1   \
--localmem=64

