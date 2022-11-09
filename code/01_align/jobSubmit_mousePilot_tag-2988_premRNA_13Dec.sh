#!/bin/bash
#$ -N jobSubmit_tag-2988_hpc_premRNA_MNT
#$ -V
#$ -pe local 8
#$ -wd /dcl01/ajaffe/data/lab/singleCell/mouse_10x/premRNA/
#$ -l mem_free=18G,h_vmem=20G,h_fsize=200G
#$ -o ./logs/jobSubmit_tag-2988_hpc_premRNA_MNT_log.txt
#$ -e ./logs/jobSubmit_tag-2988_hpc_premRNA_MNT_log.txt

/dcl01/ajaffe/data/lab/singleCell/cellranger-3.0.2/cellranger count --id=tag2988_hpc_premRNA \
		 --transcriptome=/dcl01/ajaffe/data/lab/singleCell/mouse_10x/premRNA/mm10-3.0.0_premrna \
		 --fastqs=/dcl01/ajaffe/data/lab/singleCell/mouse_10x/FASTQ/2988/ \
		 --sample=2988 \
		 --jobmode=local \
		 --localcores=8 \
		 --localmem=144
