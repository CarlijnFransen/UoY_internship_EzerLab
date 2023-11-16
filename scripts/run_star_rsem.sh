#!/bin/bash
default_dir="/shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/data/01.RawData/*"

for directory in $default_dir
do
	sample_name=$(basename $directory)
	echo $sample_name
	
	l1_file=/shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/data/01.RawData/$sample_name/*_1.fq.gz
	l2_file=/shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/data/01.RawData/$sample_name/*_2.fq.gz
	
	echo $l1_file
	echo $l2_file

	STAR \
	--genomeDir /shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/TAIR10/index_gtf \
	--runThreadN 15 \
	--readFilesIn $l1_file 	$l2_file \
	--readFilesCommand "gunzip -c" \
	--outSAMtype BAM SortedByCoordinate \
	--twopassMode Basic \
	--outFilterScoreMinOverLread 0 \
	--outFilterMatchNminOverLread 0 \
	--outFilterMatchNmin 0 \
	--outFilterMismatchNmax 1 \
	--sjdbGTFfile "/shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/TAIR10/adjusted_tair10_gtf3.gtf" \
	--quantMode TranscriptomeSAM GeneCounts \
	--outFileNamePrefix /shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/mapping_output/$sample_name/


	mkdir /shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/rsem_output/$sample_name/
	
	rsem-calculate-expression --no-bam-output --bam -p=15 --paired-end \
	"/shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/mapping_output/$sample_name/Aligned.toTranscriptome.out.bam" \
	"/shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/TAIR10/rsem_tair_reference" \
	/shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/rsem_output/$sample_name/$sample_name


done

