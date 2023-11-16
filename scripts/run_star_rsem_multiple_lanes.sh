#!/bin/bash
default_dir="/shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/data/01.RawData/*"
directory_list=("WT_10_C" "WT_14_T" "WT_23_C" "WT_40_C" "WT_47_C")
for directory in ${!directory_list[@]}
do
        sample_name=${directory_list[directory]}
        echo $sample_name

        l1_file=/shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/data/01.RawData/$sample_name/*concatenated*_1.fq.gz
        l2_file=/shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/data/01.RawData/$sample_name/*concatenated*_2.fq.gz

        STAR \
        --genomeDir /shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/TAIR10/index_gtf \
        --runThreadN 15 \
        --readFilesIn $l1_file  $l2_file \
        --readFilesCommand "gunzip -c" \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --outFilterScoreMinOverLread 0 \
        --outFilterMatchNminOverLread 0 \
        --outFilterMatchNmin 0 \
        --outFilterMismatchNmax 1 \
        --sjdbGTFfile "/shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/TAIR10/tair10_gtf_genes.gtf" \
        --quantMode TranscriptomeSAM GeneCounts \
        --outFileNamePrefix /shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/mapping_output/${sample_name}_concatenated_2/


        mkdir /shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/rsem_output/${sample_name}_concatenated_2/

        rsem-calculate-expression --no-bam-output --bam -p=15 --paired-end \
        "/shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/mapping_output/${sample_name}_concatenated_2/Aligned.toTranscriptome.out.bam" \
        "/shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/TAIR10/rsem_tair10_ref" \
        /shared/storage/biology/rsrch/de656/lab/projects/project_CarlijnFransen/rsem_output/${sample_name}_concatenated_2/$sample_name


done
