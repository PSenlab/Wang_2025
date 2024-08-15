#! /bin/bash

module load spaceranger 

spaceranger count --id ${sample_directory} \
    --transcriptome=${SPACERANGER_REF}/refdata-gex-mm10-2020-A \
    --fastqs=/path/to/sample \
    --sample=${sample_name} \
    --image=/path/to/${sample}.tif \
    --unknown-slide \
    --r1-length=26 \
    --localcores=$SLURM_CPUS_PER_TASK \
    --localmem=63
