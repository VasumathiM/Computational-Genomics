#!/bin/bash
#normal sample 
# 1. Convert SAM to BAM
samtools view -S -b Normal_S1_GCF_000001405.sam > Normal_S1.bam

# 2. Sort by read name for fixmate
samtools sort -n -o Normal_S1sort.bam Normal_S1.bam

# 3. Fix mate-pair information
samtools fixmate -m Normal_S1sort.bam Normal_S1.fixmate.bam

# 4. Sort by genomic coordinates
samtools sort -o NormalS1position.fixmat.bam Normal_S1.fixmate.bam

# 5. Mark and remove duplicates
samtools markdup -r NormalS1position.fixmat.bam NormalS1.rm.bam

# 6. Index the final BAM file
samtools index NormalS1.rm.bam

#cancer sample
# 1. Convert SAM to BAM
samtools view -S -b Tumor_S2_GCF_000001405.sam > Tumor_S2.bam

# 2. Sort by read name for fixmate
samtools sort -n -o Tumor_S2sort.bam Tumor_S2.bam

# 3. Fix mate-pair information
samtools fixmate -m Tumor_S2sort.bam Tumor_S2.fixmate.bam

# 4. Sort by genomic coordinates
samtools sort -o TumorS2position.fixmat.bam Tumor_S2.fixmate.bam

# 5. Mark and remove duplicates
samtools markdup -r TumorS2position.fixmat.bam TumorS2.rm.bam

# 6. Index the final BAM file
samtools index TumorS2.rm.bam

# ref genome

samtools faidx GCF_000001405.40_GRCh38.p14_genomic.fasta
