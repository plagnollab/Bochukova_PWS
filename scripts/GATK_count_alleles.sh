##see https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php

fasta="~/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta"
GATK=

cat /scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3/Bochukova_set3.tab | while read sample f1 f2 condition; do


    BAM=/scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3/${sample}/${sample}_unique.bam
    
    java -jar $GATK \
	-R $fasta \
	-T ASEReadCounter \
	-o test.csv \
	-I $BAM \
	-sites sites.vcf \
	-U ALLOW_N_CIGAR_READS \
	[-minDepth 10] \
	[--minMappingQuality 10] \
	[--minBaseQuality 2] \
	[-drf DuplicateRead]

done