##see https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php

fasta="/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta"
GATK=/cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar

#VCF=out_recode.vcf
VCF=data/sites.vcf

tail -n +2  /scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3/Bochukova_set3.tab | while read sample f1 f2 condition; do

    BAM=/scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3/${sample}/${sample}_unique.bam
    ls -ltrh $BAM

    java -jar $GATK \
	-R $fasta \
	-T ASEReadCounter \
	-o data/ASE/${sample}.csv \
	-I $BAM \
	-sites:VCF $VCF \
	-U ALLOW_N_CIGAR_READS

done