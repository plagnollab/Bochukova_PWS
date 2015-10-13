#$ -S /bin/bash
#$ -l h_vmem=1.8G
#$ -l tmem=1.8G
#$ -l h_rt=12:00:00
#$ -pe smp 1
#$ -R y
#$ -o cluster/out
#$ -e cluster/error
#$ -t 1-8
#$ -tc 8
#$ -cwd

##see https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php

fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta
GATK=/cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar
VCF=data/sites.vcf

##SGE_TASK_ID=1

java=/share/apps/jdk1.7.0_45/bin/java

tail -n +2  /scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3/Bochukova_set3.tab | awk '{if (NR == "'$SGE_TASK_ID'") print}' | while read sample f1 f2 condition; do

    BAM=/scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3/${sample}/${sample}_unique.bam
    ls -ltrh $BAM

    $java -Xmx1g -jar $GATK \
	-R $fasta \
	-T ASEReadCounter \
	-o data/ASE/${sample}.tab \
	-I $BAM \
	-sites:VCF $VCF \
	-minDepth 0 \
	-U ALLOW_N_CIGAR_READS
    
done