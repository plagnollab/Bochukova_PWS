find /scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3 -name *unique.bam > list_BAMs.tab

#samtools mpileup -r 1:160185788-160185788 -f /scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta -b list_BAMs.tab > chr1_160185788.tab

samtools mpileup -r 1:29529699-29529699 -f /scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta -b list_BAMs.tab > chr1_29529699.tab
