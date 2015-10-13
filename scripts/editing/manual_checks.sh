find /scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3 -name *unique.bam > list_BAMs.tab

#samtools mpileup -r 1:160185788-160185788 -f /scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta -b list_BAMs.tab > chr1_160185788.tab

#samtools mpileup -r 1:29529699-29529699 -f /scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta -b list_BAMs.tab > chr1_29529699.tab

#samtools mpileup -r 2:25042749-25042749 -f /scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta -b list_BAMs.tab > chr2_25042749.tab

#samtools mpileup -r 4:158257875-158257875 -f /scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta -b list_BAMs.tab > chr4_158257875.tab

samtools mpileup -r 2:160626593-160626593 -f /scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta -b list_BAMs.tab > chr2_160626593.tab
