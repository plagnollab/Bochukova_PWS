for sample in Control1 Control2 Control3 Control4 PWS1 PWS2 PWS3 PWS4; do
    
    for eventType in A3SS A5SS MXE RI SE; do
	summarize_miso --summarize-samples /SAN/biomed/biomed14/vyp-scratch/Bochukova/hg19_UCSC/${sample}/MISO/${eventType} MISO/summary/${eventType}_${sample}
    done

done