########## first the high level parameters about where the scripts and the bundle are located
RNASEQPIPBASE=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline  ##this probably needs to be updated, but I do not think it has to
RNASEQBUNDLE=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle  ##this probably should stay

export RNASEQPIPBASE=$RNASEQPIPBASE
export RNASEQBUNDLE=$RNASEQBUNDLE
pipeline=${RNASEQPIPBASE}/RNAseq_pipeline_v8.sh

#############
submit=yes

QC=no
starStep1a=yes
starStep1b=yes
starStep2=yes

sampleQC=no
runCufflinks=no

summary=no
prepareCounts=no
Rdeseq=no
Rdexseq=no

force=no

oFolder=/SAN/biomed/biomed14/vyp-scratch/Bochukova/hg38_kitty
iFolder=/cluster/scratch3/vyp-scratch2/Bochukova_RNASeq/fastq/set3

sh $pipeline --step0_QC ${QC} --starStep1a ${starStep1a} --starStep1b ${starStep1b} --starStep2 ${starStep2} --iFolder ${iFolder} --oFolder ${oFolder} --dataframe support/Bochukova_set3.tab --code Bochukova_set3 --prepareCounts ${prepareCounts} --Rdexseq ${Rdexseq} --Rdeseq ${Rdeseq} --submit ${submit} --summary ${summary} --species human_hg38 --force ${force} --trim_galore no





