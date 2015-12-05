########## first the high level parameters about where the scripts and the bundle are located
RNASEQPIPBASE=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline  ##this probably needs to be updated, but I do not think it has to
RNASEQBUNDLE=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle  ##this probably should stay

export RNASEQPIPBASE=$RNASEQPIPBASE
export RNASEQBUNDLE=$RNASEQBUNDLE
pipeline=${RNASEQPIPBASE}/RNAseq_pipeline_v8.sh

#############
submit=no

QC=no
starStep1a=no
starStep1b=no
starStep2=no

sampleQC=no
runCufflinks=no

summary=no
prepareCounts=no
Rdeseq=yes
Rdexseq=no

force=yes

oFolder=/scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3_hg38
iFolder=/scratch2/vyp-scratch2/ElenaBochukova

sh $pipeline --starStep1a ${starStep1a} --starStep1b ${starStep1b} --starStep2 ${starStep2} --iFolder ${iFolder} --oFolder ${oFolder} --dataframe Bochukova_PWS/support/Bochukova_set3.tab --code Bochukova_set3 --prepareCounts ${prepareCounts} --Rdexseq ${Rdexseq} --Rdeseq ${Rdeseq} --submit ${submit} --summary ${summary} --species human_hg38 --force ${force}

echo $mainscript
#qsub $mainscript
    #sh $mainscript




