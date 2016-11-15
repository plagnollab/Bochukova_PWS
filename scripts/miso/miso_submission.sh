

oFolder=/SAN/biomed/biomed14/vyp-scratch/Bochukova/hg19_UCSC
dataframe=support/Bochukova_set3.tab


step1=no  ### prepare all the scripts, but does not run them. Runs by combination of event and sample
step2=yes

if [[ $step1 == "yes" ]]; then
    
    tail -n +2 $dataframe | while read sample junk; do
	if [ ! -e ${oFolder}/${sample}/MISO ]; then mkdir ${oFolder}/${sample}/MISO; fi
	
	for event in A3SS A5SS MXE RI SE; do
	    BAM=${oFolder}/${sample}/${sample}_unique.bam
	    annotations=/SAN/biomed/biomed14/vyp-scratch/temp_reference/hg19_UCSC/hg19/pickled/${event}
	    
	    if [ ! -e ${oFolder}/${sample}/MISO/${event} ]; then mkdir ${oFolder}/${sample}/MISO/${event}; fi
	    
	    if [ -e ${oFolder}/${sample}/MISO/${event}/cluster_scripts/ ]; then rm ${oFolder}/${sample}/MISO/${event}/cluster_scripts/*; fi
	    
	    echo "
miso --run ${annotations} ${BAM} --output-dir ${oFolder}/${sample}/MISO/${event} --read-len 81 --paired-end 150 35 --use-cluster 

" > cluster/submission/MISO_${sample}_${event}.sh
	    
	    ls -ltrh cluster/submission/MISO_${sample}_${event}.sh
	    sh cluster/submission/MISO_${sample}_${event}.sh
	done
	
    done
fi

if [[ $step2 == "yes" ]]; then
    
    tail -n +2 $dataframe | while read sample junk; do
	echo $sample
	tabscript=cluster/submission/MISO_${sample}.tab
	mainscript=cluster/submission/MISO_${sample}.sh

	if [ -e $tabscript ]; then 
	    rm $tabscript
	fi

	for event in A3SS A5SS MXE RI SE; do 
	    find ${oFolder}/${sample}/MISO/${event}/cluster_scripts/ -name \*sh >> ${tabscript}
	done
	
	nfiles=`wc -l ${tabscript} |  cut -f1 -d ' '`

	ls -ltrh ${tabscript}
	echo "
#$ -S /bin/bash
#$ -l h_vmem=3.8G
#$ -l tmem=3.8G
#$ -l h_rt=12:00:00
#$ -pe smp 1
#$ -R y
#$ -o cluster/out
#$ -e cluster/error
#$ -t 1-${nfiles}
#$ -tc 50
#$ -cwd


FILE=\$(awk \"NR == \$SGE_TASK_ID\" ${tabscript} )

sh \$FILE

" > $mainscript

	ls -ltrh $mainscript
    done
fi