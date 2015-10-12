# Bochukova_PWS
PWS analysis of Elena Bochukova's samples

## Standard RNA-Seq analysis

Main submission script for RNA-Seq pipeline:
~~~~bash
ls -ltrh submit_master_set3.sh
~~~~

Support files are available here:
~~~~bash
ls -ltrh support
~~~

## Editing analysis

Download editing sites from the web:
~~~~bash
ls -ltrh download_editing_sites.sh
~~~~

Prepare the VCF file for editing analysis:
~~~~bash
ls -ltrh scripts/prepare_VCF.R
~~~~

Count edited/non-edited alleles:
~~~~~bash
ls -ltrh scripts/GATK_count_alleles.sh
~~~~

Some manual checks to confirm the GATK analysis:
~~~~bash
ls -ltrh scripts/manual_checks.sh
~~~~

