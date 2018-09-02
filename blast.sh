#blast teh output from "spacer_extractor.r " script against a combined database of mediterranei and volcanii genomes
#blast database found in "wdditional_files" in this repository


#!/bin/bash

hostname
module load blast/blast-2.6.0
yff=/scratch200/leah/Crisp2018/A
yff1=/scratch200/leah/Crisp2018
id=0.99
blastn -task blastn-short -query $yff/all_spacers_$id.txt -db $yff1/combined_medi_volc_genomes_newfromSam.fas -out $yff/blast_$id.txt -outfmt 6 -dust no -evalue 0.0001
