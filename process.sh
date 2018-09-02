
#yff - name of directry you want to work in. This script assumes:
#1. your raw fastq files are in a directory called "raw" within "yff"
#2. your raw fastq files are arranged in  sub-directories within "raw": each subdirectory contains R1 and R2 files of a single sample; the subdirectory's name is the samples name 
#Note: this script also calls on addiational files called "parameters_multisplit_2.txt" and "uc2otutab_mod.py"; these ae found in teh "additional_files" directory of this repository

#!/bin/bash
hostname
module load gcc/gcc-7.2.0
module load vsearch/vsearch-2.7.0
module load python/miniconda_python-2.7_qiime

yff=/scratch200/leah/Crisp2018/A

mkdir $yff/pear #this will hold all files created by "PEAR" program when merging per-sample forward and reverse reads
mkdir $yff/files #final merged AND quality-filtered (by PEAR) per-sample files will be put here; each file will be labelled by the sample name.


for dir in $yff/raw/*
do
if [ -d "${dir}" ] ; then
        dir_index=$( echo $dir |cut -d  '/' -f7 ) #extract the sample name from teh directory name; may need to adjust depending on strucuture and location  of sub-directory name
    cd $dir/Files

                echo "$dir"
                echo "$dir_index"
                /usr/local/bin/pear -f *R1* -r *R2* -q 20  -n 200 -o $yff/pear/"$dir_index"
        fi
  done


mv $yff/pear/*.assembled.fastq  $yff/files

# A qiime script is used to concatenate all fastq files to a single fasta file, while preserving sample identity (i.e., from which sample each seqeunce orignated)
multiple_split_libraries_fastq.py -i $yff/files --demultiplexing_method sampleid_by_file -p /scratch200/leah/Scripts/parameters_multi_split2.txt  -o $yff/spl
# Vsearch is used to dereplicate and then cluster sequences:
vsearch --derep_fulllength $yff/spl/seqs.fna --output $yff/derep.fna --minuniquesize 1 --fasta_width 0 --sizeout --uc $yff/map_derep.txt
vsearch --fasta_width 0 --sortbysize $yff/derep.fna -output $yff/derep_sorted.fna 
id=0.99 #iddef 1 ensures terminal gaps are taken into account when clustering. THE CLUSTERED OTUS WILL BE USED AS INPUT FOR THE SPACER_EXTRACTOR SCRIPT.
vsearch  -cluster_size $yff/derep.fna --centroids $yff/otus_$id.fna --uc $yff/uc.txt -iddef 1  --relabel otu_ --id $id  --sizeout --sizein --qmask none 
echo "done clustering at percent ID"
echo $id
# all sequences are then mapped back to the clusters:
vsearch  -usearch_global $yff/spl/seqs.fna -db $yff/otus_$id.fna -strand both -id $id  -uc $yff/readmap_$id.uc 
#The mapped output is converted to otu table form, called 'tabfile'. THE TABFILE WILL BE USED AS FOR QUANTIFYING HOTSPOTS OF ACQUISTION
python /scratch200/leah/drive5_files/uc2otutab_mod.py $yff/readmap_$id.uc > $yff/tabfile_$id.txt

