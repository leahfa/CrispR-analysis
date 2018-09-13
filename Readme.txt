Many of the scripts and documentation on this page are the work of Sam Marton (sam.d.marton@gmail.com).

--------------------------------
PRIMARY PIPELINE
_________________________________

To begin with, we have paired-end fastq files per array per sample. These files need to be quality filtered, merged,  dereplicated and clustered (to reduce both complexity and impact  of sequencing errors) , while retaining per-sample abundance information. New spacers are then identified and extracted, and all relevant information regarding them needs to be retrieved - this includes target location,functional annotation for that location, PAM sequence and so on. The following 5 scripts describe this process, from raw fastq to a "spacer" table with full information.

Note: clustering threshold is specified with "ID" parameter (default is 99) ; this parameter is appendixed to the name of many of the output files; thus, spacer_99 means spacers obtianed from data clustered at 99 percent ID, and so on.

The following scripts are done per array; typically, we set a directory for each array, and run all these intial steps within that directory. At a later stage, final spacer tables from all array directories are integrated to larger tables. Details on each script are also found as notes within the scripts, as well as (for some of the scripts) at the end of this page.

___Step1: process.sh _____

This script untilizes PEAR, qiime and vsearch to merge, clean, dereplicate, and cluster identical sequences at 99% identity, producing the fasta file: otus_99.fna; these OTU sequences will be screened for new spacer discovery. Additionally, all raw sequences are mapped back to the clustered OTUS, to produce an OTU per-sample abundance table (tabfile), which will be used for quantification.

___Step2: spacer_extractor.R_________ 

Input is otus_0.99.fna, outputs are fasta files of new_spacers ( designated spacers_1 for the 1st acquired spacer, spacers_2 for the 2nd acquired spacer to the same OTU, and so on;  all_spacers holds all newly acquired spacers , regardless of whether they originated in  a single, double or triple acquisition event) ,  and a text file of "original array" - names of OTUS which have the SAME number of repeats as the original array (that is, no new acquistions); this is called ori.txt

____Step 3: blast.sh _____

Blast   all_spacers.fna against a database containing Volcanii and Mediterranei genomes and all of their plasmids. Blast database is in folder "additional_file". Set e for 0.0001, and save results in tabular format. This step can also be done locally in Windows using "BioEdit"

____Step 4: After_blast.R____ 

Integrate taxonomy and abundance information per spacer.  The inputs are the tabular blast results file, the spacers fasta file,  and the corresponding table of abundance information (tabfile). The output is a text table, called "final", which gives for every unique spacer its genomic location,as well as abundance information across the samples. Note: Identical spacers acquired by different OTUS are converted  to a single entry, and their per-sample abundances are summed,  to preserve the abundance information.


 ____Step 5: final_refinement_with_annotations.R. ___
Find the PAM sequence, reverse PAM sequence,  and add functional annotation , if possible, for each spacer in the "final" table. This script reads in genome seqeunce and annotation files for all target replicons (found in "additional files" directory in this repository) , and loops over all the arrays, assuming that each array has a subdirectory of its own, which contains the "final" file for that array.  All relevant information will be extracted from the geome and annotation files to the "final" file, which will be, at the end of the process , saved as, for example:  "A_no_multiple_targets.csv".

______________________________________________________


The end product of these scripts is a summary table per array, called , for example, "A_no_multiple_targets.csv", and containing all data relevant to each unique spacer: per-sample abundance, prevalence, PAM, anntotation, location and so on. The next steps are intended to integrate spacers from ALL arrays for combined analysis:

____step6:  Integrate arrays by species.R_______

This script integrates the summay tables of all arrays of each species to a single table, titled , for example: Vol_no_multiple_targets.csv. This integrated table can then be used for downstream analyses: identifying hotspots, drawing PAMwheels and pamscapes, and so on. All this scripts  are found in the folder "R scripts for 
downstream analysis". 


________________________________________________
DETAILED DESCRIPTION OF SPACER EXTRACTION, BLAST RESULTS REFINEMENT AND PAM RETRIEVAL

   ________________________
___|SPACER EXTRACTOR (vB6)|___
   ------------------------

The purpose of this script as to filter through the supplied OTUs, find every repeat for that given array (definined as a global variable at the beginning of the script), and for OTUs which have an extra repeat, save the bases next to the new repeat as a spacer.  There is a dialogue box that appears when running the script which allows you to define whether the leader sequence occurs at the beginning or the end of the array, as well as the number of repeats at which there is an acquisition (4 by default, depending on the primer locations in the array).  
~NOTE: WHEN RUNNING THIS SCRIPT, you MUST hard-code the array you are searching in.  After all the repeats and sequences are defined, you can see variables for "ORIGINAL_SEQUENCE" and "REPEAT_SEQUENCE".  Simply type the letter of the desired array in the initialization of these variables.

Shortly after, there are variables for "REPEAT_ERRORS" and "SPACER_ERRORS".  These define how strictly to filter the results.  We have been using a "strict" filter of 2 for the first variable and 5 for the second, which are used by the bioStrings functions.  These can be adjusted to increase or decrease the number of results as desired.  This is indicated as a comment.

The program creates successive series of vectors until finally yielding the output in the desired form.  First, it converts the OTUs to a bioString type so that the search functions will work.  The first vector created I called "filter", which simply removes the OTUs for which there are NOT a sufficient number of spacers found.  If none are found, the program ends.  The remaining OTUs are condensed into a vector I called "acq", for acquisitions.  From all of these acquisitions, the next step is to compare them to the originals CRISPR arrays and remove any for which there is a match.  The whole arrays are defined at the top of the script in both the forward and reverse directions.  Excessively long (>50) or short (<10) "spacers" are also removed.  These results are copied into two parallel vectors:  The first I called "spacers", which shows the full OTU with the spacer marked by arrows as follows:

...xxxxxxxxx ->xxxxx<- xxxxxxxxx...

This is useful for visually verifying that the script is working properly.  The other vector I called "final", which simply retains the OTU and the spacer itself, in a FASTA format.  Again, if no good spacers are found, the program ends.  Finally, this vector is automatically saved to a predetermined output folder with a filename generated based on the input file (location should be edited for whatever machine you are using).  In the filename, "End" or "Start" indicates the position of the leader sequence, while the number following "out_" indicates how many acquisition events there were.  For example, out_1 is for spacers that had at least one acquisition, out_2 is for spacers in which there were at least two acquisitions, and so on.  These are found successively in a loop, and on each loop, a new generation of the "filter," "acq," "spacers," and "final" vectors are generated and numbered appropriately.  Any of these can be viewed in the command window explicitly, such as with spacers4, or with the "get" function, for example get(acqX).

#additional modification: the vecor "ori" is built durng the first iteration, and holds teh names of all  OTUS that have the origina number of repeats for that array. This is then saved as "ori.txt"

   _____________________
_
   ___________________
___|After_blast
   -------------------

This script serves two main purposes: to append the spacer sequence to each OTU for later analysis, and to remove all but the highest scoring hits for each OTU.  The first task is accomplished in a simple loop, taking data from the previous "Start" or "End" file.  For the next portion, even once the lower scores are removed, the same spacer is sometimes found in multiple locations.
[Note: we have no way to know which one is the actual Target, introducing noise into our data.  Perhaps it would be best to ignore these examples entirely, if there are not too many]
In these instances, the locations are combined so that each OTU is represented only once, using the following format (I refer to this process in the comments as "comma-combining", "colon-combining", and "semicolon-combining")-

genome: location, location; genome: location, location, location; genome: location[...]

Furthermore, the Length is recalculated (based on sStart  and sEnd) to indicate which strand the alignment is found on: positive (+) for the Top or negative (-) for the bottom.  This calculation also shortens the magnitude of the Length by 1, but the TrueLength is computed later based on the number of characters in the spacer sequence, so no information is lost.  

One complication with BLAST arises when the reported alignment does not start at the beginning of the spacer (or finish at the end of the spacer), which would be problematic when looking for the PAM.  To rectify this, a column titled "misalign" is added, which lists the distance between the reported sStart and the actual beginning of the spacer.  Luckily, this value is usually zero.

Next, all the data harvested from BLAST and refined with BLAST SCREEN is added to  the original tabfiles.  It also adds some other key information in the process: TrueLength, a flag for what I call "clipping", and a count for OTUs in which the spacer aligned to multiple locations.

The Length reported by BLAST is not always accurate, so it used only to create a new column listing which DNA strand the alignment occured on, as explained above.  Instead, TrueLength uses the spacer sequences found in the last script to count the actual length of each spacer.  Once this is created, the two lengths are compared, and if they disagree, then another column called "clipping" is listed as "YES".  This is used as a flag to indicate whether or not the full spacer length was aligned, and may or may not be used later.  Note that if misalign is >0, the clipping will always be flagged as "YES", but even if misalign =0, clipping may be "YES" OR "NO", as clipping can occur on either end or the spacer, while misalign only indicates the upstream portion.

Lastly, "Targets" is usually 1, but for OTUs that map to multiple locations in the genome, it reports the number of different places a match was found.

Rows which have no blast matches - Target is "NA"  - are removed.

Occasionally, different OTUs map to the same spacer in a given location, probably because small mutations cause the same spacer to be classified as a new OTU.  Here, any consecutive OTUs with the same location are detected, and only one is preserved.  However, the number of times that spacer was found per population (in the tabfile) is summed between the OTUs.  The OTU number itself, on the other hand, is arbitrary and only one is preserved.  All other information should be identical between the two OTUs.  
[NOTE: This could be another source of noise in our data, as two different OTUs could derive from the same spacer, but insertion or deletion could prevent us from knowing that they are actually the same.  InDels could also cause our TrueStart location to be slightly off later on, again with no way of automatically detecting it.  However, the clipping flag helps remedy this, and the incidence is low, so our data still comes out looking pretty good]

The result is saved as "final_0.99.txt"

   _____________________________
final_refinement_with_annotations.R. 


This script uses the locations of every found spacer, and searches for that location within the entirety of whatever genome the spacer came from.  To do this, the genomes are first hard-coded as string variables (another example where the file-paths will have to be modified for whatever machine is being used).  I gave all of these strings a 4-character name.  The "Location" (initially found in BLAST) is combined with "misalign" to yield the TrueStart of the spacer.  In reversed arrays, we are concerned with the opposite strand, so "trueStart" and "Strand" are updated to reflect this.  Finally, the trueStart is located within the genome of interest, also indicated by "Location", and the 10 base pairs upstream of this location are saved as the PAM (directionality of "upstream" is determined by Top vs. Bottom "Strand").  



   ___________
