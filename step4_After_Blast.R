
library(tools)			#basic library for filename manipulation
library(Biostrings)	
library(utils)	#used for reading in excel txt file as a table
library(stringr)	#simple string manipulations (specifically str_trim)

path=paste0("/scratch200/leah/Crisp2018/A/")
id=0.99
array="A"
bl.table="blast_0.99.txt"
sp.file="all_spacers_0.99.txt"
BLAST <- as.matrix(read.table(paste0(path,bl.table), strip.white = TRUE,sep="\t"))	#convert table to easier matrix
BLAST <- BLAST[, -c(3,6,11)]							#remove percentage, gaps, and e-value
BLAST <- cbind(BLAST, matrix("", nrow(BLAST), 2))		#add empty column to matrix, for spacer
colnames(BLAST) <- c("OTU", "Targets", "Length", "Errors", "qStart", "qEnd", "sStart", "sEnd", "Score", "Sequence", "misalign")
BLAST[, c(5:8)] <- str_trim(BLAST[, c(5:8)])		#strip.white doesn't carry over to matrix from table, so repeat here

#length is small in magnitude by 1, but only used to ascertain direction.  trueLength is from nchar(BLAST[, "Sequence"])
BLAST[, "Length"] <- (as.numeric(BLAST[, "sEnd"]) - as.numeric(BLAST[, "sStart"]))	#add sign to indicate direction

spacers <- scan(file = paste0(path,sp.file), what = "character", quiet = TRUE)		#read in spacers as vector
spacers <- cbind(spacers[c(TRUE, FALSE)], spacers[c(FALSE, TRUE)])	#recombine spacers as matrix, with OTU column and sequence column
colnames(spacers) <- c("OTU", "Sequence")						#give header names for spacers matrix
spacers[, "OTU"] <- substr(spacers[, "OTU"], 2, nchar(spacers[, "OTU"]))	#remove 'greater than' sign to match format of BLAST matrix

Bnum <- 1		#for BLAST number: used to index elements in BLAST
Snum <- 1		#for spacers number: used to index elements in spacers
while (Bnum <= nrow(BLAST)) { #compare OTUs. if match, add sequence and step Bnum but not Snum.  if not, step Snum but not Bnum
	if (BLAST[Bnum, "OTU"] == spacers[Snum, "OTU"]) {
		BLAST[Bnum, "Sequence"] <- spacers[Snum, "Sequence"]
		Bnum <- (Bnum+1)
	} else { Snum <- (Snum+1) }
}

#misalignment associated with FIRST location in list, (though often applies to all locations)
BLAST[, "misalign"] <- (as.numeric(BLAST[, "qStart"]) - 1)
##PREVIOUS FATAL ERROR: for top strand, did above. for bottom strand, took nchar (ie: trueLength) - qEnd
print("done matching otus and blasts")

highScore <- as.numeric(BLAST[1, "Score"])		#first score for each OTU is highest.  Reset for each new OTU
for (i in 1:(nrow(BLAST)-1)) {
	if (BLAST[i+1, "OTU"] == BLAST[i, "OTU"]) {			#if the next OTU is the same
		if (as.numeric(BLAST[i+1, "Score"]) < highScore) {	#compare the scores
			BLAST[i+1, "Score"] <- 0				#remove the lower score
			BLAST[i, "sStart"] <- paste(BLAST[i, "Targets"], BLAST[i, "sStart"], sep = ": ")	#lower Score, colon-combine and remove old Target
			BLAST[i, "Targets"] <- 1
		} else if (BLAST[i+1, "Targets"] == BLAST[i, "Targets"]) {								#for equal Scores and like Targets,
			BLAST[i+1, "sStart"] <- paste(BLAST[i+1, "sStart"], BLAST[i, "sStart"], sep = ", ")		#comma-combine the Locations 
			BLAST[i, "sStart"] <- ""												#and remove the old ones
		} else {
			BLAST[i, "sStart"] <- paste(BLAST[i, "Targets"], BLAST[i, "sStart"], sep = ": ")	#same Score but new Target, colon-combine
			BLAST[i, "Targets"] <- 1
		}
	} else {
		highScore <- as.numeric(BLAST[i+1, "Score"])	#set new highScore when finding new OTU
		BLAST[i, "sStart"] <- paste(BLAST[i, "Targets"], BLAST[i, "sStart"], sep = ": ")		#New OTU... colon-combine
		BLAST[i, "Targets"] <- 1
	}
	if ( abs(as.numeric(BLAST[i, "Length"])) < nchar(BLAST[i, "Sequence"])-4 ) { #after all scoring, simply check length. ditch segments clipped by >3
		BLAST[i, "Score"] <- 0
	}
}
BLAST[nrow(BLAST), "sStart"] <- paste(BLAST[nrow(BLAST), "Targets"], BLAST[nrow(BLAST), "sStart"], sep = ": ")	#same colon-combination as in loop,
BLAST[nrow(BLAST), "Targets"] <- 1														#to fix end-discrepancy of i+1 indexing
if ( abs(as.numeric(BLAST[nrow(BLAST), "Length"])) < nchar(BLAST[nrow(BLAST), "Sequence"])-4 ) { 			#also do score-check
	BLAST[nrow(BLAST), "Score"] <- 0
}

BLAST <- BLAST[BLAST[, "Score"] != 0,]	#remove emptied rows
BLAST <- BLAST[BLAST[, "sStart"] != "",]	#for score and location
print("done summing blast by highet score for each otu")
write.table(BLAST, file = paste0(path,id,"_Blast_refined.txt"), sep="\t",row.names = FALSE)
#file.show(paste0(path,"/R/outputs/otus_95 output_1Blast.txt"))

#STEPS
#1) run desired OTU through Spacer Extractor vB6.R, defining which array, and choosing the file, number of spacers, and leader position
#2) Run the output(s) through local blast, with combined genome database.fasta, e-value of .0001, in tabular output, default otherwise
#3) Copy results into an excel spreadsheet.  Save as .txt, changing file name from "Start" or "End" to "Table"
#4) Run new table output through BLAST screen V4.R, pressing "Okay", then choosing the Table.txt file, then the same End.txt or Start.txt file
#5) File will open in R window.  Copy to new sheet of excel, text to columns under data tab, deliminated by space.  Save as ...BLAST.txt
#6) Next run Cluster Decomposer v3.R, first choosing the new BLAST.txt file, then the associated tabfile, output shows in R info table.
#7) Copy to excel, sort Targets, remove NAs, sort Location, sort reshef columns by row1. compile all outputs per letter. combine/re-sort medi vs volc
#8) Save excel file of all sheets (_All) , and .txt copy of medi/volc to run through sum tabfile.R, (sorted.txt -> Done.txt) re-organizing list
#9) Add Done back to _All in new sheets.  use these to FIND PAM! Sort PAM by genome, THEN trueStart (save as PAM.txt)
#10) Run the new PAM.txt files through PAM wheel to find relative percentages AND wheel of each PAM (copy back to PAM excel sheet, sorted by percentage)
#11) For coverage analysis, combine Done.txt files by array.  Add new column to sum abolute wildtype abundance.  delete rows = 0. calculate relative
#... Then count singletons and total rows remaining.  use to calculate coverage.  combine all arrays into abs and rel tables

####cluster_decomposer:
clusters <- BLAST

tabfile <- read.table(paste0(path,"tabfile_0.99.txt"), header = TRUE, sep = "\t",stringsAsFactors = F)	#Leah's file, for comparison

for (i in 1:(nrow(clusters)-1)) {		#then loop again
  if (clusters[i+1, "OTU"] == clusters[i, "OTU"]) {	#when consecutive OTUs match, semi-colon combine
    clusters[i+1, "sStart"] <- paste(clusters[i+1, "sStart"], clusters[i, "sStart"], sep = "; ")
    clusters[i, "sStart"] <- ""
  } else {						#otherwise, count the number of repeats (1 + commas + semicolons)
    clusters[i, "Targets"] <- (1 + str_count(clusters[i, "sStart"], ",") + str_count(clusters[i, "sStart"], ";"))
  }
}
clusters <- clusters[clusters[, "sStart"] != "",]	#again, remove rows with emptied location

newcols <- matrix(NA, nrow(tabfile), 7)				#create 7 column matrix for appending desired information
colnames(newcols) <- c("Location", "Targets", "Sequence", "trueLength", "Strand", "misalign", "clipping")	#name the new columns
trueLength <- nchar(clusters[, "Sequence"])		#lenght of the full cluster, not just the aligned portion
clipping <- ifelse(trueLength - abs(as.numeric(clusters[, "Length"])) > 1 , "Yes", "No") #for clipping on either end, not just PAM side
clusters[, "Length"] <- ifelse(as.numeric(clusters[, "Length"]) < 0, "Bottom", "Top") #negative lengths are bottoms strand

print("done filling 'clusters' with Start and cliiping information")
write.table(clusters,paste0(path,id,"_0.99_clusters.txt"),sep="\t", row.names=FALSE)
for (j in 1:(nrow(tabfile))) {					#search every item of tabfile
  substring <- substr(tabfile[j, "OTUId"], 1, (regexpr(";", tabfile[j, "OTUId"])-1))	#store the OTU of the current tabfile element
  for (k in 1:(nrow(clusters))) {
    name=clusters[k, 1] #beacuse otu_nams in clusters show iteration
    if (substring ==name ) {				#screen through clusters until match is found
      newcols[j, c(1:3)] <- clusters[k, c("sStart", "Targets", "Sequence")]	#then addinfo to appropriate row in newcols
      newcols[j, 4] <- trueLength[k]
      newcols[j, c(5,6)] <- clusters[k, c("Length", "misalign")]
      newcols[j, 7] <- clipping[k]
      break											#after match, break to next OTU in tabfile
    }
  }
}
fin <- cbind(tabfile, newcols,stringsAsFactors=F)	#concatenate the new info to the tabfile in ONE LAST MATRIX
print("done matching clusters to tabfile")
#remove NAs:
fin.f=fin[which(is.na(fin[ ,ncol(fin)])==F),]
colnames(fin.f)=gsub(".assmebled.fastq","",fixed=T, colnames(fin.f))
#sort by location:
loc=which(colnames(fin.f)=="Location")
fin.f=fin.f[order(fin.f[ ,loc]),]
final=as.data.frame(fin.f,stringsAsFactors=F)
write.table(final, file = paste0(path,id,"_TempFinal.txt"), row.names = FALSE,sep="\t")

for (i in 1:(nrow(final)-1)) { #check successive rows for match in location
  if (final[i+1, "Location"] == final[i, "Location"]) {
    for (j in 2:(ncol(final)-7)) { #sum for each numeric column, 7 columns appended
      final[i+1, j] <- (as.numeric(final[i+1, j]) + as.numeric(final[i, j]))
    }
    final[i, "Targets"] <- 0 #after summing hits (not Targets), set previous to zero
  }
}	#only one OTU need be preserved
final <- final[final[, "Targets"] != 0,]		#remove emptied rows based on Targets
write.table(final,paste0(path,id,"_final.txt"),sep="\t",row.names=F)











