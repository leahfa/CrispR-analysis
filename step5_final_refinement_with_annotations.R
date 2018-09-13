library(stringr)
library(Biostrings)
user="Urigo10"


arrays=c("A","B","C","D","F","G","H", "I")
path1=paste0("C:/Users/",user,"/Dropbox/Crisp2018/final_tables_by_array/")
####read in genomes for to get upstream region ####  
folder <-paste0("C:/Users/",user,"/Dropbox/Crisp2018/All_replicons/")
mediFile <- paste0(folder, "mediterani complete genome.txt")
volcFile <- paste0(folder, "vol chromosome.txt")
pHM100File <- paste0(folder, "pHM100.txt")
pHM300File <- paste0(folder, "pHM300.txt")
pHM500File <- paste0(folder, "pHM500.txt")
pHV1File <- paste0(folder, "pHV1.txt")
pHV3File <- paste0(folder, "pHV3.txt")
pHV4File <- paste0(folder, "pHV4.txt")
PRL3File <- paste0(folder,"PRL3seq.txt")
medi <- paste(readLines(mediFile), collapse = "")
volc <- paste(readLines(volcFile), collapse = "")
pHM1 <- paste(readLines(pHM100File), collapse = "")	##NOTE: shorten name to preserve 4-character format
pHM3 <- paste(readLines(pHM300File), collapse = "")	##NOTE: shorten name to preserve 4-character format
pHM5 <- paste(readLines(pHM500File), collapse = "")	##NOTE: shorten name to preserve 4-character format
pHV1 <- paste(readLines(pHV1File), collapse = "")
pHV3 <- paste(readLines(pHV3File), collapse = "")
pHV4 <- paste(readLines(pHV4File), collapse = "")

for (k in 1:length(arrays)){


####read in input path and file: ####
array=arrays[k]
path=paste0(path1,arrays[k],"/")
PAM=read.table(paste0(path,"0.99_final.txt"),sep="\t",header=T,stringsAsFactors = F,check.names = F)
genome <- substr(PAM[, "Location"], 1, 4)	#save first 4 char of Targets as genome list
upstream <- trueStart <- start <- character(length = nrow(PAM))	#create new vectors to be appended
for (i in 1:nrow(PAM)) {
  if (as.numeric(PAM[i, "Targets"]) == 1) {		#if there is only one location, extract number from first space to end of location
    start[i] <- substr(PAM[i, "Location"], regexpr(" ", PAM[i, "Location"])[1] + 1, nchar(PAM[i, "Location"]))
  } else {							#for multiple locations, extract number from first space to second space
    start[i] <- substr(PAM[i, "Location"], gregexpr(" ", PAM[i, "Location"])[[1]][1] + 1, gregexpr(" ", PAM[i, "Location"])[[1]][2] - 2)
  }
}		#then use the start (just extracted) and misalign to find TrueStart
trueStart <- ifelse( (PAM[, "Strand"] == "Top"), as.numeric(start) - as.numeric(PAM[, "misalign"]), as.numeric(start) + as.numeric(PAM[, "misalign"]) )
if (array=="A"|array=="E"|array=="G") {orientation="YES"} else {orientation="NO"}
print(paste(path,"array is", arrays[k],"orientation is", orientation))
if (orientation == "YES") {
  orientation <- "Reverse"
} else if (orientation == "NO") {
  orientation <- "Forward"
} #else { error: "Must choose an orientation" }

if (orientation == "Reverse") { #For reversed arrays, traverse to position at opposite end of the spacer, and swap strand
  trueStart <- ifelse( (PAM[, "Strand"] == "Top"), trueStart + nchar(PAM[, "Sequence"]) - 1, trueStart - nchar(PAM[, "Sequence"]) + 1 )
  PAM[, "Strand"] <- ifelse( (PAM[, "Strand"] == "Top"), "Bottom", "Top")
}

for (j in 1:nrow(PAM)) {
  if (genome[j]=="gnl|"|trueStart[j]<3) {upstream[j]="No_genome"} else{
    #loop through now almost completed matrix
    if (PAM[j, "Strand"] == "Top") {	#for top strand, extract 10 leading bases. for bottom, extract trailing 10, and reverseComplement
      upstream[j] <- substr(get(genome[j]), trueStart[j] - 10, trueStart[j] - 1)
    } else if (PAM[j, "Strand"] == "Bottom") {
      upstream[j] <- substr(get(genome[j]), trueStart[j] + 1, trueStart[j] + 10)
      upstream[j] <- toString(reverseComplement(DNAString(upstream[j])))
    } #else { stop }
  }
}
PAM1=cbind(PAM,genome,trueStart,upstream)
#### get REVERSE upstream (REVERSING orienation, and FIRST extrcating again teh original truestart): ####
reverse.upstream=c()
PAM=read.table(paste0(path,"0.99_final.txt"),sep="\t",header=T,stringsAsFactors = F,check.names = F) #repeat loading of ORIGINAL TABLE,
#beacuse, for revrse arrays, we swapped strands AND CHANGED THE PAM TABLE when extrcatin the "regular" PAM
genome <- substr(PAM[, "Location"], 1, 4)	#save first 4 char of Targets as genome list
upstream <- trueStart <- start <- character(length = nrow(PAM))	#create new vectors to be appended
for (i in 1:nrow(PAM)) {
  if (as.numeric(PAM[i, "Targets"]) == 1) {		#if there is only one location, extract number from first space to end of location
    start[i] <- substr(PAM[i, "Location"], regexpr(" ", PAM[i, "Location"])[1] + 1, nchar(PAM[i, "Location"]))
  } else {							#for multiple locations, extract number from first space to second space
    start[i] <- substr(PAM[i, "Location"], gregexpr(" ", PAM[i, "Location"])[[1]][1] + 1, gregexpr(" ", PAM[i, "Location"])[[1]][2] - 2)
  }
}		#then use the start (just extracted) and misalign to find TrueStart
trueStart <- ifelse( (PAM[, "Strand"] == "Top"), as.numeric(start) - as.numeric(PAM[, "misalign"]), as.numeric(start) + as.numeric(PAM[, "misalign"]) )
if (array=="A"|array=="E"|array=="G") {orientation="NO"} else {orientation="YES"} #this is the REVERSE orientation!!!
print(paste("array is", arrays[k],"REVERSE orientation is", orientation))
if (orientation == "YES") {
  orientation <- "Reverse"
} else if (orientation == "NO") {
  orientation <- "Forward"
} #else { error: "Must choose an orientation" }

if (orientation == "Reverse") { #For reversed arrays, traverse to position at opposite end of the spacer, and swap strand
  trueStart <- ifelse( (PAM[, "Strand"] == "Top"), trueStart + nchar(PAM[, "Sequence"]) - 1, trueStart - nchar(PAM[, "Sequence"]) + 1 )
  PAM[, "Strand"] <- ifelse( (PAM[, "Strand"] == "Top"), "Bottom", "Top")
}

for (j in 1:nrow(PAM)) {
  if (genome[j]=="gnl|"|trueStart[j]<3) {upstream[j]="No_genome"} else{
    #loop through now almost completed matrix
    if (PAM[j, "Strand"] == "Top") {	#for top strand, extract 10 leading bases. for bottom, extract trailing 10, and reverseComplement
      upstream[j] <- substr(get(genome[j]), trueStart[j] - 10, trueStart[j] - 1)
    } else if (PAM[j, "Strand"] == "Bottom") {
      upstream[j] <- substr(get(genome[j]), trueStart[j] + 1, trueStart[j] + 10)
      upstream[j] <- toString(reverseComplement(DNAString(upstream[j])))
    } #else { stop }
  }
}

reverse.upstream<-upstream
####put it all toghethr and filter multiple targets:####
PAM.final <- cbind(PAM1,reverse.upstream)	#combine created info, and output desired info
PAM.f=PAM.final[which(PAM.final$Targets==1),]
multiples=PAM[which(PAM$Targets>1),]


#### read in annotaion files (SAME NAMES AS GENOMES)####
folder <- paste0("C:/Users/",user,"/Dropbox/Sam 10.1.17/new_annotation_files/")
mediFile <- paste0(folder, "medi annotation.txt")
volcFile <- paste0(folder, "volc annotation.txt")
pHM100File <- paste0(folder, "pHM100 annotation.txt")
pHM300File <- paste0(folder, "pHM300 annotation.txt")
pHM500File <- paste0(folder, "pHM500 annotation.txt")
pHV1File <- paste0(folder, "pHV1 annotation.txt")
pHV3File <- paste0(folder, "pHV3 annotation.txt")
pHV4File <- paste0(folder, "pHV4 annotation.txt")

medi <- as.matrix(read.table(mediFile, sep = "\t", quote = "\"", header = TRUE))
volc <- as.matrix(read.table(volcFile, sep = "\t", quote = "\"", header = TRUE))
pHM1 <- as.matrix(read.table(pHM100File, sep = "\t", quote = "\"", header = TRUE))	##NOTE: shorten name to 4-character format
pHM3 <- as.matrix(read.table(pHM300File, sep = "\t", quote = "\"", header = TRUE))	##NOTE: shorten name to 4-character format
pHM5 <- as.matrix(read.table(pHM500File, sep = "\t", quote = "\"", header = TRUE))	##NOTE: shorten name to 4-character format
pHV1 <- as.matrix(read.table(pHV1File, sep = "\t", quote = "\"", header = TRUE))
pHV3 <- as.matrix(read.table(pHV3File, sep = "\t", quote = "\"", header = TRUE))
pHV4 <- as.matrix(read.table(pHV4File, sep = "\t", quote = "\"", header = TRUE))
#### get annotations for PAM.f: ####
genes=PAM.f
genes=genes[order(genes$genome,genes$trueStart),]
tag <- description <- character(length = nrow(genes))		#2 new vectors for the desired annotation info

Fnum <- 1	#for "file number" used to cycle through the file of interest for each spacer
Gnum <- 1	#for "gene number" for each row of the PAM file to add the annotation info 
file <- as.character(genes[1, "genome"])		#the file of interest
while (Gnum <= nrow(genes)) {		#run loop until all spacers have been annotated
  if (genes[Gnum, "genome"] != file) {	#spacers order by genome, the trueStart
    file <- as.character(genes[Gnum, "genome"])		#When get to new genome in the list,
    Fnum <- 1					#update the file and Fnum index accordingly
  }
  #Check if the spacer trueStart falls in range of the annotated gene (removing < or > signs), or is inter-genic
  if ( as.numeric(genes[Gnum, "trueStart"]) <= as.numeric(sub(">", "", get(file)[Fnum, "geneEnd"])) ) {
    if ( as.numeric(genes[Gnum, "trueStart"]) < as.numeric(sub("<", "", get(file)[Fnum, "geneStart"])) ) {
      tag[Gnum] <- paste("Inter-genic, between", get(file)[Fnum-1, "tag"], "and", get(file)[Fnum, "tag"])
      description[Gnum] <- paste("close to", get(file)[Fnum-1, "description"], "and", get(file)[Fnum, "description"])
      Gnum <- (Gnum+1)
    } else { #if at correct range (or on Top Strand), save the tag and description info
      tag[Gnum] <- get(file)[Fnum, "tag"]
      description[Gnum] <- get(file)[Fnum, "description"]
      Gnum <- (Gnum+1)
    } #if not in range, advance to next annotated gene
  } else if ( Fnum == nrow(get(file)) ) {		#if the end of the genome is reached, store the final gene
    tag[Gnum] <- get(file)[Fnum, "tag"]		#this is only likely for the final 1-2 spacers of a genome
    description[Gnum] <- get(file)[Fnum, "description"]			
    Gnum <- (Gnum+1)
  } else { Fnum <- (Fnum+1) }	#if the correct gene isn't found yet, advance to the next one
}

res=cbind(genes,tag, description)


####add in some extras:####
z=summary(as.factor(res$genome))
jpeg(paste0(path,"Num.spacers_by_target.jpeg"),width=800,height=600)
barplot(z,ylab="Num. Unique Spacers",main=paste("array",arrays[k]))
dev.off()
z=as.data.frame(z)
write.csv(z,paste0(path,arrays[k],"_Num.unique.spacers.by.target.csv"))
ind=grep("Location",colnames(res))-1 #index number of last sample collumn (depending on number of samples)
res$sum=rowSums(res[ ,2:ind])
temp=res[ ,2:ind]
temp.b=ifelse(temp>0,1,0)
res$Prevalence=rowSums(temp.b)
#### save tables with upstream and annotations: ####
write.csv(multiples,paste0(path,arrays[k],"_multiple_targets.csv"),row.names=F)
write.csv(res,paste0(path,arrays[k],"_no_multiple_targets.csv"),row.names=F)

}












