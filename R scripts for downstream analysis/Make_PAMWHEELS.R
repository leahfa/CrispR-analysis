user="Urigo10"
library(stringr)
library(Biostrings)
library(Rcpp)
library(yaml)
library(sunburstR)		

path=paste0("C:/Users/",user,"/Dropbox/Crisp2018/integrated arrays by species/")
#path="C:/Users/Urigo10/Dropbox/Crisp2018/relooking_at_old_runs/Crisp2017/within_species/integrated arrays by species/"
tax="Vol"
dat=read.csv(paste0(path,tax,"_all_arrays_no_multiple_targets.csv"),stringsAsFactors = F)

#dat=dat[-grep("No_genome",dat$upstream),] #remove 1 spacer begining at location 1 (in vol)

####make subsets of dat####
summary(as.factor(dat$genome))
m=grep("med|pHM",dat$Location)
v=grep("vol|pHV",dat$Location)
intersect(m,v) #sanity check  -should be zero!!!
#dat=dat[which(dat$sum>1),]
all.targets=dat
all.med=dat[grep("med|pHM",dat$Location),] #med and all his plasmids
all.vol=dat[grep("vol|pHV",dat$Location),] #vol and all his plasmids
med.only=dat[grep("med",dat$Location),] #med 
vol.only=dat[grep("vol",dat$Location),]
phv4=dat[grep("pHV4",dat$Location),]
phv1=dat[grep("pHV1",dat$Location),]
phv3=dat[grep("pHV3",dat$Location),]

phm1=dat[grep("pHM1",dat$Location),]
phm3=dat[grep("pHM3",dat$Location),]
phm5=dat[grep("pHM5",dat$Location),]

temp=summary(as.factor(dat$genome))
barplot(temp)
write.csv(temp,paste0(path,tax," Num.Uniques.By.Target.csv"))#,sep="\,")
temp=tapply(dat$sum,INDEX=as.factor(dat$genome),FUN=sum)
barplot(temp)
write.csv(temp,paste0(path,tax," Num.Reads.by.target.csv"))#,sep="\,")

temp=tapply(dat$sum,INDEX=as.factor(dat$array),FUN=sum)
barplot(temp)
write.csv(temp,paste0(path,tax," Num.Reads.by.array.csv"))#,sep="\,")
targets=sets
num.spacers=c()
####make PAM percentage tables and wheel: ####
sets=c("all.targets","all.med","all.vol","med.only","vol.only","phv4","phv3","phv1","phm1","phm3","phm5")
#sets=c("vol.only","phv3","phv4","all.vol")
path2=paste0(path,"Final tables by species AND target")
for ( k in 1:length(sets)) {
set=sets[k]
PAM=eval(parse(text = set))
#print final table of dataset, showing all biological repeats
num.spacers[k]=nrow(PAM)
#change sign of spacers on bottom strand to minus: 
for (l in 1:nrow(PAM)) {
  if (PAM$Strand[l]=="Bottom") {
 
   PAM[l,2:4]=-1*PAM[l,2:4]
  #  PAM$Rep2[l]=-1*PAM$Rep2[l]
  #  PAM$Rep3[l]=-1*PAM$Rep3[l]
  }
}
write.csv(PAM, file = paste0(path2,"/species.",tax,"_targets.",set,".csv"), row.names = FALSE)
minus1 <- paste0(substr(PAM[, "upstream"], 10, 10), collapse = "")
minus2 <- paste0(substr(PAM[, "upstream"],  9,  9), collapse = "")
minus3 <- paste0(substr(PAM[, "upstream"],  8,  8), collapse = "")
minus4 <- paste0(substr(PAM[, "upstream"],  7,  7), collapse = "")

percent <- matrix(0, 4, 4)		#create an empty matrix to store the frequency of each letter, by position
dimnames(percent) = list(c("A", "C", "G", "T"), c("-4", "-3", "-2", "-1"))	#name the matrix dimensions

#each column of matrix is created by vectorized calculations on the "minusX" strings
percent["A",] <- str_count(c(minus4, minus3, minus2, minus1), "A") / nrow(PAM) * 100
percent["C",] <- str_count(c(minus4, minus3, minus2, minus1), "C") / nrow(PAM) * 100
percent["G",] <- str_count(c(minus4, minus3, minus2, minus1), "G") / nrow(PAM) * 100
percent["T",] <- str_count(c(minus4, minus3, minus2, minus1), "T") / nrow(PAM) * 100
percent									#show the table in the R console
print( rownames(percent)[apply(percent, 2, which.max)] )	#print which letters appear most, in order
print( apply(percent, 2, max) )					#display the percent (relative abundance) of these listed letters


####USE PAM in script PAMWHEEL####
seq<-PAM

DNA <- c("A", "C", "G", "T")
freq <- array(data = 0, dim = c(4, 4, 4), dimnames = list(DNA, DNA, DNA))

trips <- array(data = "", dim = c(4, 4, 4))	#"trips" for triples: all possible 3-letter combinations
for (x in 1:4) {
  trips[x,,] <- DNA[x]
}
for (y in 1:4) {
  trips[,y,] <- paste0(trips[,y,], DNA[y])
}
for (z in 1:4) {
  trips[,,z] <- paste0(trips[,,z], DNA[z])
}

for (i in 1:nrow(seq)) {		#for each PAM, navigate to appropriate cell in array, and tally by 1
  dim1 <- substr(seq[i, "upstream"], 8, 8)
  dim2 <- substr(seq[i, "upstream"], 9, 9)
  dim3 <- substr(seq[i, "upstream"], 10, 10)
  if (dim1!=""&dim2!=""&dim3!=""){ 
  freq[dim1, dim2, dim3] <- (freq[dim1, dim2, dim3] + 1)
  } else {print(paste0("missing base iteration ",i,":Location ",seq[i,"Location"]))}
}

percent <- cbind(as.vector(trips), as.vector(freq), character(length = 64))	#combine the tallies with their PAM sequences
colnames(percent) <- c("PAM", "frequency", "percentage")		#name the two columns
percent <- percent[percent[, "frequency"] != 0,]	#remove any PAMs that are not present at all
percent[, "percentage"] <- as.numeric(percent[, "frequency"]) / nrow(seq) * 100	#convert to %
percent=percent[order(as.numeric(percent[, 3]),decreasing=TRUE),]
#write.csv(percent, file = paste0(path2,"/PAM.PERCENT species.",tax,"_targets.",set,".csv"), row.names = FALSE)
print(paste(set, percent[1,3]))

wheel <- matrix("", nrow(percent), 2)		#manually create a matrix in the form of a csv file
colnames(wheel) <- c("bases", "frequency")	#first, copy the frequency values (after comma)
wheel[, "frequency"] <- percent[, "frequency"]
for (j in 1:nrow(percent)) {				#then create the hyphenated base combinations (before comma)
  wheel[j, "bases"] <- paste(substr(percent[j, "PAM"], 1, 1), substr(percent[j, "PAM"], 2, 2), substr(percent[j, "PAM"], 3, 3), sep = "-")
}

colorCode <- list(domain = c("A", "T", "G", "C"), range = c("3498DB", "F4D03F", "2ECC71", "E74C3C"))	#A-blue, T-yellow, G-Green, C-Red

#sn<-sunburst(as.data.frame(wheel), colors = colorCode,percent=T)#finally, coerce wheel matrix to data frame, and simply output as a sunburst
#sn <- htmlwidgets::prependContent(sn, htmltools::h1(paste0(tax," targeting ",sets[k]))          )
#sn <- htmlwidgets::prependContent(sn, htmltools::h3("Singeltons removed"))
#print(sn)
}
res=cbind(sets,num.spacers)
write.csv(res,paste0(path2,tax,"_filtered1_num.spacers.per.targetset"))
