user="Urigo10"
library(stringr)
library(Biostrings)
library(Rcpp)
library(yaml)
library(sunburstR)		

path=paste0("C:/Users/",user,"/Dropbox/Crisp2018/integrated arrays by species/")
#path="C:/Users/Urigo10/Dropbox/Crisp2018/relooking_at_old_runs/Crisp2017/within_species/integrated arrays by species/"
tax="Vol"
if (tax=="Vol") {
  dat=read.csv(paste0(path,tax,"_all_arrays_no_multiple_targets.csv"),stringsAsFactors = F)
  dat=dat[-c(1,25485),]
                }

if (tax=="Med")  {
dat=read.csv(paste0(path,tax,"_all_arrays_no_multiple_targets_FIXED_EXTREME_VALUES.csv"),stringsAsFactors = F)
dat$sum=rowSums(dat[,2:4])
                }




####make subsets of dat####
summary(as.factor(dat$genome))
m=grep("med|pHM",dat$Location)
v=grep("vol|pHV",dat$Location)
intersect(m,v) #sanity check  -should be zero!!!
#dat=dat[which(dat$sum>1),]
all.targets<-dat
all.med=dat[grep("med|pHM",dat$Location),] #med and all his plasmids
all.vol=dat[grep("vol|pHV",dat$Location),] #vol and all his plasmids
med.chromosome=dat[grep("med",dat$Location),] #med 
vol.chromosome=dat[grep("vol",dat$Location),]
phv4=dat[grep("pHV4",dat$Location),]
phv1=dat[grep("pHV1",dat$Location),]
phv3=dat[grep("pHV3",dat$Location),]

phm1=dat[grep("pHM1",dat$Location),]
phm3=dat[grep("pHM3",dat$Location),]
phm5=dat[grep("pHM5",dat$Location),]

num.spacers=c()

sets=c("med.chromosome","vol.chromosome","phv4","phv3","phv1","phm1","phm3","phm5")

path2=paste0(path,"Final tables by species AND target")
res=list()
for ( k in 1:length(sets)) {
set=sets[k]
####make vectors for PAM summing, abundance weighing and prevalence weighing ####
seq<-eval(parse(text = set))
pam.vec=c()
count.vec=c()
prev.vec=c()

for (i in 1:nrow(seq)){
  pam.vec[i]=substr(seq[i,"upstream"],8,10)
  count.vec[i]=seq$sum[i]
  prev.vec[i]=seq$Prevalence[i]
}


####Weighing eahc spacer by its total abundance:####
temp=as.data.frame(cbind(pam.vec,count.vec),stringsAsFactors=F)
temp$count.vec=as.numeric(temp$count.vec)
temp1=aggregate(temp[ ,2],by=list(temp$pam.vec),FUN=sum)
sum(temp1$x)==sum(seq$sum)
temp1$percentage=100*(round(temp1$x/sum(temp1$x),3))
sum(temp1$percentage)
colnames(temp1)=c("PAM","frequency","percentage")
temp1=temp1[order(temp1$frequency,decreasing = T),]
temp1$Species=rep(tax,nrow(temp1))
temp1$Target=rep(set,nrow(temp1))
percent<-temp1
print(sum(percent$frequency)==sum(seq$sum))
percent=percent[1:5,]
res[[k]]<-percent
}
final=do.call("rbind",res)
final=final[ ,c(4,5,1:3)]
write.csv(final, paste0(path,"PAMS_for_revisions/",tax,"_sum_PAM_by_target.csv"),row.names=F)

####counting each PAM from each sapcer once regardless of abundane or prevalence (equivalent to Sams original script):####

levels(as.factor(pam.vec))
df.pams=as.data.frame(summary(as.factor(pam.vec)), stringsAsFactors=F)
df.pams$PAM=rownames(df.pams)
df.pams$percentage=sapply(df.pams[ ,1],function(x) 100*(x/nrow(seq)))
colnames(df.pams)[1]="frequency"
df.pams=df.pams[,c(2,1,3)]

df.pams=df.pams[order(df.pams[ ,2],decreasing = T),]
percent<-df.pams

####Weighing eahc spacer by its PREVALENCE (and disregrarding abundance):####
temp=as.data.frame(cbind(pam.vec,prev.vec),stringsAsFactors=F)
temp$prev.vec=as.numeric(temp$prev.vec)
temp1=aggregate(temp[ ,2],by=list(temp$pam.vec),FUN=sum)
sum(temp1$x)==sum(seq$Prevalence)
temp1$percentage=100*(temp1$x/sum(temp1$x))
sum(temp1$percentage)
colnames(temp1)=c("PAM","frequency","percentage")
temp1=temp1[order(temp1$frequency,decreasing = T),]
percent<-temp1
write.csv(percent, paste0(path,"PAMS_for_revisions/temp1.csv"))

#### drawing pamwheels from 'percent' object ####

wheel <- matrix("", nrow(percent), 2)		#manually create a matrix in the form of a csv file
colnames(wheel) <- c("bases", "frequency")	#first, copy the frequency values (after comma)
wheel[, "frequency"] <- percent[, "frequency"]
for (j in 1:nrow(percent)) {				#then create the hyphenated base combinations (before comma)
  wheel[j, "bases"] <- paste(substr(percent[j, "PAM"], 1, 1), substr(percent[j, "PAM"], 2, 2), substr(percent[j, "PAM"], 3, 3), sep = "-")
}

colorCode <- list(domain = c("A", "T", "G", "C"), range = c("3498DB", "F4D03F", "2ECC71", "E74C3C"))	#A-blue, T-yellow, G-Green, C-Red

sunburst(as.data.frame(wheel), colors = colorCode,percent=T)#finally, coerce wheel matrix to data frame, and simply output as a sunburst


}

res=cbind(sets,num.spacers)
write.csv(res,paste0(path2,tax,"_filtered1_num.spacers.per.targetset"))
