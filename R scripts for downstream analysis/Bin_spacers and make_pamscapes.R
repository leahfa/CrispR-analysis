user="Urigo10"
library(Biostrings)
library(ggplot2)
library(reshape2)
####read in genomes for to get upstream region ####  
folder <-paste0("C:/Users/",user,"/Dropbox/Crisp2018/new_files_from_sam/AsFasta/")
mediFile <- paste0(folder, "mediterani complete genome.fna")
volcFile <- paste0(folder, "vol chromosome.fna")
pHM100File <- paste0(folder, "pHM100.txt")
pHM300File <- paste0(folder, "pHM300.txt")
pHM500File <- paste0(folder, "pHM500.txt")
pHV1File <- paste0(folder, "pHV1.txt")
pHV3File <- paste0(folder, "pHV3.txt")
pHV4File <- paste0(folder, "pHV4.fna")
PRL3File <- paste0(folder,"PRL3seq.txt")

#med=readDNAStringSet(mediFile)
#vol=readDNAStringSet(volcFile)
phv4=readDNAStringSet(pHV4File)
x=phv4[[1]]
####read in relevant table and add in pAM informatuin: ####
path=paste0("C:/Users/",user,"/Dropbox/Crisp2018/integrated arrays by species/Final tables by species AND target/")
filename="species.med_targets.phv4.csv"
#filename="species.med_targets.med.only_fixed_extreme_values.csv"
dat=read.csv(paste0(path,filename),stringsAsFactors = F)
#for Vol only:
#dat=dat[-1,]
####without binning spacers:####
temp=merge(dat,pamscape,by.x="trueStart",by.y="x.loc",all=T)
temp[is.na(temp)]<-0
temp->sp.not.binned
####with biining spacers same as pams:####
####choose window:####
trys=seq(from=2500,to=2550,by=2)
for (q in 1:length(trys)){
n=trys[q]#this covers almost entore genome - #1215 for med worked well (up tp 80 last bases, 1218 for vol)
round(length(x)/n)->window.size #fix manually if needed window size    - it needs to be always rounded DOWN and not up
window.size
#window.size=2458 #in this case , am rounding it up manually
r=c()
r[1]=1
for (i in 1:(n-1)){
  r[i+1]=r[i]+window.size
}
at1 <- IRanges(r, width=window.size)
at1@start[n]+at1@width[n]
print(paste(q,trys[q],(length(x)-(at1@start[n]+at1@width[n]))))
}

####use SAME range to extract PAMs and rev PAMS:####
n=2430/2
round(length(x)/n)->window.size #fix manually if needed window size    - it needs to be always rounded DOWN and not up
window.size=2427
n=length(x)/window.size #or set n according to window
#window.size=2428 #in this case , am rounding it up manually

r=c()
r[1]=1
for (i in 1:(n-1)){
  r[i+1]=r[i]+window.size
}
at1 <- IRanges(r, width=window.size)
at1@start[n]+at1@width[n]   #this is where my range ends!
nchar(x)-(at1@start[n]+at1@width[n]  ) #This amount of bases is "left over"
z=extractAt(x,at1)



#### MAKE PAMSCAPEs ####
#run twice , fro forward and revesre PAMS
PAM= "TTC" #for rev, use GAA. For med: TTC, for rev:GAA
strand="For"
res=c()
for (i in 1:n){
  q=z[[i]]
  countPattern(PAM, q,max.mismatch = 0, with.indels = FALSE, fixed = TRUE, algorithm = 'auto')->res[i]
  
}
#jpeg(paste0(path,"TAC in med genome_reverse strand.jpeg"),width=800,height=600)
#barplot(res,main=paste0("Mediterranei Chromosome, window-size= ", window.size),ylab="TAC frequency on reverse strand ",axis.lty =1)
#dev.off()
if (strand=="For") {res->pams} else {res->pams.rev}   #rev#or to pams.rev

# grab middle of each bin as a X-axis location: 



#### Start working on teh spacers file! ####

dat=dat[order(dat$trueStart),]
hist(dat$trueStart,breaks=n) #to get the same picture the long way,first switch to dat.b so we'll work on uniques:
ch=c("Rep1","Rep3","Rep2","sum")
dat.b<-dat
for (i in 1:length(ch)){
dat.b[ ,ch[i]]=ifelse(dat[ ,ch[i]]!=0,1,0)
}

hist(dat.b$trueStart,breaks=n) #to get the same picture 
#must spluit dat to Top and Bottom:
dat1=dat.b[which(dat.b$Strand=="Top"),]
dat2=dat.b[which(dat.b$Strand=="Bottom"),]
#### for total spacers use dat insteda of dat.b ####:
dat1=dat[which(dat$Strand=="Top"),]
dat2=dat[which(dat$Strand=="Bottom"),]

####bin spacers by range. first on dat1 and then on dat2:####
dat.x<-dat1  #enter dat1 OR dat 2 into dat.x

bin.counter=1
bin.name=c()
bin.loc=c()
sp.counter=1
while(sp.counter<=nrow(dat.x)) {
   while (dat.x$trueStart[sp.counter]>=at1@start[bin.counter] & dat.x$trueStart[sp.counter]<=(at1@start[bin.counter]+window.size-1))
   { bin.name[sp.counter]=paste0("Bin",bin.counter)
     bin.loc[sp.counter]=at1@start[bin.counter]+0.5*window.size
   sp.counter=sp.counter+1
   #print("exited 1st while loop")
   #print(paste0("bin counter was", bin.counter))
   
   }
  bin.counter=bin.counter+1
} 

sp.counter
bin.counter
nrow(dat.x)
dat.x=dat.x[1:(sp.counter-1),] #for rev - had to manually remove last row, spacer falls outside las bin
#IF necessay, add in manually spacers belonging to bin n+1:
#bin.name[sp.counter:nrow(dat)]="Bin2401"
#bin.loc[31806:nrow(dat)]=2947201+0.5*(length(x)-2947201)
#OR - just cut last spacers from dat, if its very few:
#dat.x=dat.x[1:16222,]
#key.loc=as.data.frame(cbind(unique(bin.name),unique(bin.loc)),stringsAsFactors=F)

dat.x$bin.name=bin.name
dat.x$bin.loc=bin.loc
temp=aggregate(dat.x[,c(2:4)],by=list(dat.x$bin.loc),sum)



####merge temp and pamscape: ####

#temp1$bin.loc=key.loc$V2[match(temp1$Group.1,key.loc$V1)]
#temp2$bin.loc=key.loc$V2[match(temp2$Group.1,key.loc$V1)]
final.F=merge(temp,pamscape.F,by.x="Group.1",by.y="x.loc",all=T)
final.F[is.na(final.F)]<-0
colnames(final.F)[1]<-"bin.loc"
final.F$bin.loc=as.numeric(final.F$bin.loc)
final.F=final.F[order(final.F$bin.loc),]
sum(rowSums(final.F[ ,2:4]))
sum(rowSums(dat1[,2:4]))

final.R=merge(temp,pamscape.R,by.x="Group.1",by.y="x.loc",all=T)
final.R[is.na(final.R)]<-0



colnames(final.R)[1]<-"bin.loc"


final.R$bin.loc=as.numeric(final.R$bin.loc)

final.R=final.R[order(final.R$bin.loc),]
sum(rowSums(final.R[ ,2:4]))
sum(rowSums(dat2[,2:4]))


write.csv(final.F,paste0(path,"species.med_target.med_TOP.Strand_",n,"_bins.csv"),row.names = F)

write.csv(final.R,paste0(path,"species.med_target.med_BOTTOM.Strand_",n,"_bins.csv"),row.names = F)
####merge the 2 tables ####

finT=merge(final.F,final.R,by="bin.loc",all=T)
#change finT "y" collumns (originating in reverse strand) to minus
#write.csv(finT,paste0(path,"temp.csv"),row.names = F)
toMinus=apply(finT[,c(6:8)],2,function(x) -1*x)
finT[ ,c(6:9)]<-toMinus
write.csv(finT,paste0(path,"Total_spacers_Med_target.phv4_BOTH.Strands_total_spacers_",round(n),"_bins.csv"),row.names = F)
#sneak a peak using ggplot:
dat.m=melt(finT,id.vars = c("bin.loc"))
ggplot(dat.m,aes(x=bin.loc,y=value,fill=variable))+
  geom_bar(stat="identity")       +
  ylim(0,100)

#### find locations of TACs:####
vmatchPattern("TAC", toString(x))->q
start=unlist(q@ends)-2
vmatchPattern("GTA", toString(x))->q
start.rev=unlist(q@ends)

#1st assign minus sign to "Bottom" strands, beacuse information will be lost in binning!
for (i in 1:nrow(dat)){
  if (dat$Strand[i]=="Bottom") {dat[i,2:4]=-1*(dat[i,2:4])}
}
