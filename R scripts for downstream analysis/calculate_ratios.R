user="Urigo10"
#path=paste0("C:/Users/",user,"/Dropbox/Crisp2018/relooking_at_old_runs/Crisp2017/final_tables_by_array/")
path=paste0("C:/Users/",user,"/Dropbox/Crisp2018/final_tables_by_array/")
#path=paste0("C:/Users/",user,"/Dropbox/Crisp2018/relooking_at_old_runs/isra3/final_tables_by_array/")
arrays=c("A","B","C","D","I","F","G","H") # add in E which has only 1 spacer by hand
num.spacers=c()
num.original.array=c()
ratio.sp.to.or=c()
reads.in.spacers=c()
reads.in.original.array=c()
ratio.reads.sp.to.reads.or=c()

for (i in 1:length(arrays)){
 
  fin=read.table(paste0(path,arrays[i],"/","0.99_final.txt"),stringsAsFactors = F,header=T)
  sp=fin$OTUId
  num.spacers[i]=length(sp)
  or=scan(paste0(path,arrays[i],"/ori_0.99.txt"),what="character")
  num.original.array[i]=length(or)
  tab=read.table(paste0(path,arrays[i],"/","tabfile_0.99.txt"),stringsAsFactors = F,header=T)
  tab$sum=rowSums(tab[ ,2:ncol(tab)])
  ratio.sp.to.or[i]=round(length(sp)/length(or),5) #num otus with real new spacers/num original array
  reads.in.spacers[i]=sum(fin[,2:(ncol(fin)-7)])
    #sum(tab$sum[which(tab$OTUId %in% sp)])
  reads.in.original.array[i]=sum(tab$sum[which(tab$OTUId %in% or)])
  ratio.reads.sp.to.reads.or[i]=round(reads.in.spacers[i]/reads.in.original.array[i],5)
  }
res=cbind(num.spacers,num.original.array,ratio.sp.to.or,reads.in.spacers,reads.in.original.array,ratio.reads.sp.to.reads.or)
res=cbind(arrays,res)

write.csv(res,paste0(path,"ratio_table.csv"),row.names=F)

#### 10.7.2018: rerun while REMOVING hits mapping to more tahn one target: ####
for (i in 1:length(arrays)){
  
  fin=read.table(paste0(path,arrays[i],"/","0.99_final.txt"),stringsAsFactors = F,header=T)
  # to identify multiple locations, count number of ":" in dat$Location:
  fin=fin[which(fin$Targets==1),]
 
  sp=fin$OTUId
  print(length(sp)==length(unique(sp)))
  num.spacers[i]=length(sp)
  or=scan(paste0(path,arrays[i],"/ori_0.99.txt"),what="character")
  num.original.array[i]=length(or)
  print(length(or)==length(unique(or)))
  tab=read.table(paste0(path,arrays[i],"/","tabfile_0.99.txt"),stringsAsFactors = F,header=T)
  tab$sum=rowSums(tab[ ,2:ncol(tab)])
  ratio.sp.to.or[i]=round(length(sp)/length(or),5) #num otus with real new spacers/num original array
  reads.in.spacers[i]=sum(fin[,2:(ncol(fin)-7)])
  #sum(tab$sum[which(tab$OTUId %in% sp)])
  reads.in.original.array[i]=sum(tab$sum[which(tab$OTUId %in% or)])
  ratio.reads.sp.to.reads.or[i]=round(reads.in.spacers[i]/reads.in.original.array[i],5)
}
res=cbind(num.spacers,num.original.array,ratio.sp.to.or,reads.in.spacers,reads.in.original.array,ratio.reads.sp.to.reads.or)
res=cbind(arrays,res)

write.csv(res,paste0(path,"ratio_table no_multiple_targets.csv"),row.names=F)
