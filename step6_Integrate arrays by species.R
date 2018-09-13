path=path1
path2=paste0("C:/Users/",user,"/Dropbox/Crisp2018/integrated arrays by species/")
#path2=paste0("C:/Users/",user,"/Dropbox/Crisp2018/relooking_at_old_runs/Crisp2017/within_species/integrated arrays by species/")
#path2=paste0("C:/Users/",user,"/Dropbox/Crisp2018/relooking_at_old_runs/isra4/vol transformed with med/integrated arrays by species/")		
set="Med" #note - "arrays" must match the "set"
arrays=c("A","B","F","G","H","I")
#arrays=c("C","D") #,"E")
as=list()
for (i in 1:length(arrays)){
  x=read.csv(paste0(path,arrays[i],"/",arrays[i],"_no_multiple_targets.csv"),stringsAsFactors = F)
  print(colnames(x))
  colnames(x)[2:ind]=c("Rep1","Rep3","Rep2") #change accroding to sample number
  x[ ,1]=paste0(arrays[i],"_",x[ ,1])
  as[[i]]=x
  
}
res=do.call("rbind",as)
res$array=substr(res$OTUId,1,1)

z=summary(as.factor(res$genome))
jpeg(paste0(path2,set,"_Num.Uniques.by.target.jpeg"),width=800,height=600)
barplot(z,ylab="Num. Unique Spacers",main=set)
dev.off()
z=as.data.frame(z)
write.csv(z,paste0(path2,set,"_Num.uniques.by.target.csv"))

z=summary(as.factor(res$array))
jpeg(paste0(path2,set," _Num.Uniques.by.array.jpeg"),width=800,height=600)
barplot(z,ylab="Num Unique Spacers",main=set)
dev.off()
z=as.data.frame(z)
write.csv(z,paste0(path2,set,"__Num.Uniques.by.array.csv"))

write.csv(res,paste0(path2,set,"_all_arrays_no_multiple_targets.csv"),row.names=F)
