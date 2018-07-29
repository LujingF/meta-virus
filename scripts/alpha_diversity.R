args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2){
	cat("##########\n")
	cat("Function: virus alph diversity analysis, first input file is a dir which contain toal species reads count result\n")
	cat("format species reads count result like this: first column is the species name(annotation information from kingdom to species, and seperated by ';', the second column is reads count)\n")
	cat("second input file is sample and group information, header is needed, first column is sample name, second column is group name\n")
	quit()
}
library(vegan)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggsignif)
library(stringr)
combine_file_2_df<-function(path){
  a=list.files(path)
  dir=paste(path,"/",a,sep="")
  n=length(dir)
  filename<-strsplit(dir[1],"/")[[1]][length(strsplit(dir[1],"/")[[1]])]
  examplename<-strsplit(filename,"_")[[1]][1]
  merge.data=read.table(file=dir[1],header=F,sep="\t")
  colnames(merge.data)<-c("species",examplename)
  for(i in 2:n){
   new.data<-read.table(file=dir[i],header=F,sep="\t")
    filename<-strsplit(dir[i],"/")[[1]][length(strsplit(dir[i],"/")[[1]])]
    examplename<-strsplit(filename,"_")[[1]][1]
    colnames(new.data)<-c("species",examplename)
    #merge.data<-cbind(merge.data,new.data)
    merge.data<-merge(merge.data,new.data,by="species",all.x=T,all.y=T)
  }
  return(merge.data)
}  

#total_df<-combine_file_2_df("counts_alph/")
#total_df<-combine_file_2_df("2_1_Annot_Percent/")
total_df<-combine_file_2_df(args[1])
total_df[is.na(total_df)]<-0
tt<-total_df
rownames(tt)<-tt[,1]
tt<-tt[,-1]
#grouptable<-read.table("sample_group.txt",header=T,sep="\t",check.names=F, comment.char="")
grouptable<-read.table(args[1],header=T,sep="\t",check.names=F, comment.char="")
colnames(grouptable)<-c("Sample","treatment")
data<-tt
group<-grouptable
paired=brewer.pal(n = 12, name = "Paired")
mycol=c("#00AED7","#FD9347","#C1E168","#319F8C","#FD6260","#FAD390",paired,colors()[c(47,83,121,562,451,97,515,518,26,490,59,70,222,13,27,33,65,104,112,258,294,367,381,428,611,634,656,178)])

reGroup=c()
for(i in 1:length(colnames(data))){
  reGroup=c(reGroup, as.character(group[which(group$Sample==colnames(data)[i]),]$treatment))
}

dataT=t(data)
set.seed(1000)

dataT2=rrarefy(dataT,min(rowSums(dataT)))
print(min(rowSums(dataT)))

shannon.wiener=diversity(dataT2, "shannon")
shannon=data.frame(values=shannon.wiener)
simpson=data.frame(values=diversity(dataT2,"simpson"))
insimpson=data.frame(values=diversity(dataT2,"inv"))
s=specnumber(dataT2)
sn=data.frame(values=s)
pielou=data.frame(values=shannon.wiener/log(s))
ca=data.frame(estimateR(dataT2)) #chao1, ACE
cat=data.frame(t(ca))
rn=rownames(shannon)
dv=data.frame(Observed=sn$values,Chao1=cat$S.chao1,
              Shannon=shannon$values,Pielou=pielou$values,
              Simpson=simpson$values,Insimpson=insimpson$values,
              g=reGroup)
dvM=melt(dv,id="g",shannon="shannon")
svg("alph_diversity.svg",width=15,height = 15)
ggplot(dvM,aes(variable,value))+stat_boxplot(aes(fill=g),geom ='errorbar', width=0.3, position=position_dodge(1.05))+
  geom_boxplot(aes(fill=g), position=position_dodge(1.05))+
  facet_wrap(~ variable, scales="free")+scale_fill_manual(values=mycol)+
  theme_calc()+labs(fill="Group")+
  theme(legend.text=element_text(size = 20),
        legend.title=element_text(size = 20),
        strip.text=element_text(size=24))
dev.off()


rearrangeGroup=data.frame(Sample=colnames(data),treatment=reGroup)
tdata=t(data)
min(rowSums(tdata))
out= rarecurve(tdata, step=20000, sample=min(rowSums(tdata)), xlab="Sample size", ylab="Species",label=F)

Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
Smax <- sapply(out, max)

svg("rare.svg",width=10,height = 8)
par(mar=c(4,4,1,8))
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Sample Size", ylab = "Species", type = "n")
groupCol=mycol[rearrangeGroup$treatment]
for (i in seq_along(out)) {
  N <- attr(out[[i]], "Subsample")
  lines(N, out[[i]],col=groupCol[i],lwd=1)
}
abline(v = min(rowSums(tdata)),lty="dashed",col=grey(0.4),lwd=1)
legend(x=xinch(7.5),y=yinch(5) ,xpd=T,inset=0.2,lty=1,col=groupCol ,horiz=FALSE,legend=c(rownames(tdata)),bty="n")
dev.off()

svg("Accum.svg",width=15,height = 12)
accum <- specaccum(tdata, method="random", permutations=100)
accumDF=data.frame(summary(accum), check.names = FALSE, stringsAsFactors = FALSE)
accumDF=accumDF[,-1]
freq=c()
for(i in 1: length(rownames(accumDF))){
  a=as.integer(str_extract_all(as.character(accumDF$Freq[i]), "[0-9]+[0-9]"))
  freq=c(freq,a)
  
}
accumDF$Freq=freq
names(accumDF)=c("Sites","Num")
#lapply(accum,summary)
ggplot(accumDF, aes(Sites,Num))+
  stat_boxplot(geom ='errorbar', width=0.2, position=position_dodge(1.05))+
  geom_boxplot(fill="#00AED7",width=0.5)+theme_calc()+labs(x="Accumulation", y="")+
  theme(axis.text.x = element_blank(),axis.title.x = element_text(size=16))
dev.off()
