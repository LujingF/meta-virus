combine_file_2_df<-function(path){
  a=list.files(path)
  dir=paste(path,"/",a,sep="")
  n=length(dir)
  filename<-strsplit(dir[1],"/")[[1]][length(strsplit(dir[1],"/")[[1]])]
  examplename<-strsplit(filename,"_")[[1]][1]
  merge.data=read.table(file=dir[1],header=F,sep="\t")
  #merge.data<-subset(merge.data,select=c("V8","V9"))
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
#total_df<-combine_file_2_df("counts_txt/")
total_df<-combine_file_2_df("../speciescounts/")
#total_df<-combine_file_2_df("Abundance/species/")
total_df[is.na(total_df)]<-0
tt<-total_df

#this is for top 10 or top 30 annotation species bar plot
rownames(total_df)<-total_df$species
total_df<-total_df[,-1]
Total_df<-total_df
for(i in 1:ncol(total_df)){
  #  Total_df[,i]<-total_df[,i]/sum(total_df[,i])*1000000
  Total_df[,i]<-total_df[,i]/sum(total_df[,i])
}
Total_df<-Total_df[order(Total_df$G5A,decreasing = T),]

top30sp<-head(rownames(Total_df),31)
species_annot<-head(Total_df,11)
for(i in 1:ncol(species_annot)){
  species_annot[11,i]<-1-sum(Total_df[1:10,i])
}
rownames(species_annot)<-c(rownames(species_annot)[1:10],"others")
#species_annot<-head(Total_df[order(Total_df$G5A,decreasing = T),],30)
species_annot$species<-rownames(species_annot)
library(reshape2)
library(ggplot2)
sp_annot<-melt(species_annot,id="species")
colnames(sp_annot)<-c("species","sample","percent")
mycol<-c('#00AED7','#C1E168','#FD9347','#319F8C', "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32","#FF00FF","#FFCC00")
shape=c(19, 17 ,15 ,18, 8, 0 ,1 ,2, 5)
#pdf("virus_sp_annot.pdf",width=24,height = 16)
pdf("virus_sp_annot.pdf")
#png("virus_sp_annot.png")
ggplot(sp_annot,aes(x=sample,y=percent,fill=species))+geom_bar(stat='identity')+scale_fill_manual(values=mycol)+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
dev.off()
###this is the end of top 10 species annotation result

###top 30 species annotation heatmap
top30sp<-rownames(Total_df)
top30sp<-top30sp[-which(top30sp=="S__undef")]
top30sp<-head(top30sp,30)
library(pheatmap)
heatmap_dat<-as.matrix(Total_df[rownames(Total_df) %in% top30sp,])
svg("top30sp_heatmap.svg")
pheatmap(heatmap_dat,scale = 'row')
dev.off()
### this is the end of top 30 species annotation heatmap

###hclust plot
library(vegan)
#distance.ex<-vegdist(t(Total_df),method='bray',na.rm=T)
distance.ex<-vegdist(t(total_df),method = 'bray',na.rm = T)
svg("virus_group_hust.svg")
#png('virus_group_hust.png')
hclust.ex<-hclust(distance.ex,method="ward.D2")
plot(hclust.ex,hang=-1)
dev.off()
###this is end of hclust

###pca analysis
library(ggrepel)
library(stringr)


da<-total_df
samplelist<-as.factor(colnames(da))
#da<-da*1000000
#da<-round(da)
data <- as.data.frame(t(da))

set.seed(1000)
newdata <- rrarefy(data, min(rowSums(data)))
data <- decostand(newdata, "hell")
tm=read.table("treatment.txt",header =T,sep="\t")
grouptable<-tm
colnames(grouptable)<-c("Sample","treatment")
MetaTable<-tt
reGroup1=c()
reGroup2=c()
for(i in 1:length(colnames(MetaTable))){
  if(is.null(tm$group)){
    reGroup1=c(reGroup1, as.character(tm[which(tm$Sample==colnames(MetaTable)[i]),]$treatment))
    #print(tm[which(tm$Sample==colnames(MetaTable)[i]),]$treatment)
    reGroup2=c(reGroup2, as.character(tm[which(tm$Sample==colnames(MetaTable)[i]),]$treatment))
  }else{
    reGroup1=c(reGroup1, as.character(tm[which(tm$Sample==colnames(MetaTable)[i]),]$treatment))
    reGroup2=c(reGroup2, as.character(tm[which(tm$Sample==colnames(MetaTable)[i]),]$site))
  }
}

groups=data.frame(Sample=rownames(data),g1=tm$treatment, g2=tm$treatment)

sampleNumber=length(tm$treatment)
groupNumber=length(levels(tm$treatment))

MetaTableTD=data
### PCA ###
svg("PCA.svg",width=8,height = 6)
par(mar=c(4,4,1,1))
pca1=rda(MetaTableTD)
pc1=c(pca1$CA$eig/sum(pca1$CA$eig))[1]*100
pc2=c(pca1$CA$eig/sum(pca1$CA$eig))[2]*100

plot(pca1,display="si",scaling=1,type="n", xlab=str_c('PC1 (',round(pc1,1),'%)'),ylab=str_c('PC2 (',round(pc2,1),'%)'),main="")
points(pca1, dis="si", scaling=1,col=mycol[groups$g1],pch=shape[groups$g2],cex=1)
legend("bottomright", legend=levels(groups$g1), col=mycol,pch=shape,bty="n",cex=0.8)

for(i in 1:sampleNumber){
  ordispider(pca1,groups$g1,dis="si",scaling=1, show.groups=levels(groups$g1)[i],col=mycol)
}

dev.off()
### NMDS ###
svg("NMDS.svg",width=8,height = 6)
par(mar=c(4,4,1,1))
nmds1=metaMDS(MetaTableTD,distance = "bray", k = 2,trymax=20)
plot(nmds1, display="si", type="n",choices = c(1, 2))
points(nmds1, dis="si", col=mycol[groups$g1],pch=shape[groups$g2],cex=1)
legend("bottomright", legend=levels(groups$g1), col=mycol,pch=shape,bty="n",cex=0.8)

for(i in 1:sampleNumber){
  #ordiellipse(nmds1,tm$TM3,dis="si",scaling=1, show.groups=levels(tm$TM3)[i],col=brewer.pal(6,"Set2")[i])
  ordispider(nmds1,groups$g1,dis="si",scaling=1, show.groups=levels(groups$g1)[i],col=mycol)
}

dev.off()

### PCoA ###
svg("PCoA.svg", width=8, height=6)
MetaTableT<-t(total_df)
MetaTableTB=vegdist(MetaTableT, method="bray")
MetaTableTBpcoA=cmdscale(MetaTableTB)
par(mar=c(4,4,1,1))
plot(MetaTableTBpcoA,xlab = "PCoA1",ylab="PCoA2")
#text(MetaTableTBpcoA, labels = row.names(MetaTableTBpcoA))
points(MetaTableTBpcoA, col=mycol[groups$g1],pch=shape[groups$g2],cex=1)
legend("bottomright", legend=levels(groups$g1), col=mycol,pch=shape,bty="n",cex=0.8)
for(i in 1:sampleNumber){
  ordispider(MetaTableTBpcoA,groups$g1,dis="si",scaling=1, show.groups=levels(groups$g1)[i],col=mycol)
}
dev.off()

###pca plot
dif_sp<-read.table("lefse_diff_species.result",header = F)
dif_sp_matrix<-as.matrix(total_df[(rownames(total_df) %in% dif_sp$V1),])
da<-dif_sp_matrix
samplelist<-as.factor(colnames(da))
#da<-da*1000000
#da<-round(da)
data <- as.data.frame(t(da))
set.seed(1000)
newdata <- rrarefy(data, min(rowSums(data)))
data <- decostand(newdata, "hell")
decorana(data)

pca = rda(data)

pc1=c(pca$CA$eig/sum(pca$CA$eig))[1]*100
pc2=c(pca$CA$eig/sum(pca$CA$eig))[2]*100
svg("diff_pca.svg")
par(mar=c(4,4,1,1))
plot(pca,display="si",scaling=1,type="n",xlim=c(-0.3, 0.8), xlab=str_c('PC1 (',round(pc1,1),'%)'),ylab=str_c('PC2 (',round(pc2,1),'%)'),main="")
#axis(1,at=seq(-0.5,0.5,0.25))  
points(pca, dis="si", scaling=1,col=mycol[groups$g1],pch=shape[groups$g2],cex=1)
legend("bottomright", legend=levels(groups$g1), col=mycol,pch=shape,bty="n",cex=0.8)

for(i in 1:sampleNumber){
  ordispider(pca,groups$g1,dis="si",scaling=1, show.groups=levels(groups$g1)[i],col=mycol)
}
dev.off()

x<-data.frame(pca$CA$u[, 1:2])
plot_dat<-cbind(grouptable,x)
tmp <- pca$CA$eig/sum(pca$CA$eig)
xlab_text <- paste("PC1 (", round(tmp[1]*100,2), "%)",sep="")
ylab_text <- paste("PC2 (", round(tmp[2]*100,2), "%)",sep="")
p <- ggplot(data=plot_dat, aes(PC1, PC2, treatment)) + geom_point(aes(color=treatment),size=2) +
  #geom_text_repel(aes(label = rownames(x))) +
  scale_color_manual(values=mycol)+
  theme_bw() +
  geom_vline(xintercept = 0,linetype='dashed',size=1)+
  geom_hline(yintercept = 0,linetype='dashed',size=1)+
  xlab(xlab_text)+
  ylab(ylab_text)+
  labs(title="PCA - P1 vs P2") +
  theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 12),legend.title=element_blank(),plot.title = element_text(face="bold",lineheight=25,hjust=0.5))
#svg("diff_pca.svg")
png("diff_pca.png")
p
dev.off()
svg("dif_sp_heatmap.svg")
pheatmap(dif_sp_matrix,scale = 'row')
dev.off()