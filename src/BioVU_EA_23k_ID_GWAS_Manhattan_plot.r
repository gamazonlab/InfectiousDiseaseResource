#load gwas results
df<-read.table('/data/g***/z***/hihost_gwas/manh_plot/pooled/sig.txt',header = T,stringsAsFactors = F)

#only keeps pheno with genome-wide sig
gwsig<-aggregate(df$P,list(df$phecode),function(x) min(x,na.rm = T))
gwsig<-gwsig[gwsig[,2]<5e-8,]
df<-df[which(df$phecode %in% gwsig[,1]),]
df<-df[!is.na(df$P),]
df$phecode=as.numeric(sub('X','',df$phecode))
df$logp=log(df$P,10)*-1

#phecode annotation
anno<-read.csv('/data/c***/z***/anno/phecode/phecode_icd9_rolled.csv',header = T,stringsAsFactors = F)
anno$phecode=anno$PheCode
anno<-anno[!duplicated(anno$PheCode),]
df<-merge(anno[,c('phecode','Phenotype')],df,by='phecode')
df[which(df$Phenotype=='Human immunodeficiency virus [HIV] disease'),'Phenotype']<-'HIV disease'


#---generate_new_pos---
df$pos<-df$BP
df$chr=df$CHR
chr_box<-as.data.frame(matrix(data=NA,nrow=22,ncol=3))
chr_df<-df[df$chr==1,]
chr_box[1,1]<-max(chr_df$pos)
chr_box[1,2]<-0
chr_box[1,3]<-chr_box[1,1]/2

for (i in 2:22){
  chr_df<-df[df$chr==i,]
  chr_box[i,1]<-max(chr_df$pos)
  chr_box[i,2]<-max(chr_df$pos)+chr_box[i-1,2]
  df[which(df$chr==i),'pos']<-df[which(df$chr==i),'pos']+chr_box[i,2]
  chr_box[i,3]<-chr_box[i,2]+chr_box[i,1]/2
}


#---------------------
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

#color
color_panel<-c(brewer.pal(9, "Set1")[1:9],brewer.pal(8, "Dark2")[c(1,3,4,6)])
phecode_char<-as.character(unique(df$Phenotype))
df$Phenotype<-as.factor(df$Phenotype)

#make a copy #for testing
df_copy=df

#---------
df<-df[df$P<1e-4,]

#label
df$label<-ifelse((df$P %in% gwsig[,2]),df$SNP,NA)
df_labeled<-df[!is.na(df$label),]
df_labeled$label<-ifelse(duplicated(df_labeled$phecode),NA,df_labeled$label)
df_unlabeled<-df[is.na(df$label),]
df_unlabeled<-df_unlabeled[sample(seq(1,nrow(df_unlabeled)),nrow(df_unlabeled),replace = F),]
df<-rbind(df_unlabeled,df_labeled)


#write.table(df,'/data/g***/z***/hihost_gwas/manh_plot/pooled/joint_13traits_1e-4.txt',quote = F,sep='\t',row.names = F)

#plot
png(paste0('/data/g***/z***/hihost_gwas/manh_plot/pooled/joint_13traits_1e-4.png'),width = 2300,height = 1600,res=300)

set.seed(1)
ggplot(df, aes(x=pos, y=logp,label=label)) + 
  scale_x_continuous(breaks=chr_box$V3, labels = c(as.character(seq(1,15)),' ','17','','19','','21','')) +  #x axis
  scale_y_continuous(breaks=c(1,2,3,5,10,20,30),trans='log10') +  #x axis
  geom_point(data = df, shape=21, color = 'black', size=1.75,aes(x = pos, y = logp,fill=Phenotype)) + #points
  scale_fill_manual(breaks=phecode_char,
                    values=color_panel)+
  labs(x = "Chromosome", y = "-log(P)", title = "") +
  theme_bw() +  #rm background
  theme(panel.grid =element_blank()) +  #rm grids
  theme(legend.position= "top",legend.text=element_text(size=9)) + 
  #theme(legend.position= "right",legend.title = element_blank()) + #set legend
  labs(fill = "ID phenotypes")+ #set legend title
  geom_hline(yintercept=log(5e-8,10)*-1, linetype="dashed", color = "black")+ 
  geom_vline(xintercept=chr_box$V2[-1], linetype="dashed", color = "ivory2")+ 
  ggtitle('') +
  guides(fill=guide_legend(nrow=5,byrow=TRUE,keyheight=0.75,title=''))+
  geom_label_repel( #non overlapped labels
    size=2.5, 
    fill=rgb(255, 255, 255, 200, maxColorValue=255),
    direction='y', #only work on Y-axis
    nudge_x=2e8, #shift to the right
    segment.alpha = 0.5,  #transparent of segment
    min.segment.length = 0.5,
    fontface='italic'
  )


dev.off()






