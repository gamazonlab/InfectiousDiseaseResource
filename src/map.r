chr = as.numeric(commandArgs(TRUE))

library(data.table)

df<-data.frame(fread('#Absolute path to trait gwas file'))
df$chr=sapply(df$variant,function(x) strsplit(x,':')[[1]][1])

df$pos1=df$pos2=sapply(df$variant,function(x) strsplit(x,":")[[1]][2])

df$a1=sapply(df$variant,function(x) strsplit(x,":")[[1]][3])
df$a2=sapply(df$variant,function(x) strsplit(x,":")[[1]][4])

df<-df[,c('chr','pos1','pos2','a1','a2')]
write.table(df,paste0('#Absolute path for generated annotation file/chr',chr,'.txt'),quote = F,sep='\t',row.names = F)


cmd=paste0('#Absolute path to annovar script/annotate_variation.pl -filter -out #Absolute path to annotation file/chr',chr,' -build hg19 -dbtype snp138 #Directory for output/chr',chr,'.txt #Directory for where annotation file is stored/annovar/humandb/')
system(cmd,wait = T)



