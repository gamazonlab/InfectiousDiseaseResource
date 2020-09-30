#First run map.r script to generate annotation file for each gwas trait

library(data.table)

df<-data.frame(fread('#Absolute path to gwas trait file'))
df$chr=sapply(df$variant,function(x) strsplit(x,":")[[1]][1])

mapped_tmp<-read.table(paste0('#Absolute path to annotation file/chr.hg19_snp138_dropped'),stringsAsFactors = F)

mapped_tmp$variant_1<-paste0(mapped_tmp$V3,':',mapped_tmp$V4,':',mapped_tmp$V6,':',mapped_tmp$V7)
mapped_tmp$variant_2<-paste0(mapped_tmp$V3,':',mapped_tmp$V4,':',mapped_tmp$V7,':',mapped_tmp$V6)

map_1<-merge(df,mapped_tmp,by.x='variant',by.y='variant_1')
map_2<-merge(df,mapped_tmp,by.x='variant',by.y='variant_2')

mapped_tmp<-rbind(map_1[,-ncol(map_1)],map_2[,-ncol(map_2)])
mapped_tmp<-mapped_tmp[,c('variant','V2','minor_allele','V6','V7','beta','se')]


mapped_tmp$major_allele<-ifelse(mapped_tmp$minor_allele==mapped_tmp$V6,mapped_tmp$V7,mapped_tmp$V6)
colnames(mapped_tmp)[2]<-'rsid'

write.table(mapped_tmp[,c('variant','rsid','minor_allele','major_allele','beta','se')],'#Absolute path to output directory for annotated trait gwas file/trait_UKBB_rsid.txt',quote = F,sep='\t',row.names = F)








