#ID GWAS IN BIOVU 23K EA

#both the genotype and phenotype data are processed by Xue.
#genotype data dir: /data/g***/z***/biovu/23k/geno/ea_chr/chr*
#binary trait dir: /data/c***/z***/data/biovu/pheno/EA23k_plink.txt
#adj cov: AGE_in_days,GENDER,PC1,PC2,PC3,PC4,PC5,ARRAY

args<-as.numeric(commandArgs(TRUE)) #528

pheno_info<-read.table('/data/g***/z***/hihost_gwas/info/N_info.txt',header = T,stringsAsFactors = F,sep='\t')
pheno_info<-pheno_info[pheno_info$N_cases>=100,] #24

#write.table(pheno_info,'/data/g***/z***/hihost_gwas/info/BioVU_23kEA_ID_GWAS_trait.txt',quote = F,sep = '\t',row.names = F)

phecode=paste0('X',pheno_info$PheCode) #'X117'

run_i=1
run_list<-list()
for (i in 1:length(phecode)){
  for (j in 1:22){
    run_list[[run_i]]<-c(i,j)
    run_i=run_i+1
  }
}

run_list<-run_list[[args]]
phecode<-phecode[run_list[1]]
chr<-run_list[2]


biovu<-read.table('/data/c***/z***/data/biovu/pheno/EA23k_plink.txt',header = T,stringsAsFactors = F)


if(phecode %in% colnames(biovu)){
  N_cases<-length(which(biovu[,phecode]==2))
  N_controls<-length(which(biovu[,phecode]==1))
  
  cmd=paste0('plink2 --bfile /data/g***/z***/biovu/23k/geno/ea_chr/chr',chr,'  --allow-no-sex --logistic hide-covar --pheno /data/c***/z***/data/biovu/pheno/EA23k_plink.txt --pheno-name ',phecode,' --geno 0.95 --hwe 0.0001 --covar /data/c***/z***/data/biovu/pheno/EA23k_plink_array.txt --covar-name AGE_in_days,GENDER,PC1,PC2,PC3,PC4,PC5,ARRAY --covar-variance-standardize --out /data/g***/z***/hihost_gwas/raw_asso/',phecode,'_chr',chr)
  
  system(cmd,wait = T)
  
  #asso_path=paste0('/data/c***/z***/projects/biovu/gwas/raw/',phecode,'_chr',chr,'.',phecode,'.glm.logistic')
  
  #asso<-read.table(asso_path,header = T,stringsAsFactors = F,comment.char = '&')
  #asso<-asso[asso$TEST=='ADD',]
  
  #write.table(asso,paste0('/data/c***/z***/projects/biovu/gwas/raw/',phecode,'_chr',chr),quote = F,sep = '\t',row.names = F)
  
  
  
}











