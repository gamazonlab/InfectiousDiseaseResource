#MR id -> complications
#ml GCC/8.2.0  OpenMPI/3.1.4 Intel/2019.1.144  IntelMPI/2018.4.274 R/3.6.0

library('MendelianRandomization')
library('ggplot2')

anno<-read.table('/data/c***/z***/anno/gencode/37/gencode.v32.GRCh37.txt',header = T,stringsAsFactors = F)

info<-read.table('/data/g***/z***/hihost_gwas/mr_id_complications/top.txt',header = T,stringsAsFactors = F,sep='\t')

info<-merge(info,anno[,c('genename','chr','left','right')],by=1)
info$chr=sub('chr','',info$chr)



p_cutoff=1e-5 

#top ranked genes  
info<-read.table('/data/g***/z***/hihost_gwas/mr_id_complications/top.txt',header = T,stringsAsFactors = F,sep='\t')
info<-merge(info,anno[,c('genename','chr','left','right')],by=1)
info$chr=sub('chr','',info$chr)

for(i in 1:nrow(info)){
  exp_phecode=paste0('X',info[i,'exp_phecode'])
  outcome_phecode=paste0('X',info[i,'outcome_phecode'])
  
  #plink clumping
  cmd=paste0('plink --bfile /data/g***/z***/biovu/23k/geno/ea --clump /data/g***/z***/hihost_gwas/asso/',exp_phecode,'.txt.gz --clump-field P --clump-p1 ',p_cutoff,' --clump-r2 0.01 --out /data/g***/z***/hihost_gwas/mr_id_complications/tmp/',exp_phecode,'_',p_cutoff)
  #system(cmd,wait = T)
  
  #load clumped snps
  clumped<-read.table(paste0('/data/g***/z***/hihost_gwas/mr_id_complications/tmp/',exp_phecode,'_',p_cutoff,'.clumped'),header = T,stringsAsFactors = F)
  
  #load exp gwas
  exp<-read.table(paste0('/data/g***/z***/hihost_gwas/asso/',exp_phecode,'.txt.gz'),header = T,stringsAsFactors = F)
  exp<-exp[(exp$SNP %in% clumped$SNP),]
  exp$beta_exp=log(exp$OR,2.718)
  exp$se_exp=exp$beta_exp/exp$Z_STAT
  exp<-exp[,c('SNP','A1','beta_exp','se_exp')]
  colnames(exp)<-c('SNP','eff_exp','beta_exp','se_exp')
  
  #load outcome gwas
  outcome<-read.table(paste0('/data/g***/z***/hihost_gwas/asso/',outcome_phecode,'.txt.gz'),header = T,stringsAsFactors = F)
  outcome$beta_outcome=log(outcome$OR,2.718)
  outcome$se_outcome=outcome$beta_outcome/outcome$Z_STAT
  outcome<-outcome[,c('SNP','A1','beta_outcome','se_outcome')]
  
  df<-merge(exp,outcome,by='SNP')
  
  #rm na
  df<-df[!is.na(df$beta_exp+df$beta_outcome),]
  
  #harmonize
  df$beta_outcome=ifelse(df$eff_exp==df$A1,df$beta_outcome,df$beta_outcome*-1)
  
  #flip the allele
  df$beta_outcome=ifelse(df$beta_exp>0,df$beta_outcome,df$beta_outcome*-1)
  df$beta_exp=abs(df$beta_exp)
  
  
  info[i,paste0('n_SNPs_allgenes_p_cutoff_',p_cutoff)]<-nrow(df)
  if(nrow(df)<3){next}
  
  #mr weighted median
  ans_median_weighted<-mr_median(mr_input(bx = df$beta_exp, bxse = df$se_exp, by = df$beta_outcome, byse = df$se_outcome),weighting = "weighted", iterations = 100)
  
  info[i,paste0('median_weighted_allgenes_p_cutoff_',p_cutoff,'_beta')]<-ans_median_weighted@Estimate
  info[i,paste0('median_weighted_allgenes_p_cutoff_',p_cutoff,'_pvalue')]<-ans_median_weighted@Pvalue
  
  #mr egger
  ans_egger<-try(mr_egger(mr_input(bx = df$beta_exp, bxse = df$se_exp, by = df$beta_outcome, byse = df$se_outcome),alpha=0.05))
  if('try-error' %in% class(ans_egger)){next}  
  info[i,paste0('egger_allgenes_p_cutoff_',p_cutoff,'_beta')]<-ans_egger@Estimate
  info[i,paste0('egger_allgenes_p_cutoff_',p_cutoff,'_pvalue')]<-ans_egger@Pvalue.Est
  
  
  #generate scatter plot
  
  png(paste0('/data/g***/z***/hihost_gwas/mr_id_complications/plot/expo_',info[i,'exp_phecode'],'_outcome_',info[i,'outcome_phecode'],'.png'),width = 1500,height = 1500,res=400)
  
  scatter_p<-ggplot(df, aes(x=beta_exp, y=beta_outcome)) +
    geom_errorbar(aes(ymin=beta_outcome-se_outcome*1.96, ymax=beta_outcome+se_outcome*1.96),size=0.05,color='ivory3') +
    geom_errorbarh(aes(xmin=beta_exp-se_exp*1.96, xmax=beta_exp+se_exp*1.96),size=0.05,color='ivory3') +
    geom_point(size=1.5,color='dodgerblue4') +
    xlab(paste0("Effect size of SNP-Exposure \n",info[i,'exp_trait'])) +
    ylab(paste0("Effect size of SNP-Outcome \n",info[i,'outcome_trait']))+
    theme_classic()
  #xlim(-0.5,3)
  print(scatter_p)
  dev.off()
  
}


write.table(info,paste0('/data/g***/z***/hihost_gwas/mr_id_complications/result/allgenes_',p_cutoff,'.txt'),quote = F,sep='\t',row.names = F)

