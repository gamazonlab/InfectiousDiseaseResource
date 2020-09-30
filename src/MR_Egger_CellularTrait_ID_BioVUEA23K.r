#MR cellular -> id
#ml GCC/8.2.0  OpenMPI/3.1.4 Intel/2019.1.144  IntelMPI/2018.4.274 R/3.6.0

args=as.numeric(commandArgs(TRUE))  #148

p_cutoff=1e-5

library('MendelianRandomization')
library('ggplot2')

#load celluar traits
file_list<-dir('/data/g***/z***/hihost_gwas/summary/')
trait_list<-sub('_Meannew7.qfam.parents.perm_rsID.txt.gz','',file_list[grep('qfam.parents.perm_rsID.txt.gz',file_list)])
trait=trait_list[args]

#mkdir
dir.create(paste0('/data/g***/z***/hihost_gwas/mr_cellular_id/tmp/',trait,'/'))
dir.create(paste0('/data/g***/z***/hihost_gwas/mr_cellular_id/tmp/',trait,'/clumping/'))
dir.create(paste0('/data/g***/z***/hihost_gwas/mr_cellular_id/tmp/',trait,'/scatter_plot_',p_cutoff,'/'))

#plink clumping
cmd=paste0('plink --bfile /data/g***/z***/biovu/23k/geno/ea --clump /data/g***/z***/hihost_gwas/summary/',trait,'_Meannew7.qfam.parents.perm_rsID.txt.gz --clump-field EMP1 --clump-p1 ',p_cutoff,' --clump-r2 0.01 --out /data/g***/z***/hihost_gwas/mr_cellular_id/tmp/',trait,'/clumping/',p_cutoff)
system(cmd,wait = T)

#load clumped snps
clumped<-read.table(paste0('/data/g***/z***/hihost_gwas/mr_cellular_id/tmp/',trait,'/clumping/',p_cutoff,'.clumped'),header = T,stringsAsFactors = F)

#load exp gwas
exp<-read.table(paste0('/data/g***/z***/hihost_gwas/summary/',trait,'_Meannew7.qfam.parents.perm_rsID.txt.gz'),header = T,stringsAsFactors = F)
exp<-exp[(exp$SNP %in% clumped$SNP),]
exp<-exp[,c('SNP','MinorA','MajorA','Beta','EMP_se')]
colnames(exp)<-c('SNP','eff_exp','ref_exp','beta_exp','se_exp')

#rm palindromic vaiants
exp<-exp[-which(paste0(exp$eff_exp,exp$ref_exp) %in% c('AT','TA','CG','GC')),]

#load outcome gwas
id_list<-sub('.txt.gz','',dir('/data/g***/z***/hihost_gwas/asso/'))

output<-data.frame(trait=id_list)

for (i in 1:length(id_list)){
  print(i)
  id=id_list[i]
  outcome<-read.table(paste0('/data/g***/z***/hihost_gwas/asso/',id,'.txt.gz'),header = T,stringsAsFactors = F)
  outcome$beta_outcome=log(outcome$OR,2.718)
  outcome$se_outcome=outcome$beta_outcome/outcome$Z_STAT
  outcome<-outcome[,c('SNP','A1','beta_outcome','se_outcome')]
  
  df<-merge(exp,outcome,by='SNP')
  
  #rm na
  df<-df[!is.na(df$beta_exp+df$beta_outcome),]
  
  #harmonize
  df$beta_outcome=ifelse(df$eff_exp==df$A1,df$beta_outcome,df$beta_outcome*-1)
  
  output[i,'n_of_SNPs_as_IV']<-nrow(df)
  
  if(nrow(df)<=3){next}
  
  #mr weighted median
  ans_median_weighted<-mr_median(mr_input(bx = df$beta_exp, bxse = df$se_exp, by = df$beta_outcome, byse = df$se_outcome),weighting = "weighted", iterations = 100)
  
  output[i,'median_weighted_beta']<-ans_median_weighted@Estimate
  output[i,'median_weighted_p']<-ans_median_weighted@Pvalue
  
  #mr egger
  ans_egger<-try(mr_egger(mr_input(bx = df$beta_exp, bxse = df$se_exp, by = df$beta_outcome, byse = df$se_outcome),alpha=0.05))
  if('try-error' %in% class(ans_egger)){next}  
  output[i,'egger_beta']<-ans_egger@Estimate
  output[i,'egger_p']<-ans_egger@Pvalue.Est
  
  # #plot
  # png(paste0('/data/g***/z***/hihost_gwas/mr_cellular_id/tmp/',trait,'/scatter_plot_',p_cutoff,'/',id,'.png'),width = 1000,height = 1000,res=150)
  # 
  # scatter_p<-ggplot(df, aes(x=beta_exp, y=beta_outcome)) +
  #   geom_point() +
  #   geom_errorbar(aes(ymin=beta_outcome-se_outcome*1.96, ymax=beta_outcome+se_outcome*1.96)) +
  #   geom_errorbarh(aes(xmin=beta_exp-se_exp*1.96, xmax=beta_exp+se_exp*1.96)) +
  #   xlab("Effect size of SNP-Exposure (Celluar trait)") +
  #   ylab("Effect size of SNP-Outcome (ID trait)")
  # print(scatter_p)
  # dev.off()
  # 
  
}

write.table(output,paste0('/data/g***/z***/hihost_gwas/mr_cellular_id/result/raw_r2_0.01/',p_cutoff,'_',trait,'.txt'),quote = F,row.names = F,sep='\t')



