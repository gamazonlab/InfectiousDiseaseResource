#ml GCC/5.4.0-2.26  OpenMPI/1.10.3 pandas/0.18.1-Python-2.7.12 numpy/1.11.1-Python-2.7.12 scipy/0.17.0-Python-2.7.12 R

args=as.numeric(commandArgs(TRUE))   #148 traits

#trait names
file_list<-dir('/data/g***/z***/hihost_gwas/summary/')
trait_list<-sub('_Meannew7.qfam.parents.perm_rsID.txt.gz','',file_list[grep('qfam.parents.perm_rsID.txt.gz',file_list)]) #
trait=trait_list[args]


#v6p
tissue_list=read.table('/data/g***/z***/data/gtex/v6p_haky/v6p.tissue',header = T,stringsAsFactors = F)
tissue_list<-tissue_list[,1]

for (i in 1:length(tissue_list)){
  tissue=tissue_list[i]
  
  cmd=paste0('python2.7 /data/c***/z***/projects/cross_tissue/metaxcan/MetaXcan/software/MetaXcan.py   --model_db_path /data/g***/z***/data/gtex/v6p_haky/GTEx-V6p-HapMap-2016-09-08/TW_',tissue,'_0.5.db   --covariance /data/g***/z***/data/gtex/v6p_haky/GTEx-V6p-HapMap-2016-09-08/TW_',tissue,'.txt.gz   --gwas_file /data/g***/z***/hihost_gwas/summary/',trait,'_Meannew7.qfam.parents.perm_rsID.txt.gz   --snp_column SNP   --effect_allele_column MinorA   --non_effect_allele_column MajorA   --beta_column Beta   --se_column EMP_se   --output_file /data/g***/z***/hihost_gwas/twas_results/v6p/',trait,'_',tissue,'.csv')
  
  system(cmd,wait = T)
  
}
