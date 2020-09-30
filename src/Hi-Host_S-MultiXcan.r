#ml GCC/6.4.0-2.28 Python/3.6.3 OpenMPI/2.1.1 numpy/1.13.1-Python-3.6.3 pandas/0.18.1-Python-3.6.3 scipy/0.19.1-Python-3.6.3  tabix/0.2.6 R
#source /data/c***/z***/tools/TIGAR/tiger_env/bin/activate
#export PYTHONPATH=$(python -c 'import sys; print(sys.path[-1])'):${PYTHONPATH}

args=as.numeric(commandArgs(TRUE))   #148 traits

#trait names
file_list<-dir('/data/g***/z***/hihost_gwas/summary/')
trait_list<-sub('_Meannew7.qfam.parents.perm_rsID.txt.gz','',file_list[grep('qfam.parents.perm_rsID.txt.gz',file_list)]) #
trait=trait_list[args]

#copy summary-stat
#one trait one folder (across all the tissues) 
setwd('/data/g***/z***/hihost_gwas/twas_results/')
dir.create(paste0('./v6p_by_trait/',trait))
cmd=(paste0('cp ./v6p/',trait,'* ./v6p_by_trait/',trait,'/'))
system(cmd,wait = T)

#run SMultixcan
cmd=paste0('python3.6 /data/c***/z***/tools/metaxcan_0.6.2/MetaXcan-master/software/SMulTiXcan.py --models_folder /data/g***/z***/data/gtex/v6p_haky/GTEx-V6p-HapMap-2016-09-08 --models_name_pattern "TW_(.*)_0.5.db" --snp_covariance /data/g***/z***/data/gtex/v6p_haky/snp_covariance_v6p.txt.gz   --metaxcan_folder /data/g***/z***/hihost_gwas/twas_results/v6p_by_trait/',trait,' --metaxcan_filter "(.*).csv" --metaxcan_file_name_parse_pattern "',trait,'_(.*).csv" --gwas_file /data/g***/z***/hihost_gwas/summary/',trait,'_Meannew7.qfam.parents.perm_rsID.txt.gz --snp_column SNP --non_effect_allele_column MajorA --effect_allele_column MinorA --beta_column Beta  --se_column EMP_se --cutoff_condition_number 30 --verbosity 7 --throw --output /data/g***/z***/hihost_gwas/smultixcan_results/v6p/',trait,'.csv')
system(cmd,wait = T)

