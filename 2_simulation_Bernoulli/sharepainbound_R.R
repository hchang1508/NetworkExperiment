source('header.R')
foso_file=paste0('/home/hc654/palmer_scratch/final_analysis_Dfs/','Stratified_nat_',expo1,expo2,'_D.Rdata')
load(foso_file)
fo=mean_pre
so=cov_pre

so_bound=so
for (i in 1:nrow(dict_group)){
  
  c=dict_group[i,1]
  index_start=dict_group[i,4]
  index_end=dict_group[i,5]    
  
 
  fo_t=fo[[c]]
  so_t=so[[c]]
  
  normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
  
  normalized_cov_t = round(normalized_cov_t,10)
  minus_ones_t=- (abs(normalized_cov_t+1)<0.001)
  
  print(dim(normalized_cov_t))
  add_on=share_pain_bound_old(minus_ones_t)
  print(paste0('max_size',(max(add_on))))
  so_bound[[c]] = normalized_cov_t+ add_on
  

}

so_bound_file=paste0('/home/hc654/palmer_scratch/final_analysis_Ds/','Stratified_nat_',expo1,expo2,'_Dbound.Rdata')
save(so_bound,file=so_bound_file)
load(so_bound_file)