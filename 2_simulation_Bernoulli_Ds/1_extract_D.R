rm(list=ls())
source('0_header.R')
input= commandArgs(trailingOnly=TRUE)

expo1=as.numeric(input[1])
expo2=as.numeric(input[2])


setwd('/home/hc654/palmer_scratch/final_analysis_Bernoulli/')

files=list.files()
#pattern=paste0('Stratified',expo1,'_',expo2,'_','FOSO*')
pattern=paste0('Bernoulli_nat_',expo1,'_',expo2,'_','FOSO*')
#pattern=paste0('Bernoulli_nat_',expo1,'_',expo2,'_','FOSO*')
temp=files[grep(pattern,files )]


print(paste0('Processing ',expo1,' ',expo2))

#load initial observation
load(temp[1])
mean_pre=output[[1]]
cov_pre=output[[2]]
n_old=10000
#
memberships=net %v% "membership"
num_components=length(unique(memberships))
membership_list=unique(memberships)


for (i in 2:length(temp)){
  print(i)
  load(temp[i]) 
  mean_new=output[[1]]
  cov_new=output[[2]]
  
  for (j in membership_list){
    n_new=10000

    output_temp=WELFORD_ONLINE_BATCH(mean_pre[[j]],cov_pre[[j]],n_old,n_new,mean_new[[j]],cov_new[[j]])
    
    mean_pre[[j]]=output_temp[[1]]
    cov_pre[[j]]=output_temp[[2]]
  }
  
  n_old=n_old + 10000
}

name=paste0('Bernoulli_',expo1,expo2,'_D.Rdata')
save(mean_pre,cov_pre,file=name)
