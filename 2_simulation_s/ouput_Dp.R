##header
##header
library(RSpectra)

working_path='/home/hc654/palmer_scratch/final_analysis_Ds/'
setwd(working_path)

expo1=1
expo2=3
for (j in 2:7){
  expo2=j
  print(j)
  foso_file=paste0('/home/hc654/palmer_scratch/final_analysis_Ds/','Stratified_nat_',expo1,expo2,'_D.Rdata')
  load(foso_file)
  fo=mean_pre
  so=cov_pre
  fo=unlist(fo)
  nrow=length(fo)
  
  data_path_D=paste0(working_path,'Stratified_nat_',expo1,expo2,'_D.csv')
  data_path_p=paste0(working_path,'Stratified_nat_',expo1,expo2,'_p.csv')
  

  D_output=matrix(0,nrow,nrow)
  
  max_component=41
  index_start=1
  
  for (c in 1:max_component){
    print(c)
    if (!is.null(dim(so[[c]])[1])){
    index_end=index_start + dim(so[[c]])[1]-1
    
    D_output[index_start:index_end,index_start:index_end]=so[[c]]
    
    index_start = index_end +1
    }
  }
  
  write.csv(fo,data_path_p)
  write.csv(D_output,data_path_D)
}
#foso_file=paste0('/home/hc654/palmer_scratch/final_analysis_outputs/',expo1,expo2,'_D.Rdata')

#

id_t = net %v% 'vertex.names'
pol_size = network.size(net)

subjects_t = (1:pol_size)[id_t %in% id_list]
memberships_t=net %v% 'membership'
#extract membership information
subjects_memberships_t=memberships_t[id_t %in% id_list]
#sort subjects according to component index
subjects_t=subjects_t[order(subjects_memberships_t)]
subject_size=length(subjects_t)

x=all_info2[match(id_t,all_info2$id),]
x=x[subjects_t,]
x=x[,c('male','age','agpop','literacy','ricearea_2010','risk_averse','disaster_prob')] 
for (i in 1:ncol(x)){
  x[,i]=x[,i]-mean(x[,i])
}
x=as.matrix(x)
data_path_x=paste0(working_path,'output/','x.csv')
write.csv(x,data_path_x)


  