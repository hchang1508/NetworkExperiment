##header
##header
library(RSpectra)

working_path='/home/hc654/palmer_scratch/Bernoulli_compareDs/'
setwd(working_path)

input= commandArgs(trailingOnly=TRUE)
expo1=as.numeric(input[1])
expo2=as.numeric(input[2])
case= as.numeric(input[3])


foso_file=paste0('/home/hc654/palmer_scratch/Bernoulli_compareDs/',case,'_Bernoulli_',expo1,expo2,'_D.Rdata')
load(foso_file)
fo=mean_pre
so=cov_pre
fo=unlist(fo)
nrow=length(fo)
  
data_path_D=paste0(working_path,case,'_Bernoulli_',expo1,expo2,'_D.csv')
data_path_p=paste0(working_path,case,'_Bernoulli_',expo1,expo2,'_p.csv')
  

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

#foso_file=paste0('/home/hc654/palmer_scratch/final_analysis_outputs/',expo1,expo2,'_D.Rdata')
