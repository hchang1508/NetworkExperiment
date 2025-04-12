library(xtable)
source('0_header.R')

input= commandArgs(trailingOnly=TRUE)
index=as.numeric(input[1])
index2=as.numeric(input[2])
expo1=1
expo2=as.numeric(index)
weight=as.numeric(index2)
#load required data
load_data_imputed_Oct14(expo1,expo2,weight) 


extract_result_imputed_Nov11(expo1,expo2,'Nov11',weight)
  
