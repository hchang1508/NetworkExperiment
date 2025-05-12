rm(list=ls())
library(xtable)
source('0_header.R')

input= commandArgs(trailingOnly=TRUE)
index=as.numeric(input[1])
index2=as.numeric(input[2])
expo1=1
expo2=as.numeric(index)
weight=as.numeric(index2)
#load required data
temp=load_data_imputed(case=4,expo1=4,expo2=6) 


extract_result_imputed(expo1,expo2,y1=temp$y1, y0=temp$y0, fo=temp$fo,so=temp$so, so_AS2=temp$so_AS2, dict_group=temp$dict_group,x=temp$x,subject_size=length(temp$y1))
  
y1=temp$y1
y0=temp$y0
fo=temp$fo
so=temp$so
so_AS2=temp$so_AS2
dict_group=temp$dict_group
x=temp$x
subject_size=length(temp$y1)