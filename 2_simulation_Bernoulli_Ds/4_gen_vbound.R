###############################################################################
#####This function pre-calculates the variance bound ##########################
###############################################################################
#clear environment
rm(list=ls())

#set working path
setwd('/home/hc654/NetworkExperiment/2_simulation_Bernoulli/')

#import functions
source('0_header.R')
input= commandArgs(trailingOnly=TRUE)
expo1=as.numeric(input[1])
expo2=as.numeric(input[2])
case=as.numeric(input[3])

##################################
#import foso information##########
##################################
foso_file=paste0('/home/hc654/palmer_scratch/Bernoulli_compareDs/',case,'_Bernoulli_',expo1,expo2,'_D.Rdata')
load(foso_file)
fo=mean_pre
so=cov_pre

##################################
#IMPORT FOSO INFORMATION##########
##################################

###basic simulation information
#network size
pol_size = network.size(net)

#id information
id_t = net %v% 'vertex.names'
degree_t=degree(net,cmode='outdegree') #warning some units have only two friends so some exposure mappings are not defined


#depending on exposure mapping some units must be dropped for identification reasons:
if (expo1 == 7 | expo1 == 12 | expo2 ==7 | expo2==12){
  
  number_dropped_t =sum( (id_t[degree_t<3] %in% id_list))
  #delete units with less than three friends
  id_list=id_list[!id_list %in% id_t[degree_t<3]] 
}else if (expo1 == 6 | expo1 == 11 | expo2 ==6 | expo2==11){
  
  number_dropped_t =sum( (id_t[degree_t<2] %in% id_list))
  #delete units with less than two friends
  id_list=id_list[!id_list %in% id_t[degree_t<2]]
  
}else{
  
  number_dropped_t =sum( (id_t[degree_t<1] %in% id_list))
  #delete units with no friends
  id_list=id_list[!id_list %in% id_t[degree_t<1]] 
  
}

######################################################
#EXTRACT SUBJECTS OF INTEREST#########################
######################################################
##divide into components####
memberships_t=net %v% 'membership'

#extract subjects of interest
subjects_t = (1:pol_size)[id_t %in% id_list]
#extract membership information
subjects_memberships_t=memberships_t[id_t %in% id_list]
#sort subjects according to component index
subjects_t=subjects_t[order(subjects_memberships_t)]


######################################################
#CREATE DICTIONARY FOR GROUP BELONGINGS###############
######################################################

group_index = cumsum(table(subjects_memberships_t))
dict_group=matrix(0,length(group_index),5)
dict_group[,1]=names(group_index)
dict_group[,2]= c(0,group_index[1:(length(group_index)-1)])+1
dict_group[,3]= group_index
dict_group[,4]=c(0,2*group_index[1:(length(group_index)-1)])+1
dict_group[,5]= group_index*2
dict_group<- matrix(as.numeric(dict_group),   ncol = ncol(dict_group))


######################################################
#CREATE VARIANCE BOUND################################
######################################################

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
  print(paste0('max_size',(max(add_on)))) #maximize value that is added on
  so_bound[[c]] = normalized_cov_t+ add_on
  

}

######################################################
#SANITY CHECK#########################################
######################################################

for (i in 1:nrow(dict_group)){
  
  c=dict_group[i,1]
  index_start=dict_group[i,4]
  index_end=dict_group[i,5]    
  
  
  fo_t=fo[[c]]
  so_t=so[[c]]
  
  normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
  
  normalized_cov_t = round(normalized_cov_t,10)
  minus_ones_t=- (abs(normalized_cov_t+1)<0.001)
  
  #this should return a 0: every entries with -1 should be zeroed out
  print(c)
  print( sum( abs( (minus_ones_t * so_bound[[c]]))>0.001))

  
}

######################################################
#CREATE OUTPUT########################################
######################################################
so_bound_file=paste0('/home/hc654/palmer_scratch/Bernoulli_compareDs/',case,'_Bernoulli_',expo1,expo2,'_Dbound.Rdata')
save(so_bound,file=so_bound_file)
