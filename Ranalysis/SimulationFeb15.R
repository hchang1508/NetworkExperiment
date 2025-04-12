##header
rm(list=ls())
library(statnet)
#set working path
working_path='/Volumes/GoogleDrive/My Drive/UnifiedDesignBasedInference/Cai(2015)Insurance/Ranalysis'
setwd(working_path)
source('functions.R')
source('functions2.R')
#load network, outcome and covar data
data_path='/Volumes/GoogleDrive/My Drive/UnifiedDesignBasedInference/Cai(2015)Insurance/Ranalysis/data/truncated_graph'
data_path=paste0(data_path,'/net.RData')
load(data_path)

##Import FO and SO matrices somewhere
components=component.dist(net,connected='weak')
memberships=components$membership
csize = components$csize

#drop cluster 12 and 29


if (expo1 == 7 | expo1 == 12 | expo2 ==7 | expo2==12){
  
  number_dropped_t =sum( (id_5[degree_5<3] %in% id_list))
  
  #delete units with less than three friends
  id_list=id_list[!id_list %in% id_5[degree_5<3]] 
  
}else if (expo1 == 6 | expo1 == 11 | expo2 ==6 | expo2==11){
  
  number_dropped_t =sum( (id_5[degree_5<2] %in% id_list))
  
  #delete units with less than two friends
  id_list=id_list[!id_list %in% id_5[degree_5<2]] 
  
}else{
  
  number_dropped_t =sum( (id_5[degree_5<1] %in% id_list))
  
  #delete units with no friends
  id_list=id_list[!id_list %in% id_5[degree_5<1]] 
  
}

comp_list=1:36
comp_list=comp_list[-c(7,8)]
net_t = get.inducedSubgraph(net,which(memberships %in% comp_list))
prob=c(0.45,0.45,0.05,0.05)
fo=list()
so=list()

components_t=component.dist(net_t,connected='weak')
memberships_t=components_t$membership
csize_t = components_t$csize
for (i in (unique(memberships_t))){
  print(i)
  net_tt=get.inducedSubgraph(net_t,which(memberships_t==i))
  foso_t=FOSO(net_tt,individual_exposure,num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list,1,2)
  fo[[i]]=foso_t[[1]]
  so[[i]]=foso_t[[2]]
}
fo_vec=unlist(fo)


#basic simulation information
pol_size = network.size(net_t)
status = 0:3
num_mappings = 12
network_subjects = net_t %v% 'vertex.names'
subjects = (1:pol_size)[network_subjects %in% id_list]
subject_size = length(subjects)

#simulation parameters
nsim=100
expo1=1
expo2=2
###read from some input files

pi1 = fo_vec[(1:length(fo_vec))%%2==1]
pi2 = fo_vec[(1:length(fo_vec))%%2==0]

##divide into components####
components_t=component.dist(net_t,connected='weak')
memberships_t=components_t$membership
csize = components_t$csize


#extract subjects of interest
subjects_t = (1:pol_size)[network_subjects %in% id_list]
#extract membership information
subjects_memberships_t=memberships_t[network_subjects %in% id_list]
#sort subjects according to component index
subjects_t=subjects_t[order(subjects_memberships_t)]

#create dictionary for group belongings
group_index = cumsum(table(subjects_memberships_t))
dict_group=matrix(0,length(group_index),4)
dict_group[,1]= c(0,group_index[1:(length(group_index)-1)])+1
dict_group[,2]= group_index
dict_group[,3]=c(0,2*group_index[1:(length(group_index)-1)])+1
dict_group[,4]= group_index*2

#contrast vectors
contrast=c(-1,1)
#prob
pi1 = pi1[subjects_t]
pi2 = pi2[subjects_t]

#method=c('HT','HA','OLS_joint','OLS_separate','WLS_joint','WLS_separate','unweighted_joint',
#         'unweighted_sepaprate','weighted_joint','weighted_separate')
method=c('HT','HA','OLS_joint','WLS_joint','unweighted_joint', 'weighted_joint')
est_sim=matrix(0,nrow=nsim,ncol=length(method)) 
var_sim=matrix(0,nrow=nsim,ncol=length(method)) 

colnames(est_sim)=method
colnames(var_sim)=method

#prepare y and x informations, in the order of group beloings
Y=Y[match(network_subjects,Y$id),]
Y=Y[subjects,] 
x=X[match(network_subjects,X$id),]
x=x[subjects,]

#demean X
#pre-treatment covarates and demeaning
x=x[,c(3,4)] #use only one covariate
for (i in 1:ncol(x)){
  x[,i]=x[,i]-mean(x[,i])
}
y=Y[,2]
x=as.matrix(x)
sample_size = matrix(0,nrow=2,ncol = nsim)
friends=matrix(FALSE,nrow=pol_size,ncol=5)
for (i in 1:pol_size){
  friend_index_t=get.neighborhood(net_t,i,'out')
  if (length(friend_index_t)==0){
    next 
  }else
    friends[i,1:length(friend_index_t)] = friend_index_t
}

for (i in 1:nsim){
  
  print(i)
  realized_assignment_t=sample(status,pol_size,replace=TRUE,prob=prob)
  exposure_t = matrix(0,nrow=pol_size,ncol=7)
  exposure_t[,1]=1:pol_size #vertex indices
  exposure_t[,2]=realized_assignment_t==0
  exposure_t[,3]=realized_assignment_t==1
  exposure_t[,4]=realized_assignment_t==2
  exposure_t[,5]=realized_assignment_t==3
  
  
  num_friend_treated_t = apply(friends,1,friend_treated,realized_assignment_t)
  exposure_t[,6]=num_friend_treated_t[1,]
  exposure_t[,7]=num_friend_treated_t[2,]  
  
  
  
  exposure_t = exposure_t[subjects,2:7]
  
  #please double check every time 
  if (is.null(dim(exposure_t))){
    state_t=individual_exposure(exposure_t)
    state_t = state_t[c(expo1,expo2)]
    y1_t = state_t[expo1] * y
    y2_t = state_t[expo2] * y
    denom_pi1_t= state_t[expo1]* 1/pi1 
    denom_pi2_t= state_t[expo2]* 1/pi2 
    pi_weight_t = denom_pi1_t + denom_pi2_t
    d1_t=state_t[expo2]
    d0_t=state_t[expo1]
  }else{
    state_t=apply(exposure_t,1,individual_exposure) 
    state_t = state_t[c(expo1,expo2),]  
    y1_t = state_t[expo1,] * y
    y2_t = state_t[expo2,] * y
    denom_pi1_t= state_t[expo1,]* 1/pi1 
    denom_pi2_t= state_t[expo2,]* 1/pi2 
    pi_weight_t = denom_pi1_t + denom_pi2_t
    d1_t=state_t[expo2,]
    d0_t=state_t[expo1,]
  }
  

  
  
  y_obs_t=y1_t+y2_t
  
  # est_sim[i,1]= HT(y2,y1,pi2,pi1,subject_size)
  # est_sim[i,2]= HA(y2,y1,pi2,pi1,denom_pi2,denom_pi1)
  # est_sim[i,3]= OLS(y_obs,d1,d0,pi_weight,x,mode='OLS_joint')
  # est_sim[i,4]= OLS(y_obs,d1,d0,pi_weight,x,mode='OLS_separate')    
  # est_sim[i,5]= OLS(y_obs,d1,d0,pi_weight,x,mode='WLS_joint')    
  # est_sim[i,6]= OLS(y_obs,d1,d0,pi_weight,x,mode='WLS_separate')    
  # est_sim[i,7]= Logit(y_obs,d1,d0,pi_weight,x,pi1,pi2,mode='unweighted_joint')
  # est_sim[i,8]= Logit(y_obs,d1,d0,pi_weight,x,pi1,pi2,mode='unweighted_separate')
  # est_sim[i,9]= Logit(y_obs,d1,d0,pi_weight,x,pi1,pi2,mode='weighted_joint')
  # est_sim[i,10]= Logit(y_obs,d1,d0,pi_weight,x,pi1,pi2,mode='weighted_separate')
  
  
  est_sim[i,1]= HT_full(y2_t,y1_t,pi2,pi1,subject_size)
  est_sim[i,2]= HA_full(y2_t,y1_t,pi2,pi1,denom_pi2_t,denom_pi1_t)
  est_sim[i,3]= OLS_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint')
  est_sim[i,4]= OLS_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint')
  est_sim[i,5]= Logit_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='unweighted_joint')
  est_sim[i,6]= Logit_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='weighted_joint')
  # 
  var_sim[i,1]= HT_var_full(y_obs_t,d1_t,d0_t,subject_size,fo,so,contrast,dict_group)
  var_sim[i,2]= HA_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,subject_size,fo,so,contrast,dict_group)
  var_sim[i,3]= OLS_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint',fo,so,contrast,dict_group)
  var_sim[i,4]= OLS_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint',fo,so,contrast,dict_group)
  var_sim[i,5]= Logit_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='unweighted_joint',fo,so,contrast,dict_group)
  var_sim[i,6]= Logit_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='weighted_joint',fo,so,contrast,dict_group)
  # 
}
apply(est_sim,2,mean)
apply(est_sim,2,var)
apply(var_sim,2,mean)
theoretical_variance_full(y,x,fo,so,contrast,dict_group)
theoretical_variance_bound_full(y,x,fo,so,contrast,dict_group)
output=list(est_sim,var_sim)
