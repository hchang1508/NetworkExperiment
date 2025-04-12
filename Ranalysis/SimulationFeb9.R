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


#First select a small cluster
components=component.dist(net,connected='weak')
memberships=components$membership
csize = components$csize
csize[10]

#select the 10th clusters
comp_10 = get.inducedSubgraph(net,which(memberships==10))
summary(comp_10,print.adj = FALSE)

#check number of experiment subjects
network.size(comp_10)  #76 units
id_10 = comp_10 %v% 'vertex.names'
sum(id_10 %in% id_list) #73 experiment subjects
degree(comp_10,cmode='outdegree') #warning some units have only two friends so some exposure mappings are not defined


#potential outcomes and sharp null
Y_10=Y[match(id_10,Y$id),]


#calculate FO; simulated estimates matches with theoretical estimates
prob=c(0.25,0.25,0.25,0.25)
pi= FO(net=comp_10, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list=id_list,expo1=1,expo2=2)
pi1 = pi[(1:length(pi))%%2==1]
pi2 = pi[(1:length(pi))%%2==0]

#################################################################
#################################################################
#################################################################
#################################################################
#HT Compare 1-2
nsim=10000
pi12= FO(net=comp_10, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list=id_list,expo1=1,expo2=2)
pi1 = pi12[(1:length(pi12))%%2==1]
pi2 = pi12[(1:length(pi12))%%2==0]
pol_size = network.size(comp_10)
status = 0:3-1
num_mappings = 12
assign_prob=prob
network_subjects = comp_10 %v% 'vertex.names'
subjects = (1:pol_size)[network_subjects %in% id_list]
expo1=1
expo2=2
subject_size = length(subjects)

#potential outcomes: sharp null
Y_10=Y_10[match(network_subjects,Y_10$id),]
Y_10=Y_10[subjects,]

#prob
pi1 = pi1[subjects]
pi2 = pi2[subjects]

est_sim = rep(NA,nsim)

friends=matrix(FALSE,nrow=pol_size,ncol=5)
for (i in 1:pol_size){
  friend_index=get.neighborhood(comp_10,i,'out')
  if (length(friend_index)==0){
    next 
  }else
    friends[i,1:length(friend_index)] = friend_index
}

for (i in 1:nsim){
  
  print(i)
  realized_assignment=sample(status,pol_size,replace=TRUE,prob=assign_prob)
  
  exposure = matrix(0,nrow=pol_size,ncol=7)
  exposure[,1]=1:pol_size #vertex indices
  exposure[,2]=realized_assignment==0
  exposure[,3]=realized_assignment==1
  exposure[,4]=realized_assignment==2
  exposure[,5]=realized_assignment==3
  

  num_friend_treated = apply(friends,1,friend_treated,realized_assignment)
  exposure[,6]= num_friend_treated[1,]
  exposure[,7]=num_friend_treated[2,]  
  
  
  exposure_temp = exposure[subjects,2:7]
  
  #please double check every time 
  
  
  state_temp=apply(exposure_temp,1,individual_exposure) 
  y1 = state_temp[expo1,] * Y_10[,2]
  y2 = state_temp[expo2,] * (Y_10[,2]+2)

  est_sim[i]= HT(y2,y1,pi2,pi1,73)
}
mean(est_sim)


#################################################################
#################################################################
#################################################################
#################################################################
#HT Compare 4-5
nsim=10000
pi45= FO(net=comp_10, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list=id_list,expo1=4,expo2=5)
pi1 = pi45[(1:length(pi56))%%2==1]
pi2 = pi45[(1:length(pi56))%%2==0]
pol_size = network.size(comp_10)
status = 0:3-1
num_mappings = 12
assign_prob=prob
network_subjects = comp_10 %v% 'vertex.names'
subjects = (1:pol_size)[network_subjects %in% id_list]
expo1=4
expo2=5
subject_size = length(subjects)

#potential outcomes: sharp null
Y_10=Y_10[match(network_subjects,Y_10$id),]
Y_10=Y_10[subjects,]

#prob
pi1 = pi1[subjects]
pi2 = pi2[subjects]

est_sim = rep(NA,nsim)

friends=matrix(FALSE,nrow=pol_size,ncol=5)
for (i in 1:pol_size){
  friend_index=get.neighborhood(comp_10,i,'out')
  if (length(friend_index)==0){
    next 
  }else
    friends[i,1:length(friend_index)] = friend_index
}

for (i in 1:nsim){
  
  print(i)
  realized_assignment=sample(status,pol_size,replace=TRUE,prob=assign_prob)
  
  exposure = matrix(0,nrow=pol_size,ncol=7)
  exposure[,1]=1:pol_size #vertex indices
  exposure[,2]=realized_assignment==0
  exposure[,3]=realized_assignment==1
  exposure[,4]=realized_assignment==2
  exposure[,5]=realized_assignment==3
  
  
  num_friend_treated = apply(friends,1,friend_treated,realized_assignment)
  exposure[,6]= num_friend_treated[1,]
  exposure[,7]=num_friend_treated[2,]  
  
  
  exposure_temp = exposure[subjects,2:7]
  
  #please double check every time 
  
  
  state_temp=apply(exposure_temp,1,individual_exposure) 
  y1 = state_temp[expo1,] * Y_10[,2]
  y2 = state_temp[expo2,] * (Y_10[,2]+2)
  
  est_sim[i]= HT(y2,y1,pi2,pi1,subject_size)
}
mean(est_sim)


#################################################################
#################################################################
#################################################################
#################################################################
#HA Compare 1-2
nsim=10000
pi12= FO(net=comp_10, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list=id_list,expo1=1,expo2=2)
pi1 = pi12[(1:length(pi12))%%2==1]
pi2 = pi12[(1:length(pi12))%%2==0]
pol_size = network.size(comp_10)
status = 0:3-1
num_mappings = 12
assign_prob=prob
network_subjects = comp_10 %v% 'vertex.names'
subjects = (1:pol_size)[network_subjects %in% id_list]
expo1=1
expo2=2
subject_size = length(subjects)

#potential outcomes: sharp null
Y_10=Y_10[match(network_subjects,Y_10$id),]
Y_10=Y_10[subjects,]

#prob
pi1 = pi1[subjects]
pi2 = pi2[subjects]

est_sim = rep(NA,nsim)

friends=matrix(FALSE,nrow=pol_size,ncol=5)
for (i in 1:pol_size){
  friend_index=get.neighborhood(comp_10,i,'out')
  if (length(friend_index)==0){
    next 
  }else
    friends[i,1:length(friend_index)] = friend_index
}

for (i in 1:nsim){
  
  print(i)
  realized_assignment=sample(status,pol_size,replace=TRUE,prob=assign_prob)
  
  exposure = matrix(0,nrow=pol_size,ncol=7)
  exposure[,1]=1:pol_size #vertex indices
  exposure[,2]=realized_assignment==0
  exposure[,3]=realized_assignment==1
  exposure[,4]=realized_assignment==2
  exposure[,5]=realized_assignment==3
  
  
  num_friend_treated = apply(friends,1,friend_treated,realized_assignment)
  exposure[,6]= num_friend_treated[1,]
  exposure[,7]=num_friend_treated[2,]  
  
  
  exposure_temp = exposure[subjects,2:7]
  
  #please double check every time 
  
  
  state_temp=apply(exposure_temp,1,individual_exposure) 
  y1 = state_temp[expo1,] * Y_10[,2]
  y2 = state_temp[expo2,] * (Y_10[,2]+2)
  denom_pi1= state_temp[expo1,]* 1/pi1 
  denom_pi2= state_temp[expo2,]* 1/pi2 
  est_sim[i]= HA(y2,y1,pi2,pi1,denom_pi2,denom_pi1)
}
mean(est_sim)

#################################################################
#################################################################
#################################################################
#################################################################
#HA Compare 4-5
nsim=10000
pi45= FO(net=comp_10, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list=id_list,expo1=4,expo2=5)
pi1 = pi45[(1:length(pi45))%%2==1]
pi2 = pi45[(1:length(pi45))%%2==0]
pol_size = network.size(comp_10)
status = 0:3-1
num_mappings = 12
assign_prob=prob
network_subjects = comp_10 %v% 'vertex.names'
subjects = (1:pol_size)[network_subjects %in% id_list]
expo1=4
expo2=5
subject_size = length(subjects)

#potential outcomes: sharp null
Y_10=Y_10[match(network_subjects,Y_10$id),]
Y_10=Y_10[subjects,]

#prob
pi1 = pi1[subjects]
pi2 = pi2[subjects]

est_sim = rep(NA,nsim)

friends=matrix(FALSE,nrow=pol_size,ncol=5)
for (i in 1:pol_size){
  friend_index=get.neighborhood(comp_10,i,'out')
  if (length(friend_index)==0){
    next 
  }else
    friends[i,1:length(friend_index)] = friend_index
}

for (i in 1:nsim){
  
  print(i)
  realized_assignment=sample(status,pol_size,replace=TRUE,prob=assign_prob)
  
  exposure = matrix(0,nrow=pol_size,ncol=7)
  exposure[,1]=1:pol_size #vertex indices
  exposure[,2]=realized_assignment==0
  exposure[,3]=realized_assignment==1
  exposure[,4]=realized_assignment==2
  exposure[,5]=realized_assignment==3
  
  
  num_friend_treated = apply(friends,1,friend_treated,realized_assignment)
  exposure[,6]= num_friend_treated[1,]
  exposure[,7]=num_friend_treated[2,]  
  
  
  exposure_temp = exposure[subjects,2:7]
  
  #please double check every time 
  
  
  state_temp=apply(exposure_temp,1,individual_exposure) 
  y1 = state_temp[expo1,] * Y_10[,2]
  y2 = state_temp[expo2,] * (Y_10[,2]+2)
  denom_pi1= state_temp[expo1,]* 1/pi1 
  denom_pi2= state_temp[expo2,]* 1/pi2 
  if (sum(state_temp[expo1,])==0){ print('empty!!')}
  if (sum(state_temp[expo2,])==0){ print('empty!!')}
  est_sim[i]= HA(y2,y1,pi2,pi1,denom_pi2,denom_pi1)
}
mean(est_sim)

#################################################################
#################################################################
comp_7=get.inducedSubgraph(net,which(memberships==7))
network.size(comp_7)  #1914 units
id_7 = comp_7 %v% 'vertex.names'
sum(id_7 %in% id_list) #1777 experiment subjects
degree(comp_7,cmode='outdegree') #warning some units have only two friends so some exposure mappings are not defined


#potential outcomes and sharp null
Y_7=Y[match(id_7,Y$id),]
X_7=X[match(id_7,X$id),]
#################################################################
#################################################################
#HA Compare 4-5, largest component
comp_7 = get.inducedSubgraph(net,which(memberships==7))
nsim=10000
pi45= FO(net=comp_7, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list=id_list,expo1=4,expo2=5)
pi1 = pi45[(1:length(pi45))%%2==1]
pi2 = pi45[(1:length(pi45))%%2==0]
pol_size = network.size(comp_7)
status = 0:3-1
num_mappings = 12
assign_prob=prob
network_subjects = comp_7 %v% 'vertex.names'
subjects = (1:pol_size)[network_subjects %in% id_list]
expo1=4
expo2=5
subject_size = length(subjects)

#potential outcomes: sharp null
Y_7=Y_7[match(network_subjects,Y_7$id),]
Y_7=Y_7[subjects,]

#prob
pi1 = pi1[subjects]
pi2 = pi2[subjects]

est_sim = rep(NA,nsim)

friends=matrix(FALSE,nrow=pol_size,ncol=5)
for (i in 1:pol_size){
  friend_index=get.neighborhood(comp_7,i,'out')
  if (length(friend_index)==0){
    next 
  }else
    friends[i,1:length(friend_index)] = friend_index
}

for (i in 1:nsim){
  
  print(i)
  realized_assignment=sample(status,pol_size,replace=TRUE,prob=assign_prob)
  
  exposure = matrix(0,nrow=pol_size,ncol=7)
  exposure[,1]=1:pol_size #vertex indices
  exposure[,2]=realized_assignment==0
  exposure[,3]=realized_assignment==1
  exposure[,4]=realized_assignment==2
  exposure[,5]=realized_assignment==3
  
  
  num_friend_treated = apply(friends,1,friend_treated,realized_assignment)
  exposure[,6]= num_friend_treated[1,]
  exposure[,7]=num_friend_treated[2,]  
  
  
  exposure_temp = exposure[subjects,2:7]
  
  #please double check every time 
  
  
  state_temp=apply(exposure_temp,1,individual_exposure) 
  y1 = state_temp[expo1,] * Y_7[,2]
  y2 = state_temp[expo2,] * (Y_7[,2]+2)
  denom_pi1= state_temp[expo1,]* 1/pi1 
  denom_pi2= state_temp[expo2,]* 1/pi2 
  if (sum(state_temp[expo1,])==0){ 
    print('empty!!') 
    break}
  if (sum(state_temp[expo2,])==0){ 
    print('empty!!') 
    break}
  est_sim[i]= HA(y2,y1,pi2,pi1,denom_pi2,denom_pi1)
}
mean(est_sim)

#################################################################
#################################################################
#################################################################
#################################################################
#OLS Compare 4-5, largest component 7
comp_7 = get.inducedSubgraph(net,which(memberships==7))
nsim=10000
pi45= FO(net=comp_7, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list=id_list,expo1=4,expo2=5)
pi1 = pi45[(1:length(pi45))%%2==1]
pi2 = pi45[(1:length(pi45))%%2==0]
pol_size = network.size(comp_7)
status = 0:3-1
num_mappings = 12
assign_prob=prob
network_subjects = comp_7 %v% 'vertex.names'
subjects = (1:pol_size)[network_subjects %in% id_list]
Y_add = degree(comp_7,cmode='outdegree')[subjects]
expo1=4
expo2=5
subject_size = length(subjects)

#potential outcomes: sharp null
Y_7=Y_7[match(network_subjects,Y_7$id),]
Y_7=Y_7[subjects,]
X_7=X_7[match(network_subjects,X_7$id),]
X_7=X_7[subjects,]
x=as.matrix(X_7)
x=x[,c(-1,-14)]
for (i in 1:ncol(x)){
  
  x[,i]=x[,i]-mean(x[,i])
}

#prob
pi1 = pi1[subjects]
pi2 = pi2[subjects]

est_sim = rep(NA,nsim)

friends=matrix(FALSE,nrow=pol_size,ncol=5)
for (i in 1:pol_size){
  friend_index=get.neighborhood(comp_7,i,'out')
  if (length(friend_index)==0){
    next 
  }else
    friends[i,1:length(friend_index)] = friend_index
}

for (i in 1:nsim){
  
  print(i)
  realized_assignment=sample(status,pol_size,replace=TRUE,prob=assign_prob)
  
  exposure = matrix(0,nrow=pol_size,ncol=7)
  exposure[,1]=1:pol_size #vertex indices
  exposure[,2]=realized_assignment==0
  exposure[,3]=realized_assignment==1
  exposure[,4]=realized_assignment==2
  exposure[,5]=realized_assignment==3
  
  
  num_friend_treated = apply(friends,1,friend_treated,realized_assignment)
  exposure[,6]= num_friend_treated[1,]
  exposure[,7]=num_friend_treated[2,]  
  
  
  exposure_temp = exposure[subjects,2:7]
  
  #please double check every time 
  
  
  state_temp=apply(exposure_temp,1,individual_exposure) 
  y1 = state_temp[expo1,] * Y_7[,2]
  y2 = state_temp[expo2,] * (Y_7[,2]+Y_add)
  y = y1 + y2 
  denom_pi1= state_temp[expo1,]* 1/pi1 
  denom_pi2= state_temp[expo2,]* 1/pi2 

  pi_weight = denom_pi1 + denom_pi2
  d= state_temp[expo2,]
  if (sum(state_temp[expo1,])==0){ 
    print('empty!!') 
    break}
  if (sum(state_temp[expo2,])==0){ 
    print('empty!!') 
    break}
  est_sim[i]= OLS(y,d,pi_weight,x,mode='WLS_joint')
}
mean(est_sim)

#################################################################
#################################################################
#################################################################
#################################################################
#Logit Compare 4-5, largest component 7
comp_7 = get.inducedSubgraph(net,which(memberships==7))
nsim=10000
pi45= FO(net=comp_7, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list=id_list,expo1=4,expo2=5)
pi1 = pi45[(1:length(pi45))%%2==1]
pi2 = pi45[(1:length(pi45))%%2==0]
pol_size = network.size(comp_7)
status = 0:3-1
num_mappings = 12
assign_prob=prob
network_subjects = comp_7 %v% 'vertex.names'
subjects = (1:pol_size)[network_subjects %in% id_list]
Y_add = degree(comp_7,cmode='outdegree')[subjects]
expo1=4
expo2=5
subject_size = length(subjects)

#potential outcomes: sharp null
Y_7=Y_7[match(network_subjects,Y_7$id),]
Y_7=Y_7[subjects,]
X_7=X_7[match(network_subjects,X_7$id),]
X_7=X_7[subjects,]
x=as.matrix(X_7)
x=x[,c(-1,-14)]
for (i in 1:ncol(x)){
  
  x[,i]=x[,i]-mean(x[,i])
}

#prob
pi1 = pi1[subjects]
pi2 = pi2[subjects]

est_sim = rep(NA,nsim)

friends=matrix(FALSE,nrow=pol_size,ncol=5)
for (i in 1:pol_size){
  friend_index=get.neighborhood(comp_7,i,'out')
  if (length(friend_index)==0){
    next 
  }else
    friends[i,1:length(friend_index)] = friend_index
}

for (i in 1:nsim){
  
  print(i)
  realized_assignment=sample(status,pol_size,replace=TRUE,prob=assign_prob)
  
  exposure = matrix(0,nrow=pol_size,ncol=7)
  exposure[,1]=1:pol_size #vertex indices
  exposure[,2]=realized_assignment==0
  exposure[,3]=realized_assignment==1
  exposure[,4]=realized_assignment==2
  exposure[,5]=realized_assignment==3
  
  
  num_friend_treated = apply(friends,1,friend_treated,realized_assignment)
  exposure[,6]= num_friend_treated[1,]
  exposure[,7]=num_friend_treated[2,]  
  
  
  exposure_temp = exposure[subjects,2:7]
  
  #please double check every time 
  
  
  state_temp=apply(exposure_temp,1,individual_exposure) 
  y1 = state_temp[expo1,] * Y_7[,2]
  y2 = state_temp[expo2,] * (Y_7[,2])
  y = y1 + y2 
  denom_pi1= state_temp[expo1,]* 1/pi1 
  denom_pi2= state_temp[expo2,]* 1/pi2 
  
  pi_weight = denom_pi1 + denom_pi2
  d1= state_temp[expo2,]
  d0=state_temp[expo1,]
  if (sum(state_temp[expo1,])==0){ 
    print('empty!!') 
    break}
  if (sum(state_temp[expo2,])==0){ 
    print('empty!!') 
    break}
  est_sim[i]= Logit(y,d1,d0,pi_weight,x,pi1,pi2,mode='weighted_joint')
}
mean(est_sim)

