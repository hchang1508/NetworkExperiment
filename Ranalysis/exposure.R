library(statnet)
library(UserNetR)
library(RColorBrewer)

#input: network matrix (directed in our case)
#the network are assumed to be a "network" object, not igraph
memberships=component.dist(net,connected=c('weak'))$membership
net = get.inducedSubgraph(net,which(memberships==11))

#basic experimental parameters
pol_size = network.size(net)
assign_prob = c(0.25,0.25,0.25,0.25)
num_mapping=12
# 0 if first+simple; 1 if first+intense; 2 if second+intense; 3 if third intense 
realized_assignment=sample(0:3,pol_size,replace=TRUE,prob=assign_prob)

#components: divide into disconnected components, 

#compute exposure: 
#Col 1: id
#Col 2-5: self status
#Col 6: #friend in first round simple
#Col 7: #friend in first round intensive
exposure = matrix(0,nrow=pol_size,ncol=7)
exposure[,1]=1:pol_size #vertex indices
exposure[,2]=realized_assignment==0
exposure[,3]=realized_assignment==1
exposure[,4]=realized_assignment==2
exposure[,5]=realized_assignment==3

for (i in 1:pol_size){
    friend_index=get.neighborhood(net,i,'out')
    exposure[i,6]= sum(realized_assignment[friend_index]==0)
    exposure[i,7]= sum(realized_assignment[friend_index]==1)    
}
table(exposure[,6])
table(exposure[,7])

individual_exposure=function(exposure){
  
  output=rep(0,12)
  if (exposure[1]==1){
    output[1]=1  #first round simple
  }else if (exposure[2]==1){
    output[2]=1 #first round intense
  }else if (exposure[3]==1 & exposure[5]==0 & exposure[6]==0){
    output[3]=1 #second round simple, no friend first round 
  }else if (exposure[3]==1 & exposure[5]>=0 & exposure[6]==0){
    output[4]=1 #second round simple, friend from simple, not friend intensive
  }else if (exposure[3]==1 & exposure[6]==1){
    output[5]=1 #second round simple, one friend intensive
  }else if (exposure[3]==1 & exposure[6]==2){
    output[6]=1 #second round simple, two friend intensive
  }else if (exposure[3]==1 & exposure[6]>2){
    output[7]=1 #second round simple, two more friend intensive
  }else if (exposure[4]==1 & exposure[5]==0 & exposure[6]==0){
    output[8]=1 #second round intensive, not friend first round
  }else if (exposure[4]==1 & exposure[5]>=0 & exposure[6]==0){
    output[9]=1 #second round intensive, friend from simple, not friend intensive
  }else if (exposure[4]==1 & exposure[6]==1){
    output[10]=1 #second round intensive, one friendintensive
  }else if (exposure[4]==1 & exposure[6]==2){
    output[11]=1 #second round intensive, two friend intensive
  }else if (exposure[4]==1 & exposure[6]>2){
    output[12]=1 #second round intensive, two more friend intensive
  }
  
  return(output)
}







###old codes####

exposure = matrix(0,nrow=pol_size,ncol=3)
exposure[,1]=1:pol_size #vertex indices
exposure[,2]=realized_assignment


for (i in 1:pol_size){
  friend_index=get.neighborhood(toy_net,i,'out')
  exposure[i,3]= sum(realized_assignment[friend_index])
  
}
table(exposure[,3])

#change to the R vector
state=c()
for (i in 1:pol_size){
  
    state=c(state,individual_exposure(exposure[i,2:3]))
}

nrep=100000
state=matrix(0,pol_size*num_mapping,nrep)
for (j in 1:nrep){
  realized_assignment=rbinom(pol_size,1,assign_prob)
  exposure = matrix(0,nrow=pol_size,ncol=3)
  exposure[,1]=1:pol_size #vertex indices
  exposure[,2]=realized_assignment
  
  for (i in 1:pol_size){
    friend_index=get.neighborhood(toy_net,i,'out')
    exposure[i,3]= any((realized_assignment[friend_index])>0)
    
  }
  
  state_temp=c()
  for (i in 1:pol_size){
    
    state_temp=c(state_temp,individual_exposure(exposure[i,2:7]))
  }
  
  state[,j]=state_temp
}
prob=apply(state,1,mean)

var_cov= cov(t(state))

second_order = diag(1/prob) %*% var_cov %*% diag(1/prob)
