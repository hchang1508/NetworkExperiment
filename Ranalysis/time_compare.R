

memberships=component.dist(net,connected = 'weak')$membership
trial_subgraph=get.inducedSubgraph(net,which(memberships==1))
#change to the R vector
nrep=c(10,50,100,500,1000,5000,10000)
nrep=50
pol_size = network.size(trial_subgraph)
assign_prob = c(0.25,0.25,0.25,0.25)
num_mapping=12

state=matrix(0,pol_size*num_mapping,nrep)
mean_pre=rep(0,pol_size*num_mapping)
cov_pre=matrix(0,pol_size*num_mapping,pol_size*num_mapping)
realized_assignments=matrix(0,nrep,pol_size)
for (j in 1:nrep){
  realized_assignments[j,]=sample(0:3,pol_size,replace=TRUE,prob=assign_prob)
}

time_start1=proc.time()

for (j in 1:nrep){
  if (j%%100==0){
    print(j)
  }
  realized_assignment=realized_assignments[j,]
  
  exposure = matrix(0,nrow=pol_size,ncol=7)
  exposure[,1]=1:pol_size #vertex indices
  exposure[,2]=realized_assignment==0
  exposure[,3]=realized_assignment==1
  exposure[,4]=realized_assignment==2
  exposure[,5]=realized_assignment==3
  
  for (i in 1:pol_size){
    friend_index=get.neighborhood(trial_subgraph,i,'out')
    exposure[i,6]= sum(realized_assignment[friend_index]==0)
    exposure[i,7]= sum(realized_assignment[friend_index]==1)    
  }
  
  state_temp=c()  
  for (i in 1:pol_size){
    state_temp=c(state_temp,individual_exposure(exposure[i,2:7]))
  }
  state[,j]=state_temp
}

prob=apply(state,1,mean)
var_cov= cov(t(state))
time_end1=proc.time()


time_start2=proc.time()
for (j in 1:nrep){
  if (j%%100==0){
    print(j)
  }
  realized_assignment=realized_assignments[j,]
  exposure = matrix(0,nrow=pol_size,ncol=7)
  exposure[,1]=1:pol_size #vertex indices
  exposure[,2]=realized_assignment==0
  exposure[,3]=realized_assignment==1
  exposure[,4]=realized_assignment==2
  exposure[,5]=realized_assignment==3
  
  for (i in 1:pol_size){
    friend_index=get.neighborhood(trial_subgraph,i,'out')
    exposure[i,6]= sum(realized_assignment[friend_index]==0)
    exposure[i,7]= sum(realized_assignment[friend_index]==1)    
  }
  
  state_temp=c()  
  for (i in 1:pol_size){
    state_temp=c(state_temp,individual_exposure(exposure[i,2:7]))
  }
  result = Welford_online(mean_pre,cov_pre,j,state_temp)
  mean_pre=result[[1]]
  cov_pre=result[[2]]
}
cov_pre=cov_pre/nrep
time_end2=proc.time()
print(time_end2-time_start2)
print(time_end1-time_start1)