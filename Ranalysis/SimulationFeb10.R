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

comp_5=get.inducedSubgraph(net,which(memberships==5))

prob=c(0.25,0.25,0.25,0.25)
foso=FOSO(net=comp_5, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list,4,5)


  
comp_5=get.inducedSubgraph(net,which(memberships==5))
network.size(comp_5)  #332 units
id_5 = comp_5 %v% 'vertex.names'
sum(id_5 %in% id_list) #1777 experiment subjects
degree(comp_5,cmode='outdegree') #warning some units have only two friends so some exposure mappings are not defined


#potential outcomes and sharp null
Y_5=Y[match(id_5,Y$id),]
X_5=X[match(id_5,X$id),]
x=as.matrix(X_5)
x=x[,c(-1,-14)]
for (i in 1:ncol(x)){
  
  x[,i]=x[,i]-mean(x[,i])
}

#prob

pol_size = network.size(comp_5)
status = 0:3-1
num_mappings = 12
assign_prob=prob
network_subjects = comp_5 %v% 'vertex.names'
subjects = (1:pol_size)[network_subjects %in% id_list]
expo1=4
expo2=5
subject_size = length(subjects)

#potential outcomes: sharp null
Y_5=Y_5[match(network_subjects,Y_5$id),]
Y_5=Y_5[subjects,]
x=x[subjects,]
inverse_pi = 1/foso[[1]]
d =foso[[2]] 
y=rep(Y_5[,2],each=2)
y=y*rep(c(-1,1),316)
var_HT = (y*inverse_pi) %*% d %*% (y*inverse_pi)/(subject_size)^2

nsim=10000
pi45= FO(net=comp_5, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list=id_list,expo1=4,expo2=5)
pi1 = pi45[(1:length(pi45))%%2==1]
pi2 = pi45[(1:length(pi45))%%2==0]
pol_size = network.size(comp_5)
status = 0:3-1
num_mappings = 12
assign_prob=prob
network_subjects = comp_5 %v% 'vertex.names'
subjects = (1:pol_size)[network_subjects %in% id_list]
expo1=4
expo2=5
subject_size = length(subjects)


#prob
pi1 = pi1[subjects]
pi2 = pi2[subjects]

est_sim = rep(NA,nsim)
sample_size = matrix(0,nrow=2,ncol = nsim)
friends=matrix(FALSE,nrow=pol_size,ncol=5)
for (i in 1:pol_size){
  friend_index=get.neighborhood(comp_5,i,'out')
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
  y1 = state_temp[expo1,] * Y_5[,2]
  y2 = state_temp[expo2,] * Y_5[,2]
  denom_pi1= state_temp[expo1,]* 1/pi1 
  denom_pi2= state_temp[expo2,]* 1/pi2 
  pi_weight = denom_pi1 + denom_pi2
  d1=state_temp[expo2,]
  d0=state_temp[expo1,]
  
  sample_size[1,i]=sum(d1)
  sample_size[2,i]=sum(d0)  
  exposure_dict[,i]=c(state_temp[expo1,],state_temp[expo2,])
  est_sim[i]= HT(y2,y1,pi2,pi1,subject_size)
  
 
}
mean(est_sim)

y=Y_5[,2]
est_sim=simulationFeb9_12(comp_5,y,x,prob,id_list)
dim(x)
