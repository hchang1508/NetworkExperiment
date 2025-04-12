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

#######################################################################
##################PREPARATIONS#########################################
#######################################################################
#First select a small cluster
components=component.dist(net,connected='weak')
memberships=components$membership
csize = components$csize

comp_5=get.inducedSubgraph(net,which(memberships==8))
network.size(comp_5)  #332 units
id_5 = comp_5 %v% 'vertex.names'
sum(id_5 %in% id_list) #1777 experiment subjects
degree(comp_5,cmode='outdegree') #warning some units have only two friends so some exposure mappings are not defined


#First Order and Second Order Matrics
prob=c(0.25,0.25,0.25,0.25)
expo1=1
expo2=2
foso=FOSO(net=comp_5, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=100000,id_list,expo1,expo2)


#potential outcomes and sharp null; extract pretreat covar data
Y_5=Y[match(id_5,Y$id),]
X_5=X[match(id_5,X$id),]
x=as.matrix(X_5)

#component statistics
pol_size = network.size(comp_5)
status = 0:3-1
num_mappings = 12
assign_prob=prob
network_subjects = comp_5 %v% 'vertex.names'
subjects = (1:pol_size)[network_subjects %in% id_list]

subject_size = length(subjects)

#potential outcomes: sharp null
Y_5=Y_5[match(network_subjects,Y_5$id),]
Y_5=Y_5[subjects,]
x=x[subjects,]
#pre-treatment covarates and demeaning
x=x[,c(3,4)] #use only one covariate
for (i in 1:ncol(x)){
  x[,i]=x[,i]-mean(x[,i])
}


########################################################################
####################Simulations$$$$$###################################
########################################################################
nsim=100000
fo= FO(net=comp_5, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list=id_list,expo1=expo1,expo2=expo2)
pi1 = fo[(1:length(fo))%%2==1]
pi2 = fo[(1:length(fo))%%2==0]
pol_size = network.size(comp_5)
status = 0:3-1
num_mappings = 12
assign_prob=prob
network_subjects = comp_5 %v% 'vertex.names'
subjects = (1:pol_size)[network_subjects %in% id_list]
subject_size = length(subjects)


#prob
pi1 = pi1[subjects]
pi2 = pi2[subjects]

est_sim = matrix(NA,nsim,15)
#est_sim=rep(0,nsim)
sample_size = matrix(0,nrow=2,ncol = nsim)
friends=matrix(FALSE,nrow=pol_size,ncol=5)
for (i in 1:pol_size){
  friend_index_t=get.neighborhood(comp_5,i,'out')
  if (length(friend_index_t)==0){
    next 
  }else
    friends[i,1:length(friend_index_t)] = friend_index_t
}

for (i in 1:nsim){
  
  print(i)
  realized_assignment_t=sample(status,pol_size,replace=TRUE,prob=assign_prob)
  exposure_t = matrix(0,nrow=pol_size,ncol=7)
  exposure_t[,1]=1:pol_size #vertex indices
  exposure_t[,2]=realized_assignment_t==0
  exposure_t[,3]=realized_assignment_t==1
  exposure_t[,4]=realized_assignment_t==2
  exposure_t[,5]=realized_assignment_t==3
  
  
  num_friend_treated_t = apply(friends,1,friend_treated,realized_assignment_t)
  exposure_t[,6]= num_friend_treated_t[1,]
  exposure_t[,7]=num_friend_treated_t[2,]  
  
  
  exposure_t = exposure_t[subjects,2:7]
  
  #please double check every time 
  state_t=apply(exposure_t,1,individual_exposure) 
  y1_t = state_t[expo1,] * Y_5[,2]
  y2_t = state_t[expo2,] * Y_5[,2]
  denom_pi1_t= state_t[expo1,]* 1/pi1 
  denom_pi2_t= state_t[expo2,]* 1/pi2 
  pi_weight_t = denom_pi1_t + denom_pi2_t
  d1_t=state_t[expo2,]
  d0_t=state_t[expo1,]
  
  y_obs_t=y1_t+y2_t
  sample_size[1,i]=sum(d1_t)
  sample_size[2,i]=sum(d0_t)  
  est_sim[i,]= Logit(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='weighted_joint')
  
  #est_sim[i,]= OLS(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint')
  #est_sim[i]=HT(y2,y1,pi2,pi1,subject_size)
}
mean(est_sim)
var(est_sim[,1])
#reg=lm(y~c1+c2+c1:x1+c2:x1-1)
#reg=glm(y~c1+c2+x1-1,family='quasibinomial',weights=pi12)
#reg=glm(y~c1+c2+c1:x1+c2:x1-1,family='quasibinomial')

#temp=reg$coefficients
#cbind(apply(est_sim,2,mean),c(0,temp))
#apply(est_sim,2,mean)-c(0,temp)

fo=foso[[1]]
so=foso[[2]]
normalized_cov = diag(1/fo) %*% so %*% diag(1/fo)
y=rep(Y_5[,2],each=2)
x_predict = x[rep(1:nrow(x),each=2),]
d1 = rep(c(0,1),length(y)/2)
d0 = 1-d1
contrast=c(-1,1)
contrast_OLS = c(contrast,rep(0,ncol(x)))
contrast_vec=rep(contrast,length(y)/2)
x_OLS_joint = data.frame(cbind(d0,d1,x_predict))
reg_OLS_joint =  lm(y~d1+d0+x_predict-1,weights=fo)
z_OLS_joint=  diag(y-predict(reg_OLS_joint,x_OLS_joint)) %*% diag(fo) %*%as.matrix(x_OLS_joint) %*% solve( t(as.matrix(x_OLS_joint)) %*%diag(fo)%*%as.matrix(x_OLS_joint))
z_OLS_joint = z_OLS_joint %*% contrast_OLS
var_OLS_joint = t(z_OLS_joint) %*% normalized_cov %*% z_OLS_joint


###Logit Weighted#####
data_t = cbind(y,d0,d1,x_predict)
data_t= data.frame(data_t)
x_Logit_joint_w = data.frame(cbind(d0,d1,x_predict))
reg_Logit_joint_w =  glm(y~d0+d1+x_predict-1,family='quasibinomial',data=data_t)
z_Logit_joint_w = y-predict(reg_Logit_joint_w,x_Logit_joint_w,type='response')
z_Logit_joint_w = z_Logit_joint_w * contrast_vec
var_Logit_joint_w = t(z_Logit_joint_w) %*% normalized_cov %*% z_Logit_joint_w 
var_Logit_joint_w /(subject_size)^2


###Logit Unweighted####
data_t = cbind(y,d0,d1,x_predict)
data_t= data.frame(data_t)
x_Logit_joint_u = data.frame(cbind(d0,d1,x_predict))
colnames(x_Logit_joint_u) = names(data_t)[2:length(data_t)]
reg_Logit_joint_u =  glm(y~d0+d1+x_predict-1,family='quasibinomial',weights=fo,data=data_t)
z_Logit_joint_u = y-predict(reg_Logit_joint_u,x_Logit_joint_u,type='response')
z_Logit_joint_u = z_Logit_joint_u * contrast_vec
var_Logit_joint_u = t(z_Logit_joint_u) %*% normalized_cov %*% z_Logit_joint_u 
var_Logit_joint_u /(sample_size)^2



