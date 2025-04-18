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

comp_5=get.inducedSubgraph(net,which(memberships==5))
network.size(comp_5)  #332 units
id_5 = comp_5 %v% 'vertex.names'
sum(id_5 %in% id_list) #1777 experiment subjects
degree_5=degree(comp_5,cmode='outdegree') #warning some units have only two friends so some exposure mappings are not defined


#First Order and Second Order Matrics
prob=c(0.25,0.25,0.25,0.25)
expo1=4
expo2=5

#depending on exposure mapping some units must be dropped for identification reasons:

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

prompt_t = paste0(number_dropped_t,' units dropped.')
print(prompt_t)

foso=FOSO(net=comp_5, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=100000,id_list,expo1,expo2)
fo=foso[[1]]
so=foso[[2]]
#potential outcomes and sharp null; extract pretreat covar data
Y_5=Y[match(id_5,Y$id),]
X_5=X[match(id_5,X$id),]
x=as.matrix(X_5)
contrast=c(-1,1)
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
y=Y_5[,2]

#computing  bound for estimation
normalized_cov=diag(1/fo) %*% so %*% diag(1/fo)
cov_bound = AS_bound(normalized_cov)
p = so + (fo) %*% t(fo) 
p=round(p,10)

p_zeros=(p==0)
cov_bound_zeros = (cov_bound==0)
min(cov_bound_zeros-p_zeros) #should expect 0

p_denomiator = p + cov_bound_zeros
cov_bound_HT= cov_bound / p_denomiator
########################################################################
####################Simulations$$$$$###################################
########################################################################
sim_value=simulationFeb13(comp_5,y,x,prob,id_list,expo1,expo2,cov_bound_HT)
est_sim=sim_value[[1]]
var_sim=sim_value[[2]]
simulated_value_var = apply(var_sim,2,mean)
simulated_value_est = apply(est_sim,2,var)

theoreitcal_value=theoretical_variance(y,x,fo,so,subject_size,contrast)
theoretical_variance_bound=theoretical_variance_bound(y,x,fo,so,subject_size,contrast)
#compare=cbind(theoreitcal_value,simulated_value)
compare=cbind(theoreitcal_value,theoretical_variance_bound)
#save(est_sim,compare,file='SimFeb11.Rdata')
