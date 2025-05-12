##header
source('0_header.R')

#index for parallelizing
input= commandArgs(trailingOnly=TRUE)
index=as.numeric(input[1])
expo1=as.numeric(input[2])
expo2=as.numeric(input[3])
prob_cases=as.numeric(input[4])

contrast=c(-1,1)

#set random seed
set.seed(index)

case_to_do=read.csv(paste0('/home/hc654/NetworkExperiment/2_simulation_Bernoulli_Ds/Sim_',prob_cases,'_',expo1,'_',expo2,'.csv'))[,2]
index=case_to_do[index]
#simulation parameters
nsim=1

time_start=Sys.time()


#initial membership informations
#drop cluster 7, 12 and 29, exposure mapping not defined
#comp_list=1:31
#comp_list=comp_list[-c(7,12,29)] #no exposure mapping for these two  
#comp_list=comp_list[-7]
#net_t = get.inducedSubgraph(net,which(memberships %in% comp_list))
if (prob_cases==1){
  prob=rep(0.25,4)
}else if ( prob_cases ==2){
  prob=c(1/5,1/5,3/10,3/10)
}else if ( prob_cases ==3){
  prob=c(1/6,1/6,4/12,4/12)
}else if (prob_cases==4){
  prob=c(1/7,1/7,5/14,5/14)
}else if (prob_cases==5){
  prob=c(1/8,1/8,6/16,6/16)
}else if (prob_cases==6){
  prob=c(1/9,2/9,4/12,4/12)
}


#Load treatment probability and FO design matrices
foso_file=paste0('/home/hc654/palmer_scratch/final_analysis_Ds/','Stratified_nat_',expo1,expo2,'_D.Rdata')
load(foso_file)
fo=mean_pre
so=cov_pre


#Load variance bound matrices 
so_bound_file=paste0('/home/hc654/palmer_scratch/final_analysis_Ds/','Stratified_nat_',expo1,expo2,'_Dbound.Rdata')
load(so_bound_file)
so_AS2=so_bound


#unfold fo for later calculations
fo_vec=unlist(fo)


###basic simulation information
#network size
pol_size = network.size(net)

#treatment status vectors/num of mappings
status = 1:4
num_mappings = 12

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

###read from some input files
pi0 = fo_vec[(1:length(fo_vec))%%2==1]
pi1 = fo_vec[(1:length(fo_vec))%%2==0]

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
dict_group[,2]= c(0,group_index[(1:(length(group_index))-1)])+1
dict_group[,3]= group_index
dict_group[,4]=c(0,2*group_index[(1:(length(group_index))-1)])+1
dict_group[,5]= group_index*2
dict_group<- matrix(as.numeric(dict_group),   ncol = ncol(dict_group))

###########################################################################################
#SANITY CHECK: MAKE SURE THE ORDERING FOR FO AND SO MATRICES ARE CONSISTENT################
###########################################################################################

print('Checking inconsisntey in id orderings of foso probabilities (no error msg is good)')


compare1=compute_FOSO_all_components_sanity_check(net,expo1,expo2 ,id_list,option='Bernoulli')
for (i in 1:nrow(dict_group)){
  
  
  #network index
  m = dict_group[i,1]
  
  #for those subjects who are in group m, what are their id?
  compare2_temp = id_t[subjects_t[dict_group[i,2]:dict_group[i,3]]]
  
  compare1_temp = compare1[[m]]  
  
  temp = (compare1_temp)==(compare2_temp)
  #if there is any inequality !temp will contain a true and this prompts an error msg
  if (any(!temp)==TRUE ){
    print(paste0('Inconsistency in network component m:', m))
  }
  
}

#################################################################
####SORT OUTCOMES################################################
#################################################################
#the outcomes are sorted according to the orderings of subjects_t

Y=matrix(0,0,3)
for (i in 1:nrow(dict_group)){
  
  #network index
  m = dict_group[i,1]
  
  #for those subjects who are in group m, what are their id?
  id_temp = id_t[subjects_t[dict_group[i,2]:dict_group[i,3]]]
  
  
  Y_temp = Y_impute[match(id_temp,Y_impute[,1]),c(1,expo1+1,expo2+1)]
  
  Y = rbind(Y,Y_temp)
  
}


#Y=Y[subjects_t,]
#Y=Y_all[,c(1,2,3)]
#weights_random = runif(subject_size)
y0=Y[,2] #control
y1=Y[,3] #treated
#y0=y0_opt * (weights_random <= weight) + y0_imp * (weights_random > weight) 
#y1=y1_opt * (weights_random <= weight) + y1_imp * (weights_random > weight) 

############################################################################
####SANITY CHECK AGAIN: MAKE THE OUTCOMES AND FOSO ARE ALIGNED##############
############################################################################

print('Checking nnconsisntey in id orderings of outcomes (no error msg is good)')

for (i in 1:nrow(dict_group)){
  
  
  #network index
  m = dict_group[i,1]
  
  #for those subjects who are in group m, what are their id?
  compare3_temp =Y[(dict_group[i,2]:dict_group[i,3]),1]
  
  
  compare1_temp = compare1[[m]]  
  temp = (compare1_temp)==(compare3_temp)
  
  #if there is any inequality !temp will contain a true and this prompts an error msg
  if (any(!temp)==TRUE ){
    
    print(paste0('Inconsistency in network component m:', m))
  }
  
}



#################################################################
#######COVARIATES################################################
#################################################################

X=all_info2[match(id_t,all_info2$id),]
X=X[,c('id','male','age','agpop','literacy','ricearea_2010','risk_averse','disaster_prob')] 
x=c()
for (i in 1:nrow(dict_group)){
  
  #network index
  m = dict_group[i,1]
  
  #for those subjects who are in group m, what are their id?
  id_temp = id_t[subjects_t[dict_group[i,2]:dict_group[i,3]]]
  
  
  X_temp = X[match(id_temp,X[,1]),]
  
  x = rbind(x,X_temp)
  
}


for (i in 2:ncol(x)){
  x[,i]=(x[,i]-mean(x[,i]))/sd(x[,i])
}
x=as.matrix(x)

############################################################################
####SANITY CHECK AGAIN: MAKE THE OUTCOMES AND COVARATES ALIGNED##############
############################################################################

print('Checking nnconsisntey in id orderings of x (no error msg is good)')
temp = (Y[,1]==x[,1])

if (any(!temp)==TRUE ){
  
  print('Inconsistency in between outcomes and covariates')
}

x=x[,2:ncol(x)]

#################################################################
#######FRIENDSHIP################################################
#################################################################

#friendship information
friends=matrix(FALSE,nrow=pol_size,ncol=5)
for (i in 1:pol_size){
  friend_index_t=get.neighborhood(net,i,'out')
  if (length(friend_index_t)==0){
    next 
  }else
    friends[i,1:length(friend_index_t)] = friend_index_t
}


#################################################################
#######REMOVE UNITS##############################################
#################################################################


temp = (1:subject_size) [which( (pi0<1e-4)| (pi1<1e-4)) ] #some units who are not assigned to the exposures


#y1
#y0
#fo
#so
#





#################################################################
#######SIMULATION################################################
#################################################################


for (i in 1:nsim){
  
  
  result=SIM_ONERUN_IMPUTED_AUG20_AS2(status,pol_size,prob,y1,y0,pi1,pi0,index,expo1,expo2,option='Bernoulli',x,subjects_t=subjects_t,so_AS2=so_AS2)
  
  est_sim=result[[1]]
  var_sim=result[[2]]
  var_sim2=result[[3]]
    
}


time_end=Sys.time()

print(time_start-time_end)
output=list(est_sim,var_sim,var_sim2)

print(output)

file_name=paste0('/home/hc654/palmer_scratch/final_analysis_Ds/simulation_output/',expo1,'_',expo2,'_Sim_',index,'.Rdata')
#file_name2=paste0('/home/hc654/palmer_scratch/final_analysis_simulation/',weight,'_',expo1,expo2,'_SimSampleSize_',index,'.Rdata')
save(output,file=file_name)

