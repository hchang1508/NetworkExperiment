# This file contain functions used in simulations:
# 
# A. Welford_online: To calculate first order and second order assignment probabilities
# 
# B. individual_exposure: To calculate individual_exposure information
# 
# C. FOSO: To calculate first order and second order probability matrix
#
# D. Target Probability: Analytically compute first order probabilities
#


WELFORD_ONLINE=function(mean_pre,cov_pre,n,new_data){
  
  mean_new = (n-1)*mean_pre + new_data
  mean_new = mean_new / n
  diff= mean_pre-mean_new
  cov_new = (n-1)*cov_pre + (new_data-mean_pre) %*% t(new_data-mean_new)
  cov_new=cov_new/n
  
  return(list(mean_new,cov_new))
}

WELFORD_ONLINE_BATCH=function(mean_pre,cov_pre,n_old,n_new,mean_new,cov_new){
  
  n = n_old + n_new
  mean_new = n_old /n * mean_pre + n_new / n * mean_new

  diff=mean_new-mean_pre
  cov_new=n_old/n * cov_pre+ n_new/n * cov_new+ (diff)%*% t(diff) *(n_old/n)*(n_new/n)

  return(list(mean_new,cov_new))
}

Welford_online_mean=function(mean_pre,n,new_data){
  
  mean_new = (n-1)*mean_pre + new_data
  mean_new = mean_new / n
  return(mean_new)
}

INDIVIDUAL_EXPOSURE=function(exposure){
  
  output=rep(0,12)
  if (exposure[1]==1){
    output[1]=1  #first round simple
  }else if (exposure[2]==1){
    output[2]=1 #first round intense
  }else if (exposure[3]==1 & exposure[5]==0 & exposure[6]==0){
    output[3]=1 #second round simple, no friend first round 
  }else if (exposure[3]==1 & exposure[5]>=0 & exposure[6]==0){
    output[4]=1 #second round simple, friend from simple, no friend intensive
  }else if (exposure[3]==1 & exposure[6]==1){
    output[5]=1 #second round simple, one friend intensive
  }else if (exposure[3]==1 & exposure[6]==2){
    output[6]=1 #second round simple, two friend intensive
  }else if (exposure[3]==1 & exposure[6]>2){
    output[7]=1 #second round simple, two more friend intensive
  }else if (exposure[4]==1 & exposure[5]==0 & exposure[6]==0){
    output[8]=1 #second round intensive, no friend first round
  }else if (exposure[4]==1 & exposure[5]>=0 & exposure[6]==0){
    output[9]=1 #second round intensive, friend from simple, not friend intensive
  }else if (exposure[4]==1 & exposure[6]==1){
    output[10]=1 #second round intensive, one friendi ntensive
  }else if (exposure[4]==1 & exposure[6]==2){
    output[11]=1 #second round intensive, two friend intensive
  }else if (exposure[4]==1 & exposure[6]>2){
    output[12]=1 #second round intensive, two more friend intensive
  }
  
  return(output)
}

individual_exposure_factorial=function(exposure){
  
  output=rep(0,24)
  
  if (exposure[7]==0){ #default
  if (exposure[1]==1){
    output[1]=1  #first round simple
  }else if (exposure[2]==1){
    output[2]=1 #first round intense
  }else if (exposure[3]==1 & exposure[5]==0 & exposure[6]==0){
    output[3]=1 #second round simple, no friend first round 
  }else if (exposure[3]==1 & exposure[5]>=0 & exposure[6]==0){
    output[4]=1 #second round simple, friend from simple, no friend intensive
  }else if (exposure[3]==1 & exposure[6]==1){
    output[5]=1 #second round simple, one friend intensive
  }else if (exposure[3]==1 & exposure[6]==2){
    output[6]=1 #second round simple, two friend intensive
  }else if (exposure[3]==1 & exposure[6]>2){
    output[7]=1 #second round simple, two more friend intensive
  }else if (exposure[4]==1 & exposure[5]==0 & exposure[6]==0){
    output[8]=1 #second round intensive, no friend first round
  }else if (exposure[4]==1 & exposure[5]>=0 & exposure[6]==0){
    output[9]=1 #second round intensive, friend from simple, not friend intensive
  }else if (exposure[4]==1 & exposure[6]==1){
    output[10]=1 #second round intensive, one friendi ntensive
  }else if (exposure[4]==1 & exposure[6]==2){
    output[11]=1 #second round intensive, two friend intensive
  }else if (exposure[4]==1 & exposure[6]>2){
    output[12]=1 #second round intensive, two more friend intensive
  }
  }
  
  if (exposure[7]==1){
    if (exposure[1]==1){
      output[13]=1  #first round simple
    }else if (exposure[2]==1){
      output[14]=1 #first round intense
    }else if (exposure[3]==1 & exposure[5]==0 & exposure[6]==0){
      output[15]=1 #second round simple, no friend first round 
    }else if (exposure[3]==1 & exposure[5]>=0 & exposure[6]==0){
      output[16]=1 #second round simple, friend from simple, no friend intensive
    }else if (exposure[3]==1 & exposure[6]==1){
      output[17]=1 #second round simple, one friend intensive
    }else if (exposure[3]==1 & exposure[6]==2){
      output[18]=1 #second round simple, two friend intensive
    }else if (exposure[3]==1 & exposure[6]>2){
      output[19]=1 #second round simple, two more friend intensive
    }else if (exposure[4]==1 & exposure[5]==0 & exposure[6]==0){
      output[20]=1 #second round intensive, no friend first round
    }else if (exposure[4]==1 & exposure[5]>=0 & exposure[6]==0){
      output[21]=1 #second round intensive, friend from simple, not friend intensive
    }else if (exposure[4]==1 & exposure[6]==1){
      output[22]=1 #second round intensive, one friendi ntensive
    }else if (exposure[4]==1 & exposure[6]==2){
      output[23]=1 #second round intensive, two friend intensive
    }else if (exposure[4]==1 & exposure[6]>2){
      output[24]=1 #second round intensive, two more friend intensive
    }
  } 
  return(output)
}
FOSO_nppl = function(net, num_status,prob,num_mappings,nrep,id_list,expo1,expo2,option='Bernoulli'){
  
  #This function depends on the following functions:
  # 1. WELFORD_ONLINE: for online estimation of mean vector and covariance matrix.
  # 2. ASSIGNMENT: treatment assignment for the village-level stratified design.
  # 3. FRIEND_TREATED: compute number of friends treated.
  # 4. INDIVIDUAL_EXPOSURE: caculate individuals exposure mapping.
  # 5. 
  
  
  
  # This function calculate First Order and Second Order Assignment Probabilities 
  # for a network.
  # 
  # Input:
  #     1. net: a network object. 
  #     2. FUN: user-defined function to calculate exposure mappings
  #     3. num_status: number of treatment status
  #        --> Please check the detail of the code to make sure the output of the 
  #        num_treatments matches with the input of the FUN 
  #     4. probability of each treatment status. Currently we are assuing a Bernoulli
  #     5. number of mappings: number of distinct exposure mappings
  #     6. nrep: number of reptitions
  #   7&8. expo1, expo2: a pair of exposures that is of interest
  
  
  #calculate network size
  pol_size = network.size(net)
  status = 1:num_status 
  num_mappings = num_mappings
  assign_prob=prob
  pairs = 2 # a pair of exposures
  #track progress
  threshold = nrep/10
  
  
  #define experient targets
  network_subjects = net %v% 'vertex.names'
  subjects = (1:pol_size)[network_subjects %in% id_list]
  
  subject_size = length(subjects)
  #output 
  result=matrix(0,2,1)  
  #extract directed friend nomination
  friends=matrix(FALSE,nrow=pol_size,ncol=5)
  for (i in 1:pol_size){
    friend_index=get.neighborhood(net,i,'out')
    if (length(friend_index)==0){
      next 
    }else
      friends[i,1:length(friend_index)] = friend_index
  }
  
  for (j in 1:nrep){
    if (j%%threshold==0){
      prompt_sim=paste0(round(j/threshold),"% completed")
      print(prompt_sim)
    }
    if (option=='Bernoulli'){
      realized_assignment=sample(status,pol_size,replace=TRUE,prob=assign_prob)
    }else if (option=='Stratified'){
      
      temp=cbind(net %v% "group",net %v% "assignment")
      realized_assignment=ASSIGNMENT(temp)
    }
    
    #calculate exposures
    exposure = matrix(0,nrow=pol_size,ncol=7)
    exposure[,1]=1:pol_size #vertex indices
    exposure[,2]=realized_assignment==1
    exposure[,3]=realized_assignment==2
    exposure[,4]=realized_assignment==3
    exposure[,5]=realized_assignment==4
    num_friend_treated = apply(friends,1,FRIEND_TREATED,realized_assignment)
    exposure[,6]=num_friend_treated[1,]
    exposure[,7]=num_friend_treated[2,]  
    
    
    exposure_temp = exposure[subjects,2:7]
    
    #please double check every time 
    
    if (is.null(dim(exposure_temp))){
      state_temp=individual_exposure(exposure_temp)
      state_temp = state_temp[c(expo1,expo2)]
    }else{
      state_temp=apply(exposure_temp,1,INDIVIDUAL_EXPOSURE) 
      state_temp = state_temp[c(expo1,expo2),]      
    }
    
    nppl=apply(state_temp,1,sum)
    #please double check every time 

    #store results and update FOSO

    result = (result * (j-1)+ nppl )/j
    
  }
  print(result)
  return(result)
}

FOSO = function(net, num_status,prob,num_mappings,nrep,id_list,expo1,expo2,option='Bernoulli'){
  

  # This function calculate First Order and Second Order Assignment Probabilities 
  # for a network.
  # 
  # Input:
  #     1. net: a network object. 
  #     2. FUN: user-defined function to calculate exposure mappings
  #     3. num_status: number of treatment status
  #        --> Please check the detail of the code to make sure the output of the 
  #        num_treatments matches with the input of the FUN 
  #     4. probability of each treatment status. Currently we are assuing a Bernoulli
  #     5. number of mappings: number of distinct exposure mappings
  #     6. nrep: number of reptitions
  #   7&8. expo1, expo2: a pair of exposures that is of interest
  
  
  #calculate network size
  pol_size = network.size(net)
  status = 1:num_status #number of treatment status

  num_mappings = num_mappings
  
  assign_prob=prob
  
  pairs = 2 # a pair of exposures
  
  #track progress
  threshold = nrep/10
  
  
  #define experient targets
  network_subjects = net %v% 'vertex.names'
  subjects = (1:pol_size)[network_subjects %in% id_list]
  subject_size = length(subjects)
  
  #output 
  mean_pre=rep(0,subject_size*pairs )
  cov_pre=matrix(0,subject_size*pairs ,subject_size*pairs )  
  
  #extract directed friend nomination
  friends=matrix(FALSE,nrow=pol_size,ncol=5)
  for (i in 1:pol_size){
    friend_index=get.neighborhood(net,i,'out')
    if (length(friend_index)==0){
      next 
    }else
      friends[i,1:length(friend_index)] = friend_index
  }
  
  for (j in 1:nrep){

    if (j%%threshold==0){
      prompt_sim=paste0(round(j/threshold),"% completed")
      print(prompt_sim)
    }
    
    if (option=='Bernoulli'){
      realized_assignment=sample(status,pol_size,replace=TRUE,prob=assign_prob)
    }else if (option=='Stratified'){
      
      temp=cbind(net %v% "group",net %v% "assignment")
      realized_assignment=ASSIGNMENT(temp)
    }

    ##########################################
    ###Calculate  Exposures###################
    ##########################################

    exposure = matrix(0,nrow=pol_size,ncol=7)
    exposure[,1]=1:pol_size #vertex indices
    exposure[,2]=realized_assignment==1 #first round simple
    exposure[,3]=realized_assignment==2 #first round intense
    exposure[,4]=realized_assignment==3 #second round simple
    exposure[,5]=realized_assignment==4 #second round intense

    num_friend_treated = apply(friends,1,FRIEND_TREATED,realized_assignment)
    exposure[,6]=num_friend_treated[1,] #number of friends treated in first round
    exposure[,7]=num_friend_treated[2,] #number of friends treated in second round
    
    
    exposure_temp = exposure[subjects,2:7]
    
    #please double check every time 
    
    if (is.null(dim(exposure_temp))){
      state_temp=individual_exposure(exposure_temp)
      state_temp = state_temp[c(expo1,expo2)]
    }else{
      state_temp=apply(exposure_temp,1,INDIVIDUAL_EXPOSURE) 
      state_temp = state_temp[c(expo1,expo2),]      
    }
    
    #please double check every time 
    state_temp = as.vector(state_temp) 
    
    #store results and update FOSO
    result = WELFORD_ONLINE(mean_pre,cov_pre,j,state_temp)
    
    mean_pre=result[[1]]
    cov_pre=result[[2]]

  }

  return(result)
  
}
FOSO_sanity_check = function(net, num_status,prob,num_mappings,nrep,id_list,expo1,expo2,option='Bernoulli'){
  
  
  # This function returns the list (order) of subjects 
  
  #calculate network size
  pol_size = network.size(net)

  #define experient targets
  network_subjects = net %v% 'vertex.names'
  subjects = (1:pol_size)[network_subjects %in% id_list]
  subject_size = length(subjects)
  
  #output the list of subjects in the network
  output_list = network_subjects[subjects]
  
  return(output_list)
  
}

#true treatment probability under Bernoulli asignment  
target_prob = function(prob,num){
  
  second_round_prob = prob[3] + prob[4]
  output=rep(0,12)
  output[1]=prob[1]
  output[2]=prob[2]
  output[3]=prob[3] * dbinom(num,num,second_round_prob)
  output[4]=prob[3] * (1-pbinom(0,num,prob[1])) * dbinom(0,num,prob[2])
  output[5]=prob[3] * dbinom(1,num,prob[2])
  output[6]=prob[3] * dbinom(2,num,prob[2])
  output[7]=prob[3] * (1-pbinom(2,num,prob[2]))
  output[8]=prob[4] * dbinom(num,num,second_round_prob)
  output[9]=prob[4] * (1-pbinom(0,num,prob[1])) * dbinom(0,num,prob[2])
  output[10]=prob[4] * dbinom(1,num,prob[2])
  output[11]=prob[4] * dbinom(2,num,prob[2])
  output[12]=prob[4] * (1-pbinom(2,num,prob[2]))
  
  return(output)
}

FO = function(net, FUN, num_status,prob,num_mappings,nrep,id_list,expo1,expo2){
  
  
  # This function calculate First Order and Second Order Assignment Probabilities 
  # for a network.
  # 
  # Input:
  #     1. net: a network object. 
  #     2. FUN: user-defined function to calculate exposure mappings
  #     3. num_status: number of treatment status
  #        --> Please check the detail of the code to make sure the output of the 
  #        num_treatments matches with the input of the FUN 
  #     4. probability of each treatment status. Currently we are assuing a Bernoulli
  #     5. number of mappings: number of distinct exposure mappings
  #     6. nrep: number of reptitions
  #   
  
  
  #calculate network size
  pol_size = network.size(net)
  status = 1:num_status - 1
  num_mappings = num_mappings
  assign_prob=prob
  pairs = 2 # a pair of exposures
  #track progress
  threshold = round(nrep/100,0)
  
  
  #define experient targets
  network_subjects = net %v% 'vertex.names'
  subjects = (1:pol_size)[network_subjects %in% id_list]
  
  subject_size = length(subjects)
  #output 
  mean_pre=rep(0,subject_size*pairs )
  
  #neigh
  friends=matrix(FALSE,nrow=pol_size,ncol=5)
  for (i in 1:pol_size){
    friend_index=get.neighborhood(net,i,'out')
    if (length(friend_index)==0){
      next 
    }else
      friends[i,1:length(friend_index)] = friend_index
    }

  
  for (j in 1:nrep){
    if (j%%threshold==0){
      prompt_sim=paste0(round(j/threshold),"% completed")
      print(prompt_sim)
    }
    
    realized_assignment=sample(status,pol_size,replace=TRUE,prob=assign_prob)

    #calculate exposures
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
    state_temp = state_temp[c(expo1,expo2),]
    #please double check every time 
    state_temp = as.vector(state_temp) 
 
    #store results and update FOSO

    result = Welford_online_mean(mean_pre,j,state_temp)
    time.end=proc.time()

    mean_pre=result

  }
  
  return(result)
}

FRIEND_TREATED = function(friend_index,realized_assignment){
  
  
  num_friend_treated =c( sum(realized_assignment[friend_index]==1),sum(realized_assignment[friend_index]==2) )

  return(num_friend_treated)
}

ASSIGNMENT=function(temp){
  
  
  #column names
  colnames(temp)=c('group','assignment')

  #list of group_index
  group_index=unique(temp[,'group'])
  
  #within each group, completely randomized
  for (i in group_index){
    
    #units in this group
    ind=temp[,'group']==i
    
    #number of units in this group
    nobs=sum(ind)
    
    temp[ind,'assignment']=sample(temp[ind,'assignment'],nobs)
    
  }

  return(temp[,'assignment'])
  
}
