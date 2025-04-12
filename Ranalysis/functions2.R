# This file contain functions used in simulations:
# 
# A. HT Estimator
# 
# B. HA Estimaator
# 



HT = function(epsilon1,epsilon0,pi1,pi0,sample_size){
  
  est = sum(epsilon1 * (1/pi1)) - sum(epsilon0 * (1/pi0))
  
  est = est/sample_size
  
  return(est)
  
  

}

HA =function(epsilon1,epsilon0,pi1,pi0,denom_pi1,denom_pi0){
  
  est = sum(epsilon1 * (1/pi1))/sum(denom_pi1) - sum(epsilon0 * (1/pi0))/sum(denom_pi0)
  
  return(est)
  
}

OLS=function(y,d1,d0,pi,x,mode){
  

  if(mode=='OLS_joint'){
    weights = as.numeric((pi>0))
    reg=lm(y~d1+d0-1+x,weights=weights)
    
    avg_effect1=reg$coefficients['d1']
    avg_effect0=reg$coefficients['d0']
    coef=reg$coefficients
  }else if (mode=='OLS_separate'){
    
    weights = as.numeric((pi>0))
    reg=lm(y~d1+d0-1+d1:x+d0:x,weights=weights) 
    avg_effect1=reg$coefficients['d1']
    avg_effect0=reg$coefficients['d0']
    coef=reg$coefficients
  }else if (mode=='WLS_joint'){
    weights = pi
    reg=lm(y~d1+d0-1+x,weights=weights)
    avg_effect1=reg$coefficients['d1']
    avg_effect0=reg$coefficients['d0']
    coef=reg$coefficients
  }else if (mode=='WLS_separate'){
    weights = pi
    reg=lm(y~d1+d0-1+d1:x+d0:x,weights=weights) 
    avg_effect1=reg$coefficients['d1']
    avg_effect0=reg$coefficients['d0']
    coef=reg$coefficients
  }
  
  return(avg_effect1-avg_effect0)
  #return(c(avg_effect1-avg_effect0,coef))
  
}

Logit=function(y,d1,d0,weights,x,pi1,pi0,mode){
  

  if (mode=='unweighted_joint'){

    data=cbind(y,d0,d1,x)
    data=data.frame(data)
    x1=cbind(rep(0,nrow(x)),rep(1,nrow(x)),x)
    x0=cbind(rep(1,nrow(x)),rep(0,nrow(x)),x)
    x1=data.frame(x1)
    x0=data.frame(x0)
    colnames(x1)=colnames(data)[2:length(data)]
    colnames(x0)=colnames(data)[2:length(data)]
    
    weights = as.numeric((weights>0))
    reg = glm(y~d0+d1+x-1,family='quasibinomial',weights=weights,data=data)
    coef=reg$coefficients
    y_p1 = predict(reg,x1,type='response')
    y_p0 = predict(reg,x0,type='response')
    
    avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
    avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)
  }else if (mode=='unweighted_separate'){
    
    data1=data.frame(cbind(y,x))
    data0=data.frame(cbind(y,x))
    x_predict=data.frame(x)
    colnames(x_predict)=colnames(data1)[2:length(data1)]
    reg1 = glm(y~x,family='quasibinomial',weights=d1,data=data1)
    reg0 = glm(y~x,family='quasibinomial',weights=d0,data=data0)
    
    y_p1 = predict(reg1,x_predict,type='response')
    y_p0 = predict(reg0,x_predict,type='response')    
    coef=c(reg1$coefficients,reg0$coefficients)
    avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
    avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)
  }else if (mode=='weighted_joint'){
    
    data=cbind(y,d0,d1,x)
    data=data.frame(data)
    x1=cbind(rep(0,nrow(x)),rep(1,nrow(x)),x)
    x0=cbind(rep(1,nrow(x)),rep(0,nrow(x)),x)
    x1=data.frame(x1)
    x0=data.frame(x0)
    colnames(x1)=colnames(data)[2:length(data)]
    colnames(x0)=colnames(data)[2:length(data)]

    weights = weights
    reg = glm(y~d0+d1+x-1,family='quasibinomial',weights=weights)
    
    y_p1 = predict(reg,x1,type='response')
    y_p0 = predict(reg,x0,type='response')
    coef=reg$coefficients 
    avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
    avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)
  }else if (mode=='weighted_separate'){
    data1=data.frame(cbind(y,x))
    data0=data.frame(cbind(y,x))
    x_predict=data.frame(x)
    colnames(x_predict)=colnames(data1)[2:length(data1)]
    reg1 = glm(y~x,family='quasibinomial',weights=d1*weights,data=data1)
    reg0 = glm(y~x,family='quasibinomial',weights=d0*weights,data=data0)
    
    y_p1 = predict(reg1,x_predict,type='response')
    y_p0 = predict(reg0,x_predict,type='response')    
    
    avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
    avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)
    coef=c(reg1$coefficients,reg0$coefficients)
  }
  return(avg_effect1-avg_effect0)
  #return(c(avg_effect1-avg_effect0,coef))
}


simulationFeb10=function(net,y,x,prob,id_list,expo1,expo2){
  
  
  pol_size = network.size(net)
  status = 0:3-1
  num_mappings = 12
  assign_prob=prob
  network_subjects = net %v% 'vertex.names'
  subjects = (1:pol_size)[network_subjects %in% id_list]

  subject_size = length(subjects)
  

  nsim=100000
  expo1=expo1
  expo2=expo2
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
  
  #method=c('HT','HA','OLS_joint','OLS_separate','WLS_joint','WLS_separate','unweighted_joint',
  #         'unweighted_sepaprate','weighted_joint','weighted_separate')
  method=c('HT','HA','OLS_joint','WLS_joint','unweighted_joint', 'weighted_joint')
  est_sim=matrix(0,nrow=nsim,ncol=length(method)) 
  colnames(est_sim)=method
  
  #est_sim=rep(0,nsim)
  sample_size = matrix(0,nrow=2,ncol = nsim)
  friends=matrix(FALSE,nrow=pol_size,ncol=5)
  for (i in 1:pol_size){
    friend_index_t=get.neighborhood(net,i,'out')
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
    exposure_t[,6]=num_friend_treated_t[1,]
    exposure_t[,7]=num_friend_treated_t[2,]  
    
    
    exposure_t = exposure_t[subjects,2:7]
    
    #please double check every time 
    state_t=apply(exposure_t,1,individual_exposure) 
    y1_t = state_t[expo1,] * y
    y2_t = state_t[expo2,] * y
    denom_pi1_t= state_t[expo1,]* 1/pi1 
    denom_pi2_t= state_t[expo2,]* 1/pi2 
    pi_weight_t = denom_pi1_t + denom_pi2_t
    d1_t=state_t[expo2,]
    d0_t=state_t[expo1,]
    
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
   
    est_sim[i,1]= HT(y2_t,y1_t,pi2,pi1,subject_size)
    est_sim[i,2]= HA(y2_t,y1_t,pi2,pi1,denom_pi2_t,denom_pi1_t)
    est_sim[i,3]= OLS(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint')
    est_sim[i,4]= OLS(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint')    
    est_sim[i,5]= Logit(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='unweighted_joint')
    est_sim[i,6]= Logit(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='weighted_joint')

  }

  return(est_sim)
}

theoretical_variance=function(y,x,fo,so,sample_size,contrast){
  
  y=rep(y,each=2)
  inverse_prob=1/fo
  normalized_cov=diag(1/fo) %*% so %*% diag(1/fo)

  
  y_length = length(y)
  x_predict = x[rep(1:nrow(x),each=2),]
  d1 = rep(c(0,1),y_length/2)
  d0 = 1-d1
  
  contrast_OLS = c(contrast,rep(0,ncol(x)))
  contrast_vec=rep(contrast,length(y)/2)

  ####HT#####
  z_HT = y * contrast_vec
  var_HT = z_HT %*% normalized_cov %*% z_HT
  var_HT = var_HT/(sample_size)^2

  ####HA#####    
  x_HA = data.frame(cbind(d0,d1))
  reg_HA =  lm(y~d1+d0-1)
  z_HA= y-predict(reg_HA,x_HA)
  z_HA = contrast_vec * z_HA
  var_HA = z_HA %*% normalized_cov %*% z_HA
  var_HA = var_HA/(sample_size)^2
  
  ####OLS_joint####
  x_OLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_OLS_joint =  lm(y~d1+d0+x_predict-1,weights=fo)
  z_OLS_joint=  diag(y-predict(reg_OLS_joint,x_OLS_joint)) %*% diag(fo) %*%as.matrix(x_OLS_joint) %*% solve( t(as.matrix(x_OLS_joint)) %*%diag(fo)%*%as.matrix(x_OLS_joint))
  z_OLS_joint = z_OLS_joint %*% contrast_OLS
  var_OLS_joint = t(z_OLS_joint) %*% normalized_cov %*% z_OLS_joint
  
  ####WLS_joint####
  x_WLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_WLS_joint =  lm(y~d1+d0+x_predict-1)
  z_WLS_joint=  diag(y-predict(reg_WLS_joint,x_WLS_joint)) %*% as.matrix(x_WLS_joint) %*% solve( t(as.matrix(x_WLS_joint)) %*%as.matrix(x_WLS_joint))
  z_WLS_joint = z_WLS_joint %*% contrast_OLS
  var_WLS_joint = t(z_WLS_joint) %*% normalized_cov %*% z_WLS_joint

  ###Logit_joint_Unweighted###
  x_Logit_joint_u = data.frame(cbind(d0,d1,x_predict))
  reg_Logit_joint_u =  glm(y~d0+d1+x_predict-1,family='quasibinomial',weights=fo)
  z_Logit_joint_u = y-predict(reg_Logit_joint_u,x_Logit_joint_u,type='response')
  z_Logit_joint_u = z_Logit_joint_u * contrast_vec
  var_Logit_joint_u = t(z_Logit_joint_u) %*% normalized_cov %*% z_Logit_joint_u 
  var_Logit_joint_u=var_Logit_joint_u /(sample_size)^2
  
  ###Logit_joint_Weighted###
  x_Logit_joint_w = data.frame(cbind(d0,d1,x_predict))
  reg_Logit_joint_w =  glm(y~d0+d1+x_predict-1,family='quasibinomial')
  z_Logit_joint_w = y-predict(reg_Logit_joint_w,x_Logit_joint_w,type='response')
  z_Logit_joint_w = z_Logit_joint_w * contrast_vec
  
  var_Logit_joint_w = t(z_Logit_joint_w) %*% normalized_cov %*% z_Logit_joint_w 
  var_Logit_joint_w=var_Logit_joint_w /(sample_size)^2
  
  output=c(var_HT,var_HA,var_OLS_joint,var_WLS_joint,var_Logit_joint_u,var_Logit_joint_w)
  
  return(output)
}

theoretical_variance_bound=function(y,x,fo,so,sample_size,contrast){
  
  
  normalized_cov=diag(1/fo) %*% so %*% diag(1/fo)
  cov_bound = AS_bound(normalized_cov)
  
  y=rep(y,each=2)
  
  
  y_length = length(y)
  x_predict = x[rep(1:nrow(x),each=2),]
  d1 = rep(c(0,1),y_length/2)
  d0 = 1-d1
  
  contrast_OLS = c(contrast,rep(0,ncol(x)))
  contrast_vec=rep(contrast,length(y)/2)
  
  ####HT#####
  z_HT = y * contrast_vec
  var_HT = z_HT %*% cov_bound %*% z_HT
  var_HT = var_HT/(sample_size)^2
  
  ####HA#####    
  x_HA = data.frame(cbind(d0,d1))
  reg_HA =  lm(y~d1+d0-1)
  z_HA= y-predict(reg_HA,x_HA)
  z_HA = contrast_vec * z_HA
  var_HA = z_HA %*% cov_bound %*% z_HA
  var_HA = var_HA/(sample_size)^2
  
  ####OLS_joint####
  x_OLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_OLS_joint =  lm(y~d1+d0+x_predict-1,weights=fo)
  z_OLS_joint=  diag(y-predict(reg_OLS_joint,x_OLS_joint)) %*% diag(fo) %*%as.matrix(x_OLS_joint) %*% solve( t(as.matrix(x_OLS_joint)) %*%diag(fo)%*%as.matrix(x_OLS_joint))
  z_OLS_joint = z_OLS_joint %*% contrast_OLS
  var_OLS_joint = t(z_OLS_joint) %*% cov_bound %*% z_OLS_joint
  
  ####WLS_joint####
  x_WLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_WLS_joint =  lm(y~d1+d0+x_predict-1)
  z_WLS_joint=  diag(y-predict(reg_WLS_joint,x_WLS_joint)) %*% as.matrix(x_WLS_joint) %*% solve( t(as.matrix(x_WLS_joint)) %*%as.matrix(x_WLS_joint))
  z_WLS_joint = z_WLS_joint %*% contrast_OLS
  var_WLS_joint = t(z_WLS_joint) %*% cov_bound %*% z_WLS_joint
  
  ###Logit_joint_Unweighted###
  x_Logit_joint_u = data.frame(cbind(d0,d1,x_predict))
  reg_Logit_joint_u =  glm(y~d0+d1+x_predict-1,family='quasibinomial',weights=fo)
  z_Logit_joint_u = y-predict(reg_Logit_joint_u,x_Logit_joint_u,type='response')
  z_Logit_joint_u = z_Logit_joint_u * contrast_vec
  var_Logit_joint_u = t(z_Logit_joint_u) %*% cov_bound %*% z_Logit_joint_u 
  var_Logit_joint_u=var_Logit_joint_u /(sample_size)^2
  
  ###Logit_joint_Weighted###
  x_Logit_joint_w = data.frame(cbind(d0,d1,x_predict))
  reg_Logit_joint_w =  glm(y~d0+d1+x_predict-1,family='quasibinomial')
  z_Logit_joint_w = y-predict(reg_Logit_joint_w,x_Logit_joint_w,type='response')
  z_Logit_joint_w = z_Logit_joint_w * contrast_vec
  
  var_Logit_joint_w = t(z_Logit_joint_w) %*% cov_bound %*% z_Logit_joint_w 
  var_Logit_joint_w=var_Logit_joint_w /(sample_size)^2
  
  output=c(var_HT,var_HA,var_OLS_joint,var_WLS_joint,var_Logit_joint_u,var_Logit_joint_w)
  
  return(output)
}


AS_bound=function(normalized_cov){
  
  #rounds to 10th digit, this is arbitrary
  normalized_cov = round(normalized_cov,10)
  minus_ones=normalized_cov==-1
  add_on = apply(minus_ones,2,sum)
  
  
  #AS formula
  AS_bound = normalized_cov + minus_ones + diag(add_on)
  return(AS_bound)
  
  
}

HT_variance=function(y,sample_size,cov_bound,contrast){
  
  #make contrast estimator
  contrast_vec=rep(contrast,length(y)/2)
  
  #computing variance bound, accounting for contrasts
  var_est = 1/(sample_size)^2 * (contrast_vec*y) %*% cov_bound %*% (y*contrast_vec)
  return(var_est)
  
}

HT_var=function(y,d1,d0,sample_size,cov_bound,contrast){
 
  
   #make a vector of y's
   output=rep(y,each=2)
   
   #interaction of treatment indicator and y
   output[(1:length(output))%%2==0]=output[(1:length(output))%%2==0]*d1
   output[(1:length(output))%%2==1]=output[(1:length(output))%%2==1]*d0
   
   #computing variance
   var_HT=HT_variance(output,sample_size,cov_bound,contrast)
   
   return(var_HT)
   
  
  
}

HA_var=function(y,d1,d0,pi,sample_size,cov_bound,contrast){
  
  #create matrix for regression
  data=cbind(y,d0,d1)
  data=data.frame(data)
  weights = pi
  output=rep(y,each=2)

  #for formatting(lm+glm)
  x1=cbind(rep(0,sample_size),rep(1,sample_size))
  x0=cbind(rep(1,sample_size),rep(0,sample_size))
  x1=data.frame(x1)
  x0=data.frame(x0)
  colnames(x1)=colnames(data)[2:length(data)]
  colnames(x0)=colnames(data)[2:length(data)]
  
  reg=lm(y~d0+d1-1,weights=weights,data=data)
  y_p1 = predict(reg,x1)
  y_p0 = predict(reg,x0)
  
  fitted_1 = (y-y_p1) * d1
  fitted_0 = (y-y_p0) * d0
  
  output[(1:length(output))%%2==1]=fitted_1
  output[(1:length(output))%%2==0]=fitted_0 
  
  var_HA=HT_variance(output,sample_size,cov_bound,contrast)
  return(var_HA)
  
}

Logit_var=function(y,d1,d0,weights,x,pi1,pi0,mode,sample_size,cov_bound,contrast){
  
  
  if (mode=='unweighted_joint'){
    
    data=cbind(y,d0,d1,x)
    data=data.frame(data)
    x1=cbind(rep(0,nrow(x)),rep(1,nrow(x)),x)
    x0=cbind(rep(1,nrow(x)),rep(0,nrow(x)),x)
    x1=data.frame(x1)
    x0=data.frame(x0)
    colnames(x1)=colnames(data)[2:length(data)]
    colnames(x0)=colnames(data)[2:length(data)]
    
    weights = as.numeric((weights>0))
    reg = glm(y~d0+d1+x-1,family='quasibinomial',weights=weights,data=data)
 
    y_p1 = predict(reg,x1,type='response')
    y_p0 = predict(reg,x0,type='response')
    
    fitted_1 = (y-y_p1) * d1
    fitted_0 = (y-y_p0) * d0
  }else if (mode=='unweighted_separate'){
    
    data1=data.frame(cbind(y,x))
    data0=data.frame(cbind(y,x))
    x_predict=data.frame(x)
    colnames(x_predict)=colnames(data1)[2:length(data1)]
    reg1 = glm(y~x,family='quasibinomial',weights=d1,data=data1)
    reg0 = glm(y~x,family='quasibinomial',weights=d0,data=data0)
    
    y_p1 = predict(reg1,x_predict,type='response')
    y_p0 = predict(reg0,x_predict,type='response')    
    
    fitted_1 = (y-y_p1) * d1
    fitted_0 = (y-y_p0) * d0
  }else if (mode=='weighted_joint'){
    
    data=cbind(y,d0,d1,x)
    data=data.frame(data)
    x1=cbind(rep(0,nrow(x)),rep(1,nrow(x)),x)
    x0=cbind(rep(1,nrow(x)),rep(0,nrow(x)),x)
    x1=data.frame(x1)
    x0=data.frame(x0)
    colnames(x1)=colnames(data)[2:length(data)]
    colnames(x0)=colnames(data)[2:length(data)]
    
    weights = weights
    reg = glm(y~d0+d1+x-1,family='quasibinomial',weights=weights)
    
    y_p1 = predict(reg,x1,type='response')
    y_p0 = predict(reg,x0,type='response')
    
    
    #fitted value for variance calculation
    fitted_1 = (y-y_p1) * d1
    fitted_0 = (y-y_p0) * d0
  }else if (mode=='weighted_separate'){
    data1=data.frame(cbind(y,x))
    data0=data.frame(cbind(y,x))
    x_predict=data.frame(x)
    weights = weights
    
    colnames(x_predict)=colnames(data1)[2:length(data1)]
    reg1 = glm(y~x,family='quasibinomial',weights=d1*weights,data=data1)
    reg0 = glm(y~x,family='quasibinomial',weights=d0*weights,data=data0)
    
    y_p1 = predict(reg1,x_predict,type='response')
    y_p0 = predict(reg0,x_predict,type='response')    
    
    #fitted value for variance calculation
    fitted_1 = (y-y_p1) * d1
    fitted_0 = (y-y_p0) * d0
  }
  output=rep(y,each=2)
  output[(1:length(output))%%2==1]=fitted_0
  output[(1:length(output))%%2==0]=fitted_1  
  
  var_Logit=HT_variance(output,sample_size,cov_bound,contrast)
  return(var_Logit)
  #return(c(avg_effect1-avg_effect0,coef))
}

OLS_var=function(y,d1,d0,pi,x,mode,cov_bound,contrast){
  
  output=rep(y,each=2)

  contrast_OLS = c(contrast,rep(0,ncol(x)))
  

  if(mode=='OLS_joint'){
    
    data=cbind(y,d0,d1,x)
    data=data.frame(data)
    weights = as.numeric((pi>0))

    x1=cbind(rep(0,nrow(x)),rep(1,nrow(x)),x)
    x0=cbind(rep(1,nrow(x)),rep(0,nrow(x)),x)
    x1=data.frame(x1)
    x0=data.frame(x0)
    colnames(x1)=colnames(data)[2:length(data)]
    colnames(x0)=colnames(data)[2:length(data)]
    
    weights = as.numeric((weights>0))
    reg=lm(y~d0+d1+x-1,weights=weights,data=data)
    y_p1 = predict(reg,x1,type='response')
    y_p0 = predict(reg,x0,type='response')
    
    fitted_1 = (y-y_p1) * d1
    fitted_0 = (y-y_p0) * d0
 
    c1 = rep(c(0,1),length(output)/2)
    c0 = 1-c1
    x_predict = x[rep(1:nrow(x),each=2),]
    x_var = data.frame(cbind(c0,c1,x_predict))   
    output[(1:length(output))%%2==0]=fitted_1
    output[(1:length(output))%%2==1]=fitted_0  
    
    z_OLS_joint_hat=  diag(output) %*% diag(fo) %*%as.matrix(x_var) %*% solve( t(as.matrix(x_var)) %*%diag(fo)%*%as.matrix(x_var))
    z_OLS_joint_hat = z_OLS_joint_hat %*% contrast_OLS
    var_output =   t(z_OLS_joint_hat) %*%   cov_bound  %*% z_OLS_joint_hat
  }else if (mode=='OLS_separate'){
    
    data1=data.frame(cbind(y,x))
    data0=data.frame(cbind(y,x))
    x_predict=data.frame(x)
    colnames(x_predict)=colnames(data1)[2:length(data1)]
    reg1 = lm(y~x,weights=d1,data=data1)
    reg0 = lm(y~x,weights=d0,data=data0)
    
    y_p1 = predict(reg1,x_predict,type='response')
    y_p0 = predict(reg0,x_predict,type='response')    
    
    fitted_1 = (y-y_p1) * d1
    fitted_0 = (y-y_p0) * d0
 
    output[(1:length(output))%%2==1]=fitted_0
    output[(1:length(output))%%2==0]=fitted_1  
    
    z_OLS_separate_hat=  output %*% diag(fo) %*%as.matrix(x_var) %*% solve( t(as.matrix(x_var)) %*%diag(fo)%*%as.matrix(x_var))
    z_OLS_separatt_hat = z_separate_joint_hat %*% contrast_OLS
    var_output =   z_separate_joint_hat %*%   cov_bound  %*% z_separate_joint_hat
  }else if (mode=='WLS_joint'){
    data=cbind(y,d0,d1,x)
    data=data.frame(data)
    weights = as.numeric((pi>0))
    
    x1=cbind(rep(0,nrow(x)),rep(1,nrow(x)),x)
    x0=cbind(rep(1,nrow(x)),rep(0,nrow(x)),x)
    x1=data.frame(x1)
    x0=data.frame(x0)
    colnames(x1)=colnames(data)[2:length(data)]
    colnames(x0)=colnames(data)[2:length(data)]
    
    reg=lm(y~d0+d1+x-1,weights=weights,data=data)
    y_p1 = predict(reg,x1,type='response')
    y_p0 = predict(reg,x0,type='response')
    
    fitted_1 = (y-y_p1) * d1
    fitted_0 = (y-y_p0) * d0
    
    
    c1 = rep(c(0,1),length(output)/2)
    c0 = 1-c1
    x_predict = x[rep(1:nrow(x),each=2),]
    x_var = data.frame(cbind(c0,c1,x_predict))   
    output[(1:length(output))%%2==0]=fitted_1
    output[(1:length(output))%%2==1]=fitted_0  
    
    
    z_WLS_joint_hat=  diag(output)  %*%as.matrix(x_var) %*% solve( t(as.matrix(x_var)) %*% as.matrix(x_var))
    z_WLS_joint_hat = z_WLS_joint_hat %*% contrast_OLS
    var_output =   t(z_WLS_joint_hat) %*%   cov_bound  %*% z_WLS_joint_hat
  }else if (mode=='WLS_separate'){
    data1=data.frame(cbind(y,x))
    data0=data.frame(cbind(y,x))
    x_predict=data.frame(x)
    colnames(x_predict)=colnames(data1)[2:length(data1)]
    reg1 = lm(y~x,weights=d1*weights,data=data1)
    reg0 = lm(y~x,weights=d0*weights,data=data0)
    
    y_p1 = predict(reg1,x_predict,type='response')
    y_p0 = predict(reg0,x_predict,type='response')    
    
    fitted_1 = (y-y_p1) * d1
    fitted_0 = (y-y_p0) * d0
    
    z_WLS_separate_hat=  output  %*%as.matrix(x_var) %*% solve( t(as.matrix(x_var)) %*% as.matrix(x_var))
    z_WLS_separate_hat = z_WLS_separate_hat %*% contrast_OLS
    var_output =   z_WLS_separate_hat %*%   cov_bound  %*% z_WLS_separate_hat
  }

  

  return(var_output)
  #return(c(avg_effect1-avg_effect0,coef))
  
}

simulationFeb11=function(net,y,x,prob,id_list,expo1,expo2){
  
  
  pol_size = network.size(net)
  status = 0:3-1
  num_mappings = 12
  assign_prob=prob
  network_subjects = net %v% 'vertex.names'
  subjects = (1:pol_size)[network_subjects %in% id_list]
  
  subject_size = length(subjects)
  
  
  nsim=100000
  expo1=expo1
  expo2=expo2
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
  
  #method=c('HT','HA','OLS_joint','OLS_separate','WLS_joint','WLS_separate','unweighted_joint',
  #         'unweighted_sepaprate','weighted_joint','weighted_separate')
  method=c('HT','HA','OLS_joint','WLS_joint','unweighted_joint', 'weighted_joint')
  est_sim=matrix(0,nrow=nsim,ncol=length(method)) 
  var_sim=matrix(0,nrow=nsim,ncol=length(method)) 
  
  colnames(est_sim)=method
  
  #est_sim=rep(0,nsim)
  sample_size = matrix(0,nrow=2,ncol = nsim)
  friends=matrix(FALSE,nrow=pol_size,ncol=5)
  for (i in 1:pol_size){
    friend_index_t=get.neighborhood(net,i,'out')
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
    exposure_t[,6]=num_friend_treated_t[1,]
    exposure_t[,7]=num_friend_treated_t[2,]  
    
    
    exposure_t = exposure_t[subjects,2:7]
    
    #please double check every time 
    state_t=apply(exposure_t,1,individual_exposure) 
    y1_t = state_t[expo1,] * y
    y2_t = state_t[expo2,] * y
    denom_pi1_t= state_t[expo1,]* 1/pi1 
    denom_pi2_t= state_t[expo2,]* 1/pi2 
    pi_weight_t = denom_pi1_t + denom_pi2_t
    d1_t=state_t[expo2,]
    d0_t=state_t[expo1,]


    
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
    
    est_sim[i,1]= HT(y2_t,y1_t,pi2,pi1,subject_size)
    est_sim[i,2]= HA(y2_t,y1_t,pi2,pi1,denom_pi2_t,denom_pi1_t)
    est_sim[i,3]= OLS(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint')
    est_sim[i,4]= OLS(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint')    
    est_sim[i,5]= Logit(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='unweighted_joint')
    est_sim[i,6]= Logit(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='weighted_joint')
    
    
    var_est[i,1]=HT(y,subject_size,cov_bound)
    
    
    
  }
  
  return(est_sim)
}

simulationFeb13=function(net,y,x,prob,id_list,expo1,expo2,cov_bound_HT){
  
  
  pol_size = network.size(net)
  status = 0:3-1
  num_mappings = 12
  assign_prob=prob
  network_subjects = net %v% 'vertex.names'
  subjects = (1:pol_size)[network_subjects %in% id_list]
  subject_size = length(subjects)
  
  
  nsim=10000
  expo1=expo1
  expo2=expo2
  fo= FO(net=net, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list=id_list,expo1=expo1,expo2=expo2)
  pi1 = fo[(1:length(fo))%%2==1]
  pi2 = fo[(1:length(fo))%%2==0]
  pol_size = network.size(net)
  status = 0:3-1
  num_mappings = 12
  assign_prob=prob
  network_subjects = net %v% 'vertex.names'
  subjects = (1:pol_size)[network_subjects %in% id_list]
  subject_size = length(subjects)
  contrast=c(-1,1)
  
  #prob
  pi1 = pi1[subjects]
  pi2 = pi2[subjects]
  
  #method=c('HT','HA','OLS_joint','OLS_separate','WLS_joint','WLS_separate','unweighted_joint',
  #         'unweighted_sepaprate','weighted_joint','weighted_separate')
  method=c('HT','HA','OLS_joint','WLS_joint','unweighted_joint', 'weighted_joint')
  est_sim=matrix(0,nrow=nsim,ncol=length(method)) 
  var_sim=matrix(0,nrow=nsim,ncol=length(method)) 
  
  colnames(est_sim)=method
  
  #est_sim=rep(0,nsim)
  sample_size = matrix(0,nrow=2,ncol = nsim)
  friends=matrix(FALSE,nrow=pol_size,ncol=5)
  for (i in 1:pol_size){
    friend_index_t=get.neighborhood(net,i,'out')
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
    exposure_t[,6]=num_friend_treated_t[1,]
    exposure_t[,7]=num_friend_treated_t[2,]  
    
    
    exposure_t = exposure_t[subjects,2:7]
    
    #please double check every time 
    state_t=apply(exposure_t,1,individual_exposure) 
    y1_t = state_t[expo1,] * y
    y2_t = state_t[expo2,] * y
    denom_pi1_t= state_t[expo1,]* 1/pi1 
    denom_pi2_t= state_t[expo2,]* 1/pi2 
    pi_weight_t = denom_pi1_t + denom_pi2_t
    d1_t=state_t[expo2,]
    d0_t=state_t[expo1,]
    
    
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
    
    
    est_sim[i,1]= HT(y2_t,y1_t,pi2,pi1,subject_size)
    est_sim[i,2]= HA(y2_t,y1_t,pi2,pi1,denom_pi2_t,denom_pi1_t)
    est_sim[i,3]= OLS(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint')
    est_sim[i,4]= OLS(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint')    
    est_sim[i,5]= Logit(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='unweighted_joint')
    est_sim[i,6]= Logit(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='weighted_joint')
    

    var_sim[i,1]= HT_var(y_obs_t,d1_t,d0_t,subject_size,cov_bound_HT,contrast)
    var_sim[i,2]= HA_var(y_obs_t,d1_t,d0_t,pi_weight_t,subject_size,cov_bound_HT,contrast)
    var_sim[i,3]= OLS_var(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint',cov_bound_HT,contrast)
    var_sim[i,4]= OLS_var(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint',cov_bound_HT,contrast)
    var_sim[i,5]= Logit_var(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='unweighted_joint',subject_size,cov_bound_HT,contrast)
    var_sim[i,6]= Logit_var(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='weighted_joint',subject_size,cov_bound_HT,contrast)
    
  }
  output=list(est_sim,var_sim)
  return(output)
}

simulationFeb15=function(net,y,x,prob,id_list,expo1,expo2,cov_bound_HT){
  
  
  pol_size = network.size(net)
  status = 0:3-1
  num_mappings = 12
  assign_prob=prob
  network_subjects = net %v% 'vertex.names'
  subjects = (1:pol_size)[network_subjects %in% id_list]
  subject_size = length(subjects)
  
  
  nsim=10000
  expo1=expo1
  expo2=expo2
  
  ###read from some input files
  fo= FO(net=net, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=100,id_list=id_list,expo1=expo1,expo2=expo2)

  
  pi1 = fo[(1:length(fo))%%2==1]
  pi2 = fo[(1:length(fo))%%2==0]
  pol_size = network.size(net)
  status = 0:3-1
  num_mappings = 12
  assign_prob=prob
  network_subjects = comp_5 %v% 'vertex.names'
  subjects = (1:pol_size)[network_subjects %in% id_list]
  subject_size = length(subjects)
  contrast=c(-1,1)
  
  #prob
  pi1 = pi1[subjects]
  pi2 = pi2[subjects]
  
  #method=c('HT','HA','OLS_joint','OLS_separate','WLS_joint','WLS_separate','unweighted_joint',
  #         'unweighted_sepaprate','weighted_joint','weighted_separate')
  method=c('HT','HA','OLS_joint','WLS_joint','unweighted_joint', 'weighted_joint')
  est_sim=matrix(0,nrow=nsim,ncol=length(method)) 
  var_sim=matrix(0,nrow=nsim,ncol=length(method)) 
  
  colnames(est_sim)=method
  
  #est_sim=rep(0,nsim)
  sample_size = matrix(0,nrow=2,ncol = nsim)
  friends=matrix(FALSE,nrow=pol_size,ncol=5)
  for (i in 1:pol_size){
    friend_index_t=get.neighborhood(net,i,'out')
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
    exposure_t[,6]=num_friend_treated_t[1,]
    exposure_t[,7]=num_friend_treated_t[2,]  
    
    
    exposure_t = exposure_t[subjects,2:7]
    
    #please double check every time 
    state_t=apply(exposure_t,1,individual_exposure) 
    y1_t = state_t[expo1,] * y
    y2_t = state_t[expo2,] * y
    denom_pi1_t= state_t[expo1,]* 1/pi1 
    denom_pi2_t= state_t[expo2,]* 1/pi2 
    pi_weight_t = denom_pi1_t + denom_pi2_t
    d1_t=state_t[expo2,]
    d0_t=state_t[expo1,]
    
    
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
    
    
    est_sim[i,1]= HT(y2_t,y1_t,pi2,pi1,subject_size)
    est_sim[i,2]= HA(y2_t,y1_t,pi2,pi1,denom_pi2_t,denom_pi1_t)
    est_sim[i,3]= OLS(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint')
    est_sim[i,4]= OLS(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint')    
    est_sim[i,5]= Logit(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='unweighted_joint')
    est_sim[i,6]= Logit(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='weighted_joint')
    
    
    var_sim[i,1]= HT_var(y_obs_t,d1_t,d0_t,subject_size,cov_bound_HT,contrast)
    var_sim[i,2]= HA_var(y_obs_t,d1_t,d0_t,pi_weight_t,subject_size,cov_bound_HT,contrast)
    var_sim[i,3]= OLS_var(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint',cov_bound_HT,contrast)
    var_sim[i,4]= OLS_var(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint',cov_bound_HT,contrast)
    var_sim[i,5]= Logit_var(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='unweighted_joint',subject_size,cov_bound_HT,contrast)
    var_sim[i,6]= Logit_var(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi1,pi2,mode='weighted_joint',subject_size,cov_bound_HT,contrast)
    
  }
  output=list(est_sim,var_sim)
  return(output)
}

