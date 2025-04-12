# This file contain functions used in simulations:
# 
# A. HT Estimator
# 
# B. HA Estimaator
# 



HT_full = function(epsilon1,epsilon0,pi1,pi0,sample_size){
  
  est = sum(epsilon1 * (1/pi1)) - sum(epsilon0 * (1/pi0))
  
  est = est/sample_size
  
  return(est)
  
  
  
}

HA_full =function(epsilon1,epsilon0,pi1,pi0,denom_pi1,denom_pi0){
  
  est = sum(epsilon1 * (1/pi1))/sum(denom_pi1) - sum(epsilon0 * (1/pi0))/sum(denom_pi0)
  
  return(est)
  
}

OLS_full=function(y,d1,d0,pi,x,mode){
  
  
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

Logit_full=function(y,d1,d0,weights,x,pi0,pi1,mode){
  
  
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


theoretical_variance_full=function(y,x,fo,so,contrast,dict_group){
  
  
  subject_size=length(y)
  y=rep(y,each=2)
  fo_vec=unlist(fo)
  y_length = length(y)
  x_predict = x[rep(1:nrow(x),each=2),]
  d1 = rep(c(0,1),y_length/2)
  d0 = 1-d1
  
  contrast_OLS = c(contrast,rep(0,ncol(x)))
  contrast_vec=rep(contrast,length(y)/2)
  
  ####HT#####
  z_HT = y * contrast_vec
  
  ####HA#####    
  x_HA = data.frame(cbind(d0,d1))
  reg_HA =  lm(y~d1+d0-1)
  z_HA= y-predict(reg_HA,x_HA)
  z_HA = contrast_vec * z_HA
  
  
  ####OLS_joint####
  x_OLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_OLS_joint =  lm(y~d1+d0+x_predict-1,weights=fo_vec)
  temp_OLS =  (y-predict(reg_OLS_joint,x_OLS_joint)) * fo_vec
  z_OLS_joint=  temp_OLS * as.matrix(x_OLS_joint) %*% solve( t(as.matrix(x_OLS_joint)) %*%diag(fo_vec)%*%as.matrix(x_OLS_joint)*(1/subject_size))
  z_OLS_joint = z_OLS_joint %*% contrast_OLS
  
  ####WLS_joint####
  x_WLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_WLS_joint =  lm(y~d1+d0+x_predict-1)
  temp_WLS = y-predict(reg_WLS_joint,x_WLS_joint)
  z_WLS_joint=  temp_WLS* as.matrix(x_WLS_joint) %*% solve( t(as.matrix(x_WLS_joint)) %*%as.matrix(x_WLS_joint)* (1/subject_size))
  z_WLS_joint = z_WLS_joint %*% contrast_OLS
  
  ###Logit_joint_Unweighted###
  x_Logit_joint_u = data.frame(cbind(d0,d1,x_predict))
  reg_Logit_joint_u =  glm(y~d0+d1+x_predict-1,family='quasibinomial',weights=fo_vec)
  z_Logit_joint_u = y-predict(reg_Logit_joint_u,x_Logit_joint_u,type='response')
  z_Logit_joint_u = z_Logit_joint_u * contrast_vec
  
  
  ###Logit_joint_Weighted###
  x_Logit_joint_w = data.frame(cbind(d0,d1,x_predict))
  reg_Logit_joint_w =  glm(y~d0+d1+x_predict-1,family='quasibinomial')
  z_Logit_joint_w = y-predict(reg_Logit_joint_w,x_Logit_joint_w,type='response')
  z_Logit_joint_w = z_Logit_joint_w * contrast_vec
  
  output=rep(0,6)
  names(output)=c('var_HT','var_HA','var_OLS_joint','var_WLS_joint','var_Logit_joint_u','var_Logit_joint_w')
  for (i in 1:nrow(dict_group)){
    
    index_start=dict_group[i,3]
    index_end=dict_group[i,4]    
    fo_t=fo[[i]]
    so_t=so[[i]]
    normalized_cov=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    
    z_HT_t = z_HT[index_start:index_end]
    z_HA_t = z_HA[index_start:index_end]
    z_OLS_joint_t = z_OLS_joint[index_start:index_end]
    z_WLS_joint_t = z_WLS_joint[index_start:index_end]
    z_Logit_joint_u_t =z_Logit_joint_u[index_start:index_end]
    z_Logit_joint_w_t = z_Logit_joint_w[index_start:index_end]
    
    output[1]=output[1] +  t(z_HT_t) %*% normalized_cov %*% z_HT_t
    output[2]=output[2] +  t(z_HA_t) %*% normalized_cov %*% z_HA_t
    output[3]=output[3] +  t(z_OLS_joint_t) %*% normalized_cov %*% z_OLS_joint_t
    output[4]=output[4] +  t(z_WLS_joint_t) %*% normalized_cov %*% z_WLS_joint_t
    output[5]=output[5] +  t(z_Logit_joint_u_t) %*% normalized_cov %*% z_Logit_joint_u_t
    output[6]=output[6] +  t(z_Logit_joint_w_t) %*% normalized_cov %*%  z_Logit_joint_w_t
    
    
  }
  
  output=output/(subject_size^2)
  
  return(output)
}

theoretical_variance_bound_full=function(y,x,fo,so,contrast,dict_group){
  
  
  subject_size=length(y)
  y=rep(y,each=2)
  fo_vec=unlist(fo)
  y_length = length(y)
  x_predict = x[rep(1:nrow(x),each=2),]
  d1 = rep(c(0,1),y_length/2)
  d0 = 1-d1
  
  contrast_OLS = c(contrast,rep(0,ncol(x)))
  contrast_vec=rep(contrast,length(y)/2)
  
  ####HT#####
  z_HT = y * contrast_vec
  
  ####HA#####    
  x_HA = data.frame(cbind(d0,d1))
  reg_HA =  lm(y~d1+d0-1)
  z_HA= y-predict(reg_HA,x_HA)
  z_HA = contrast_vec * z_HA
  
  ####OLS_joint####
  x_OLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_OLS_joint =  lm(y~d1+d0+x_predict-1,weights=fo_vec)
  temp_OLS =  (y-predict(reg_OLS_joint,x_OLS_joint)) * fo_vec
  z_OLS_joint=  temp_OLS * as.matrix(x_OLS_joint) %*% solve( t(as.matrix(x_OLS_joint)) %*%diag(fo_vec)%*%as.matrix(x_OLS_joint)*(1/subject_size))
  z_OLS_joint = z_OLS_joint %*% contrast_OLS
  
  ####WLS_joint####
  x_WLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_WLS_joint =  lm(y~d1+d0+x_predict-1)
  temp_WLS = y-predict(reg_WLS_joint,x_WLS_joint)
  z_WLS_joint=  temp_WLS* as.matrix(x_WLS_joint) %*% solve( t(as.matrix(x_WLS_joint)) %*%as.matrix(x_WLS_joint)* (1/subject_size))
  z_WLS_joint = z_WLS_joint %*% contrast_OLS
  
  ###Logit_joint_Unweighted###
  x_Logit_joint_u = data.frame(cbind(d0,d1,x_predict))
  reg_Logit_joint_u =  glm(y~d0+d1+x_predict-1,family='quasibinomial',weights=fo_vec)
  z_Logit_joint_u = y-predict(reg_Logit_joint_u,x_Logit_joint_u,type='response')
  z_Logit_joint_u = z_Logit_joint_u * contrast_vec
  
  
  ###Logit_joint_Weighted###
  x_Logit_joint_w = data.frame(cbind(d0,d1,x_predict))
  reg_Logit_joint_w =  glm(y~d0+d1+x_predict-1,family='quasibinomial')
  z_Logit_joint_w = y-predict(reg_Logit_joint_w,x_Logit_joint_w,type='response')
  z_Logit_joint_w = z_Logit_joint_w * contrast_vec
  
  output=rep(0,6)
  names(output)=c('var_HT','var_HA','var_OLS_joint','var_WLS_joint','var_Logit_joint_u','var_Logit_joint_w')
  for (i in 1:nrow(dict_group)){
    
    index_start=dict_group[i,3]
    index_end=dict_group[i,4]    
    fo_t=fo[[i]]
    so_t=so[[i]]
    normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    
    normalized_cov_t=round(normalized_cov_t,10)
    cov_bound_t = AS_bound(normalized_cov_t)
    
    
    z_HT_t = z_HT[index_start:index_end]
    z_HA_t = z_HA[index_start:index_end]
    z_OLS_joint_t = z_OLS_joint[index_start:index_end]
    z_WLS_joint_t = z_WLS_joint[index_start:index_end]
    z_Logit_joint_u_t =z_Logit_joint_u[index_start:index_end]
    z_Logit_joint_w_t = z_Logit_joint_w[index_start:index_end]
    
    output[1]=output[1] +  t(z_HT_t) %*% cov_bound_t %*% z_HT_t
    output[2]=output[2] +  t(z_HA_t) %*% cov_bound_t %*% z_HA_t
    output[3]=output[3] +  t(z_OLS_joint_t) %*% cov_bound_t %*% z_OLS_joint_t
    output[4]=output[4] +  t(z_WLS_joint_t) %*% cov_bound_t %*% z_WLS_joint_t
    output[5]=output[5] +  t(z_Logit_joint_u_t) %*% cov_bound_t %*% z_Logit_joint_u_t
    output[6]=output[6] +  t(z_Logit_joint_w_t) %*% cov_bound_t %*%  z_Logit_joint_w_t
  }
  
  output=output/(subject_size^2)
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

HT_variance_full=function(y,sample_size,cov_bound,contrast){
  
  #make contrast estimator
  contrast_vec=rep(contrast,length(y)/2)
  
  #computing variance bound, accounting for contrasts
  var_est = 1/(sample_size)^2 * (contrast_vec*y) %*% cov_bound %*% (y*contrast_vec)
  return(var_est)
  
}

HT_var_full=function(y,d1,d0,sample_size,fo,so,contrast,dict_group){
  
  
  
  #computing variance
  num_component=nrow(dict_group)
  subject_size=length(y)
  
  var_est=0
  for (i in 1:num_component){
    
    index_start=dict_group[i,1]
    index_end=dict_group[i,2]
    
    #select needed component
    y_t = y[index_start:index_end]  
    
    d1_t=d1[index_start:index_end]  
    d0_t=d0[index_start:index_end]  
    output_t=rep(y_t,each=2)    
    output_t[(1:length(output_t))%%2==0]= (output_t[(1:length(output_t))%%2==0]*d1_t)
    output_t[(1:length(output_t))%%2==1]=output_t[(1:length(output_t))%%2==1]*d0_t   
    
    contrast_vec=rep(contrast,length(output_t)/2)
    
    fo_t=fo[[i]]
    so_t=so[[i]]
    
    normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    
    cov_bound_t = AS_bound(normalized_cov_t)
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    p_zeros_t=(p_t==0)
    cov_bound_zeros_t = (cov_bound_t==0)
    #print('Check if denomiator matrix is correct')
    if (min(cov_bound_zeros_t-p_zeros_t)<0){
      print('The denomiator matrix is wrong')
    }
    #print(min(cov_bound_zeros_t-p_zeros_t)) #should expect 0
    print('HT')
    print(i)
    
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    var_est = var_est + (contrast_vec*output_t) %*% cov_bound_HT_t %*% (output_t*contrast_vec)
    
  }
  
  var_est=var_est/(subject_size)^2
  return(var_est)
  
  
  
}

HA_var_full=function(y,d1,d0,pi,sample_size,fo,so,contrast,dict_group){
  
  num_component=nrow(dict_group)
  subject_size=length(y)
  var_est=0
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
  
  for (i in 1:num_component){
    
    print('HA')
    print(i)
    index_start=dict_group[i,1]
    index_end=dict_group[i,2]
    
    #select needed component
    y_t = y[index_start:index_end]  
    
    output_t=rep(y_t,each=2)    
    output_t[(1:length(output_t))%%2==0]=fitted_1[index_start:index_end]
    output_t[(1:length(output_t))%%2==1]=fitted_0[index_start:index_end]   
    
    contrast_vec=rep(contrast,length(output_t)/2)
    
    fo_t=fo[[i]]
    so_t=so[[i]]
    
    normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    
    cov_bound_t = AS_bound(normalized_cov_t)
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    p_zeros_t=(p_t==0)
    cov_bound_zeros_t = (cov_bound_t==0)
    if (min(cov_bound_zeros_t-p_zeros_t)<0){
      print('The denomiator matrix is wrong')
    }
    
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    contrast_vec=rep(contrast,length(output_t)/2)
    
    var_est = var_est + (contrast_vec*output_t) %*% cov_bound_HT_t %*% (output_t*contrast_vec)
    
  }
  
  var_est=var_est/(subject_size)^2
  
  return(var_est)
  
}

Logit_var_full=function(y,d1,d0,weights,x,mode,fo,so,contrast,dict_group){
  
  
  num_component=nrow(dict_group)
  subject_size=length(y)
  var_est=0
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
  
  for (i in 1:num_component){
    
    index_start=dict_group[i,1]
    index_end=dict_group[i,2]
    
    #select needed component
    y_t = y[index_start:index_end]  
    
    output_t=rep(y_t,each=2)    
    output_t[(1:length(output_t))%%2==0]=fitted_1[index_start:index_end]
    output_t[(1:length(output_t))%%2==1]=fitted_0[index_start:index_end]   
    
    contrast_vec=rep(contrast,length(output_t)/2)
    
    print('Logit')
    print(i)
    fo_t=fo[[i]]
    so_t=so[[i]]
    
    normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    
    cov_bound_t = AS_bound(normalized_cov_t)
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    p_zeros_t=(p_t==0)
    cov_bound_zeros_t = (cov_bound_t==0)
    if (min(cov_bound_zeros_t-p_zeros_t)<0){
      print('The denomiator matrix is wrong')
    }
    
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    contrast_vec=rep(contrast,length(output_t)/2)
    
    var_est = var_est + (contrast_vec*output_t) %*% cov_bound_HT_t %*% (output_t*contrast_vec)
    
  }
  
  var_est=var_est/(subject_size)^2
  
  return(var_est)
  #return(c(avg_effect1-avg_effect0,coef))
}

OLS_var_full=function(y,d1,d0,pi,x,mode,fo,so,contrast,dict_group){
  
  num_component=nrow(dict_group)
  subject_size=length(y)
  output=rep(y,each=2)
  contrast_OLS = c(contrast,rep(0,ncol(x)))
  fo_vec=unlist(fo)
  var_est=0
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
    
    print('in OLS')
    output = output*fo_vec
    #z_OLS_joint_hat=  diag(output) %*% diag(fo_vec) %*%as.matrix(x_var) %*% solve( t(as.matrix(x_var)) %*%diag(fo_vec)%*%as.matrix(x_var)*1/subject_size)
    z_OLS_joint_hat=  output*as.matrix(x_var) %*% solve( t(as.matrix(x_var)) %*%diag(fo_vec)%*%as.matrix(x_var)*1/subject_size)
    z_OLS_joint_hat = z_OLS_joint_hat %*% contrast_OLS
    print('finished')
    for (i in 1:num_component){
      
      print('OLS')
      print(i)
      index_start=dict_group[i,3]
      index_end=dict_group[i,4]
      
      #select needed component
      y_t = y[index_start:index_end]  
      
      output_t=rep(y_t,each=2)    
      output_t=z_OLS_joint_hat[index_start:index_end]
      
      fo_t=fo[[i]]
      so_t=so[[i]]
      
      normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
      
      cov_bound_t = AS_bound(normalized_cov_t)
      p_t = so_t + (fo_t) %*% t(fo_t) 
      p_t=round(p_t,10)
      p_zeros_t=(p_t==0)
      cov_bound_zeros_t = (cov_bound_t==0)
      if (min(cov_bound_zeros_t-p_zeros_t)<0){
        print('The denomiator matrix is wrong')
      }
      
      p_denomiator_t = p_t + cov_bound_zeros_t
      cov_bound_HT_t= cov_bound_t / p_denomiator_t
      
      var_est = var_est + (output_t) %*% cov_bound_HT_t %*% (output_t)
      
    }
    
    var_est = var_est/(subject_size^2)
    
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
    
    
    z_WLS_joint_hat=  output*as.matrix(x_var) %*% solve( t(as.matrix(x_var)) %*% as.matrix(x_var)*1/subject_size)
    z_WLS_joint_hat = z_WLS_joint_hat %*% contrast_OLS
    
    for (i in 1:num_component){
      
      print('WLS')
      print(i)
      index_start=dict_group[i,3]
      index_end=dict_group[i,4]
      
      #select needed component
      y_t = y[index_start:index_end]  
      
      output_t=rep(y_t,each=2)    
      output_t=z_WLS_joint_hat[index_start:index_end]
      
      fo_t=fo[[i]]
      so_t=so[[i]]
      
      normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
      
      cov_bound_t = AS_bound(normalized_cov_t)
      p_t = so_t + (fo_t) %*% t(fo_t) 
      p_t=round(p_t,10)
      p_zeros_t=(p_t==0)
      cov_bound_zeros_t = (cov_bound_t==0)
      if (min(cov_bound_zeros_t-p_zeros_t)<0){
        print('The denomiator matrix is wrong')
      }
      
      p_denomiator_t = p_t + cov_bound_zeros_t
      cov_bound_HT_t= cov_bound_t / p_denomiator_t
      
      var_est = var_est + (output_t) %*% cov_bound_HT_t %*% (output_t)
      
    }
    
    var_est = var_est/(subject_size^2)
    
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
  
  
  
  return(var_est)
  #return(c(avg_effect1-avg_effect0,coef))
  
}

