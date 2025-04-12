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

OLS_full=function(y,d1,d0,weights,x,pi0,pi1,mode){
  if (mode=='OLS_joint'){
    
    data=cbind(y,d0,d1,x)
    data=data.frame(data)
    x1=cbind(rep(0,nrow(x)),rep(1,nrow(x)),x)
    x0=cbind(rep(1,nrow(x)),rep(0,nrow(x)),x)
    x1=data.frame(x1)
    x0=data.frame(x0)
    colnames(x1)=colnames(data)[2:length(data)]
    colnames(x0)=colnames(data)[2:length(data)]
    
    weights = as.numeric((weights>0))
    reg = lm(y~d0+d1+x-1,weights=weights,data=data)
    coef=reg$coefficients
    y_p1 = predict(reg,x1,type='response')
    y_p0 = predict(reg,x0,type='response')
    avg_effect1=reg$coefficients['d1']
    avg_effect0=reg$coefficients['d0']
  }else if (mode=='WLS_joint'){
    
    data=cbind(y,d0,d1,x)
    data=data.frame(data)
    x1=cbind(rep(0,nrow(x)),rep(1,nrow(x)),x)
    x0=cbind(rep(1,nrow(x)),rep(0,nrow(x)),x)
    x1=data.frame(x1)
    x0=data.frame(x0)
    colnames(x1)=colnames(data)[2:length(data)]
    colnames(x0)=colnames(data)[2:length(data)]
    
    weights = weights
    reg = lm(y~d0+d1+x-1,weights=weights)
    
    y_p1 = predict(reg,x1,type='response')
    y_p0 = predict(reg,x0,type='response')
    coef=reg$coefficients 
    avg_effect1=reg$coefficients['d1']
    avg_effect0=reg$coefficients['d0']
  }
  return(list(avg_effect1-avg_effect0,y_p1,y_p0))
  #return(c(avg_effect1-avg_effect0,coef))
  
}

logistic=function(u){
  
  return(exp(u)/(1+exp(u)))
}
Optimal_Reg = function(y1,y0,x){
  
  
  subject_size=length(y1)
  
  #create the outcome vector
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=-y0 #make sure the sign is correct: it should be negative
  y_temp[2*(1:subject_size)]=y1
  
  treated_ind=rep(c(0,1),subject_size)
  control_ind=rep(c(-1,0),subject_size)
  
  x_temp=matrix(0,nrow=subject_size*2,ncol(x))
  x_temp[2*(1:subject_size)-1,]=-x
  x_temp[2*(1:subject_size),]=x
  
  x_temp=cbind(control_ind,treated_ind,x_temp)
  
  
  XDX=matrix(0,nrow=ncol(x)+2,ncol=ncol(x)+2)
  XDY=rep(0,ncol(x)+2)
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    #normalized_cov=so_t
    
    
     XDX = XDX + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% x_temp[index_start:index_end,]/subject_size
    
    XDY = XDY + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% y_temp[index_start:index_end]/subject_size
    
    
  }
  #solve for optimal coefficients
  optimal_b = solve(XDX,XDY)
  print(optimal_b)
  #optimal predictions
  optimal_pred = x_temp %*% optimal_b
  
  #z_optimal
  z_optimal = y_temp-optimal_pred
  
  return(z_optimal)
  
}
Noharm_Reg = function(y1,y0,x){
  
  
  subject_size=length(y1)
  
  #create the outcome vector
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=-y0 #make sure the sign is correct: it should be negative
  y_temp[2*(1:subject_size)]=y1
  
  
  x_temp=matrix(0,nrow=subject_size*2,1)
  x_temp[2*(1:subject_size)-1,]=-x[2*(1:subject_size)-1]
  x_temp[2*(1:subject_size),]=x[2*(1:subject_size)]
  
  
  
  XDX=matrix(0,nrow=1,ncol=1)
  XDY=rep(0,1)
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    #normalized_cov=so_t
    
    
    XDX = XDX + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% x_temp[index_start:index_end,]/subject_size
    
    XDY = XDY + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% y_temp[index_start:index_end]/subject_size
    
    
  }
  #solve for optimal coefficients
  optimal_b = solve(XDX,XDY)
  
  #optimal predictions
  optimal_pred = x_temp %*% optimal_b
  
  #z_optimal
  z_optimal = y_temp-optimal_pred
  
  return(z_optimal)
  
}
Noharm_Optimal_Reg = function(y1,y0,x){
  
  
  subject_size=length(y1)
  
  #create the outcome vector
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=-y0 #make sure the sign is correct: it should be negative
  y_temp[2*(1:subject_size)]=y1
  
  
  x_temp=matrix(0,nrow=subject_size*2,3)
  x_temp[2*(1:subject_size)-1,1]=-1
  x_temp[2*(1:subject_size),2]=1
  x_temp[2*(1:subject_size)-1,3]=-x[2*(1:subject_size)-1]
  x_temp[2*(1:subject_size),3]=x[2*(1:subject_size)]
  
  XDX=matrix(0,nrow=3,ncol=3)
  XDY=rep(0,3)
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    #normalized_cov=so_t
    
    
    XDX = XDX + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% x_temp[index_start:index_end,]/subject_size
    
    XDY = XDY + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% y_temp[index_start:index_end]/subject_size
    
    
  }
  
  #solve for optimal coefficients
  optimal_b = solve(XDX,XDY)
  print('Optimal Coefficient is')
  print(optimal_b)
  #optimal predictions
  optimal_pred = x_temp %*% optimal_b
  
  #z_optimal
  z_optimal = y_temp-optimal_pred
  
  return(z_optimal)
  
}

Optimal_Reg_Est = function(y,d1,d0,x,pi0,pi1){
  
  
  subject_size=length(y)
  #create the outcome vector
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=-y*d0/pi0  #make sure the sign is correct: it should be negative
  y_temp[2*(1:subject_size)]=y*d1/pi1 
  
  treated_ind=rep(c(0,1),subject_size)
  control_ind=rep(c(-1,0),subject_size)
  
  x_temp=matrix(0,nrow=subject_size*2,ncol(x))
  x_temp[2*(1:subject_size)-1,]=-x
  x_temp[2*(1:subject_size),]=x
  x_temp=cbind(control_ind,treated_ind,x_temp)
  
  x_pred=matrix(0,nrow=subject_size*2,ncol(x))
  x_pred[2*(1:subject_size)-1,]=x
  x_pred[2*(1:subject_size),]=x  
  x_pred=cbind(-control_ind,treated_ind,x_pred)
  
  
  ###fixed denominator
  
  XDX2=matrix(0,nrow=ncol(x)+2,ncol=ncol(x)+2)
  XDY2=rep(0,ncol(x)+2)
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]] 
    normalized_cov=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    #normalized_cov=so_t
    #extract variance bound
    # cov_bound_t = so_AS2[[c]]
    # 
    # #second order probabilities
    # p_t = so_t + (fo_t) %*% t(fo_t) 
    # p_t=round(p_t,10)
    # 
    # cov_bound_zeros_t = (abs(cov_bound_t)<1e-2)
    # 
    # p_denomiator_t = p_t + cov_bound_zeros_t
    # cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    XDX2 = XDX2 + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% x_temp[index_start:index_end,]/subject_size
    
    XDY2 = XDY2 + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% y_temp[index_start:index_end]/subject_size
    
    
  }
  
  XDX2=XDX2

  #solve for optimal coefficients
  optimal_b2 = solve(XDX2,XDY2)
  print(optimal_b2)
  #optimal predictions
  optimal_pred2 = x_pred %*% optimal_b2

  y_p1_2=optimal_pred2[2*(1:subject_size)]
  y_p0_2=optimal_pred2[2*(1:subject_size)-1]
  
  #DR
  avg_effect1_2 = mean(y_p1_2) +  mean((y-y_p1_2) * d1/pi1)
  avg_effect0_2 = mean(y_p0_2) +  mean((y-y_p0_2) * d0/pi0)    
  
  effect_est2=avg_effect1_2-avg_effect0_2
  
  
  return(list(effect_est2,optimal_b2,y_p1_2,y_p0_2,optimal_b2))
  
}

Logit_full=function(y,d1,d0,weights,x,pi0,pi1,mode){
  
  print(mode)
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
    print(coef)
    avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
    avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)
  }
  return(list(avg_effect1-avg_effect0,y_p1,y_p0))
  #return(c(avg_effect1-avg_effect0,coef))
}

theoretical_variance_full=function(y1,y0,x,fo,so,contrast,dict_group){
  
  
  subject_size=length(y1)
  
  #create the potential outcome vectors
  y=rep(0,subject_size*2)
  y[2*(1:subject_size)-1]=y0 
  y[2*(1:subject_size)]=y1
  
  fo_vec=unlist(fo)
  y_length = length(y)
  x_predict = x[rep(1:nrow(x),each=2),]
  x_predict=as.matrix(x_predict)
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
  reg_OLS_joint =  lm(y~d0+d1+x_predict-1,weights=fo_vec)
  temp_OLS =  (y-predict(reg_OLS_joint,x_OLS_joint)) * fo_vec
  z_OLS_joint=  temp_OLS * as.matrix(x_OLS_joint) %*% solve( t(as.matrix(x_OLS_joint)) %*%diag(fo_vec)%*%as.matrix(x_OLS_joint)*(1/subject_size))
  z_OLS_joint = z_OLS_joint %*% contrast_OLS
  
  ####WLS_joint####
  x_WLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_WLS_joint =  lm(y~d0+d1+x_predict-1)
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
  print(reg_Logit_joint_w$coefficients)
  z_Logit_joint_w = y-predict(reg_Logit_joint_w,x_Logit_joint_w,type='response')
  z_Logit_joint_w = z_Logit_joint_w * contrast_vec
  
  ####OLS_DR####
  x_OLS_DR = data.frame(cbind(d0,d1,x_predict))
  reg_OLS_DR =  lm(y~d0+d1+x_predict-1,weights=fo_vec)
  temp_OLS_DR =  y-predict(reg_OLS_DR,x_OLS_DR)
  z_OLS_DR = temp_OLS_DR *  contrast_vec
  
  ####Optimal Coef####
  z_optimal = Optimal_Reg(y1,y0,x)
  
  
  
  ####Optimal Coef
  z_optimal_Logit_MS = Optimal_Logit_B(y1,y0,x,Maxiter=1000,Optim=F)
  
  
  print('noharm')
  ###No_harm_WLS####
  
  WLS_pred=predict(reg_WLS_joint,x_WLS_joint)
  
  z_Noharm_WLS=Noharm_Reg(y1,y0,  WLS_pred)
  
  print('noharm')
  
  ##No_harm_Wlogit
  WLogit_pred=  predict(reg_Logit_joint_w,x_Logit_joint_w,type='response')
  z_Noharm_WLogit=Noharm_Reg(y1,y0,WLogit_pred)  
  
  ####Optimal Coef
  #z_optimal_Logit_OS = Optimal_Logit_A(y1,y0,x,Maxiter=10,Optim=F)
  
  WOLS_pred=  predict(reg_OLS_DR,x_OLS_DR,type='response')
  z_Noharm_WOLS=Noharm_Reg(y1,y0,WOLS_pred)  
  
  output=rep(0,13)
  
  names(output)=c('var_HT','var_HA','var_OLS_joint','var_WLS_joint','var_Logit_joint_u','var_Logit_joint_w','var_OLS_DR','var_Linear_Optimal','var_Logit_Optimal_MS','no_harm_wls','no_harm_logit','var_Logit_Optimal_OS','no_harm_OLS')
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    #normalized_cov=so_t
    z_HT_t = z_HT[index_start:index_end]
    z_HA_t = z_HA[index_start:index_end]
    z_OLS_joint_t = z_OLS_joint[index_start:index_end]
    z_WLS_joint_t = z_WLS_joint[index_start:index_end]
    z_Logit_joint_u_t =z_Logit_joint_u[index_start:index_end]
    z_Logit_joint_w_t = z_Logit_joint_w[index_start:index_end]
    z_OLS_DR_t=z_OLS_DR[index_start:index_end]
    z_optimal_t = z_optimal[index_start:index_end]
    #z_optimal_Logit_OS_t = z_optimal_Logit_OS[index_start:index_end]
    z_optimal_Logit_MS_t = z_optimal_Logit_MS[index_start:index_end]
    z_Noharm_WLS_t=z_Noharm_WLS[index_start:index_end]
    z_Noharm_WLogit_t=z_Noharm_WLogit[index_start:index_end]  
    z_Noharm_WOLS_t=z_Noharm_WOLS[index_start:index_end]      
    
    output[1]=output[1] +  t(z_HT_t) %*% normalized_cov %*% z_HT_t/(subject_size^2)
    output[2]=output[2] +  t(z_HA_t) %*% normalized_cov %*% z_HA_t/(subject_size^2)
    output[3]=output[3] +  t(z_OLS_joint_t) %*% normalized_cov %*% z_OLS_joint_t/(subject_size^2)
    output[4]=output[4] +  t(z_WLS_joint_t) %*% normalized_cov %*% z_WLS_joint_t/(subject_size^2)
    output[5]=output[5] +  t(z_Logit_joint_u_t) %*% normalized_cov %*% z_Logit_joint_u_t/(subject_size^2)
    output[6]=output[6] +  t(z_Logit_joint_w_t) %*% normalized_cov %*%  z_Logit_joint_w_t/(subject_size^2)
    output[7]=output[7] +  t(z_OLS_DR_t) %*% normalized_cov %*%  z_OLS_DR_t/(subject_size^2)
    output[8]=output[8] +  t(z_optimal_t) %*% normalized_cov %*%  z_optimal_t/(subject_size^2)
    #output[9]=output[9] +  t(z_optimal_Logit_OS_t) %*% normalized_cov %*%  z_optimal_Logit_OS_t/(subject_size^2)
    output[10]=output[10] +  t(z_optimal_Logit_MS_t) %*% normalized_cov %*%  z_optimal_Logit_MS_t/(subject_size^2)
    output[11]=output[11] +  t(z_Noharm_WLS_t) %*% normalized_cov %*%  z_Noharm_WLS_t/(subject_size^2)
    output[12]=output[12] +  t(z_Noharm_WLogit_t) %*% normalized_cov %*%  z_Noharm_WLogit_t/(subject_size^2)    
    output[13]=output[13] +  t(z_Noharm_WOLS_t) %*% normalized_cov %*%  z_Noharm_WOLS_t/(subject_size^2)    
    
  }
  
  return(output)
}
theoretical_variance_full2=function(y1,y0,x,fo,so,contrast,dict_group){
  
  
  subject_size=length(y1)
  
  #create the potential outcome vectors
  y=rep(0,subject_size*2)
  y[2*(1:subject_size)-1]=y0 
  y[2*(1:subject_size)]=y1
  
  fo_vec=unlist(fo)
  y_length = length(y)
  x_predict = x[rep(1:nrow(x),each=2),]
  x_predict=as.matrix(x_predict)
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
  reg_OLS_joint =  lm(y~d0+d1+x_predict-1,weights=fo_vec)
  OLS_pred=predict(reg_OLS_joint,x_OLS_joint)
  temp_OLS =  (y-predict(reg_OLS_joint,x_OLS_joint)) * fo_vec
  z_OLS_joint=  temp_OLS * as.matrix(x_OLS_joint) %*% solve( t(as.matrix(x_OLS_joint)) %*%diag(fo_vec)%*%as.matrix(x_OLS_joint)*(1/subject_size))
  z_OLS_joint = z_OLS_joint %*% contrast_OLS
  
  ####WLS_joint####
  x_WLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_WLS_joint =  lm(y~d0+d1+x_predict-1)
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
  print(reg_Logit_joint_w$coefficients)
  z_Logit_joint_w = y-predict(reg_Logit_joint_w,x_Logit_joint_w,type='response')
  z_Logit_joint_w = z_Logit_joint_w * contrast_vec
  
  ####OLS_DR####
  x_OLS_DR = data.frame(cbind(d0,d1,x_predict))
  reg_OLS_DR =  lm(y~d0+d1+x_predict-1,weights=fo_vec)
  temp_OLS_DR =  y-predict(reg_OLS_DR,x_OLS_DR)
  z_OLS_DR = temp_OLS_DR *  contrast_vec
  
  ####Optimal Coef####
  z_optimal = Optimal_Reg(y1,y0,x)
  
  
  
  ####Optimal Coef
  # z_optimal_Logit_MS = Optimal_Logit_B(y1,y0,x,Maxiter=1000,Optim=F)
  
  
  print('noharm')
  ###No_harm_WLS####
  
  WLS_pred=predict(reg_WLS_joint,x_WLS_joint)
  
  z_Noharm_WLS=Noharm_Reg(y1,y0,  WLS_pred)
  
  print('noharm')
  
  ##No_harm_Wlogit
  ULogit_pred=  predict(reg_Logit_joint_u,x_Logit_joint_u,type='response')
  z_Noharm_WLogit=Noharm_Reg(y1,y0,ULogit_pred)  
  
  ####Optimal Coef
  #z_optimal_Logit_OS = Optimal_Logit_A(y1,y0,x,Maxiter=10,Optim=F)
  
  z_Noharm_WLS_Optimal=Noharm_Optimal_Reg(y1,y0,  OLS_pred)
  z_Noharm_Logit_Optimal=Noharm_Optimal_Reg(y1,y0,  ULogit_pred)
  
  output=rep(0,13)
  
  names(output)=c('var_HT','var_HA','var_OLS_joint','var_WLS_joint','var_Logit_joint_u','var_Logit_joint_w','var_OLS_DR','var_Linear_Optimal','var_Logit_Optimal_MS','no_harm_wls','no_harm_logit','var_Logit_Optimal_OS','no_harm_optimal_WLS')
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    #normalized_cov=so_t
    z_HT_t = z_HT[index_start:index_end]
    z_HA_t = z_HA[index_start:index_end]
    z_OLS_joint_t = z_OLS_joint[index_start:index_end]
    z_WLS_joint_t = z_WLS_joint[index_start:index_end]
    z_Logit_joint_u_t =z_Logit_joint_u[index_start:index_end]
    z_Logit_joint_w_t = z_Logit_joint_w[index_start:index_end]
    z_OLS_DR_t=z_OLS_DR[index_start:index_end]
    z_optimal_t = z_optimal[index_start:index_end]
    #z_optimal_Logit_OS_t = z_optimal_Logit_OS[index_start:index_end]
    #z_optimal_Logit_MS_t = z_optimal_Logit_MS[index_start:index_end]
    z_Noharm_WLS_t=z_Noharm_WLS[index_start:index_end]
    z_Noharm_WLogit_t=z_Noharm_WLogit[index_start:index_end]  
    z_Noharm_WLS_Optimal_t=z_Noharm_WLS_Optimal[index_start:index_end]      
    z_Noharm_Logit_t = z_Noharm_Logit_Optimal[index_start:index_end]   
    output[1]=output[1] +  t(z_HT_t) %*% normalized_cov %*% z_HT_t/(subject_size^2)
    output[2]=output[2] +  t(z_HA_t) %*% normalized_cov %*% z_HA_t/(subject_size^2)
    output[3]=output[3] +  t(z_OLS_joint_t) %*% normalized_cov %*% z_OLS_joint_t/(subject_size^2)
    output[4]=output[4] +  t(z_WLS_joint_t) %*% normalized_cov %*% z_WLS_joint_t/(subject_size^2)
    output[5]=output[5] +  t(z_Logit_joint_u_t) %*% normalized_cov %*% z_Logit_joint_u_t/(subject_size^2)
    output[6]=output[6] +  t(z_Logit_joint_w_t) %*% normalized_cov %*%  z_Logit_joint_w_t/(subject_size^2)
    output[7]=output[7] +  t(z_OLS_DR_t) %*% normalized_cov %*%  z_OLS_DR_t/(subject_size^2)
    output[8]=output[8] +  t(z_optimal_t) %*% normalized_cov %*%  z_optimal_t/(subject_size^2)
    #output[10]=output[10] +  t(z_optimal_Logit_MS_t) %*% normalized_cov %*%  z_optimal_Logit_MS_t/(subject_size^2)
    output[9]=output[9] +  t(z_Noharm_WLS_t) %*% normalized_cov %*%  z_Noharm_WLS_t/(subject_size^2)
    output[10]=output[10] +  t(z_Noharm_WLogit_t) %*% normalized_cov %*%  z_Noharm_WLogit_t/(subject_size^2)    
    output[11]=output[11] +  t(z_Noharm_WLS_Optimal_t) %*% normalized_cov %*%  z_Noharm_WLS_Optimal_t/(subject_size^2)    
    output[12]=output[12] +  t(z_Noharm_Logit_t) %*% normalized_cov %*%  z_Noharm_Logit_t/(subject_size^2)
    
  }
  
  return(output)
}
theoretical_variance_bound_full2_AS2=function(y1,y0,x,fo,so,contrast,dict_group){
  
  
  subject_size=length(y1)
  
  #create the potential outcome vectors
  y=rep(0,subject_size*2)
  y[2*(1:subject_size)-1]=y0 
  y[2*(1:subject_size)]=y1
  
  fo_vec=unlist(fo)
  y_length = length(y)
  x_predict = x[rep(1:nrow(x),each=2),]
  x_predict=as.matrix(x_predict)
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
  reg_OLS_joint =  lm(y~d0+d1+x_predict-1,weights=fo_vec)
  temp_OLS =  (y-predict(reg_OLS_joint,x_OLS_joint)) * fo_vec
  z_OLS_joint=  temp_OLS * as.matrix(x_OLS_joint) %*% solve( t(as.matrix(x_OLS_joint)) %*%diag(fo_vec)%*%as.matrix(x_OLS_joint)*(1/subject_size))
  z_OLS_joint = z_OLS_joint %*% contrast_OLS
  
  ####WLS_joint####
  x_WLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_WLS_joint =  lm(y~d0+d1+x_predict-1)
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
  print(reg_Logit_joint_w$coefficients)
  z_Logit_joint_w = y-predict(reg_Logit_joint_w,x_Logit_joint_w,type='response')
  z_Logit_joint_w = z_Logit_joint_w * contrast_vec
  
  ####OLS_DR####
  x_OLS_DR = data.frame(cbind(d0,d1,x_predict))
  reg_OLS_DR =  lm(y~d0+d1+x_predict-1,weights=fo_vec)
  temp_OLS_DR =  y-predict(reg_OLS_DR,x_OLS_DR)
  z_OLS_DR = temp_OLS_DR *  contrast_vec
  
  ####Optimal Coef####
  z_optimal = Optimal_Reg(y1,y0,x)
  
  
  
  ####Optimal Coef
  # z_optimal_Logit_MS = Optimal_Logit_B(y1,y0,x,Maxiter=1000,Optim=F)
  
  
  print('noharm')
  ###No_harm_WLS####
  
  WLS_pred=predict(reg_WLS_joint,x_WLS_joint)
  
  z_Noharm_WLS=Noharm_Reg(y1,y0,  WLS_pred)
  
  print('noharm')
  
  ##No_harm_Wlogit
  ULogit_pred=  predict(reg_Logit_joint_u,x_Logit_joint_u,type='response')
  z_Noharm_Logit=Noharm_Reg(y1,y0,ULogit_pred)  
  
  ####Optimal Coef
  #z_optimal_Logit_OS = Optimal_Logit_A(y1,y0,x,Maxiter=10,Optim=F)
  
  z_Noharm_WLS_Optimal=Noharm_Optimal_Reg(y1,y0,  WLS_pred)
  z_Noharm_Logit_Optimal=Noharm_Optimal_Reg(y1,y0,  ULogit_pred)
  
  output=rep(0,13)
  
  names(output)=c('var_HT','var_HA','var_OLS_joint','var_WLS_joint','var_Logit_joint_u','var_Logit_joint_w','var_OLS_DR','var_Linear_Optimal','var_Logit_Optimal_MS','no_harm_wls','no_harm_logit','var_Logit_Optimal_OS','no_harm_optimal_WLS')
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    
    normalized_cov_t=round(normalized_cov_t,10)
    cov_bound_t = so_AS2[[c]]
    #normalized_cov=so_t
    z_HT_t = z_HT[index_start:index_end]
    z_HA_t = z_HA[index_start:index_end]
    z_OLS_joint_t = z_OLS_joint[index_start:index_end]
    z_WLS_joint_t = z_WLS_joint[index_start:index_end]
    z_Logit_joint_u_t =z_Logit_joint_u[index_start:index_end]
    z_Logit_joint_w_t = z_Logit_joint_w[index_start:index_end]
    z_OLS_DR_t=z_OLS_DR[index_start:index_end]
    z_optimal_t = z_optimal[index_start:index_end]
    #z_optimal_Logit_OS_t = z_optimal_Logit_OS[index_start:index_end]
    #z_optimal_Logit_MS_t = z_optimal_Logit_MS[index_start:index_end]
    z_Noharm_WLS_t=z_Noharm_WLS[index_start:index_end]
    z_Noharm_Logit_t=z_Noharm_Logit[index_start:index_end]  
    z_Noharm_WLS_Optimal_t=z_Noharm_WLS_Optimal[index_start:index_end]      
    z_Noharm_Logit_Optimal_t = z_Noharm_Logit_Optimal[index_start:index_end]   
    output[1]=output[1] +  t(z_HT_t) %*% cov_bound_t %*% z_HT_t/(subject_size^2)
    output[2]=output[2] +  t(z_HA_t) %*% cov_bound_t %*% z_HA_t/(subject_size^2)
    output[3]=output[3] +  t(z_OLS_joint_t) %*% cov_bound_t %*% z_OLS_joint_t/(subject_size^2)
    output[4]=output[4] +  t(z_WLS_joint_t) %*% cov_bound_t %*% z_WLS_joint_t/(subject_size^2)
    output[5]=output[5] +  t(z_Logit_joint_u_t) %*% cov_bound_t %*% z_Logit_joint_u_t/(subject_size^2)
    output[6]=output[6] +  t(z_Logit_joint_w_t) %*% cov_bound_t %*%  z_Logit_joint_w_t/(subject_size^2)
    output[7]=output[7] +  t(z_OLS_DR_t) %*% cov_bound_t %*%  z_OLS_DR_t/(subject_size^2)
    output[8]=output[8] +  t(z_optimal_t) %*% cov_bound_t %*%  z_optimal_t/(subject_size^2)
    #output[10]=output[10] +  t(z_optimal_Logit_MS_t) %*% normalized_cov %*%  z_optimal_Logit_MS_t/(subject_size^2)
    output[9]=output[9] +  t(z_Noharm_WLS_t) %*% cov_bound_t %*%  z_Noharm_WLS_t/(subject_size^2)
    output[10]=output[10] +  t(z_Noharm_Logit_t) %*% cov_bound_t %*%  z_Noharm_Logit_t/(subject_size^2)    
    output[11]=output[11] +  t(z_Noharm_WLS_Optimal_t) %*% cov_bound_t %*%  z_Noharm_WLS_Optimal_t/(subject_size^2)    
    output[12]=output[12] +  t(z_Noharm_Logit_Optimal_t) %*% cov_bound_t %*%  z_Noharm_Logit_Optimal_t/(subject_size^2)
    
  }
  
  return(output)
}

theoretical_variance_bound_full_AS2=function(y1,y0,x,fo,so,contrast,dict_group){
  
  
  subject_size=length(y1)
  
  #create the potential outcome vectors
  y=rep(0,subject_size*2)
  y[2*(1:subject_size)-1]=y0 
  y[2*(1:subject_size)]=y1
  
  fo_vec=unlist(fo)
  y_length = length(y)
  x_predict = x[rep(1:nrow(x),each=2),]
  d1 = rep(c(0,1),y_length/2)
  d0 = 1-d1
  
  contrast_OLS = c(contrast,rep(0,ncol(x)))
  contrast_vec= rep(contrast,length(y)/2)
  
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
  
  ####OLS_DR####
  x_OLS_DR = data.frame(cbind(d0,d1,x_predict))
  reg_OLS_DR =  lm(y~d0+d1+x_predict-1,weights=fo_vec)
  temp_OLS_DR =  y-predict(reg_OLS_DR,x_OLS_DR,type='response')
  z_OLS_DR = temp_OLS_DR *  contrast_vec
  
  
  
  
  
  ####Optimal Coef####
  z_optimal = Optimal_Reg(y1,y0,x)
  
  ####Optimal Coef
#  z_optimal_Logit_OS = Optimal_Logit_A(y1,y0,x,Maxiter=1,Optim=T)
  
  ####Optimal Coef
  # z_optimal_Logit_MS = Optimal_Logit_B(y1,y0,x,Maxiter=10000,Optim=F)
  
  
  
  ###No_harm_WLS####
  
   WLS_pred=predict(reg_WLS_joint,x_WLS_joint)

   
   z_Noharm_WLS=Noharm_Reg(y1,y0,  WLS_pred)
   
   
   ##No_harm_Wlogit
   WLogit_pred=  predict(reg_Logit_joint_w,x_Logit_joint_w,type='response')
   

   z_Noharm_WLogit=Noharm_Reg(y1,y0,WLogit_pred)  
  
  
  ####Optimal Coef
  #z_optimal_Logit_OS = Optimal_Logit_A(y1,y0,x,Maxiter=10,Optim=F)
  
  WOLS_pred=  predict(reg_OLS_DR,x_OLS_DR,type='response')
  z_Noharm_WOLS=Noharm_Reg(y1,y0,WOLS_pred)  
  
  output=rep(0,13)
  
  names(output)=c('var_HT','var_HA','var_OLS_joint','var_WLS_joint','var_Logit_joint_u','var_Logit_joint_w','var_OLS_DR','var_Linear_Optimal','var_Logit_Optimal_MS','no_harm_wls','no_harm_logit','var_Logit_Optimal_OS','no_harm_OLS')
  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    
    normalized_cov_t=round(normalized_cov_t,10)
    cov_bound_t = so_AS2[[c]]
    
    
    z_HT_t = z_HT[index_start:index_end]
    z_HA_t = z_HA[index_start:index_end]
    z_OLS_joint_t = z_OLS_joint[index_start:index_end]
    z_WLS_joint_t = z_WLS_joint[index_start:index_end]
    z_Logit_joint_u_t =z_Logit_joint_u[index_start:index_end]
    z_Logit_joint_w_t = z_Logit_joint_w[index_start:index_end]
    z_OLS_DR_t=z_OLS_DR[index_start:index_end]
    z_optimal_t = z_optimal[index_start:index_end]
    #z_optimal_Logit_OS_t = z_optimal_Logit_OS[index_start:index_end]
    #z_optimal_Logit_MS_t = z_optimal_Logit_MS[index_start:index_end]
    z_Noharm_WLS_t=z_Noharm_WLS[index_start:index_end]
    z_Noharm_WLogit_t=z_Noharm_WLogit[index_start:index_end]  
    z_Noharm_WOLS_t=z_Noharm_WOLS[index_start:index_end]      
    
    
    output[1]=output[1] +  t(z_HT_t) %*% cov_bound_t %*% z_HT_t/(subject_size^2)
    output[2]=output[2] +  t(z_HA_t) %*% cov_bound_t %*% z_HA_t/(subject_size^2)
    output[3]=output[3] +  t(z_OLS_joint_t) %*% cov_bound_t %*% z_OLS_joint_t/(subject_size^2)
    output[4]=output[4] +  t(z_WLS_joint_t) %*% cov_bound_t %*% z_WLS_joint_t/(subject_size^2)
    output[5]=output[5] +  t(z_Logit_joint_u_t) %*% cov_bound_t %*% z_Logit_joint_u_t/(subject_size^2)
    output[6]=output[6] +  t(z_Logit_joint_w_t) %*% cov_bound_t %*%  z_Logit_joint_w_t/(subject_size^2)
    output[7]=output[7] +  t(z_OLS_DR_t) %*% cov_bound_t %*%  z_OLS_DR_t/(subject_size^2)
    output[8]=output[8] +  t(z_optimal_t) %*% cov_bound_t %*%  z_optimal_t/(subject_size^2)
    #output[9]=output[9] +  t(z_optimal_Logit_OS_t) %*% cov_bound_t %*%  z_optimal_Logit_OS_t/(subject_size^2)
    #output[10]=output[10] +  t(z_optimal_Logit_MS_t) %*% cov_bound_t %*%  z_optimal_Logit_MS_t/(subject_size^2)
    output[11]=output[11] +  t(z_Noharm_WLS_t) %*% cov_bound_t %*%  z_Noharm_WLS_t/(subject_size^2)
    output[12]=output[12] +  t(z_Noharm_WLogit_t) %*% cov_bound_t %*%  z_Noharm_WLogit_t/(subject_size^2)    
    output[13]=output[13] +  t(z_Noharm_WOLS_t) %*% cov_bound_t %*%  z_Noharm_WOLS_t/(subject_size^2)    
    
  }
  
  return(output)
}

AS_bound=function(normalized_cov){
  
  #rounds to 10th digit, this is arbitrary
  normalized_cov = round(normalized_cov,10)
  minus_ones= (abs(normalized_cov+1)<0.001)
  add_on = apply(minus_ones,2,sum)
  
  
  #AS formula
  AS_bound = normalized_cov - minus_ones * normalized_cov+ diag(add_on)
  return(AS_bound)
  
  
}
share_pain_bound_old=function(minus_ones_t){
  
  
  bound_matrix=matrix(0,nrow(minus_ones_t),nrow(minus_ones_t))
  i=0
  while (sum( abs(minus_ones_t+1)<1e-3)>0){
    
    if (i%%100==0){
      print(paste0(i, ' iterations have passed' ))
      print(sum( abs(minus_ones_t+1)<1e-5))
      
    }
    #the center guy
    degree=-apply(minus_ones_t,1,sum)
    # degree=degree * (degree>0) + nrow(minus_ones_t) * (degree==0)
    
    # -1's
    index_temp=which( (minus_ones_t[which.max(degree),]+1)<1e-3)
    
    #
    index_temp=c(which.max(degree),index_temp)
    
    
    bound_matrix[index_temp,index_temp] =bound_matrix[index_temp,index_temp]+1
    minus_ones_t[index_temp,index_temp]=0
    
    i=i+1
  }
  
  return(bound_matrix)
}

share_pain_bound=function(minus_ones_t){
  
  
  bound_matrix=matrix(0,nrow(minus_ones_t),nrow(minus_ones_t))
  i=0
  while (sum( abs(minus_ones_t+1)<1e-3)>0){
    
    if (i%%100==0){
      print(paste0(i, ' iterations have passed' ))
      print(sum( abs(minus_ones_t+1)<1e-5))
      
    }
    #the center guy
    degree=-apply(minus_ones_t,1,sum)
    # degree=degree * (degree>0) + nrow(minus_ones_t) * (degree==0)
    
    # -1's
    index_temp=which( (minus_ones_t[which.max(degree),]+1)<1e-3)
    
    #
    index_temp1=c(which.max(degree),index_temp)
    
    # how many -1s are we filling?
    eff1 = sum(abs(minus_ones_t[index_temp1,index_temp1])) /  length(index_temp1)^2
    
    index_temp2=c(index_temp,which.max(degree))
    index_temp2=unique( (index_temp2-1)%/%2+1)
    
    index_temp2=c(index_temp2*2-1,index_temp2*2)
    eff2 = sum(abs(minus_ones_t[index_temp2,index_temp2])) /  length(index_temp2)^2
    print(paste0('eff1= ',eff1))
    print(paste0('eff2= ',eff2))
    if (eff2>eff1){
      index_temp=index_temp2
    }else{
      index_temp=index_temp1
    }
    
    bound_matrix[index_temp,index_temp] =bound_matrix[index_temp,index_temp]+1
    minus_ones_t[index_temp,index_temp]=0
    
    i=i+1
  }
  
  return(bound_matrix)
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
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,2]
    index_end=dict_group[i,3]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    
    #select needed component
    y_t = y[index_start:index_end]  
    
    d1_t=d1[index_start:index_end]  
    d0_t=d0[index_start:index_end]  
    output_t=rep(y_t,each=2)    
    output_t[(1:length(output_t))%%2==0]=output_t[(1:length(output_t))%%2==0]*d1_t
    output_t[(1:length(output_t))%%2==1]=output_t[(1:length(output_t))%%2==1]*d0_t   
    
    contrast_vec=rep(contrast,length(output_t)/2)
    
    
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
    #print('HT')
    #print(i)
    
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    var_est = var_est + (contrast_vec*output_t) %*% cov_bound_HT_t %*% (output_t*contrast_vec)
    
  }
  
  var_est=var_est/(subject_size)^2
  return(var_est)
  
  
  
}
HT_var_full_AS2=function(y,d1,d0,sample_size,fo,so,contrast,dict_group){
  
  
  
  #computing variance
  num_component=nrow(dict_group)
  subject_size=length(y)
  
  var_est=0
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,2]
    index_end=dict_group[i,3]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    
    #select needed component
    y_t = y[index_start:index_end]  
    
    d1_t=d1[index_start:index_end]  
    d0_t=d0[index_start:index_end]  
    output_t=rep(y_t,each=2)    
    output_t[(1:length(output_t))%%2==0]=output_t[(1:length(output_t))%%2==0]*d1_t
    output_t[(1:length(output_t))%%2==1]=output_t[(1:length(output_t))%%2==1]*d0_t   
    
    contrast_vec=rep(contrast,length(output_t)/2)
    
    
    
    cov_bound_t = so_AS2[[c]]
    
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    p_zeros_t=(p_t==0)
    cov_bound_zeros_t = (abs(cov_bound_t)<1e-2)
    #print('Check if denomiator matrix is correct')
    if (min(cov_bound_zeros_t-p_zeros_t)<0){
      print('The denomiator matrix is wrong')
    }
    #print(min(cov_bound_zeros_t-p_zeros_t)) #should expect 0
    #print('HT')
    #print(i)
    
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
  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,2]
    index_end=dict_group[i,3]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    
    #select needed component
    y_t = y[index_start:index_end]  
    
    output_t=rep(y_t,each=2)    
    output_t[(1:length(output_t))%%2==0]=fitted_1[index_start:index_end]
    output_t[(1:length(output_t))%%2==1]=fitted_0[index_start:index_end]   
    
    contrast_vec=rep(contrast,length(output_t)/2)
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
HA_var_full_AS2=function(y,d1,d0,pi,sample_size,fo,so,contrast,dict_group){
  
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
  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,2]
    index_end=dict_group[i,3]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    
    #select needed component
    y_t = y[index_start:index_end]  
    
    output_t=rep(y_t,each=2)    
    output_t[(1:length(output_t))%%2==0]=fitted_1[index_start:index_end]
    output_t[(1:length(output_t))%%2==1]=fitted_0[index_start:index_end]   
    
    contrast_vec=rep(contrast,length(output_t)/2)
    
    cov_bound_t = so_AS2[[c]]
    
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    p_zeros_t=(p_t==0)
    cov_bound_zeros_t = (abs(cov_bound_t)<1e-2)
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
    #print('finished')
    for (i in 1:nrow(dict_group)){
      
      c=dict_group[i,1]
      #print(c)
      index_start=dict_group[i,4]
      index_end=dict_group[i,5]    
      
      #select needed component
      
      output_t=z_OLS_joint_hat[index_start:index_end]
      #print(length(output_t))      
      fo_t=fo[[c]]
      so_t=so[[c]]
      
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
      #print(dim(cov_bound_HT_t))
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
    
    for (i in 1:nrow(dict_group)){
      
      c=dict_group[i,1]
      index_start=dict_group[i,4]
      index_end=dict_group[i,5]    
      
      #select needed component
      
      output_t=z_WLS_joint_hat[index_start:index_end]
      
      fo_t=fo[[c]]
      so_t=so[[c]]
      
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
  }
  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,2]
    index_end=dict_group[i,3]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    
    #select needed component
    #y_t = y[index_start:index_end]  
    
    
    output_t=rep(0,(index_end-index_start+1)*2)
    output_t[(1:length(output_t))%%2==0]=fitted_1[index_start:index_end]
    output_t[(1:length(output_t))%%2==1]=fitted_0[index_start:index_end]   
    #print(length(output_t))
    print(output_t[1:10])
    contrast_vec=rep(contrast,length(output_t)/2)
    
    #print('Logit')
    #print(i)
    
    
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
    #print(dim(cov_bound_HT_t))    
    contrast_vec=rep(contrast,length(output_t)/2)
    
    var_est = var_est + (contrast_vec*output_t) %*% cov_bound_HT_t %*% (output_t*contrast_vec)/(subject_size)^2
    
  }
  
  
  return(var_est)
  #return(c(avg_effect1-avg_effect0,coef))
}

Logit_var_full_AS2=function(y,d1,d0,weights,x,mode,fo,so,contrast,dict_group){
  
  
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
  }
  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,2]
    index_end=dict_group[i,3]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    
    #select needed component
    #y_t = y[index_start:index_end]  
    
    
    output_t=rep(0,(index_end-index_start+1)*2)
    output_t[(1:length(output_t))%%2==0]=fitted_1[index_start:index_end]
    output_t[(1:length(output_t))%%2==1]=fitted_0[index_start:index_end]   
    #print(length(output_t))
    print(output_t[1:10])
    contrast_vec=rep(contrast,length(output_t)/2)
    
    #print('Logit')
    #print(i)
    
    
    
    cov_bound_t = so_AS2[[c]]
    
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,100)
    
    p_zeros_t=(p_t==0)
    cov_bound_zeros_t = (abs(cov_bound_t)<1e-2)
    if (min(cov_bound_zeros_t-p_zeros_t)<0){
      print('The denomiator matrix is wrong')
    }
    
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    #print(dim(cov_bound_HT_t))    
    contrast_vec=rep(contrast,length(output_t)/2)
    
    var_est = var_est + (contrast_vec*output_t) %*% cov_bound_HT_t %*% (output_t*contrast_vec)/(subject_size)^2
    
  }
  
  
  return(var_est)
  #return(c(avg_effect1-avg_effect0,coef))
}

OLS_var_full_AS2=function(y,d1,d0,weights,x,mode,fo,so,contrast,dict_group){
  
  num_component=nrow(dict_group)
  subject_size=length(y)
  output=rep(y,each=2)
  contrast_OLS = c(contrast,rep(0,ncol(x)))
  fo_vec=unlist(fo)
  var_est=0
  if(mode=='OLS_joint'){
    
    data=cbind(y,d0,d1,x)
    data=data.frame(data)
    
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
    #print('finished')
    for (i in 1:nrow(dict_group)){
      
      c=dict_group[i,1]
      #print(c)
      index_start=dict_group[i,4]
      index_end=dict_group[i,5]    
      
      #select needed component
      
      output_t=z_OLS_joint_hat[index_start:index_end]
      #print(length(output_t))      
      fo_t=fo[[c]]
      so_t=so[[c]]
      
      cov_bound_t = so_AS2[[c]]
      
      p_t = so_t + (fo_t) %*% t(fo_t) 
      p_t=round(p_t,10)
      p_zeros_t=(p_t==0)
      cov_bound_zeros_t = (abs(cov_bound_t)<1e-2)
      if (min(cov_bound_zeros_t-p_zeros_t)<0){
        print('The denomiator matrix is wrong')
      }
      
      p_denomiator_t = p_t + cov_bound_zeros_t
      cov_bound_HT_t= cov_bound_t / p_denomiator_t
      #print(dim(cov_bound_HT_t))
      var_est = var_est + (output_t) %*% cov_bound_HT_t %*% (output_t)
      
    }
    
    var_est = var_est/(subject_size^2)
    
  }else if (mode=='WLS_joint'){
    data=cbind(y,d0,d1,x)
    data=data.frame(data)
    weights = weights
    
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
    
    for (i in 1:nrow(dict_group)){
      
      c=dict_group[i,1]
      index_start=dict_group[i,4]
      index_end=dict_group[i,5]    
      
      #select needed component
      
      output_t=z_WLS_joint_hat[index_start:index_end]
      
      fo_t=fo[[c]]
      so_t=so[[c]]
      
      
      cov_bound_t = so_AS2[[c]]
      p_t = so_t + (fo_t) %*% t(fo_t) 
      p_t=round(p_t,10)
      p_zeros_t=(p_t==0)
      cov_bound_zeros_t = (cov_bound_t==0)
      if (min(cov_bound_zeros_t-p_zeros_t)<0){
        print('The denomiator matrix is wrong')
      }
      
      p_denomiator_t = p_t + cov_bound_zeros_t
      cov_bound_HT_t= cov_bound_t / p_denomiator_t
      
      var_est = var_est + (output_t) %*% cov_bound_HT_t %*% (output_t)/(subject_size^2)
      
    }
    
  }
  
  
  
  return(var_est)
  #return(c(avg_effect1-avg_effect0,coef))
  
}

Optimal_Reg_var_full_AS2=function(y,d1,d0,x,optimal_b){
  
  
  num_component=nrow(dict_group)
  subject_size=length(y)
  var_est=0
  
  #create the outcome vector
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=y*d0 
  y_temp[2*(1:subject_size)]=y*d1
  
  treated_ind=rep(c(0,1),subject_size)
  control_ind=rep(c(1,0),subject_size)
  
  x_temp=matrix(0,nrow=subject_size*2,ncol(x))
  x_temp[2*(1:subject_size)-1,]=x
  x_temp[2*(1:subject_size),]=x
  
  x_temp=cbind(control_ind,treated_ind,x_temp)
  optimal_pred= x_temp %*% optimal_b
  optimal_pred = sapply(optimal_pred,function(x){ return( min(max(0,x),1) )})
  
  z_est=  y_temp - optimal_pred
  z_est[2*(1:subject_size)]=z_est[2*(1:subject_size)] * d1
  z_est[2*(1:subject_size)-1]=-z_est[2*(1:subject_size)-1] * d0
  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    
    z_est_t=z_est[index_start:index_end]
    
    #print('Optimal_Reg')
    #print(i)
    fo_t=fo[[c]]
    so_t=so[[c]]
    
    
    cov_bound_t = so_AS2[[c]]
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    p_zeros_t=(round(p_t,3)==0)
    cov_bound_zeros_t = (abs(cov_bound_t)<1e-2)
    
    if (max(cov_bound_zeros_t-p_zeros_t)<0){
      print('The denomiator matrix is wrong')
    }
    
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    
    var_est = var_est + z_est_t %*% cov_bound_HT_t %*% z_est_t/(subject_size)^2
  }
  
  return(var_est)
  #return(c(avg_effect1-avg_effect0,coef))
}

Optimal_Reg_var_full=function(y,d1,d0,x,optimal_b){
  
  
  num_component=nrow(dict_group)
  subject_size=length(y)
  var_est=0
  
  #create the outcome vector
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=y*d0 #make sure the sign is correct: it should be negative
  y_temp[2*(1:subject_size)]=y*d1
  
  treated_ind=rep(c(0,1),subject_size)
  control_ind=rep(c(1,0),subject_size)
  
  x_temp=matrix(0,nrow=subject_size*2,ncol(x))
  x_temp[2*(1:subject_size)-1,]=x
  x_temp[2*(1:subject_size),]=x
  
  x_temp=cbind(control_ind,treated_ind,x_temp)
  
  z_est=  y_temp - x_temp %*% optimal_b
  z_est[2*(1:subject_size)]=z_est[2*(1:subject_size)] * d1
  z_est[2*(1:subject_size)-1]=-z_est[2*(1:subject_size)-1] * d0
  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    
    z_est_t=z_est[index_start:index_end]
    
    #print('Optimal_Reg')
    #print(i)
    fo_t=fo[[c]]
    so_t=so[[c]]
    
    normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    
    cov_bound_t = AS_bound(normalized_cov_t)
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    p_zeros_t=(p_t==0)
    cov_bound_zeros_t = (abs(cov_bound_t)<1e-2)
    if (min(cov_bound_zeros_t-p_zeros_t)<0){
      print('The denomiator matrix is wrong')
    }
    
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    
    var_est = var_est + z_est_t %*% cov_bound_HT_t %*% z_est_t/(subject_size)^2
  }
  
  return(var_est)
  #return(c(avg_effect1-avg_effect0,coef))
}
