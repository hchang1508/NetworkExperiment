# This file contain functions used in simulations:
# 
# A. HT Estimator
# 
# B. HA Estimator
# 
# C. OLS Estimator
# 
# D. 
# 

###############################################################################
########################HT Estimator###########################################
###############################################################################
HT_full = function(epsilon1,epsilon0,pi1,pi0,sample_size){
  
  est = sum(epsilon1 * (1/pi1)) - sum(epsilon0 * (1/pi0))
  
  est = est/sample_size
  
  return(est)
  
  
  
}


###############################################################################
########################HA Estimator###########################################
###############################################################################

HA_full =function(epsilon1,epsilon0,pi1,pi0,denom_pi1,denom_pi0){
  
  est = sum(epsilon1 * (1/pi1))/sum(denom_pi1) - sum(epsilon0 * (1/pi0))/sum(denom_pi0)
  
  return(est)
  
}

###############################################################################
########################OLS Estimator##########################################
###############################################################################

OLS_full=function(y,d1,d0,weights,x,pi0,pi1,mode){
  
  if (mode=='OLS_joint'){
    
    data=cbind(y,d0,d1,x)
    data=data.frame(data)
    
    #weights
    weights = as.numeric((weights>(0+1e-6)))
    
    #OLS estimator
    reg = lm(y~d0+d1+x-1,weights=weights,data=data)
    coef=reg$coefficients
    
    #create covariates for prediction
    x1=cbind(rep(0,nrow(x)),rep(1,nrow(x)),x)
    x0=cbind(rep(1,nrow(x)),rep(0,nrow(x)),x)
    x1=data.frame(x1)
    x0=data.frame(x0)
    colnames(x1)=colnames(data)[2:length(data)]
    colnames(x0)=colnames(data)[2:length(data)]
    
    #predicted outcomes 
    y_p1 = predict(reg,x1,type='response')
    y_p0 = predict(reg,x0,type='response')
    
    #estimate average effects
    avg_effect1=reg$coefficients['d1']
    avg_effect0=reg$coefficients['d0']
    
  }else if (mode=='WLS_joint'){
    
    data=cbind(y,d0,d1,x)
    data=data.frame(data)
    
    #weights
    weights = weights
    
    #regression
    reg = lm(y~d0+d1+x-1,weights=weights)
    
    #create covariates for prediction
    x1=cbind(rep(0,nrow(x)),rep(1,nrow(x)),x)
    x0=cbind(rep(1,nrow(x)),rep(0,nrow(x)),x)
    x1=data.frame(x1)
    x0=data.frame(x0)
    colnames(x1)=colnames(data)[2:length(data)]
    colnames(x0)=colnames(data)[2:length(data)]
    
    #weighted least squaress
    y_p1 = predict(reg,x1,type='response')
    y_p0 = predict(reg,x0,type='response')
    
    #estimate average effects
    coef=reg$coefficients 
    avg_effect1=reg$coefficients['d1']
    avg_effect0=reg$coefficients['d0']
    
  }
  
  return(list(avg_effect1-avg_effect0,y_p1,y_p0,coef))

}

###############################################################################
########################Logisitc Imputation####################################
###############################################################################

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
    
    weights = as.numeric((weights>(0+1e-6)))
    reg = glm(y~d0+d1+x-1,family='binomial',weights=weights,data=data)
    
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
    reg = glm(y~d0+d1+x-1,family='binomial',weights=weights)
    
    y_p1 = predict(reg,x1,type='response')
    y_p0 = predict(reg,x0,type='response')
    coef=reg$coefficients 
    print(coef)
    avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
    avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)
  }
  return(list(avg_effect1-avg_effect0,y_p1,y_p0,coef))
}

