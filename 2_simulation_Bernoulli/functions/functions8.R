NO_HARM=function(y,d1,d0,pi0,pi1,y_p1,y_p0){
  
  
  print(y_p1[1:10])
  print(y_p0[1:10])
  
  
  
  
  #######
  #Step 1: Compute the correction factor
  #######
  
  #import prediction
  subject_size=length(y)
  
  #this the optimal coefficient colde
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=-y*d0/pi0 #make sure the sign is correct: it should be negative
  y_temp[2*(1:subject_size)]=y*d1/pi1
  
  
  x_temp=matrix(0,nrow=subject_size*2,1)
  x_temp[2*(1:subject_size)-1,]=-y_p0 #imputed
  x_temp[2*(1:subject_size),]=y_p1
  #  x_temp=cbind(control_ind,treated_ind,x_temp)
  
  x_pred=matrix(0,nrow=subject_size*2,1)
  x_pred[2*(1:subject_size)-1,]=y_p0
  x_pred[2*(1:subject_size),]=y_p1 
  #  x_pred=cbind(-control_ind,treated_ind,x_pred)
  
  
  XDX2=matrix(0,nrow=1,ncol=1)
  XDY2=rep(0,1)
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]] 
    
    normalized_cov=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    #normalized_cov=so_t
    
    XDX2 = XDX2 + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% x_temp[index_start:index_end,]/subject_size
    
    XDY2 = XDY2 + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% y_temp[index_start:index_end]/subject_size
    
  }
  
  #XDX2=XDX2+diag(dim(XDX2)[1]) * 1/(subject_size)^(1/3)
  
  optimal_b2 = solve(XDX2,XDY2)
  
  #optimal predictions
  optimal_pred2 = x_pred %*% optimal_b2
  
  y_p1_2=optimal_pred2[2*(1:subject_size)]
  y_p0_2=optimal_pred2[2*(1:subject_size)-1]
  
  #DR
  avg_effect1_2 = mean(y_p1_2) +  mean((y-y_p1_2) * d1/pi1)
  avg_effect0_2 = mean(y_p0_2) +  mean((y-y_p0_2) * d0/pi0)    
  
  effect_est2=avg_effect1_2-avg_effect0_2
  
  
  #optimal_b2 = 1
  #######
  #Step 2: Compute the implied variance
  #######
  
  subject_size=length(y)
  var_est=0
  
  #create the outcome vector
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=y*d0 
  y_temp[2*(1:subject_size)]=y*d1
  
  
  z_est=  y_temp - x_pred %*% optimal_b2
  z_est[2*(1:subject_size)]=z_est[2*(1:subject_size)] * d1
  z_est[2*(1:subject_size)-1]=-z_est[2*(1:subject_size)-1] * d0
  
  
  print(optimal_b2)
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    
    z_est_t=z_est[index_start:index_end]
    #print(z_est_t[1:10])
    #print('Optimal_Reg')
    #print(i)
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
    
    # print(dim(cov_bound_HT_t))
    #  print(length(z_est_t))
    var_est = var_est + z_est_t %*% cov_bound_HT_t %*% z_est_t/(subject_size)^2
  }
  
  return(list(effect_est2,var_est))
  
}

NO_HARM_AS2=function(y,d1,d0,pi0,pi1,y_p1,y_p0){
  
  
  print(y_p1[1:10])
  print(y_p0[1:10])
  
  
  
  
  #######
  #Step 1: Compute the correction factor
  #######
  
  #import prediction
  subject_size=length(y)
  
  #this the optimal coefficient colde
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=-y*d0/pi0 #make sure the sign is correct: it should be negative
  y_temp[2*(1:subject_size)]=y*d1/pi1
  
  
  x_temp=matrix(0,nrow=subject_size*2,1)
  x_temp[2*(1:subject_size)-1,]=-y_p0 #imputed
  x_temp[2*(1:subject_size),]=y_p1
  #  x_temp=cbind(control_ind,treated_ind,x_temp)
  
  x_pred=matrix(0,nrow=subject_size*2,1)
  x_pred[2*(1:subject_size)-1,]=y_p0
  x_pred[2*(1:subject_size),]=y_p1 
  #  x_pred=cbind(-control_ind,treated_ind,x_pred)
  
  
  XDX2=matrix(0,nrow=1,ncol=1)
  XDY2=rep(0,1)
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]] 
    
    normalized_cov=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    #normalized_cov=so_t
    
    XDX2 = XDX2 + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% x_temp[index_start:index_end,]/subject_size
    
    XDY2 = XDY2 + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% y_temp[index_start:index_end]/subject_size
    
  }
  
  #XDX2=XDX2+diag(dim(XDX2)[1]) * 1/(subject_size)^(1/3)
  
  optimal_b2 = solve(XDX2,XDY2)
  
  #optimal predictions
  optimal_pred2 = x_pred %*% optimal_b2
  
  y_p1_2=optimal_pred2[2*(1:subject_size)]
  y_p0_2=optimal_pred2[2*(1:subject_size)-1]
  
  #DR
  avg_effect1_2 = mean(y_p1_2) +  mean((y-y_p1_2) * d1/pi1)
  avg_effect0_2 = mean(y_p0_2) +  mean((y-y_p0_2) * d0/pi0)    
  
  effect_est2=avg_effect1_2-avg_effect0_2
  
  
  #optimal_b2 = 1
  #######
  #Step 2: Compute the implied variance
  #######
  
  subject_size=length(y)
  var_est=0
  
  #create the outcome vector
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=y*d0 
  y_temp[2*(1:subject_size)]=y*d1
  
  
  z_est=  y_temp - x_pred %*% optimal_b2
  z_est[2*(1:subject_size)]=z_est[2*(1:subject_size)] * d1
  z_est[2*(1:subject_size)-1]=-z_est[2*(1:subject_size)-1] * d0
  
  
  print(optimal_b2)
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    
    z_est_t=z_est[index_start:index_end]
    #print(z_est_t[1:10])
    #print('Optimal_Reg')
    #print(i)
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
    
    # print(dim(cov_bound_HT_t))
    #  print(length(z_est_t))
    var_est = var_est + z_est_t %*% cov_bound_HT_t %*% z_est_t/(subject_size)^2
  }
  
  return(list(effect_est2,var_est))
  
}
NO_HARM2=function(y,d1,d0,pi0,pi1,y_p1,y_p0){
  
  print('no harm2')
  print(y_p1[1:10])
  print(y_p0[1:10])
  
  
  
  
  #######
  #Step 1: Compute the correction factor
  #######
  
  #import prediction
  subject_size=length(y)
  
  #this the optimal coefficient colde
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=-y*d0/pi0 #make sure the sign is correct: it should be negative
  y_temp[2*(1:subject_size)]=y*d1/pi1
  
  
  x_temp=matrix(0,nrow=subject_size*2,2)
  x_temp[2*(1:subject_size)-1,1]=-y_p0 #imputed
  x_temp[2*(1:subject_size),2]=y_p1
  #  x_temp=cbind(control_ind,treated_ind,x_temp)
  
  x_pred=matrix(0,nrow=subject_size*2,2)
  x_pred[2*(1:subject_size)-1,1]=y_p0
  x_pred[2*(1:subject_size),2]=y_p1 
  #  x_pred=cbind(-control_ind,treated_ind,x_pred)
  
  
  XDX2=matrix(0,nrow=2,ncol=2)
  XDY2=rep(0,2)
  XDY3=rep(0,2)
  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]] 
    
    normalized_cov=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    #normalized_cov=so_t
    
    XDX2 = XDX2 + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% x_temp[index_start:index_end,]/subject_size
    
    XDY2 = XDY2 + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% y_temp[index_start:index_end]/subject_size
    
  }
  
  #XDX2=XDX2+diag(dim(XDX2)[1]) * 1/(subject_size)^(1/3)
  
  optimal_b2 = solve(XDX2,XDY2)
  
  #optimal predictions
  optimal_pred2 = x_pred %*% optimal_b2
  
  y_p1_2=optimal_pred2[2*(1:subject_size)]
  y_p0_2=optimal_pred2[2*(1:subject_size)-1]
  
  #DR
  avg_effect1_2 = mean(y_p1_2) +  mean((y-y_p1_2) * d1/pi1)
  avg_effect0_2 = mean(y_p0_2) +  mean((y-y_p0_2) * d0/pi0)    
  
  effect_est2=avg_effect1_2-avg_effect0_2
  
  
  #optimal_b2 = 1
  #######
  #Step 2: Compute the implied variance
  #######
  
  subject_size=length(y)
  var_est=0
  
  #create the outcome vector
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=y*d0 
  y_temp[2*(1:subject_size)]=y*d1
  
  
  z_est=  y_temp - x_pred %*% optimal_b2
  z_est[2*(1:subject_size)]=z_est[2*(1:subject_size)] * d1
  z_est[2*(1:subject_size)-1]=-z_est[2*(1:subject_size)-1] * d0
  
  
  print(optimal_b2)
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    
    z_est_t=z_est[index_start:index_end]
    #print(z_est_t[1:10])
    #print('Optimal_Reg')
    #print(i)
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
    
    # print(dim(cov_bound_HT_t))
    #  print(length(z_est_t))
    var_est = var_est + z_est_t %*% cov_bound_HT_t %*% z_est_t/(subject_size)^2
  }
  
  return(list(effect_est2,var_est))
  
}

NO_HARM2_AS2=function(y,d1,d0,pi0,pi1,y_p1,y_p0){
  
  print('no harmAS2')
  
  print(y_p1[1:10])
  print(y_p0[1:10])
  
  
  
  
  #######
  #Step 1: Compute the correction factor
  #######
  
  #import prediction
  subject_size=length(y)
  
  #this the optimal coefficient colde
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=-y*d0/pi0 #make sure the sign is correct: it should be negative
  y_temp[2*(1:subject_size)]=y*d1/pi1
  
  
  x_temp=matrix(0,nrow=subject_size*2,3)
  x_temp[2*(1:subject_size)-1,1]=-1 #imputed
  x_temp[2*(1:subject_size),2]=1
  x_temp[2*(1:subject_size)-1,3]=-y_p0 #imputed
  x_temp[2*(1:subject_size),3]=y_p1
  #  x_temp=cbind(control_ind,treated_ind,x_temp)
  
  
  x_pred=matrix(0,nrow=subject_size*2,3)
  x_pred[2*(1:subject_size)-1,1]=1 #imputed
  x_pred[2*(1:subject_size),2]=1
  x_pred[2*(1:subject_size)-1,3]=y_p0
  x_pred[2*(1:subject_size),3]=y_p1
  
  #  x_pred=cbind(-control_ind,treated_ind,x_pred)
  
  
  XDX2=matrix(0,nrow=3,ncol=3)
  XDY2=rep(0,3)
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]] 
    
    normalized_cov=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    #normalized_cov=so_t
    
    XDX2 = XDX2 + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% x_temp[index_start:index_end,]/subject_size
    
    XDY2 = XDY2 + t(x_temp[index_start:index_end,]) %*% normalized_cov %*% y_temp[index_start:index_end]/subject_size
    
  }
  
  #XDX2=XDX2+diag(dim(XDX2)[1]) * 1/(subject_size)^(1/3)
  
  optimal_b2 = solve(XDX2,XDY2)
  print('The Coefficient Is')
  print(optimal_b2)
  #optimal predictions
  optimal_pred2 = x_pred %*% optimal_b2
  
  y_p1_2=optimal_pred2[2*(1:subject_size)]
  y_p0_2=optimal_pred2[2*(1:subject_size)-1]
  
  #DR
  avg_effect1_2 = mean(y_p1_2) +  mean((y-y_p1_2) * d1/pi1)
  avg_effect0_2 = mean(y_p0_2) +  mean((y-y_p0_2) * d0/pi0)    
  
  effect_est2=avg_effect1_2-avg_effect0_2
  
  
  #optimal_b2 = 1
  #######
  #Step 2: Compute the implied variance
  #######
  
  subject_size=length(y)
  var_est=0
  
  #create the outcome vector
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=y*d0 
  y_temp[2*(1:subject_size)]=y*d1
  
  
  z_est=  y_temp - x_pred %*% optimal_b2
  z_est[2*(1:subject_size)]=z_est[2*(1:subject_size)] * d1
  z_est[2*(1:subject_size)-1]=-z_est[2*(1:subject_size)-1] * d0
  
  
  print(optimal_b2)
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    
    z_est_t=z_est[index_start:index_end]
    #print(z_est_t[1:10])
    #print('Optimal_Reg')
    #print(i)
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
    
    # print(dim(cov_bound_HT_t))
    #  print(length(z_est_t))
    var_est = var_est + z_est_t %*% cov_bound_HT_t %*% z_est_t/(subject_size)^2
  }
  
  return(list(effect_est2,var_est))
  
}
