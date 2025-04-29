#####################################################################
#############Variance Bound Estimation###############################
#####################################################################

######################################
####Generate Variance Bounds##########
######################################
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
    
    #rows who has joint-observaility problem with the row with largest unobserved-centraility
    index_temp=which( (minus_ones_t[which.max(degree),]+1)<1e-3)
    
    #
    index_temp=c(which.max(degree),index_temp)
    
    
    bound_matrix[index_temp,index_temp] =bound_matrix[index_temp,index_temp]+1
    minus_ones_t[index_temp,index_temp]=0 #for those that are taken care of assigned 0
    
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

###############################################
#####HT variances #############################
###############################################
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
    
    ##################################
    ####Generate Variance Bound#######
    ##################################
    
    normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    
    normalized_cov_t = round(normalized_cov_t,10)
    minus_ones_t=- (abs(normalized_cov_t+1)<0.001)
    
    add_on=share_pain_bound_old(minus_ones_t)
    cov_bound_t = normalized_cov_t+ add_on
    
    
    ###################################
    
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    #want to avoid division by zero problem    
    p_zeros_t=as.matrix((abs(p_t - 0))<1e-3)
    cov_bound_zeros_t = as.matrix((abs(cov_bound_t - 0))<1e-3)
    
    #p_zeros_t * (1-cov_bound_zeros_t) are entries where the denominaotr is zero but the numerator is nonzero
    #this is inconsistent because p_zeros_t indicates joint nonobservaility
    #this sometimes happens due to simulation errors when calculating FOSO probabilities
    #if this happens too often, it means something is wrong and we return a error message
    if ( sum( p_zeros_t * (1-cov_bound_zeros_t))/length(p_zeros_t)>0.01){
      print(c)
      print(sum( p_zeros_t * (1-cov_bound_zeros_t)))
      print('The denomiator matrix may be wrong')
    }
    
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
    
    #want to avoid division by zero problem    
    p_zeros_t=as.matrix((abs(p_t - 0))<1e-3)
    cov_bound_zeros_t = as.matrix((abs(cov_bound_t - 0))<1e-3)
    
    #p_zeros_t * (1-cov_bound_zeros_t) are entries where the denominaotr is zero but the numerator is nonzero
    #this is inconsistent because p_zeros_t indicates joint nonobservaility
    #this sometimes happens due to simulation errors when calculating FOSO probabilities
    #if this happens too often, it means something is wrong and we return a error message
    if ( sum( p_zeros_t * (1-cov_bound_zeros_t))/length(p_zeros_t)>0.01){
      print(c)
      print(sum( p_zeros_t * (1-cov_bound_zeros_t)))
      print('The denomiator matrix may be wrong')
    }

    
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    var_est = var_est + (contrast_vec*output_t) %*% cov_bound_HT_t %*% (output_t*contrast_vec)
    
  }
  
  var_est=var_est/(subject_size)^2
  return(var_est)
  
  
  
}

###############################################
#####HA variances #############################
###############################################
HA_var_full=function(y,d1,d0,pi,sample_size,fo,so,contrast,dict_group){
  
  num_component=nrow(dict_group)
  subject_size=length(y)
  var_est=0
  var_est2=0
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

    ##################################
    ####Generate Variance Bound#######
    ##################################
    
    normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    
    normalized_cov_t = round(normalized_cov_t,10)
    minus_ones_t=- (abs(normalized_cov_t+1)<0.001)
    
    add_on=share_pain_bound_old(minus_ones_t)
    cov_bound_t = normalized_cov_t+ add_on
    
    
    ###################################
    p_t = so_t + (fo_t) %*% t(fo_t) 
    
    p_t=round(p_t,10)
    
    #want to avoid division by zero problem    
    p_zeros_t=as.matrix((abs(p_t - 0))<1e-3)
    cov_bound_zeros_t = as.matrix((abs(cov_bound_t - 0))<1e-3)
    
    #p_zeros_t * (1-cov_bound_zeros_t) are entries where the denominaotr is zero but the numerator is nonzero
    #this is inconsistent because p_zeros_t indicates joint nonobservaility
    #this sometimes happens due to simulation errors when calculating FOSO probabilities
    #if this happens too often, it means something is wrong and we return a error message
    if ( sum( p_zeros_t * (1-cov_bound_zeros_t))/length(p_zeros_t)>0.01){
      print(c)
      print(sum( p_zeros_t * (1-cov_bound_zeros_t)))
      print('The denomiator matrix may be wrong')
    }
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    contrast_vec=rep(contrast,length(output_t)/2)
    
    var_est = var_est + (contrast_vec*output_t) %*% cov_bound_HT_t %*% (output_t*contrast_vec)
    
    ###################################################
    
    cov_bound_t2 = so_AS2[[c]]
    
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    p_zeros_t=(p_t==0)
    cov_bound_zeros_t2 = (abs(cov_bound_t)<1e-2)
    if (min(cov_bound_zeros_t-p_zeros_t)<0){
      print('The denomiator matrix is wrong')
    }
    
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    contrast_vec=rep(contrast,length(output_t)/2)
    
    var_est2 = var_est2 + (contrast_vec*output_t) %*% cov_bound_HT_t %*% (output_t*contrast_vec)
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
    
    #want to avoid division by zero problem    
    p_zeros_t=as.matrix((abs(p_t - 0))<1e-3)
    cov_bound_zeros_t = as.matrix((abs(cov_bound_t - 0))<1e-3)
    
    #p_zeros_t * (1-cov_bound_zeros_t) are entries where the denominaotr is zero but the numerator is nonzero
    #this is inconsistent because p_zeros_t indicates joint nonobservaility
    #this sometimes happens due to simulation errors when calculating FOSO probabilities
    #if this happens too often, it means something is wrong and we return a error message
    if ( sum( p_zeros_t * (1-cov_bound_zeros_t))/length(p_zeros_t)>0.01){
      print(c)
      print(sum( p_zeros_t * (1-cov_bound_zeros_t)))
      print('The denomiator matrix may be wrong')
    }
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    contrast_vec=rep(contrast,length(output_t)/2)
    
    var_est = var_est + (contrast_vec*output_t) %*% cov_bound_HT_t %*% (output_t*contrast_vec)
    
  }
  
  var_est=var_est/(subject_size)^2
  
  return(var_est)
  
}

###############################################
####OLS variances #############################
###############################################
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
      
      ##################################
      ####Generate Variance Bound#######
      ##################################
      
      normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
      
      normalized_cov_t = round(normalized_cov_t,10)
      minus_ones_t=- (abs(normalized_cov_t+1)<0.001)
      
      add_on=share_pain_bound_old(minus_ones_t)
      cov_bound_t = normalized_cov_t+ add_on
      
      
      ###################################

      p_t = so_t + (fo_t) %*% t(fo_t) 
      p_t=round(p_t,10)
      
      #want to avoid division by zero problem    
      p_zeros_t=as.matrix((abs(p_t - 0))<1e-3)
      cov_bound_zeros_t = as.matrix((abs(cov_bound_t - 0))<1e-3)
      
      #p_zeros_t * (1-cov_bound_zeros_t) are entries where the denominaotr is zero but the numerator is nonzero
      #this is inconsistent because p_zeros_t indicates joint nonobservaility
      #this sometimes happens due to simulation errors when calculating FOSO probabilities
      #if this happens too often, it means something is wrong and we return a error message
      if ( sum( p_zeros_t * (1-cov_bound_zeros_t))/length(p_zeros_t)>0.01){
        print(c)
        print(sum( p_zeros_t * (1-cov_bound_zeros_t)))
        print('The denomiator matrix may be wrong')
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
      
      ##################################
      ####Generate Variance Bound#######
      ##################################
      
      normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
      
      normalized_cov_t = round(normalized_cov_t,10)
      minus_ones_t=- (abs(normalized_cov_t+1)<0.001)
      
      add_on=share_pain_bound_old(minus_ones_t)
      cov_bound_t = normalized_cov_t+ add_on
      
      
      ###################################
      p_t = so_t + (fo_t) %*% t(fo_t) 
      p_t=round(p_t,10)
      
      #want to avoid division by zero problem    
      p_zeros_t=as.matrix((abs(p_t - 0))<1e-3)
      cov_bound_zeros_t = as.matrix((abs(cov_bound_t - 0))<1e-3)
      
      #p_zeros_t * (1-cov_bound_zeros_t) are entries where the denominaotr is zero but the numerator is nonzero
      #this is inconsistent because p_zeros_t indicates joint nonobservaility
      #this sometimes happens due to simulation errors when calculating FOSO probabilities
      #if this happens too often, it means something is wrong and we return a error message
      if ( sum( p_zeros_t * (1-cov_bound_zeros_t))/length(p_zeros_t)>0.01){
        print(c)
        print(sum( p_zeros_t * (1-cov_bound_zeros_t)))
        print('The denomiator matrix may be wrong')
      }
      
      p_denomiator_t = p_t + cov_bound_zeros_t
      cov_bound_HT_t= cov_bound_t / p_denomiator_t
      
      var_est = var_est + (output_t) %*% cov_bound_HT_t %*% (output_t)
      
    }
    
    var_est = var_est/(subject_size^2)
    
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
      
      #want to avoid division by zero problem    
      p_zeros_t=as.matrix((abs(p_t - 0))<1e-3)
      cov_bound_zeros_t = as.matrix((abs(cov_bound_t - 0))<1e-3)
      
      #p_zeros_t * (1-cov_bound_zeros_t) are entries where the denominaotr is zero but the numerator is nonzero
      #this is inconsistent because p_zeros_t indicates joint nonobservaility
      #this sometimes happens due to simulation errors when calculating FOSO probabilities
      #if this happens too often, it means something is wrong and we return a error message
      if ( sum( p_zeros_t * (1-cov_bound_zeros_t))/length(p_zeros_t)>0.01){
        print(c)
        print(sum( p_zeros_t * (1-cov_bound_zeros_t)))
        print('The denomiator matrix may be wrong')
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
      
      #want to avoid division by zero problem    
      p_zeros_t=as.matrix((abs(p_t - 0))<1e-3)
      cov_bound_zeros_t = as.matrix((abs(cov_bound_t - 0))<1e-3)
      
      #p_zeros_t * (1-cov_bound_zeros_t) are entries where the denominaotr is zero but the numerator is nonzero
      #this is inconsistent because p_zeros_t indicates joint nonobservaility
      #this sometimes happens due to simulation errors when calculating FOSO probabilities
      #if this happens too often, it means something is wrong and we return a error message
      if ( sum( p_zeros_t * (1-cov_bound_zeros_t))/length(p_zeros_t)>0.01){
        print(c)
        print(sum( p_zeros_t * (1-cov_bound_zeros_t)))
        print('The denomiator matrix may be wrong')
      }
      p_denomiator_t = p_t + cov_bound_zeros_t
      cov_bound_HT_t= cov_bound_t / p_denomiator_t
      
      var_est = var_est + (output_t) %*% cov_bound_HT_t %*% (output_t)/(subject_size^2)
      
    }
    
  }
  
  
  
  return(var_est)
  #return(c(avg_effect1-avg_effect0,coef))
  
}
###############################################
####Logit variances ###########################
###############################################
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
    #print(output_t[1:10])
    contrast_vec=rep(contrast,length(output_t)/2)
    
    #print('Logit')
    #print(i)
    
    
    ##################################
    ####Generate Variance Bound#######
    ##################################
    
    normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    
    normalized_cov_t = round(normalized_cov_t,10)
    minus_ones_t=- (abs(normalized_cov_t+1)<0.001)
    
    add_on=share_pain_bound_old(minus_ones_t)
    cov_bound_t = normalized_cov_t+ add_on
    
    
    ###################################
    
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    #want to avoid division by zero problem    
    p_zeros_t=as.matrix((abs(p_t - 0))<1e-3)
    cov_bound_zeros_t = as.matrix((abs(cov_bound_t - 0))<1e-3)
    
    #p_zeros_t * (1-cov_bound_zeros_t) are entries where the denominaotr is zero but the numerator is nonzero
    #this is inconsistent because p_zeros_t indicates joint nonobservaility
    #this sometimes happens due to simulation errors when calculating FOSO probabilities
    #if this happens too often, it means something is wrong and we return a error message
    if ( sum( p_zeros_t * (1-cov_bound_zeros_t))/length(p_zeros_t)>0.01){
      print(c)
      print(sum( p_zeros_t * (1-cov_bound_zeros_t)))
      print('The denomiator matrix may be wrong')
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
    reg = glm(y~d0+d1+x-1,family='binomial',weights=weights,data=data)
    
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
    reg = glm(y~d0+d1+x-1,family='binomial',weights=weights)
    
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
    #print(output_t[1:10])
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

###############################################
####Optimal Regression variances ##############
###############################################

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
    
    #want to avoid division by zero problem    
    p_zeros_t=as.matrix((abs(p_t - 0))<1e-3)
    cov_bound_zeros_t = as.matrix((abs(cov_bound_t - 0))<1e-3)
    
    #p_zeros_t * (1-cov_bound_zeros_t) are entries where the denominaotr is zero but the numerator is nonzero
    #this is inconsistent because p_zeros_t indicates joint nonobservaility
    #this sometimes happens due to simulation errors when calculating FOSO probabilities
    #if this happens too often, it means something is wrong and we return a error message
    if ( sum( p_zeros_t * (1-cov_bound_zeros_t))/length(p_zeros_t)>0.01){
      print(c)
      print(sum( p_zeros_t * (1-cov_bound_zeros_t)))
      print('The denomiator matrix may be wrong')
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
    
    ##################################
    ####Generate Variance Bound#######
    ##################################
    
    normalized_cov_t=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    
    normalized_cov_t = round(normalized_cov_t,10)
    minus_ones_t=- (abs(normalized_cov_t+1)<0.001)
    
    add_on=share_pain_bound_old(minus_ones_t)
    cov_bound_t = normalized_cov_t+ add_on
    
    
    ###################################
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    #want to avoid division by zero problem    
    p_zeros_t=as.matrix((abs(p_t - 0))<1e-3)
    cov_bound_zeros_t = as.matrix((abs(cov_bound_t - 0))<1e-3)
    
    #p_zeros_t * (1-cov_bound_zeros_t) are entries where the denominaotr is zero but the numerator is nonzero
    #this is inconsistent because p_zeros_t indicates joint nonobservaility
    #this sometimes happens due to simulation errors when calculating FOSO probabilities
    #if this happens too often, it means something is wrong and we return a error message
    if ( sum( p_zeros_t * (1-cov_bound_zeros_t))/length(p_zeros_t)>0.01){
      print(c)
      print(sum( p_zeros_t * (1-cov_bound_zeros_t)))
      print('The denomiator matrix may be wrong')
    }
    
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    
    var_est = var_est + z_est_t %*% cov_bound_HT_t %*% z_est_t/(subject_size)^2
  }
  
  return(var_est)
  #return(c(avg_effect1-avg_effect0,coef))
}
