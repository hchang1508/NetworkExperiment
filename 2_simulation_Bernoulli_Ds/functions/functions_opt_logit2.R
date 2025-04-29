##############################################################
# This Script Contains Scripts for Optimal Logit Adjustments##
##############################################################
##############################################################
##############################################################
Optimal_Logit_Est_M = function(y,d1,d0,x,pi0,pi1,Maxiter=1,Optim=FALSE,Option='BFGS'){
  
  
  #Iterated Weighted Least Square
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x) 
  control_covar = cbind(control_ind,x)
  
  #optimal coefficients
  optimal_b=rnorm(ncol(treated_covar))
  
  #If Optim==TRUE, search for local optimum
  if (Optim==TRUE ){
    print('Optimizing')
    #print(optimal_b)
    optimal_b = optim(par=optimal_b,fn=Optimal_Logit_M_Criterion,gr=Optimal_Logit_Gradient,method=Option,
                      y1=y1_temp,y0=y0_temp,d1=d1,d0=d0,pi1=pi1,pi0=pi0,treated_covar=treated_covar,control_covar=control_covar,
                      control=c(maxit=1000,trace=1))   
    optimal_b=optimal_b$par
    
  }
  
  treat_util = treated_covar %*% optimal_b
  control_util = control_covar %*% optimal_b
  #predicition
  y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  #DR
  avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
  avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)    
  
  effect_est=avg_effect1-avg_effect0
  
  
  return(list(effect_est,optimal_b))
}

Optimal_Logit_M_Criterion=function(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar){
  
  subject_size=length(y1)
  treat_util = treated_covar %*% beta
  control_util = control_covar %*% beta
  
  
  treat_impute = sapply(treat_util,logistic)
  control_impute = sapply(control_util,logistic)
  
  #intercept1=(1-sapply(treated_covar[,2] *beta[2],logistic))*sapply(treated_covar[,2] *beta[2],logistic)
  
  #intercept2=(1-sapply( control_covar[,1] *beta[1],logistic))*sapply( control_covar[,1] *beta[1],logistic)
  #normalization_impute = (intercept1 * intercept2)
  
  #for (i in 3:9){
  #  normalization_impute =  normalization_impute* rep(logistic(beta[i]),subject_size) * rep(1-logistic(beta[i]),subject_size) *2^2
  #}
  #print(max(normalization_impute))
  
  #print(max(normalization_impute))
  normalization_impute=1
  ######c'(y-f)
  residual1=rep(0,2*length(y1)) 
  residual1[2*(1:subject_size)-1]=-y0 * d0/pi0 
  residual1[2*(1:subject_size)]=y1 * d1/pi1 
  
  residual2=rep(0,2*length(y1)) 
  residual2[2*(1:subject_size)-1]= -control_impute
  residual2[2*(1:subject_size)]=+ treat_impute
  
  M_criterion1=0
  M_criterion2=0
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]      
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=sweep(so_t,1,1/fo_t,'*')
    normalized_cov=sweep(normalized_cov,2,1/fo_t,'*')
    
    
    M_criterion1= M_criterion1 + t(residual1[index_start:index_end])%*% normalized_cov %*% residual2[index_start:index_end]/subject_size
    M_criterion2= M_criterion2 + t(residual2[index_start:index_end]) %*% normalized_cov %*% residual2[index_start:index_end]/subject_size
    
  }  
  M_criterion= -2* M_criterion1 + M_criterion2
  return(M_criterion)
  
}



Optimal_Logit_Gradient=function(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar){
  
  subject_size=length(y1)
  treat_util = treated_covar %*% beta
  control_util = control_covar %*% beta
  
  
  treat_impute = sapply(treat_util,logistic)
  control_impute = sapply(control_util,logistic)
  

  ######c'(y-f)
  residual1=rep(0,2*length(y1)) 
  residual1[2*(1:subject_size)-1]=-y0 * d0/pi0 
  residual1[2*(1:subject_size)]=y1 * d1/pi1 
  
  residual2=rep(0,2*length(y1)) 
  residual2[2*(1:subject_size)-1]=-control_impute
  residual2[2*(1:subject_size)]= treat_impute
  
  Nablaf_treat = sweep(treated_covar,1,treat_impute*(1-treat_impute)/normalization_impute,'*')
  Nablaf_control=sweep(control_covar,1,-control_impute*(1-control_impute)/normalization_impute,'*')
  
  
  Nablaf=matrix(0,nrow=2*subject_size,ncol=ncol(Nablaf_treat))
  Nablaf[2*(1:subject_size)-1,]=Nablaf_control
  Nablaf[2*(1:subject_size),]=Nablaf_treat
  
  M_criterion1=rep(0,1,ncol(treated_covar))
  M_criterion2=rep(0,1,ncol(treated_covar))  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]      
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=sweep(so_t,1,1/fo_t,'*')
    normalized_cov=sweep(normalized_cov,2,1/fo_t,'*')
    
    
    M_criterion1= M_criterion1 + residual1[index_start:index_end] %*% normalized_cov %*% Nablaf[index_start:index_end,]/subject_size
    M_criterion2= M_criterion2 + residual2[index_start:index_end] %*% normalized_cov %*% Nablaf[index_start:index_end,]/subject_size
    
  }  
  M_criterion=-2*(M_criterion1) +2* M_criterion2

  return(M_criterion)
  
}


######################################
########Variance Bound Estimator######
######################################
Optimal_Logit_var_full_AS2=function(y,d1,d0,x,optimal_b){
  
  
  num_component=nrow(dict_group)
  subject_size=length(y)
  var_est=0
  
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x)
  control_covar = cbind(control_ind,x)
  
  treat_util = treated_covar %*% optimal_b
  control_util = control_covar %*% optimal_b
  
  #predicition
  y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  z_est=rep(0,subject_size*2)
  
  z_est[2*(1:subject_size)]=  (y1_temp - y_p1) *d1
  z_est[2*(1:subject_size)-1]=  -(y0_temp - y_p0) *d0
  
  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    
    z_est_t=z_est[index_start:index_end]
    
    #print('Optimal_Logistic')
    
    
    
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
    
    
    var_est = var_est + z_est_t %*% cov_bound_HT_t %*% z_est_t/(subject_size)^2
  }
  
  return(var_est)
  #return(c(avg_effect1-avg_effect0,coef))
}

Optimal_Logit_var_full=function(y,d1,d0,x,optimal_b){
  
  
  num_component=nrow(dict_group)
  subject_size=length(y)
  var_est=0
  
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x)
  control_covar = cbind(control_ind,x)
  
  treat_util = treated_covar %*% optimal_b
  control_util = control_covar %*% optimal_b
  
  #predicition
  y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  z_est=rep(0,subject_size*2)
  
  z_est[2*(1:subject_size)]=  (y1_temp - y_p1) *d1
  z_est[2*(1:subject_size)-1]=  -(y0_temp - y_p0) *d0
  
  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    
    z_est_t=z_est[index_start:index_end]
    
    #print('Optimal_Logistic')
    
    
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
    
    p_zeros_t=(p_t==0)
    cov_bound_zeros_t = (cov_bound_t==0)
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


########################################################
########Theoretical Optimal Coefficient Estimation######
########################################################
Optimal_Logit_M=function(y1,y0,x,Maxiter=1,Optim=TRUE){
  
  subject_size=length(y1)
  
  #create the outcome vector
  y0_temp=y0
  y1_temp=y1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x)
  control_covar = cbind(control_ind,x)
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x)
  control_covar = cbind(control_ind,x)
  
  #optimal coefficients
  beta_init=rnorm(ncol(treated_covar))
  if (Optim==TRUE){
    print('Optimizing')
    optimal_b = optim(par=beta_init,fn=Optimal_Logit_M_Population_Criterion,gr=Optimal_Logit_Population_Gradient,method='BFGS',
                      y1=y1_temp,y0=y0_temp,treated_covar=treated_covar,control_covar=control_covar,control=c(trace=1,maxit=500))
    
    optimal_b=optimal_b$par
  }
  treat_util = treated_covar %*% optimal_b
  control_util = control_covar %*% optimal_b
  
  treat_impute = sapply(treat_util,logistic)
  control_impute = sapply(control_util,logistic)
  
  y1_res= y1-treat_impute
  y0_res= y0-control_impute
  
  z_Optim_Logistic= rep(0,2*subject_size)
  
  z_Optim_Logistic[2*(1:subject_size)]=y1_res
  z_Optim_Logistic[2*(1:subject_size)-1]=-y0_res
  
  return(z_Optim_Logistic)
}


Optimal_Logit_M_Population_Criterion=function(beta,y1,y0,treated_covar,control_covar){
  
  subject_size=length(y1)
  treat_util = treated_covar %*% beta
  control_util = control_covar %*% beta
  
  
  treat_impute = sapply(treat_util,logistic)
  control_impute = sapply(control_util,logistic)
  

  ######c'(y-f)
  residual1=rep(0,2*length(y1)) 
  residual1[2*(1:subject_size)-1]=-y0 
  residual1[2*(1:subject_size)]=y1
  
  residual2=rep(0,2*length(y1)) 
  residual2[2*(1:subject_size)-1]= -control_impute
  residual2[2*(1:subject_size)]=+ treat_impute
  
  M_criterion1=0
  M_criterion2=0
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]      
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=sweep(so_t,1,1/fo_t,'*')
    normalized_cov=sweep(normalized_cov,2,1/fo_t,'*')
    
    
    M_criterion1= M_criterion1 + t(residual1[index_start:index_end])%*% normalized_cov %*% residual2[index_start:index_end]/subject_size
    M_criterion2= M_criterion2 + t(residual2[index_start:index_end]) %*% normalized_cov %*% residual2[index_start:index_end]/subject_size
    
  }  
  
  M_criterion= -2* M_criterion1 + M_criterion2
  
  return(M_criterion)
  
}

Optimal_Logit_Population_Gradient=function(beta,y1,y0,treated_covar,control_covar){
  
  subject_size=length(y1)
  treat_util = treated_covar %*% beta
  control_util = control_covar %*% beta
  
  
  treat_impute = sapply(treat_util,logistic)
  control_impute = sapply(control_util,logistic)
  
  
  ######c'(y-f)
  residual1=rep(0,2*length(y1)) 
  residual1[2*(1:subject_size)-1]=-y0  
  residual1[2*(1:subject_size)]=y1 
  
  residual2=rep(0,2*length(y1)) 
  residual2[2*(1:subject_size)-1]=-control_impute
  residual2[2*(1:subject_size)]= treat_impute
  
  Nablaf_treat = sweep(treated_covar,1,treat_impute*(1-treat_impute),'*')
  Nablaf_control=sweep(control_covar,1,-control_impute*(1-control_impute),'*')
  
  
  Nablaf=matrix(0,nrow=2*subject_size,ncol=ncol(Nablaf_treat))
  Nablaf[2*(1:subject_size)-1,]=Nablaf_control
  Nablaf[2*(1:subject_size),]=Nablaf_treat
  
  M_criterion1=rep(0,1,ncol(treated_covar))
  M_criterion2=rep(0,1,ncol(treated_covar))  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]      
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=sweep(so_t,1,1/fo_t,'*')
    normalized_cov=sweep(normalized_cov,2,1/fo_t,'*')
    
    
    M_criterion1= M_criterion1 + residual1[index_start:index_end] %*% normalized_cov %*% Nablaf[index_start:index_end,]/subject_size
    M_criterion2= M_criterion2 + residual2[index_start:index_end] %*% normalized_cov %*% Nablaf[index_start:index_end,]/subject_size
    
  }  
  M_criterion=-2*(M_criterion1) +2* M_criterion2
  
  return(M_criterion)
  
}

####Testing cases########
#beta=rnorm(9)

#print('First Coordinate')
#beta1= beta
#beta1[1]=beta1[1]+0.005
#beta2= beta
#beta2[1]=beta2[1]-0.005

#(Optimal_Logit_M_Criterion(beta1,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar) - Optimal_Logit_M_Criterion(beta2,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar))/(2*0.005) 
#Optimal_Logit_Gradient(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar )

#print('Second Coordinate')
#beta1= beta
#beta1[2]=beta1[2]+0.005
#beta2= beta
#beta2[2]=beta2[2]-0.005

#(Optimal_Logit_M_Criterion(beta1,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar) - Optimal_Logit_M_Criterion(beta2,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar))/(2*0.005) 
#Optimal_Logit_Gradient(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar )


#print('Third Coordinate')
#beta1= beta
#beta1[3]=beta1[3]+0.005
#beta2= beta
#beta2[3]=beta2[3]-0.005

#(Optimal_Logit_M_Criterion(beta1,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar) - Optimal_Logit_M_Criterion(beta2,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar))/(2*0.005) 
#Optimal_Logit_Gradient(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar )


#print('Fourth Coordinate')
#beta1= beta
#beta1[4]=beta1[4]+0.005
#beta2= beta
#beta2[4]=beta2[4]-0.005

#(Optimal_Logit_M_Criterion(beta1,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar) - Optimal_Logit_M_Criterion(beta2,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar))/(2*0.005) 
#Optimal_Logit_Gradient(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar )

#print('Fifth Coordinate')
#beta1= beta
#beta1[5]=beta1[5]+0.005
#beta2= beta
#beta2[5]=beta2[5]-0.005

#(Optimal_Logit_M_Criterion(beta1,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar) - Optimal_Logit_M_Criterion(beta2,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar))/(2*0.005) 
#Optimal_Logit_Gradient(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar )

