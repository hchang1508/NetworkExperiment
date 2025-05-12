##############################################################
# This Script Contains Scripts for Optimal Logit Adjustments##
##############################################################
##############################################################
##############################################################
Optimal_Logit_Est_A = function(y,d1,d0,x,pi0,pi1,Maxiter=1,Optim=FALSE,Option='BFGS'){
  
  
  #Iterated Weighted Least Square
  
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x) 
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  print(dim(treated_covar))
  iter=0
  tol=1e-4
  diff=1
  beta_old=rnorm(ncol(treated_covar)) 
  #beta_old=rep(0,ncol(treated_covar))
  # print(beta_old)
  while (diff>tol & iter<Maxiter){
    #   
    beta_new=Solve_Logit_OLS_Est(y,d1,d0,pi0,pi1,treated_covar,control_covar,beta_old,X)
    #   
    beta_new=as.matrix(beta_new)
    #   
    #    print(beta_new)
    #  #    print((beta_new-beta_old)^2)
    diff= sum((beta_new-beta_old)^2)
    #   
    beta_old = 0.5 * beta_new + 0.5* beta_old
    #   
    iter=iter+1
  }
  print(beta_old)
  ind=0
  if (iter==Maxiter){
    ind=1
  }
  #optimal coefficients
  optimal_b=rep(0,ncol(treated_covar))
  
  #If Optim==TRUE, search for local optimum
  if (Optim==TRUE ){
    print('Optimizing')
    print(optimal_b)
    optimal_b = optim(par=rep(0,ncol(treated_covar)),fn=Optimal_Logit_Z_Criterion,gr=gradient_sample,method=Option,
                      y1=y1_temp,y0=y0_temp,d1=d1,d0=d0,pi1=pi1,pi0=pi0,treated_covar=treated_covar,control_covar=control_covar,X=X,XX=XX,
                      control=c(maxit=1000,trace=1))   
    optimal_b=optimal_b$par
    
  }
  print(optimal_b)
  treat_util = treated_covar %*% optimal_b
  control_util = control_covar %*% optimal_b
  
  #predicition
  y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  #DR
  avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
  avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)    
  
  effect_est=avg_effect1-avg_effect0
  
  
  return(list(effect_est,optimal_b,ind))
}

Optimal_Logit_Est_B = function(y,d1,d0,x,pi0,pi1,Maxiter=1,Optim=FALSE,Option='BFGS'){
  
  
  #Iterated Weighted Least Square
  print('In B')
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x) 
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  print(dim(treated_covar))
  iter=0
  tol=1e-4
  diff=1
  beta_old=rep(0,ncol(treated_covar))
  
  while (diff>tol & iter<Maxiter){
    #
    beta_new=Solve_Logit_OLS_Est(y,d1,d0,pi0,pi1,treated_covar,control_covar,beta_old,X)
    #
    beta_new=as.matrix(beta_new)
    #
    #    print(beta_new)
    #    print((beta_new-beta_old)^2)
    diff= sqrt(sum((beta_new-beta_old)^2))
    #
    beta_old= beta_new + beta_old
    print(diff)
    iter=iter+1
  }
  print(paste0('iter is',iter))
  print('nonlinear 2')
  print(beta_old)
  print(iter)
  
  optimal_b = beta_old
  print(optimal_b)
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
Optimal_Logit_Est_C = function(y,d1,d0,x,pi0,pi1,Maxiter=1,Optim=FALSE,Option='BFGS',optimal_b_fixed){
  
  
  #local search with multiple random starts
  print('in C')
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x) 
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  
  
  # #random starts
  # Max_out_iter=1000
  # terminate=0
  # out_iter=1
  # while (terminate==0 & out_iter<Max_out_iter){
  # 
  #   print(paste0('In outer loop', out_iter))
  #   terminate=1
  #   iter=1
  #   tol=1e-4
  #   diff=1
  #   if (out_iter==1){
  #     beta_old=rep(0,9)
  #   }else{
  #     beta_old=rnorm(9)*0.1
  #     
  #   }
  #   weight=0.1
  #   foc=10
  #   print(beta_old)
  #   loss_new=10
  #   ind=0
  #   while (foc>tol & iter<Maxiter & diff>tol){
  #     
  #     loss_old=Optimal_Logit_Z_Criterion(beta_old,y1_temp,y0_temp,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
  #     ind=0
  #     weight=0.1
  #     iter=0
  #     while (ind==0){
  #       
  #       updates=update(beta_old,weight,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
  #     
  #       beta_new=updates[[1]]
  #       foc=sum(updates[[2]]^2)
  #     
  #       beta_new=as.matrix(beta_new)
  #       diff= sqrt(sum((beta_new-beta_old)^2))
  #       print('iter')
  #       print(iter)
  #       print('foc')
  #       print(foc)
  #       print(foc>tol)
  #       print('diff')
  #       print(diff)
  #       print(beta_new)
  #       loss_new=Optimal_Logit_Z_Criterion(beta_new,y1_temp,y0_temp,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
  # 
  #       if (loss_new>loss_old){
  #         weight = weight/2
  #         print('weights halved')
  #         print('backtracking')
  #       
  #       }else{
  #         ind =1
  #       }
  #       iter=iter+1
  #     }
  #       if (max(abs(beta_old)>(3+out_iter*0.1))){
  #         print('Exceed search buondary')
  #         terminate=0
  #         break
  #     }
  #       beta_old=beta_new
  #       loss_old=loss_new 
  #   }
  # 
  #     out_iter = out_iter + 1
  # }
  # 
  # if (out_iter==Max_out_iter){
  #   stop()
  # }  
  # if (max(abs(beta_old)>10)){
  #   stop()
  # }
  # if (ind==0){
  #   stop()
  # }
  # if (iter==Maxiter){
  #   stop()
  # }
  # print(paste0('iter is',iter))
  # print('newton')
  # print(beta_old)
  # 
  # optimal_b=beta_old
  
  beta_old = rnorm(ncol(treated_covar))
  min_val=Optimal_Logit_Z_Criterion(beta_old,y1_temp,y0_temp,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
  #If Optim==TRUE, search for local optimum
  if (Optim==TRUE & min_val>1e-4){
    print('Optimizing')
    
    optimal_b = optim(par= beta_old,fn=Optimal_Logit_Z_Criterion,gr=gradient_sample,method=Option,
                      y1=y1_temp,y0=y0_temp,d1=d1,d0=d0,pi1=pi1,pi0=pi0,treated_covar=treated_covar,control_covar=control_covar,X=X,XX=XX,
                      control=list(maxit=500,trace=1))   
    optimal_b=optimal_b$par
    
  }
  print(optimal_b)
  treat_util = treated_covar %*% optimal_b
  control_util = control_covar %*% optimal_b
  
  #predicition
  y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  #DR
  avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
  avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)    
  
  effect_est=avg_effect1-avg_effect0
  
  
  return(list(effect_est,optimal_b,ind))
}

Optimal_Logit_Est_D = function(y,d1,d0,x,pi0,pi1,Maxiter=1,Optim=FALSE,Option='BFGS'){
  
  
  #Iterated Weighted Least Square
  
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x) 
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  print(dim(treated_covar))
  iter=0
  tol=1e-4
  diff=1
  beta_old=rep(0,ncol(treated_covar))
  
  print(beta_old)
  while (diff>tol & iter<Maxiter){
    #
    beta_new=Solve_Logit_OLS_Est(y,d1,d0,pi0,pi1,treated_covar,control_covar,beta_old,X)
    #
    beta_new=as.matrix(beta_new)
    #
    #    print(beta_new)
    #    print((beta_new-beta_old)^2)
    diff= sum((beta_new-beta_old)^2)
    #
    beta_old = 0.5* beta_new + 0.5* beta_old
    #
    iter=iter+1
  }
  #
  ind=0
  if (iter==Maxiter){
    ind=1
  }
  print(iter)
  optimal_b=rnorm(ncol(treated_covar))
  print(Optimal_Logit_Z_Criterion(optimal_b,y1_temp,y0_temp,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX))
  print(Optimal_Logit_Z_Criterion(rep(0,ncol(treated_covar)),y1_temp,y0_temp,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)) 
  #If Optim==TRUE, search for local optimum
  if (Optim==TRUE ){
    print('Optimizing')
    print(optimal_b)
    print(Option)
    optimal_b = optim(par= optimal_b,fn=Optimal_Logit_Z_Criterion,gr=gradient_sample,method=Option,
                      y1=y1_temp,y0=y0_temp,d1=d1,d0=d0,pi1=pi1,pi0=pi0,treated_covar=treated_covar,control_covar=control_covar,X=X,XX=XX,
                      control=c(maxit=1000,trace=1,abstol=1e-4,gamma=1.01))   
    optimal_b=optimal_b$par
    
  }
  print(optimal_b)
  treat_util = treated_covar %*% optimal_b
  control_util = control_covar %*% optimal_b
  
  #predicition
  y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  #DR
  avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
  avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)    
  
  effect_est=avg_effect1-avg_effect0
  
  
  return(list(effect_est,optimal_b,ind))
}
Optimal_Logit_Est_E = function(y,d1,d0,x,pi0,pi1,Maxiter=1,Optim=FALSE,Option='BFGS',optimal_b_fixed){
  
  
  #local search with multiple random starts
  print('in C')
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x) 
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  
  #If Optim==TRUE, search for local optimum
  if (Optim==TRUE ){
    print('Optimizing')
    
    optimal_b = optimx(par= rep(0,9),fn=Optimal_Logit_Z_Criterion,gr=gradient_sample,lower=-3,upper=3,method='spg',
                      y1=y1_temp,y0=y0_temp,d1=d1,d0=d0,pi1=pi1,pi0=pi0,treated_covar=treated_covar,control_covar=control_covar,X=X,XX=XX,
                      control=list(trace=1))   
    optimal_b=optimal_b$par
    
  }
  print(optimal_b)
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

Optimal_Bound_Logit_A = function(y,d1,d0,x,pi0,pi1,Maxiter=1,Optim=FALSE,Option='BFGS'){
  
  
  #Iterated Weighted Least Square
  print('in A Bound')
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x) 
  control_covar = cbind(control_ind,x)
  
  beta_old=rnorm(9)
  #If Optim==TRUE, search for local optimum
  if (Optim==TRUE ){
    print('Optimizing')
    
    optimal_b = optim(par= beta_old,fn=Optimal_Logit_M_Criterion_Bound_EST,gr=gradient_sample_Bound,method=Option,
                      y1=y1_temp,y0=y0_temp,d1=d1,d0=d0,pi1=pi1,pi0=pi0,treated_covar=treated_covar,control_covar=control_covar,X=X,XX=XX,
                      control=list(maxit=500,trace=1))   
    optimal_b=optimal_b$par
    
  }
  print(optimal_b)
  treat_util = treated_covar %*% optimal_b
  control_util = control_covar %*% optimal_b
  
  #predicition
  y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  #DR
  avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
  avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)    
  
  effect_est=avg_effect1-avg_effect0
  
  
  return(list(effect_est,optimal_b,ind))
}
Optimal_Bound_Logit_B = function(y,d1,d0,x,pi0,pi1,Maxiter=1,Optim=FALSE,Option='BFGS',optimal_b_fixed){
  
  
  #Iterated Weighted Least Square
  print('in B bound')
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x) 
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  
  Max_out_iter=1000
  terminate=0
  out_iter=1
  while (terminate==0 & out_iter<Max_out_iter){
    
    print(paste0('In outer loop', out_iter))
    terminate=1
    iter=1
    tol=1e-4
    diff=1
    if (out_iter==1){
      beta_old=rep(0,9)
    }else{
      beta_old=rnorm(9)*0.1
      
    }
    weight=1
    foc=10
    print(beta_old)
    loss_new=10
    ind=0
    while (foc>tol & iter<Maxiter & diff>tol){
      loss_old=Optimal_Logit_M_Criterion_Bound_EST(beta_old,y1_temp,y0_temp,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
      ind=0
      weight=0.1
      
      while (ind==0){
        print(weight)
        updates=update_bound(beta_old,weight,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
        
        beta_new=updates[[1]]
        foc=sum(updates[[2]]^2)
        
        beta_new=as.numeric(beta_new)
        diff= sqrt(sum((beta_new-beta_old)^2))
        print('iter')
        print(iter)
        print('foc')
        print(foc)
        loss_new=Optimal_Logit_M_Criterion_Bound_EST(beta_new,y1_temp,y0_temp,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
        
        if (loss_new>loss_old){
          weight = weight/2
          print('weights halved')
          print('backtracking')
          print(paste0('loss_new',loss_new))
          print(paste0('loss_old',loss_old))
          
          
        }else{
          ind =1
        }
        if (max(abs(beta_old)>2)){
          print('Exceed search buondary')
          terminate=0
          break
        }
      }
      iter=iter+1
      
      if (max(abs(beta_old)>2)){
        print('Exceed search buondary')
        terminate=0
        break
      }
      beta_old=beta_new
      loss_old=loss_new 
    }
    
    out_iter = out_iter + 1
  }
  
  if (out_iter==Max_out_iter){
    stop()
  }  
  if (max(abs(beta_old)>1.5)){
    stop()
  }
  #
  ind=0
  if (iter==Maxiter){
    stop()
  }
  print(paste0('iter is',iter))
  print('newton')
  print(beta_old)
  
  optimal_b=beta_old
  
  min_val=Optimal_Logit_M_Criterion_Bound_EST(beta_old,y1_temp,y0_temp,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
  #If Optim==TRUE, search for local optimum
  if (Optim==TRUE & min_val>1e-4){
    print('Optimizing')
    
    optimal_b = optim(par= beta_old,fn=Optimal_Logit_Z_Criterion,gr=gradient_sample,method=Option,
                      y1=y1_temp,y0=y0_temp,d1=d1,d0=d0,pi1=pi1,pi0=pi0,treated_covar=treated_covar,control_covar=control_covar,X=X,XX=XX,
                      control=list(maxit=500,trace=1))   
    optimal_b=optimal_b$par
    
  }
  print(optimal_b)
  treat_util = treated_covar %*% optimal_b
  control_util = control_covar %*% optimal_b
  
  #predicition
  y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  #DR
  avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
  avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)    
  
  effect_est=avg_effect1-avg_effect0
  
  
  return(list(effect_est,optimal_b,ind))
}
Optimal_Bound_Logit_C = function(y,d1,d0,x,pi0,pi1,Maxiter=1,Optim=FALSE,Option='BFGS',optimal_b_fixed){
  
  
  #Iterated Weighted Least Square
  print('in B bound')
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x) 
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  
  Max_out_iter=1000
  terminate=0
  out_iter=1
  while (terminate==0 & out_iter<Max_out_iter){
    
    print(paste0('In outer loop', out_iter))
    terminate=1
    iter=1
    tol=1e-4
    diff=1
    if (out_iter==1){
      beta_old=as.numeric(optimal_b_fixed)
    }else{
      beta_old=rnorm(9)*0.1
      
    }
    foc=10
    print(beta_old)
    loss_new=10
    ind=0
    coord = 0
    while (foc>tol & iter<Maxiter & diff>tol){
      
      loss_old=Optimal_Logit_M_Criterion_Bound_EST(beta_old,y1_temp,y0_temp,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
      ind=0
      weight=0.1
      coord = (coord ) %% 9 +1
       
      #output new beta_old and loss foc
      output=Inner_Optimal_Bound_Logit_C(y1_temp,y0_temp,d1,d0,pi1,pi0,beta_old,weighting,coord,treated_covar,control_covar,X,XX)
      
      beta_old = output[[1]]
      loss_new=Optimal_Logit_M_Criterion_Bound_EST(beta_old,y1_temp,y0_temp,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
      
      diff = loss_old-loss_new
      print('beta_old')
      print(beta_old)
      print('difference is')
      print(diff)
      
      updates=update_bound(beta_old,weight,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
      print('foc')
      foc=sum(updates[[2]]^2)
      print(foc)
      out_iter = out_iter + 1
  }
  
  if (out_iter==Max_out_iter){
    stop()
  }  
  if (max(abs(beta_old)>1.5)){
    stop()
  }
  #
  ind=0
  if (iter==Maxiter){
    stop()
  }
  print(paste0('iter is',iter))
  print('newton')
  print(beta_old)
  
  optimal_b=beta_old

  print(optimal_b)
  treat_util = treated_covar %*% optimal_b
  control_util = control_covar %*% optimal_b
  
  #predicition
  y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  #DR
  avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
  avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)    
  
  effect_est=avg_effect1-avg_effect0
  
  
  return(list(effect_est,optimal_b,ind))
}

}
Optimal_Bound_Logit_D = function(y,d1,d0,x,pi0,pi1,Maxiter=1,Optim=FALSE,Option='BFGS'){
  
  
  #Iterated Weighted Least Square
  
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x) 
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  print(dim(treated_covar))
  iter=0
  tol=1e-4
  diff=1
  beta_old=rnorm(ncol(treated_covar)) 
  #beta_old=rep(0,ncol(treated_covar))
  # print(beta_old)
  while (diff>tol & iter<Maxiter){
    #   
    beta_new=Solve_Logit_OLS_Est_Bound(y,d1,d0,pi0,pi1,treated_covar,control_covar,beta_old,X)
    #   
    beta_new=as.matrix(beta_new)
    #
    
    print('beta_old')
    print(beta_old)
    print('beta_new')
    print(beta_new)
    #  #    print((beta_new-beta_old)^2)
    diff= sum((beta_new-beta_old)^2)
    #   
    beta_old = 0.1* beta_new + beta_old
    #   
    iter=iter+1
  }
  print(beta_old)
  ind=0
  if (iter==Maxiter){
    ind=1
  }
  #optimal coefficients
  optimal_b=rep(0,ncol(treated_covar))
  
  #If Optim==TRUE, search for local optimum
  if (Optim==TRUE ){
    print('Optimizing')
    print(optimal_b)
    optimal_b = optim(par=rep(0,ncol(treated_covar)),fn=Optimal_Logit_Z_Criterion,gr=gradient_sample,method=Option,
                      y1=y1_temp,y0=y0_temp,d1=d1,d0=d0,pi1=pi1,pi0=pi0,treated_covar=treated_covar,control_covar=control_covar,X=X,XX=XX,
                      control=c(maxit=1000,trace=1))   
    optimal_b=optimal_b$par
    
  }
  print(optimal_b)
  treat_util = treated_covar %*% optimal_b
  control_util = control_covar %*% optimal_b
  
  #predicition
  y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  #DR
  avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
  avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)    
  
  effect_est=avg_effect1-avg_effect0
  
  
  return(list(effect_est,optimal_b,ind))
}

Inner_Optimal_Bound_Logit_C=function(y1,y0,d1,d0,pi1,pi0,beta_old,weighting,coord,treated_covar,control_covar,X,XX){
  
  foc=10 #track first order condition
  tol=1e-4 #tolerance
  Maxiter=1000 #maximum iteration
  iter=1 #iteration counter
  diff=10 #loss differences
  ind=0
  weight=0.3
  print(paste0('Optimizing ',coord, 'th coordinate'))
  while (foc>tol & iter<Maxiter & diff>tol){
    ind=0
    while (ind==0){
      print(beta_old)
      loss_old= Optimal_Logit_M_Criterion_Bound_EST(beta_old,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
      updates=update_bound_coord(beta_old,weight,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX,coord)
    
      beta_new=updates[[1]]
      foc=sum(updates[[2]]^2)
    
      beta_new=as.numeric(beta_new)
      diff= sqrt(sum((beta_new-beta_old)^2))
      print('iter')
      print(iter)
      print('foc')
      print(foc)
      loss_new=Optimal_Logit_M_Criterion_Bound_EST(beta_new,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
    
      if (loss_new>loss_old){
        weight = weight/2
        print('weights halved')
        print('backtracking')
        print(paste0('loss_new ',loss_new))
        print(paste0('loss_old ',loss_old))
      }else{
        ind =1
    }
      if (max(abs(beta_old)>10)){
        print('Exceed search buondary')
        terminate=0
        break
    }
  }
  iter=iter+1
  print(iter)
  beta_old=beta_new
  loss_old=loss_new 
  }
  print('Terminating')
  print(paste0('foc: ',foc))
  
  print(paste0('iter: ',iter))
  
  print(paste0('diff: ',diff))
  return(list(beta_old,loss_old))
}

update_inner=function(beta,weighting,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX){
  Jtemp=J(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
  #Jtemp[j,]
  Fc_temp=Fc(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
  Fc_temp=t(Fc_temp)
  foc=2*Jtemp %*% Fc_temp
  
  return(foc)
}

update=function(beta,weighting,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX){
  
  j=sample(1:9,1)
  diff= rep(0,9) 
  diff[j]=0.005
  #beta=rnorm(9)
  
  #(Fc(beta+diff,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)-Fc(beta-diff,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX))/(0.005*2)
  
  #Jtemp[j,]

  foc=update_inner(beta,weighting,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
  print('Numerical Derivative')
  print((Optimal_Logit_Z_Criterion(beta+diff,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)-Optimal_Logit_Z_Criterion(beta-diff,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX))/(0.005*2))
  print('Formula Derivative')
  
  print(foc[j])
  if (sum(foc^2)<1e-4){
  #second order condition
  H=matrix(0,9,9)
  for (j in 1:9){
    diff= rep(0,9)
    diff[j]=0.005
    H[j,]= (update_inner(beta+diff,weighting,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)-update_inner(beta-diff,weighting,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX))/0.005*2
  }
  print('Hessian')
  print(eigen(H)$values)}
  
  beta_new = beta -foc * weighting
  
  return(list(beta_new,foc))
  
}

Solve_Logit_OLS_Est=function(y,d1,d0,pi0,pi1,treated_covar,control_covar,beta0,X){
  
  
  treat_util = treated_covar %*% beta0
  control_util = control_covar %*% beta0
  
  normalization = treated_covar[,1] *beta0[1] +  treated_covar[,2] *beta0[2]
  
  normalization_weight= 1
    
    #predicition
    y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  #first step regression  
  p=rep(0,2*length(y))
  p[2*(1:subject_size)-1]=y_p1
  p[2*(1:subject_size)]=y_p0
  
  d=rep(0,2*length(d1))
  d[2*(1:subject_size)-1]=d1
  d[2*(1:subject_size)]=d0
  
  #create covariate matrix
  X_right=sweep(X,1,p*(1-p),'*')
  X_left=sweep(X,1,p*(1-p),'*')
  
  #create the outcome vector
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=-(y-y_p0)*d0/pi0 #make sure the sign is correct: it should be negative
  y_temp[2*(1:subject_size)]=(y-y_p1)*d1/pi1
  
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
    
    
    XDX = XDX + t(X_left[index_start:index_end,]) %*% normalized_cov %*% X_right[index_start:index_end,]/subject_size
    
    XDY = XDY + t(X_left[index_start:index_end,]) %*% normalized_cov %*% y_temp[index_start:index_end]/subject_size
    
    
  }
  
  XDX = XDX
  optimal_b = solve(XDX,XDY)
  
  return(optimal_b)
}

Optimal_Logit_M_Criterion=function(beta,y1,y0,treated_covar,control_covar,X,XX){
  
  
  subject_size=length(y1)
  treat_util = treated_covar %*% beta
  control_util = control_covar %*% beta
  
  treat_impute = sapply(treat_util,logistic)
  control_impute = sapply(control_util,logistic)
  
  
  residual=rep(0,subject_size*2)
  residual[2*(1:subject_size)]=y1-treat_impute
  residual[2*(1:subject_size)-1]=-(y0-control_impute)
  
  criterion=0
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=sweep(so_t,1,1/fo_t,'*')
    normalized_cov=sweep(normalized_cov,2,1/fo_t,'*')
    
    criterion=criterion + residual[index_start:index_end] %*% normalized_cov %*% residual[index_start:index_end]/subject_size
    
  }  
  
  return(criterion)
  
}

Optimal_Logit_M_Criterion_Bound=function(beta,y1,y0,treated_covar,control_covar,X,XX){
  
  
  subject_size=length(y1)
  treat_util = treated_covar %*% beta
  control_util = control_covar %*% beta
  
  treat_impute = sapply(treat_util,logistic)
  control_impute = sapply(control_util,logistic)
  
  
  residual=rep(0,subject_size*2)
  residual[2*(1:subject_size)]=y1-treat_impute
  residual[2*(1:subject_size)-1]=-(y0-control_impute)
  
  criterion=0
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    #extract variance bound
    cov_bound_t = so_AS2[[c]]
    
    #second order probabilities
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    cov_bound_zeros_t = (abs(cov_bound_t)<1e-2)
    
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    
    criterion=criterion + residual[index_start:index_end] %*% cov_bound_HT_t %*% residual[index_start:index_end]/subject_size
    
  }  
  
  return(criterion)
  
}


##########Auxiliary Functions##############################
J=function(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX){
  
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
  residual=rep(0,2*length(y1)) 
  residual[2*(1:subject_size)-1]=-y0 * d0/pi0 + control_impute
  residual[2*(1:subject_size)]=y1 * d1/pi1 - treat_impute
  p=rep(0,2*length(y1))
  p[2*(1:subject_size)-1]=control_impute
  p[2*(1:subject_size)]=treat_impute
  
  
  Nablaf_treat = sweep(treated_covar,1,-treat_impute*(1-treat_impute)/normalization_impute,'*')
  Nablaf_control=sweep(control_covar,1,control_impute*(1-control_impute)/normalization_impute,'*')
  
  
  
  Nablaf=matrix(0,nrow=2*subject_size,ncol=ncol(Nablaf_treat))
  Nablaf[2*(1:subject_size)-1,]=Nablaf_control
  Nablaf[2*(1:subject_size),]=Nablaf_treat
  
  Z_criterion1=rep(0,ncol(treated_covar),ncol(treated_covar))
  Z_criterion2=rep(0,1,ncol(treated_covar)^2)  
  
  XX = sweep(XX,1,p*(1-p)*(1-2*p),'*')
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]      
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=sweep(so_t,1,1/fo_t,'*')
    normalized_cov=sweep(normalized_cov,2,1/fo_t,'*')
    
    
    Z_criterion1= Z_criterion1 + t(Nablaf[index_start:index_end,]) %*% normalized_cov %*% Nablaf[index_start:index_end,]/subject_size
    Z_criterion2= Z_criterion2 + residual[index_start:index_end] %*% normalized_cov %*% XX[index_start:index_end,]/subject_size
    
  }  
  Z_criterion2 = matrix(Z_criterion2, ncol(X),ncol(X))
  Z_criterion2 = t(matrix(Z_criterion2, ncol(X),ncol(X)))
  Z_criterion=-Z_criterion1 + Z_criterion2
  
  return( Z_criterion)
  
}
J_bound=function(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX){
  
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
  residual=rep(0,2*length(y1)) 
  residual[2*(1:subject_size)-1]=-y0 * d0/pi0 + control_impute
  residual[2*(1:subject_size)]=y1 * d1/pi1 - treat_impute
  p=rep(0,2*length(y1))
  p[2*(1:subject_size)-1]=control_impute
  p[2*(1:subject_size)]=treat_impute
  
  
  Nablaf_treat = sweep(treated_covar,1,-treat_impute*(1-treat_impute)/normalization_impute,'*')
  Nablaf_control=sweep(control_covar,1,control_impute*(1-control_impute)/normalization_impute,'*')
  
  
  
  Nablaf=matrix(0,nrow=2*subject_size,ncol=ncol(Nablaf_treat))
  Nablaf[2*(1:subject_size)-1,]=Nablaf_control
  Nablaf[2*(1:subject_size),]=Nablaf_treat
  
  Z_criterion1=rep(0,1,ncol(treated_covar))
  Z_criterion2=rep(0,1,ncol(treated_covar))  
  
  XX = sweep(XX,1,p*(1-p)*(1-2*p),'*')
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]      
    fo_t=fo[[c]]
    so_t=so[[c]]
   #extract variance bound
    cov_bound_t = so_AS2[[c]]
    
    #second order probabilities
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    cov_bound_zeros_t = (abs(cov_bound_t)<1e-2)
    
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    Z_criterion1= Z_criterion1 + t(Nablaf[index_start:index_end,]) %*% cov_bound_HT_t %*% Nablaf[index_start:index_end,]/subject_size
    Z_criterion2= Z_criterion2 + residual[index_start:index_end] %*% cov_bound_HT_t %*% XX[index_start:index_end,]/subject_size
    
  }  
  Z_criterion2 = matrix(Z_criterion2, ncol(X),ncol(X))
  Z_criterion2 = t(matrix(Z_criterion2, ncol(X),ncol(X)))
  Z_criterion=-Z_criterion1 + Z_criterion2
  
  return( -2*Z_criterion)
  
}
J_pop=function(beta,y1,y0,pi1,pi0,treated_covar,control_covar,X,XX){
  
  subject_size=length(y1)
  treat_util = treated_covar %*% beta
  control_util = control_covar %*% beta
  
  treat_impute = sapply(treat_util,logistic)
  control_impute = sapply(control_util,logistic)
  
  ######c'(y-f)
  residual=rep(0,2*length(y1)) 
  residual[2*(1:subject_size)-1]=-y0  + control_impute
  residual[2*(1:subject_size)]=y1 - treat_impute
  p=rep(0,2*length(y1))
  p[2*(1:subject_size)-1]=control_impute
  p[2*(1:subject_size)]=treat_impute
  
  Nablaf_treat = sweep(treated_covar,1,treat_impute*(1-treat_impute),'*')
  Nablaf_control=sweep(control_covar,1,-control_impute*(1-control_impute),'*')
  
  Nablaf=matrix(0,nrow=2*subject_size,ncol=ncol(Nablaf_treat))
  Nablaf[2*(1:subject_size)-1,]=Nablaf_control
  Nablaf[2*(1:subject_size),]=Nablaf_treat
  
  Z_criterion1=rep(0,ncol(treated_covar),ncol(treated_covar))
  Z_criterion2=rep(0,1,ncol(treated_covar)^2)  
  
  XX = sweep(XX,1,p*(1-p)*(1-2*p),'*')
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]      
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=sweep(so_t,1,1/fo_t,'*')
    normalized_cov=sweep(normalized_cov,2,1/fo_t,'*')
    
    
    Z_criterion1= Z_criterion1 + t(Nablaf[index_start:index_end,]) %*% normalized_cov %*% Nablaf[index_start:index_end,]/subject_size
    Z_criterion2= Z_criterion2 + residual[index_start:index_end] %*% normalized_cov %*% XX[index_start:index_end,]/subject_size
    
  }  
  Z_criterion2 = matrix(Z_criterion2, ncol(X),ncol(X))

  Z_criterion= Z_criterion1 - Z_criterion2
  
  return(2* Z_criterion)
  
}

Fc=function(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX){
  
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
  residual2[2*(1:subject_size)-1]=control_impute
  residual2[2*(1:subject_size)]=- treat_impute
  
  Nablaf_treat = sweep(treated_covar,1,-treat_impute*(1-treat_impute)/normalization_impute,'*')
  Nablaf_control=sweep(control_covar,1,control_impute*(1-control_impute)/normalization_impute,'*')
  
  
  
  Nablaf=matrix(0,nrow=2*subject_size,ncol=ncol(Nablaf_treat))
  Nablaf[2*(1:subject_size)-1,]=Nablaf_control
  Nablaf[2*(1:subject_size),]=Nablaf_treat
  
  Z_criterion1=rep(0,1,ncol(treated_covar))
  Z_criterion2=rep(0,1,ncol(treated_covar))  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]      
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=sweep(so_t,1,1/fo_t,'*')
    normalized_cov=sweep(normalized_cov,2,1/fo_t,'*')
    
    
    Z_criterion1= Z_criterion1 + residual1[index_start:index_end] %*% normalized_cov %*% Nablaf[index_start:index_end,]/subject_size
    Z_criterion2= Z_criterion2 + residual2[index_start:index_end] %*% normalized_cov %*% Nablaf[index_start:index_end,]/subject_size
    
  }  
  Z_criterion=Z_criterion1 + Z_criterion2
  
  return(-Z_criterion)
  
}

Fc_pop=function(beta,y1,y0,pi1,pi0,treated_covar,control_covar,X,XX){
  
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
  residual2[2*(1:subject_size)-1]=control_impute
  residual2[2*(1:subject_size)]=- treat_impute
  
  Nablaf_treat = sweep(treated_covar,1,-treat_impute*(1-treat_impute),'*')
  Nablaf_control=sweep(control_covar,1,control_impute*(1-control_impute),'*')
  
  
  Nablaf=matrix(0,nrow=2*subject_size,ncol=ncol(Nablaf_treat))
  Nablaf[2*(1:subject_size)-1,]=Nablaf_control
  Nablaf[2*(1:subject_size),]=Nablaf_treat
  
  Z_criterion1=rep(0,1,ncol(treated_covar))
  Z_criterion2=rep(0,1,ncol(treated_covar))  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]      
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=sweep(so_t,1,1/fo_t,'*')
    normalized_cov=sweep(normalized_cov,2,1/fo_t,'*')
    
    
    Z_criterion1= Z_criterion1 + residual1[index_start:index_end] %*% normalized_cov %*% Nablaf[index_start:index_end,]/subject_size
    Z_criterion2= Z_criterion2 + residual2[index_start:index_end] %*% normalized_cov %*% Nablaf[index_start:index_end,]/subject_size
    
  }  
  Z_criterion=Z_criterion1 + Z_criterion2
  
  return(2*Z_criterion)
  
}

update_pop=function(beta,weighting,y1,y0,pi1,pi0,treated_covar,control_covar,X,XX){
  
  Jtemp=J_pop(beta,y1,y0,pi1,pi0,treated_covar,control_covar,X,XX)
  Fc_temp=Fc_pop(beta,y1,y0,pi1,pi0,treated_covar,control_covar,X,XX)
  Fc_temp=t(Fc_temp)
  foc=solve(Jtemp,Fc_temp)
  #foc=Fc_temp

  ###Test cases
 # beta=rnorm(9)
  j=sample(1:9,1)
  deltah=rep(0,9)
  deltah[j]=0.0005
  print('Numerical Derivative')
  print((Optimal_Logit_M_Criterion(beta+deltah,y1,y0,treated_covar,control_covar,X,XX)-Optimal_Logit_M_Criterion(beta-deltah,y1,y0,treated_covar,control_covar,X,XX))/(2*0.0005))
  print('Analytical Derivative')
  Fc_temp=Fc_pop(beta,y1,y0,pi1,pi0,treated_covar,control_covar,X,XX)
  print(Fc_temp[j])

  #print('Numerical Hessian')
  #print((Fc_pop(beta+deltah,y1,y0,pi1,pi0,treated_covar,control_covar,X,XX)-Fc_pop(beta-deltah,y1,y0,pi1,pi0,treated_covar,control_covar,X,XX))/(2*0.0005))
  #print('Analytical Hessian')
  
  #Jtemp=J_pop(beta,y1,y0,pi1,pi0,treated_covar,control_covar,X,XX)
  #print(Jtemp[j,])
  ###########
  print('beta')
  print(beta)
  print('foc')
  print(foc*weighting)
  beta_new = beta - foc * weighting
  print('beta_new')
  print(beta_new)
  return(list(beta_new,foc))
  
}

update_bound=function(beta,weighting,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX){
  

  J=J_bound(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
  descent= as.numeric(gradient_sample_Bound(y1,y0,d1,d0,pi1,pi0,beta,treated_covar,control_covar) )
 
  foc=solve(J,descent)
  
  beta_new = beta -foc * weighting
  print('beta_old')
  print(beta)
  print('Newton Direction')
  print(foc*weighting)
  print('beta_new')
  print(beta_new)
  return(list(beta_new,foc))
  

}

update_bound_coord=function(beta,weighting,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX,coord){
  
 # J=J_bound(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
  descent= as.numeric(gradient_sample_Bound(y1,y0,d1,d0,pi1,pi0,beta,treated_covar,control_covar) )
  
 # foc=solve(J,descent)
  beta_new = beta
  beta_new[coord] = beta[coord] -descent[coord] * weighting
  print('beta_old')
  print(beta)
  print('Newton Direction')
  print(descent*weighting)
  print('beta_new')
  print(beta_new)
  return(list(beta_new,descent[coord]))
  
  
}

Optimal_Logit_Z_Criterion=function(beta,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX){
  
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
  residual2[2*(1:subject_size)-1]=control_impute
  residual2[2*(1:subject_size)]=- treat_impute
  
  Nablaf_treat = sweep(treated_covar,1,-treat_impute*(1-treat_impute)/normalization_impute,'*')
  Nablaf_control=sweep(control_covar,1,control_impute*(1-control_impute)/normalization_impute,'*')
  
  
  Nablaf=matrix(0,nrow=2*subject_size,ncol=ncol(Nablaf_treat))
  Nablaf[2*(1:subject_size)-1,]=Nablaf_control
  Nablaf[2*(1:subject_size),]=Nablaf_treat
  
  Z_criterion1=rep(0,1,ncol(treated_covar))
  Z_criterion2=rep(0,1,ncol(treated_covar))  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]      
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=sweep(so_t,1,1/fo_t,'*')
    normalized_cov=sweep(normalized_cov,2,1/fo_t,'*')
    
    
    Z_criterion1= Z_criterion1 + residual1[index_start:index_end] %*% normalized_cov %*% Nablaf[index_start:index_end,]/subject_size
    Z_criterion2= Z_criterion2 + residual2[index_start:index_end] %*% normalized_cov %*% Nablaf[index_start:index_end,]/subject_size
    
  }  
  Z_criterion=Z_criterion1 + Z_criterion2
  criterion=sum(Z_criterion^2) + 1/subject_size*sum(beta^2) #add a penalty
  print(criterion)
  return(criterion)
  
}

Optimal_Logit_Est = function(y,d1,d0,x,pi0,pi1,beta_initial){
  
  ### local search given an intial guess
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x)
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  start_time=proc.time()
  # optimal_b = optim(par=beta_initial,fn=Optimal_Logit_Z_Criterion,gr=gradient_sample,method='Nelder-Mead',
  #                   y1=y1_temp,y0=y0_temp,d1=d1,d0=d0,pi1=pi1,pi0=pi0,treated_covar=treated_covar,control_covar=control_covar,X=X,XX=XX,
  #                   control=c(maxit=1000,trace=1))
  optimal_b = optim(par=beta_initial,fn=Optimal_Logit_Z_Criterion,gr=gradient_sample,method='BFGS',
                    y1=y1_temp,y0=y0_temp,d1=d1,d0=d0,pi1=pi1,pi0=pi0,treated_covar=treated_covar,control_covar=control_covar,X=X,XX=XX,
                    control=c(maxit=1000,trace=1))
  
  end_time=proc.time()
  
  #optimal coefficients
  optimal_b=optimal_b$par
  
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

Solve_Logit_OLS=function(y1,y0,treated_covar,control_covar,beta0,X){
  
  
  treat_util = treated_covar %*% beta0
  control_util = control_covar %*% beta0
  
  #predicition
  y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  #first step regression  
  p=rep(0,2*length(y1))
  p[2*(1:subject_size)-1]=y_p1
  p[2*(1:subject_size)]=y_p0
  
  #create covariate matrix
  X_right=sweep(X,1,p*(1-p),'*')
  X_left=sweep(X,1,p*(1-p),'*')
  
  #create the outcome vector
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=-(y0-y_p0) #make sure the sign is correct: it should be negative
  y_temp[2*(1:subject_size)]=(y1-y_p1)
  
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
    
    
    XDX = XDX + t(X_left[index_start:index_end,]) %*% normalized_cov %*% X_right[index_start:index_end,]/subject_size
    
    XDY = XDY + t(X_left[index_start:index_end,]) %*% normalized_cov %*% y_temp[index_start:index_end]/subject_size
    
    
  }
  
  optimal_b = solve(XDX,XDY)
  
  return(optimal_b)
}

Optimal_Logit=function(y1,y0,x){
  
  subject_size=length(y1)
  
  #create the outcome vector
  y0_temp=y0
  y1_temp=y1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x)
  control_covar = cbind(control_ind,x)
  
  X = second_order_matrix(treated_covar,control_covar )[[1]]
  optimal_b = optim(par=rep(0,ncol(treated_covar)),fn=Optimal_Logit_M_Criterion,gr=gradient_ppl,method='CG',
                    y1=y1_temp,y0=y0_temp,treated_covar=treated_covar,control_covar=control_covar,X=X,control=c(trace=1,maxit=250))
  
  beta=optimal_b$par
  
  treat_util = treated_covar %*% beta
  control_util = control_covar %*% beta
  
  treat_impute = sapply(treat_util,logistic)
  control_impute = sapply(control_util,logistic)
  
  y1_res= y1-treat_impute
  y0_res= y0-control_impute
  
  z_Optim_Logistic= rep(0,2*subject_size)
  
  z_Optim_Logistic[2*(1:subject_size)]=y1_res
  z_Optim_Logistic[2*(1:subject_size)-1]=-y0_res
  
  return(z_Optim_Logistic)
}

Optimal_Logit_2=function(y1,y0,x){
  
  subject_size=length(y1)
  
  #create the outcome vector
  y0_temp=y0
  y1_temp=y1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x)
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  
  
  #first step regression  
  p=rep(0,2*length(y1))
  p[2*(1:subject_size)-1]=0.5
  p[2*(1:subject_size)]=0.5
  
  
  #create covariate matrix
  X_right=sweep(X,1,1,'*')
  X_left=sweep(X,1,p*(1-p),'*')
  
  #create the outcome vector
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=-(y0-0.5)#make sure the sign is correct: it should be negative
  y_temp[2*(1:subject_size)]=(y1-0.5)
  
  XDX=matrix(0,nrow=ncol(x)+2,ncol=ncol(x)+2)
  XDY=rep(0,ncol(x)+2)
  for (i in 1:nrow(dict_group)){
    
    index_start=dict_group[i,3]
    index_end=dict_group[i,4]    
    fo_t=fo[[i]]
    so_t=so[[i]]
    normalized_cov=diag(1/fo_t) %*% so_t %*% diag(1/fo_t)
    #normalized_cov=so_t
    
    
    XDX = XDX + 0.25*t(X_left[index_start:index_end,]) %*% normalized_cov %*% X_right[index_start:index_end,]/subject_size
    
    XDY = XDY + t(X_left[index_start:index_end,]) %*% normalized_cov %*% y_temp[index_start:index_end]/subject_size
    
    
  }
  
  optimal_b = solve(XDX,XDY)
  
  optimal_b = optim(par=optimal_b,fn=Optimal_Logit_M_Criterion,gr=gradient_ppl,method='Nelder-Mead',
                    y1=y1_temp,y0=y0_temp,treated_covar=treated_covar,control_covar=control_covar,X=X,control=c(trace=1,maxit=500))
  
  beta=optimal_b$par
  
  treat_util = treated_covar %*% beta
  control_util = control_covar %*% beta
  
  treat_impute = sapply(treat_util,logistic)
  control_impute = sapply(control_util,logistic)
  
  y1_res= y1-treat_impute
  y0_res= y0-control_impute
  
  z_Optim_Logistic= rep(0,2*subject_size)
  
  z_Optim_Logistic[2*(1:subject_size)]=y1_res
  z_Optim_Logistic[2*(1:subject_size)-1]=-y0_res
  
  return(z_Optim_Logistic)
}

Optimal_Logit_A=function(y1,y0,x,Maxiter=1,Optim=FALSE){
  
  subject_size=length(y1)
  
  #create the outcome vector
  y0_temp=y0
  y1_temp=y1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x)
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x)
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  
  iter=0
  tol=1e-4
  diff=1
  beta_old=rep(0,ncol(treated_covar))
  
  while (diff>tol & iter<Maxiter){
    print(iter)
    beta_new=Solve_Logit_OLS(y1,y0,treated_covar,control_covar,beta_old,X)
    
    beta_new=as.matrix(beta_new)
    
    print(beta_new)
    print((beta_new-beta_old)^2)
    diff= sum((beta_new-beta_old)^2)
    
    beta_old = 0.5 * beta_new + 0.5* beta_old
    
    iter = iter+1
  }
  
  ind=0
  if (iter==Maxiter){
    ind=1
  }
  
  #optimal coefficients
  optimal_b=beta_old
  if (Optim==TRUE){
    print('Optimizing')
    optimal_b = optim(par=optimal_b,fn=Optimal_Logit_M_Criterion,gr=gradient_ppl,method='BFGS',
                      y1=y1_temp,y0=y0_temp,treated_covar=treated_covar,control_covar=control_covar,X=X,control=c(trace=1,maxit=1000))
    
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
Optimal_Logit_B=function(y1,y0,x,Maxiter=1,Optim=FALSE){
  
  subject_size=length(y1)
  
  #create the outcome vector
  y0_temp=y0
  y1_temp=y1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x)
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x)
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  
  iter=0
  tol=1e-4
  diff=1
  beta_old=rep(0,9)
  weight=0.1
  #random starts
  Max_out_iter=1000
  terminate=0
  out_iter=1
  while (terminate==0 & out_iter<Max_out_iter){
    
    print(paste0('In outer loop', out_iter))
    terminate=1
    iter=1
    tol=1e-4
    diff=1
    if (out_iter==1){
      beta_old=rep(0,9)
    }else{
      beta_old=rnorm(9) * (0.1+out_iter * 0.1)
      
    }
    weight=0.1
    foc=10
    print(beta_old)
    loss_new=10
    ind=0
  
     while ( iter<Maxiter & sum(foc^2)>tol){
      loss_old=Optimal_Logit_M_Criterion(beta_old,y1,y0,treated_covar,control_covar,X,XX)
      print(iter)
      weight=0.1
      ind=0
      iter=1
     ###########Backtracking block
        while (ind==0 ){
        updates=update_pop(beta_old,weight,y1,y0,pi1,pi0,treated_covar,control_covar,X,XX)
        beta_new=updates[[1]]
        foc=updates[[2]]
         #update new_loss
        loss_new=Optimal_Logit_M_Criterion(beta_new,y1,y0,treated_covar,control_covar,X,XX)
        beta_new=as.matrix(beta_new)
        print(iter)
        print('beta_new')
        print(beta_new)
        print('foc')
        print(sum(foc^2))
        print(foc)
        #print((beta_new-beta_old)^2)
        print('difference')
        diff= sqrt(sum((beta_new-beta_old)^2))
        print(diff)
        if (loss_new>loss_old){
          weight=weight/2
          print('backtracking')
          print('loss_old')
          print(loss_old)          
          
          print('loss_new')
          print(loss_new)          
          
        }else{
          ind=1
        }
        iter = iter+1
        
    
    }
    #################################
      if (max(abs(beta_old)>5)){
        print('Exceed search buondary')
        terminate=0
        break
      }
      beta_old = beta_new
      
     }
    out_iter = out_iter + 1

  }
  

  print('final result')

  #optimal coefficients
  optimal_b=beta_old
  if (Optim==TRUE){
    print('Optimizing')
    optimal_b = optim(par=optimal_b,fn=Optimal_Logit_M_Criterion,gr=gradient_ppl,method='Nelder-Mead',
                      y1=y1_temp,y0=y0_temp,treated_covar=treated_covar,control_covar=control_covar,X=X,control=c(trace=1,maxit=500))
    
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

logistic=function(u){
  
  return(1/(1+exp(-1*u)))
}

second_order_matrix = function(treated_covar,control_covar){
  
  
  ####
  
  #This is auxiliary function for calculating the gradient of the criterion function
  
  ####
  
  subject_size=nrow(treated_covar)
  X_temp =matrix(0,nrow=subject_size*2,ncol=ncol(treated_covar))  
  X_temp2=matrix(0,nrow=subject_size*2,ncol=ncol(treated_covar))  
  X_temp[2*(1:subject_size)-1,]=-control_covar
  X_temp[2*(1:subject_size),]=treated_covar
  X_temp2[2*(1:subject_size)-1,]=control_covar
  X_temp2[2*(1:subject_size),]=treated_covar
  XX_temp=matrix(0,nrow=subject_size*2,ncol=ncol(X_temp)^2)
  
  counter=1
  
  for (i in 1:ncol(X_temp)){
    
    for (j in 1:ncol(X_temp)){      
      
      XX_temp[,counter] = X_temp2[,j]*X_temp[,i]
      counter = counter + 1}
  }  
  return(list(X_temp,XX_temp))
}

gradient_sample = function(y1,y0,d1,d0,pi1,pi0,X,XX,beta,treated_covar,control_covar){
  gradient=rep(0,length(beta))
  for (i in 1:length(beta)){
    
    point_temp1=beta
    point_temp1[i]=beta[i]+0.005
    
    point_temp2=beta
    point_temp2[i]=beta[i]-0.005
    
    temp_up=Optimal_Logit_Z_Criterion(point_temp1,y1,y0,d1,d0,pi1,pi0,treated_covar ,control_covar,X,XX)
    temp_low=Optimal_Logit_Z_Criterion(point_temp2,y1,y0,d1,d0,pi1,pi0,treated_covar ,control_covar,X,XX)
    
    
    gradient[i]=(temp_up-temp_low)/0.01
    
  }
  return(gradient)
}
gradient_sample_Bound_num = function(y1,y0,d1,d0,pi1,pi0,X,XX,beta,treated_covar,control_covar){
  gradient=rep(0,length(beta))
  for (i in 1:length(beta)){
    point_temp1=beta
    point_temp1[i]=beta[i]+0.005
    
    point_temp2=beta
    point_temp2[i]=beta[i]-0.005
    
    temp_up=Optimal_Logit_M_Criterion_Bound_EST(point_temp1,y1,y0,d1,d0,pi1,pi0,treated_covar ,control_covar,X,XX)
    temp_low=Optimal_Logit_M_Criterion_Bound_EST(point_temp2,y1,y0,d1,d0,pi1,pi0,treated_covar ,control_covar,X,XX)
    
    
    gradient[i]=(temp_up-temp_low)/0.01
    
  }
  return(gradient)
}

gradient_sample_Bound = function(y1,y0,d1,d0,pi1,pi0,beta,treated_covar,control_covar){
  
  gradient=rep(0,length(beta))
  
  subject_size=length(y1)
  treat_util = treated_covar %*% beta
  control_util = control_covar %*% beta
  
  treat_impute = sapply(treat_util,logistic)
  control_impute = sapply(control_util,logistic)
  
  ######c'y
  residual1=rep(0,2*length(y1)) 
  residual1[2*(1:subject_size)-1]=-y0 * d0/pi0 
  residual1[2*(1:subject_size)]=y1 * d1/pi1 
  
  #####c'f
  residual2=rep(0,2*length(y1)) 
  residual2[2*(1:subject_size)-1]=-control_impute
  residual2[2*(1:subject_size)]= treat_impute

  Nablaf_treat = sweep(treated_covar,1,treat_impute*(1-treat_impute),'*')
  Nablaf_control=sweep(control_covar,1,-control_impute*(1-control_impute),'*')
  
  
  Nablaf=matrix(0,nrow=2*subject_size,ncol=ncol(Nablaf_treat))
  Nablaf[2*(1:subject_size)-1,]=Nablaf_control
  Nablaf[2*(1:subject_size),]=Nablaf_treat
  
  M_criterion1=0
  M_criterion2=0

  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    fo_t=fo[[c]]
    so_t=so[[c]]
    
    #extract variance bound
    cov_bound_t = so_AS2[[c]]
    
    #second order probabilities
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    cov_bound_zeros_t = (abs(cov_bound_t)<1e-2)
    
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    M_criterion1= M_criterion1 - 2 *residual1[index_start:index_end] %*% cov_bound_HT_t %*% Nablaf[index_start:index_end,]/subject_size
    M_criterion2= M_criterion2 + 2* residual2[index_start:index_end] %*% cov_bound_HT_t %*% Nablaf[index_start:index_end,]/subject_size
  }  
  gradient=M_criterion1 + M_criterion2
  return(gradient)
}
gradient_sample_deprecated = function(y1,y0,d1,d0,pi1,pi0,X,XX,beta,treated_covar,control_covar){
  
  #output wrong gradiant
  
  #code up the matrix first:
  subject_size=length(y1)
  
  
  ######probabilities
  treat_util = treated_covar %*% beta
  control_util = control_covar %*% beta
  treat_impute = sapply(treat_util,logistic)
  control_impute = sapply(control_util,logistic)
  p=rep(0,2*length(y1))
  p[2*(1:subject_size)-1]=control_impute
  p[2*(1:subject_size)]=treat_impute
  
  ######c'(y-f)
  residual1=y1-treat_impute
  residual0=y0-control_impute
  residual=rep(0,2*length(y1))
  residual[2*(1:subject_size)-1]=-residual0 * d0/pi0
  residual[2*(1:subject_size)]=residual1 * d1/pi1
  
  
  ##### prepare ingredients
  XX = sweep(XX,1,-p*(1-p)*(1-2*p),'*')
  Rdivpi=rep(0,2*length(d1))
  Rdivpi[2*(1:subject_size)-1]=d0/pi0
  Rdivpi[2*(1:subject_size)]=d1/pi1  
  
  X_right = sweep(X,1,p*(1-p)*Rdivpi,'*')
  X_left =sweep(X,1,p*(1-p),'*')
  
  term1 = rep(0,ncol(treated_covar)*ncol(treated_covar))
  term2 = matrix(0,ncol(treated_covar),ncol(treated_covar))
  term3 = rep(0,ncol(treated_covar))
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]      
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=sweep(so_t,1,1/fo_t,'*')
    normalized_cov=sweep(normalized_cov,2,1/fo_t,'*')
    
    term1= term1 + as.numeric(residual[index_start:index_end] %*% normalized_cov %*%  XX[index_start:index_end,])/subject_size
    
    term2 =term2 + t(X_left[index_start:index_end,]) %*% normalized_cov %*% X_right[index_start:index_end,]/subject_size
    
    term3 =  term3 + as.numeric(residual[index_start:index_end] %*% normalized_cov %*%  X_left[index_start:index_end,])/subject_size
    
    
    
  }    
  term1 = matrix(term1, ncol(X),ncol(X))
  gradient = -term3 %*% term1 - term3 %*% term2
  
  gradient = 2 * gradient
  return(gradient)
  
}

gradient_ppl = function(y1,y0,X,beta,treated_covar,control_covar){
  
  #code up the matrix first:
  subject_size=length(y1)
  
  
  ######probabilities
  treat_util = treated_covar %*% beta
  control_util = control_covar %*% beta
  treat_impute = sapply(treat_util,logistic)
  control_impute = sapply(control_util,logistic)
  p=rep(0,2*length(y1))
  p[2*(1:subject_size)-1]=treat_impute
  p[2*(1:subject_size)]=control_impute
  
  ######c'(y-f)
  residual1=y1-treat_impute
  residual0=y0-control_impute
  residual=rep(0,2*length(y1))
  residual[2*(1:subject_size)-1]=-residual0 
  residual[2*(1:subject_size)]=residual1 
  
  #the name f does not matter
  X_left =sweep(X,1,p*(1-p),'*')
  
  term1 = rep(0,ncol(treated_covar))
  
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]      
    fo_t=fo[[c]]
    so_t=so[[c]]
    normalized_cov=sweep(so_t,1,1/fo_t,'*')
    normalized_cov=sweep(normalized_cov,2,1/fo_t,'*')
    
    term1= term1 + as.numeric(residual[index_start:index_end] %*% normalized_cov %*%  X_left[index_start:index_end,])/subject_size
    
    
  }    
  term1=-term1*2
  return(term1)
}

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
    
    
    var_est = var_est + z_est_t %*% cov_bound_HT_t %*% z_est_t/(subject_size)^2
  }
  
  return(var_est)
  #return(c(avg_effect1-avg_effect0,coef))
}

Optimal_Logit_Est_3 = function(y,d1,d0,x,pi0,pi1){
  
  
  #Iterated Weighted Least Square
  
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x)
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  
  Maxiter=100
  iter=0
  tol=1e-4
  diff=1
  beta_old=rep(0,ncol(treated_covar))
  
  while (diff>tol & iter<Maxiter){
    
    beta_new=Solve_Logit_OLS_Est(y,d1,d0,pi0,pi1,treated_covar,control_covar,beta_old,X)
    
    beta_new=as.matrix(beta_new)
    
    print(beta_new)
    print((beta_new-beta_old)^2)
    diff= sum((beta_new-beta_old)^2)
    
    beta_old = 0.5 * beta_new + 0.5* beta_old
    
  }
  
  ind=0
  if (iter==Maxiter){
    ind=1
  }
  
  #optimal coefficients
  optimal_b=beta_old
  
  treat_util = treated_covar %*% optimal_b
  control_util = control_covar %*% optimal_b
  
  #predicition
  y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  #DR
  avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
  avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)    
  
  effect_est=avg_effect1-avg_effect0
  
  
  return(list(effect_est,optimal_b,ind))
}

Optimal_Logit_Est_4 = function(y,d1,d0,x,pi0,pi1,Maxiter=100){
  
  
  #Iterated Weighted Least Square followed by local search
  
  subject_size=length(y)
  
  #create the outcome vector
  y0_temp=y*d0 
  y1_temp=y*d1
  
  treated_ind=cbind(rep(0,subject_size),rep(1,subject_size))
  control_ind=cbind(rep(1,subject_size),rep(0,subject_size))
  
  treated_covar = cbind(treated_ind,x)
  control_covar = cbind(control_ind,x)
  
  auxiliary_matrices = second_order_matrix(treated_covar,control_covar)
  X=auxiliary_matrices[[1]]
  XX=auxiliary_matrices[[2]]
  
  iter=0
  tol=1e-4
  diff=1
  beta_old=rep(0,ncol(treated_covar))
  
  while (diff>tol & iter<Maxiter){
    
    beta_new=Solve_Logit_OLS_Est(y,d1,d0,pi0,pi1,treated_covar,control_covar,beta_old,X)
    
    beta_new=as.matrix(beta_new)
    
    print(beta_new)
    print((beta_new-beta_old)^2)
    diff= sum((beta_new-beta_old)^2)
    
    beta_old = 0.5 * beta_new + 0.5* beta_old
    
  }
  
  ind=0
  if (iter==Maxiter){
    ind=1
  }
  
  #optimal coefficients
  optimal_b=beta_old
  
  optimal_b = optim(par=optimal_b,fn=Optimal_Logit_Z_Criterion,gr=gradient_sample,method='Nelder-Mead',
                    y1=y1_temp,y0=y0_temp,d1=d1,d0=d0,pi1=pi1,pi0=pi0,treated_covar=treated_covar,control_covar=control_covar,X=X,XX=XX,
                    control=c(maxit=1000,trace=1))
  
  optimal_b=optimal_b$par
  treat_util = treated_covar %*% optimal_b
  control_util = control_covar %*% optimal_b
  
  #predicition
  y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  #DR
  avg_effect1 = mean(y_p1) +  mean((y-y_p1) * d1/pi1)
  avg_effect0 = mean(y_p0) +  mean((y-y_p0) * d0/pi0)    
  
  effect_est=avg_effect1-avg_effect0
  
  
  return(list(effect_est,optimal_b,ind))
}


Solve_Logit_OLS_Est_Bound=function(y,d1,d0,pi0,pi1,treated_covar,control_covar,beta0,X){
  
  
  treat_util = treated_covar %*% beta0
  control_util = control_covar %*% beta0
  

    
    #predicition
    y_p1 = sapply(treat_util,logistic)
  y_p0 = sapply(control_util,logistic)
  
  #first step regression  
  p=rep(0,2*length(y))
  p[2*(1:subject_size)-1]=y_p1
  p[2*(1:subject_size)]=y_p0
  
  d=rep(0,2*length(d1))
  d[2*(1:subject_size)-1]=d1
  d[2*(1:subject_size)]=d0
  
  #create covariate matrix
  X_right=sweep(X,1,p*(1-p),'*')
  X_left=sweep(X,1,p*(1-p),'*')
  
  #create the outcome vector
  y_temp=rep(0,subject_size*2)
  y_temp[2*(1:subject_size)-1]=-(y)*d0/pi0 +y_p0#make sure the sign is correct: it should be negative
  y_temp[2*(1:subject_size)]=(y)*d1/pi1 -y_p1
  
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
    #extract variance bound
    cov_bound_t = so_AS2[[c]]
    
    #second order probabilities
    p_t = so_t + (fo_t) %*% t(fo_t) 
    p_t=round(p_t,10)
    
    cov_bound_zeros_t = (abs(cov_bound_t)<1e-2)
    
    p_denomiator_t = p_t + cov_bound_zeros_t
    cov_bound_HT_t= cov_bound_t / p_denomiator_t
    
    XDX = XDX + t(X_left[index_start:index_end,]) %*% cov_bound_HT_t %*% X_right[index_start:index_end,]/subject_size
    
    XDY = XDY + t(X_left[index_start:index_end,]) %*% cov_bound_HT_t %*% y_temp[index_start:index_end]/subject_size
    
    
  }
  
  XDX = XDX
  optimal_b = solve(XDX,XDY)
  
  return(optimal_b)
}

# grid=expand.grid(-5:5,-5:5)  
# 
# landscape = function(x){
#   
#   result=Optimal_Logit_Z_Criterion(x,y1,y0,d1,d0,pi1,pi0,treated_covar,control_covar,X,XX)
#   print(x)
#   return(result) 
# }
# 
# temp=apply(grid,1,landscape)
# temp=matrix(temp,11,11)
# write.csv(temp,'result.csv')
