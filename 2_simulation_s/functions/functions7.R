SIM_ONERUN_IMPUTED_AUG20_AS2=function(status,pol_size,prob,y1,y0,pi1,pi0,Y_all,index,expo1,expo2,option,x){
  
  print('here')
  est_sim=rep(0,11)
  var_sim=rep(0,11) 
  var_sim2=rep(0,11) 
  
  #generate random assignments
  if (option=='Bernoulli'){
    realized_assignment=sample(status,pol_size,replace=TRUE,prob=prob)
  }else if (option=='Stratified'){
    
    temp=cbind(net %v% "group",net %v% "assignment")
    realized_assignment=ASSIGNMENT(temp)
    
  }
  
  ############################################
  ##########calculate exposures###############
  ############################################
  
  exposure = matrix(0,nrow=pol_size,ncol=7)
  exposure[,1]=1:pol_size #vertex indices
  exposure[,2]=realized_assignment==1
  exposure[,3]=realized_assignment==2
  exposure[,4]=realized_assignment==3
  exposure[,5]=realized_assignment==4
  num_friend_treated = apply(friends,1,FRIEND_TREATED,realized_assignment)
  exposure[,6]=num_friend_treated[1,]
  exposure[,7]=num_friend_treated[2,]  
  
  
  exposure_t = exposure[subjects_t,2:7]
  
  ###############################
  #please double check every time 
  if (is.null(dim(exposure_t))){
    state_t=INDIVIDUAL_EXPOSURE(exposure_t)
    state_t = state_t[c(expo1,expo2)]
    y0_t = state_t[expo1] * y0
    y1_t = state_t[expo2] * y1
    denom_pi0_t= state_t[expo1]* 1/pi0 
    denom_pi1_t= state_t[expo2]* 1/pi1 
    pi_weight_t = denom_pi0_t + denom_pi1_t
    d1_t=state_t[expo2]
    d0_t=state_t[expo1]
  }else{
    state_t=apply(exposure_t,1,INDIVIDUAL_EXPOSURE) 
    #y_observed=apply(t(state_t) * Y_all[,1:12],1,sum)
    
    state_t = state_t[c(expo1,expo2),]
    
    y0_t = state_t[1,] * y0
    y1_t = state_t[2,] * y1
    denom_pi0_t= state_t[1,]* 1/pi0 
    denom_pi1_t= state_t[2,]* 1/pi1 
    pi_weight_t = denom_pi0_t + denom_pi1_t
    d1_t=state_t[2,]
    d0_t=state_t[1,]
  }
  
  print(sum(is.na((pi1))))
  
  print(sum(is.na(sum(pi0))))
  y_obs_t=y0_t+y1_t
  
  #HT, HA, OLS and WLS
  est_sim[1]= HT_full(y1_t,y0_t,pi1,pi0,subject_size)
  est_sim[2]= HA_full(y1_t,y0_t,pi1,pi0,denom_pi1_t,denom_pi0_t)
  print('here')  
  OLS_result = OLS_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint')
  WLS_result = OLS_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint')  
  est_sim[3]= OLS_result[[1]]
  est_sim[4]= WLS_result[[1]]
  
  WLS_p1=WLS_result[[2]]
  WLS_p0=WLS_result[[3]]  
  
  OLS_p1=OLS_result[[2]]
  OLS_p0=OLS_result[[3]]    
  
  print('NO_HARM_Linear')  
  WLS_No_HARM=NO_HARM_AS2(y_obs_t,d1_t,d0_t,pi0,pi1,WLS_p1,WLS_p0)
  
  est_sim[[5]]=WLS_No_HARM[[1]]
  
  #Logit
  #Logit_u = Logit_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi0,pi1,mode='unweighted_joint')
  Logit_w = Logit_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi0,pi1,mode='unweighted_joint')
  est_sim[6]= Logit_w[[1]]
  #  est_sim[6]= Logit_w[[2]]
  Logit_p1=Logit_w[[2]]
  Logit_p0=Logit_w[[3]]  
  
  print('NO_HARM_Logit')  
  Logit_No_HARM=NO_HARM_AS2(y_obs_t,d1_t,d0_t,pi0,pi1,Logit_p1,Logit_p0)
  
  est_sim[[7]]=Logit_No_HARM[[1]]
  
  print('Optimal_Linear')
  print('here')
  ##Optimal Linear Adjustment
  result_temp=Optimal_Reg_Est(y_obs_t,d1_t,d0_t,x,pi0,pi1)
  
  est_sim[8]= result_temp[[1]]
  optimal_b_fixed =result_temp[[2]]  
  
  
  ##Optimal Linear Adjustment
  result_temp=Optimal_Reg_Est(y_obs_t,d1_t,d0_t,x,pi0,pi1)
  
  est_sim[8]= result_temp[[1]]
  optimal_b_fixed =result_temp[[2]]    
  
  print('Optimal Linear Generated Regressor')
  result_optimal_linear_genereated=NO_HARM2_AS2(y_obs_t,d1_t,d0_t,pi0,pi1,OLS_p1,OLS_p0)
  est_sim[9]=result_optimal_linear_genereated[[1]]
  
  print('Optimal Logit Generated Regressor')
  result_optimal_logit_genereated=NO_HARM2_AS2(y_obs_t,d1_t,d0_t,pi0,pi1,Logit_p1,Logit_p0)
  est_sim[10]=result_optimal_logit_genereated[[1]]

  #result_temp_Logistic_Bound=Optimal_Logit_Est_C(y_obs_t,d1_t,d0_t,x,pi0,pi1,1000,F,'BFGS')
  #est_sim[11]= result_temp_Logistic_Bound[[1]]
  #optimal_b_Logit_Bound =result_temp_Logistic_Bound[[2]] 
  
  
  #result_temp_Logistic_MS=Optimal_Logit_Est_C(y_obs_t,d1_t,d0_t,x,pi0,pi1,10000,F,'Nelder-Mead',optimal_b_fixed)
  #est_sim[9]= result_temp_Logistic_MS[[1]]
  #optimal_b_Logit_C =result_temp_Logistic_MS[[2]]       
  
  #result_temp_Logistic_mS=Optimal_Logit_Est_B(y_obs_t,d1_t,d0_t,x,pi0,pi1,5,T,'Nelder-Mead')
  #est_sim[10]= result_temp_Logistic_mS[[1]]
  #optimal_b_Logit_B =result_temp_Logistic_mS[[2]]       


  
  
  #variance bound estimation
  var_sim[1]= HT_var_full_AS2(y_obs_t,d1_t,d0_t,subject_size,fo,so,contrast,dict_group)
  var_sim[2]= HA_var_full_AS2(y_obs_t,d1_t,d0_t,pi_weight_t,subject_size,fo,so,contrast,dict_group)
  print('OLS')  
  var_sim[3]= OLS_var_full_AS2(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint',fo,so,contrast,dict_group)
  print('WLS')
  var_sim[4]= OLS_var_full_AS2(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint',fo,so,contrast,dict_group)
  
  var_sim[5]=  WLS_No_HARM[[2]]
  print('logit')
  var_sim[6]= Logit_var_full_AS2(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='weighted_joint',fo,so,contrast,dict_group)
  var_sim[7]=  Logit_No_HARM[[2]]  
  print('Optimal_Linear')  
  
  var_sim[8]= Optimal_Reg_var_full_AS2(y_obs_t,d1_t,d0_t,x,optimal_b_fixed)  
  
  #print('Optimal_Logit') 
  var_sim[9]= result_optimal_linear_genereated[[2]]

  #var_sim[9]=Optimal_Logit_var_full_AS2(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_C) 
  var_sim[10]=result_optimal_logit_genereated[[2]]
  
  #variance bound estimation
  var_sim2[1]= HT_var_full(y_obs_t,d1_t,d0_t,subject_size,fo,so,contrast,dict_group)
  var_sim2[2]= HA_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,subject_size,fo,so,contrast,dict_group)
  print('OLS')  
  var_sim2[3]= OLS_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint',fo,so,contrast,dict_group)
  print('WLS')
  var_sim2[4]= OLS_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint',fo,so,contrast,dict_group)
  
  print('logit')
  var_sim2[6]= Logit_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='weighted_joint',fo,so,contrast,dict_group)
  print('Optimal_Linear')  
  
  var_sim2[8]= Optimal_Reg_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_fixed)  
  
  #print('Optimal_Logit') 
  #var_sim2[9]= Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_MS)  
 # var_sim2[9]=Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_C) 
  #var_sim2[10]=Optimal_Logit_var_full_AS2(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_B) 
  
  #var_sim2[11]=Optimal_Reg_var_full_finite_AS3(y_obs_t,d1_t,d0_t,x,optimal_b_fixed) 
  
  print('NO_HARM_Linear')  
  WLS_No_HARM_AS=NO_HARM(y_obs_t,d1_t,d0_t,pi0,pi1,WLS_p1,WLS_p0)
  
  var_sim2[5]=  WLS_No_HARM_AS[[2]]
  
  print('NO_HARM_Logit')  
  Logit_No_HARM_AS=NO_HARM(y_obs_t,d1_t,d0_t,pi0,pi1,Logit_p1,Logit_p0)
  var_sim2[7]=  Logit_No_HARM_AS[[2]]  
  
  
  return(list(est_sim,var_sim,var_sim2))
}
SIM_ONERUN_IMPUTED_Oct21_AS2=function(status,pol_size,prob,y1,y0,pi1,pi0,Y_all,index,expo1,expo2,option){
  
  est_sim=rep(0,11)
  var_sim=rep(0,11) 
  var_sim2=rep(0,11)
  if (option=='Bernoulli'){
    realized_assignment=sample(status,pol_size,replace=TRUE,prob=prob)
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
  
  
  exposure_t = exposure[subjects_t,2:7]
  
  #please double check every time 
  if (is.null(dim(exposure_t))){
    state_t=individual_exposure(exposure_t)
    state_t = state_t[c(expo1,expo2)]
    y0_t = state_t[expo1] * y0
    y1_t = state_t[expo2] * y1
    denom_pi0_t= state_t[expo1]* 1/pi0 
    denom_pi1_t= state_t[expo2]* 1/pi1 
    pi_weight_t = denom_pi0_t + denom_pi1_t
    d1_t=state_t[expo2]
    d0_t=state_t[expo1]
  }else{
    state_t=apply(exposure_t,1,INDIVIDUAL_EXPOSURE) 
    y_observed=apply(t(state_t) * Y_all[,1:12],1,sum)
    
    state_t = state_t[c(expo1,expo2),]
    
    y0_t = state_t[1,] * y0
    y1_t = state_t[2,] * y1
    denom_pi0_t= state_t[1,]* 1/pi0 
    denom_pi1_t= state_t[2,]* 1/pi1 
    pi_weight_t = denom_pi0_t + denom_pi1_t
    d1_t=state_t[2,]
    d0_t=state_t[1,]
  }
  
  
  y_obs_t=y0_t+y1_t
  
  #HT, HA, OLS and WLS
  est_sim[1]= HT_full(y1_t,y0_t,pi1,pi0,subject_size)
  est_sim[2]= HA_full(y1_t,y0_t,pi1,pi0,denom_pi1_t,denom_pi0_t)
  print('here')  
  OLS_result = OLS_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint')
  WLS_result = OLS_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint')  
  est_sim[3]= OLS_result[[1]]
  est_sim[4]= WLS_result[[1]]
  
  WLS_p1=WLS_result[[2]]
  WLS_p0=WLS_result[[3]]  
  print('NO_HARM_Linear')  
  WLS_No_HARM=NO_HARM_AS2(y_obs_t,d1_t,d0_t,pi0,pi1,WLS_p1,WLS_p0)
  
  est_sim[[5]]=WLS_No_HARM[[1]]
  
  #Logit
  #Logit_u = Logit_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi0,pi1,mode='unweighted_joint')
  Logit_w = Logit_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi0,pi1,mode='weighted_joint')
  est_sim[6]= Logit_w[[1]]
  #  est_sim[6]= Logit_w[[2]]
  Logit_p1=Logit_w[[2]]
  Logit_p0=Logit_w[[3]]  
  
  print('NO_HARM_Logit')  
  Logit_No_HARM=NO_HARM_AS2(y_obs_t,d1_t,d0_t,pi0,pi1,Logit_p1,Logit_p0)
  
  est_sim[[7]]=Logit_No_HARM[[1]]
  
  print('Optimal_Linear')
  ##Optimal Linear Adjustment
  result_temp=Optimal_Reg_Est(y_obs_t,d1_t,d0_t,x,pi0,pi1)
  
  est_sim[8]= result_temp[[1]]
  optimal_b_fixed =result_temp[[2]]  
  optimal_p1=result_temp[[3]]
  optimal_p0=result_temp[[4]]    
  
  print('Optimal_Logit')
  result_temp_Logistic_MS=Optimal_Logit_Est_A(y_obs_t,d1_t,d0_t,x,pi0,pi1,100,F)
  est_sim[9]= result_temp_Logistic_MS[[1]]
  optimal_b_Logit_MS =result_temp_Logistic_MS[[2]]       
  
  result_temp_Logistic_mS=Optimal_Logit_Est_A(y_obs_t,d1_t,d0_t,x,pi0,pi1,10,F)
  est_sim[10]= result_temp_Logistic_mS[[1]]
  optimal_b_Logit_mS =result_temp_Logistic_mS[[2]]       
  
  print('NO_HARM_Optimal_Linear')
  ##Optimal Linear Adjustment
  result_temp_NHOL=NO_HARM_AS2(y_obs_t,d1_t,d0_t,pi0,pi1,optimal_p1,optimal_p0)
  
  est_sim[11]= result_temp_NHOL[[1]]
  
  #variance bound estimation
  var_sim[1]= HT_var_full_AS2(y_obs_t,d1_t,d0_t,subject_size,fo,so,contrast,dict_group)
  var_sim[2]= HA_var_full_AS2(y_obs_t,d1_t,d0_t,pi_weight_t,subject_size,fo,so,contrast,dict_group)
  print('OLS')  
  var_sim[3]= OLS_var_full_AS2(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint',fo,so,contrast,dict_group)
  print('WLS')
  var_sim[4]= OLS_var_full_AS2(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint',fo,so,contrast,dict_group)
  
  var_sim[5]=  WLS_No_HARM[[2]]
  print('logit')
  var_sim[6]= Logit_var_full_AS2(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='weighted_joint',fo,so,contrast,dict_group)
  var_sim[7]=  Logit_No_HARM[[2]]  
  print('Optimal_Linear')  
  
  var_sim[8]= Optimal_Reg_var_full_AS2(y_obs_t,d1_t,d0_t,x,optimal_b_fixed)  
  
  print('Optimal_Logit') 
  var_sim[9]= Optimal_Logit_var_full_AS2(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_MS)  
  var_sim[10]=Optimal_Logit_var_full_AS2(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_mS) 
  var_sim[11]=Optimal_Reg_var_full_finite_AS3(y_obs_t,d1_t,d0_t,x,optimal_b_fixed) 
  
  #variance bound estimation
  var_sim2[1]= HT_var_full(y_obs_t,d1_t,d0_t,subject_size,fo,so,contrast,dict_group)
  var_sim2[2]= HA_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,subject_size,fo,so,contrast,dict_group)
  print('OLS')  
  var_sim2[3]= OLS_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint',fo,so,contrast,dict_group)
  print('WLS')
  var_sim2[4]= OLS_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint',fo,so,contrast,dict_group)
  
  print('logit')
  var_sim2[6]= Logit_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='weighted_joint',fo,so,contrast,dict_group)
  print('Optimal_Linear')  
  
  var_sim2[8]= Optimal_Reg_var_full_AS2(y_obs_t,d1_t,d0_t,x,optimal_b_fixed)  
  
  print('Optimal_Logit') 
  var_sim2[9]= Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_MS)  
  var_sim2[10]=Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_mS) 
  #var_sim2[11]=Optimal_Reg_var_full_finite_AS3(y_obs_t,d1_t,d0_t,x,optimal_b_fixed) 
  
  print('NO_HARM_Linear')  
  WLS_No_HARM_AS=NO_HARM(y_obs_t,d1_t,d0_t,pi0,pi1,WLS_p1,WLS_p0)
  
  var_sim2[5]=  WLS_No_HARM_AS[[2]]
  
  print('NO_HARM_Logit')  
  Logit_No_HARM_AS=NO_HARM(y_obs_t,d1_t,d0_t,pi0,pi1,Logit_p1,Logit_p0)
  var_sim2[7]=  Logit_No_HARM_AS[[2]]  
  
  
  return(list(est_sim,var_sim,var_sim2))
}

sim_onerun_sharpnull=function(status,pol_size,prob,y1,y0,pi1,pi0,Y_all,index,expo1,expo2){
  
  est_sim=rep(0,15)
  var_sim=rep(0,15) 
  realized_assignment_t=sample(status,pol_size,replace=TRUE,prob=prob)
  exposure_t = matrix(0,nrow=pol_size,ncol=7)
  exposure_t[,1]=1:pol_size #vertex indices
  exposure_t[,2]=realized_assignment_t==0
  exposure_t[,3]=realized_assignment_t==1
  exposure_t[,4]=realized_assignment_t==2
  exposure_t[,5]=realized_assignment_t==3
  
  num_friend_treated_t = apply(friends,1,friend_treated,realized_assignment_t)
  exposure_t[,6]=num_friend_treated_t[1,]
  exposure_t[,7]=num_friend_treated_t[2,]  
  
  exposure_t = exposure_t[subjects_t,2:7] #compute exposure
  
  #please double check every time 
  if (is.null(dim(exposure_t))){
    state_t=individual_exposure(exposure_t)
    state_t = state_t[c(expo1,expo2)]
    y0_t = state_t[expo1] * y0
    y1_t = state_t[expo2] * y1
    denom_pi0_t= state_t[expo1]* 1/pi0 
    denom_pi1_t= state_t[expo2]* 1/pi1 
    pi_weight_t = denom_pi0_t + denom_pi1_t
    d1_t=state_t[expo2]
    d0_t=state_t[expo1]
  }else{
    state_t=apply(exposure_t,1,individual_exposure) 
    y_observed=apply(t(state_t) * Y_all[,2:13],1,sum)
    
    write.csv(y_observed,paste0('/home/hc654/palmer_scratch/Unified_Complete_Simulation_Sharpnull/',expo1,expo2,'_',index,'_Sim_est.csv'))
    state_t = state_t[c(expo1,expo2),]
    
    y0_t = state_t[1,] * y0
    y1_t = state_t[2,] * y1
    denom_pi0_t= state_t[1,]* 1/pi0 
    denom_pi1_t= state_t[2,]* 1/pi1 
    pi_weight_t = denom_pi0_t + denom_pi1_t
    d1_t=state_t[2,]
    d0_t=state_t[1,]
  }
  
  sample_size[1,i]=sum(d1_t)
  sample_size[2,i]=sum(d0_t)  
  y_obs_t=y0_t+y1_t
  
  # est_sim[1]= HT(y2,y1,pi2,pi1,subject_size)
  # est_sim[2]= HA(y2,y1,pi2,pi1,denom_pi2,denom_pi1)
  # est_sim[3]= OLS(y_obs,d1,d0,pi_weight,x,mode='OLS_joint')
  # est_sim[4]= OLS(y_obs,d1,d0,pi_weight,x,mode='OLS_separate')    
  # est_sim[5]= OLS(y_obs,d1,d0,pi_weight,x,mode='WLS_joint')    
  # est_sim[6]= OLS(y_obs,d1,d0,pi_weight,x,mode='WLS_separate')    
  # est_sim[7]= Logit(y_obs,d1,d0,pi_weight,x,pi1,pi2,mode='unweighted_joint')
  # est_sim[8]= Logit(y_obs,d1,d0,pi_weight,x,pi1,pi2,mode='unweighted_separate')
  # est_sim[9]= Logit(y_obs,d1,d0,pi_weight,x,pi1,pi2,mode='weighted_joint')
  # est_sim[10]= Logit(y_obs,d1,d0,pi_weight,x,pi1,pi2,mode='weighted_separate')
  
  
  #HT, HA, OLS and WLS
  est_sim[1]= HT_full(y1_t,y0_t,pi1,pi0,subject_size)
  est_sim[2]= HA_full(y1_t,y0_t,pi1,pi0,denom_pi1_t,denom_pi0_t)
  
  OLS_result = OLS_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint')
  WLS_result = OLS_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint')  
  est_sim[3]= OLS_result[[2]]
  est_sim[4]= WLS_result[[2]]
  
  
  #DR
  est_sim[7]= OLS_DR(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi0,pi1,mode='unweighted_joint')
  
  #Logit
  Logit_u = Logit_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi0,pi1,mode='unweighted_joint')
  Logit_w = Logit_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi0,pi1,mode='weighted_joint')
  est_sim[5]= Logit_u[[2]]
  est_sim[6]= Logit_w[[2]]
  
  
  ##Optimal Linear Adjustment
  result_temp=Optimal_Reg_Est(y_obs_t,d1_t,d0_t,x,pi0,pi1)
  
  ##Using a variable denominator
  est_sim[8]= result_temp[[1]]
  optimal_b =result_temp[[2]]  
  
  ##Using the fixed demoniator
  est_sim[9]= result_temp[[3]]
  optimal_b_fixed =result_temp[[4]]  
  
  
  ##Optimal Logistic Adjustment, using unweighted logit regression as initial values
  result_temp_Logistic=Optimal_Logit_Est(y_obs_t,d1_t,d0_t,x,pi0,pi1,Logit_u[[1]])
  est_sim[10]= result_temp_Logistic[[1]]
  optimal_b_Logit =result_temp_Logistic[[2]]  
  
  ##Optimal Logistic Adjustment, using the one step procedure to generate initial values
  result_temp_Logistic_OS=Optimal_Logit_Est_A(y_obs_t,d1_t,d0_t,x,pi0,pi1,1,T)
  est_sim[11]= result_temp_Logistic_OS[[1]]
  optimal_b_Logit_OS =result_temp_Logistic_OS[[2]]   
  
  ##Optimal Logistic Adjustment, using multiple step procedure to generate initial values
  result_temp_Logistic_MS=Optimal_Logit_Est_A(y_obs_t,d1_t,d0_t,x,pi0,pi1,100,F)
  est_sim[12]= result_temp_Logistic_MS[[1]]
  optimal_b_Logit_MS =result_temp_Logistic_MS[[2]]       
  optimization_ind[i]=result_temp_Logistic_MS[[3]] 
  
  ##Optimal Logistic Adjustment, using the multiple step procedure to generate initial values, followed by local search
  result_temp_Logistic_MSO=Optimal_Logit_Est_A(y_obs_t,d1_t,d0_t,x,pi0,pi1,100,T)
  est_sim[13]= result_temp_Logistic_MSO[[1]]
  optimal_b_Logit_MSO =result_temp_Logistic_MSO[[2]]       
  
  
  ##Optimal Logistic Adjustment, using multiple step procedure to generate initial values
  result_temp_Logistic_mS=Optimal_Logit_Est_A(y_obs_t,d1_t,d0_t,x,pi0,pi1,10,F)
  est_sim[14]= result_temp_Logistic_mS[[1]]
  optimal_b_Logit_mS =result_temp_Logistic_mS[[2]]       
  
  ##Optimal Logistic Adjustment, using the multiple step procedure to generate initial values, followed by local search
  #result_temp_Logistic_mSO=Optimal_Logit_Est_A(y_obs_t,d1_t,d0_t,x,pi0,pi1,10,T)
  #est_sim[15]= result_temp_Logistic_mSO[[1]]
  #optimal_b_Logit_mSO =result_temp_Logistic_mSO[[2]]   
  
  coef_sim_OLS=OLS_result[[1]]
  coef_sim_WLS=WLS_result[[1]]
  coef_sim_Logit[i,]=Logit_u[[1]]
  coef_sim_WLogit[i,]=Logit_w[[1]]
  coef_sim_OptimalOLS[i,]=optimal_b_fixed
  coef_sim_OptimalLogit[i,]=result_temp_Logistic[[2]]
  coef_sim_OptimalLogitOS[i,]=result_temp_Logistic_OS[[2]] 
  coef_sim_OptimalLogitMS[i,]=result_temp_Logistic_MS[[2]]   
  coef_sim_OptimalLogitMSO[i,]=result_temp_Logistic_MSO[[2]]    
  coef_sim_OptimalLogitmS[i,]=result_temp_Logistic_mS[[2]]   
  coef_sim_OptimalLogitmSO[i,]=result_temp_Logistic_mSO[[2]]  
  #variance bound estimation
  var_sim[1]= HT_var_full(y_obs_t,d1_t,d0_t,subject_size,fo,so,contrast,dict_group)
  var_sim[2]= HA_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,subject_size,fo,so,contrast,dict_group)
  var_sim[3]= OLS_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint',fo,so,contrast,dict_group)
  var_sim[4]= OLS_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint',fo,so,contrast,dict_group)
  var_sim[5]= Logit_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='unweighted_joint',fo,so,contrast,dict_group)
  var_sim[6]= Logit_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='weighted_joint',fo,so,contrast,dict_group)
  var_sim[7]= OLS_DR_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='unweighted_joint',fo,so,contrast,dict_group)
  var_sim[8]= Optimal_Reg_var_full(y_obs_t,d1_t,d0_t,x,optimal_b)
  var_sim[9]= Optimal_Reg_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_fixed)  
  var_sim[10]=Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit)  
  var_sim[11]=Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_OS)  
  var_sim[12]=Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_MS)  
  var_sim[13]=Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_MSO)    
  var_sim[14]=Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_mS)  
  var_sim[15]=Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_mSO)
  
  return(list(est_sim,var_sim))
}

sim_onerun_imputed_rerand=function(status,pol_size,prob,y1,y0,pi1,pi0,expo1,expo2){
  
  est_sim=rep(0,6)
  var_sim=rep(0,6) 
  realized_assignment_t=sample(status,pol_size,replace=TRUE,prob=prob)
  exposure_t = matrix(0,nrow=pol_size,ncol=7)
  exposure_t[,1]=1:pol_size #vertex indices
  exposure_t[,2]=realized_assignment_t==0
  exposure_t[,3]=realized_assignment_t==1
  exposure_t[,4]=realized_assignment_t==2
  exposure_t[,5]=realized_assignment_t==3
  
  num_friend_treated_t = apply(friends,1,friend_treated,realized_assignment_t)
  exposure_t[,6]=num_friend_treated_t[1,]
  exposure_t[,7]=num_friend_treated_t[2,]  
  
  exposure_t = exposure_t[subjects_t,2:7] #compute exposure
  
  #please double check every time 
  if (is.null(dim(exposure_t))){
    state_t=individual_exposure(exposure_t)
    state_t = state_t[c(expo1,expo2)]
    y0_t = state_t[expo1] * y0
    y1_t = state_t[expo2] * y1
    denom_pi0_t= state_t[expo1]* 1/pi0 
    denom_pi1_t= state_t[expo2]* 1/pi1 
    pi_weight_t = denom_pi0_t + denom_pi1_t
    d1_t=state_t[expo2]
    d0_t=state_t[expo1]
  }else{
    state_t=apply(exposure_t,1,individual_exposure) 
    #y_observed=apply(t(state_t) * Y_all[,2:13],1,sum)
    
    #write.csv(y_observed,paste0('/home/hc654/palmer_scratch/Unified_Complete_Simulation_Imputed_Aug9/',expo1,expo2,'_',index,'_Sim_est.csv'))
    state_t = state_t[c(expo1,expo2),]
    
    y0_t = state_t[1,] * y0
    y1_t = state_t[2,] * y1
    denom_pi0_t= state_t[1,]* 1/pi0 
    denom_pi1_t= state_t[2,]* 1/pi1 
    pi_weight_t = denom_pi0_t + denom_pi1_t
    d1_t=state_t[2,]
    d0_t=state_t[1,]
  }
  
  sample_size[1,i]=sum(d1_t)
  sample_size[2,i]=sum(d0_t)  
  y_obs_t=y0_t+y1_t
  
  
  #HT, HA, OLS and WLS
  est_sim[1]= HT_full(y1_t,y0_t,pi1,pi0,subject_size)
  est_sim[2]= HA_full(y1_t,y0_t,pi1,pi0,denom_pi1_t,denom_pi0_t)
  
  OLS_result = OLS_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint')
  WLS_result = OLS_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint')  
  #est_sim[3]= OLS_result[[2]]
  est_sim[3]= WLS_result[[2]]
  
  
  #DR
  #est_sim[7]= OLS_DR(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi0,pi1,mode='unweighted_joint')
  
  #Logit
  #Logit_u = Logit_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi0,pi1,mode='unweighted_joint')
  Logit_w = Logit_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi0,pi1,mode='weighted_joint')
  #est_sim[5]= Logit_u[[2]]
  est_sim[4]= Logit_w[[2]]
  
  
  ##Optimal Linear Adjustment
  time.start=Sys.time()
  result_temp=Optimal_Reg_Est(y_obs_t,d1_t,d0_t,x,pi0,pi1)
  print('Optimal_OLS')
  time.end=Sys.time()
  print(time.start-time.end)  
  ##Using a variable denominator
  #est_sim[8]= result_temp[[1]]
  #optimal_b =result_temp[[2]]  
  
  ##Using the fixed demoniator
  est_sim[5]= result_temp[[3]]
  optimal_b_fixed =result_temp[[4]]  
  
  
  ##Optimal Logistic Adjustment, using unweighted logit regression as initial values
  #  result_temp_Logistic=Optimal_Logit_Est(y_obs_t,d1_t,d0_t,x,pi0,pi1,Logit_u[[1]])
  #  est_sim[10]= result_temp_Logistic[[1]]
  #  optimal_b_Logit =result_temp_Logistic[[2]]  
  
  ##Optimal Logistic Adjustment, using the one step procedure to generate initial values
  #  result_temp_Logistic_OS=Optimal_Logit_Est_A(y_obs_t,d1_t,d0_t,x,pi0,pi1,1,T)
  #  est_sim[11]= result_temp_Logistic_OS[[1]]
  #  optimal_b_Logit_OS =result_temp_Logistic_OS[[2]]   
  
  ##Optimal Logistic Adjustment, using multiple step procedure to generate initial values
  time.start=Sys.time()
  result_temp_Logistic_MS=Optimal_Logit_Est_C(y_obs_t,d1_t,d0_t,x,pi0,pi1,100)
  est_sim[6]= result_temp_Logistic_MS[[1]]
  optimal_b_Logit_MS =result_temp_Logistic_MS[[2]]       
  optimization_ind[i]=result_temp_Logistic_MS[[3]] 
  
  #est_sim[11]= result_temp_Logistic_MS[[4]]
  #optimal_b_Logit_MSO =result_temp_Logistic_MS[[5]] 
  
  time.end=Sys.time()
  print('MS')
  print(time.start-time.end)
  
  ##Optimal Logistic Adjustment, using multiple step procedure to generate initial values
  #  result_temp_Logistic_mS=Optimal_Logit_Est_A(y_obs_t,d1_t,d0_t,x,pi0,pi1,10,F)
  #  est_sim[12]= result_temp_Logistic_mS[[1]]
  #  optimal_b_Logit_mS =result_temp_Logistic_mS[[2]]       
  
  ##Optimal Logistic Adjustment, using the multiple step procedure to generate initial values, followed by local search
  #  result_temp_Logistic_mSO=Optimal_Logit_Est_A(y_obs_t,d1_t,d0_t,x,pi0,pi1,10,T)
  #  est_sim[15]= result_temp_Logistic_mSO[[1]]
  #  optimal_b_Logit_mSO =result_temp_Logistic_mSO[[2]]   
  
  #variance bound estimation
  time.start=Sys.time()
  var_sim[1]= HT_var_full(y_obs_t,d1_t,d0_t,subject_size,fo,so,contrast,dict_group)
  print('HT_var')
  time.end=Sys.time()
  print(time.start-time.end)   
  
  var_sim[2]= HA_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,subject_size,fo,so,contrast,dict_group)
  #var_sim[3]= OLS_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint',fo,so,contrast,dict_group)
  var_sim[3]= OLS_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint',fo,so,contrast,dict_group)
  #var_sim[5]= Logit_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='unweighted_joint',fo,so,contrast,dict_group)
  var_sim[4]= Logit_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='weighted_joint',fo,so,contrast,dict_group)
  #var_sim[7]= OLS_DR_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='unweighted_joint',fo,so,contrast,dict_group)
  
  time.start=Sys.time()
  #var_sim[8]= Optimal_Reg_var_full(y_obs_t,d1_t,d0_t,x,optimal_b)
  var_sim[5]= Optimal_Reg_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_fixed)  
  print('ols_var')
  time.end=Sys.time()
  print(time.start-time.end) 
  #  var_sim[10]=Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit)  
  #  var_sim[11]=Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_OS)  
  time.start=Sys.time()
  var_sim[6]=Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_MS)  
  #var_sim[11]=Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_MSO)    
  print('Logit_var')
  time.end=Sys.time()
  print(time.start-time.end)    #var_sim[12]=Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_mS)  
  #  var_sim[15]=Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_mSO)
  
  return(list(est_sim,var_sim))
}
