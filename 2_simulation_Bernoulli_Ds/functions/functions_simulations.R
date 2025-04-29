SIM_ONERUN_IMPUTED_AUG20_AS2=function(status,pol_size,prob,y1,y0,pi1,pi0,index,expo1,expo2,option,x,subjects_t,so_AS2){
  
  print('Simulation')
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
  
  ##########################################################################
  ####################Observed Y############################################
  ##########################################################################
  
  y_obs_t=y0_t+y1_t

  subject_size=length(y_obs_t)
  ##########################################################################
  ####################Estimators############################################
  ##########################################################################
  
  
  #HT, HA, OLS and WLS
  print('HT')
  est_sim[1]= HT_full(y1_t,y0_t,pi1,pi0,subject_size)
  print('HA')
  est_sim[2]= HA_full(y1_t,y0_t,pi1,pi0,denom_pi1_t,denom_pi0_t)
  print('OLS')  
  OLS_result = OLS_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint')
  print('WLS')
  WLS_result = OLS_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint')  
  est_sim[3]= OLS_result[[1]]
  est_sim[4]= WLS_result[[1]]
  
  #####################################
  ######Imputed values#################
  #####################################
  WLS_p1=WLS_result[[2]]
  WLS_p0=WLS_result[[3]]  
  
  OLS_p1=OLS_result[[2]]
  OLS_p0=OLS_result[[3]]    
  
  
  #####################################
  ######Imputed values#################
  #####################################
  
  print('NO_HARM_Linear')  
  WLS_No_HARM=NO_HARM_AS2(y_obs_t,d1_t,d0_t,pi0,pi1,OLS_p1,OLS_p0,so_AS2)
  
  est_sim[[5]]=WLS_No_HARM[[1]]
  
  #####################################
  ######LOGIT##########################
  #####################################
  print('Logit')
  Logit_u = Logit_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,pi0,pi1,mode='unweighted_joint')
  est_sim[6]= Logit_u[[1]]
  Logit_p1=Logit_u[[2]]
  Logit_p0=Logit_u[[3]]  
  
  print('NO_HARM_Logit')  
  Logit_No_HARM=NO_HARM_AS2(y_obs_t,d1_t,d0_t,pi0,pi1,Logit_p1,Logit_p0,so_AS2)
  
  est_sim[[7]]=Logit_No_HARM[[1]]
  
  ##############################################
  ######OPTIMAL LINEAR##########################
  ##############################################  
  
  print('Optimal Linear')
  ##Optimal Linear Adjustment
  result_temp=Optimal_Reg_Est(y_obs_t,d1_t,d0_t,x,pi0,pi1)
  
  est_sim[8]= result_temp[[1]]
  optimal_b =result_temp[[2]]  
  
  ###############################################
  ######OPTIMAL LOGIT###########################
  ##############################################  
  print('Optimal Logit')
  result_temp_Logistic_MS=Optimal_Logit_Est_M(y_obs_t,d1_t,d0_t,x,pi0,pi1,10000,T,'BFGS')
  est_sim[9]= result_temp_Logistic_MS[[1]]
  optimal_b_Logit_MS =result_temp_Logistic_MS[[2]]       
  
  ##############################################
  ######OPTIMAL Imputed#########################
  ##############################################  
  
  
  print('Optimal Linear Generated Regressor')
  result_optimal_linear_genereated=NO_HARM2_AS2(y_obs_t,d1_t,d0_t,pi0,pi1,OLS_p1,OLS_p0,so_AS2)
  est_sim[10]=result_optimal_linear_genereated[[1]]
  
  print('Optimal Logit Generated Regressor')
  result_optimal_logit_genereated=NO_HARM2_AS2(y_obs_t,d1_t,d0_t,pi0,pi1,Logit_p1,Logit_p0,so_AS2)
  est_sim[11]=result_optimal_logit_genereated[[1]]


  ##########################################################################
  ####################Variance Bound Estimation#############################
  ##########################################################################
  
  print('Variance Bound Estimation')
  #variance bound estimation
  var_sim[1]= HT_var_full_AS2(y_obs_t,d1_t,d0_t,subject_size,fo,so,contrast,dict_group)
  var_sim[2]= HA_var_full_AS2(y_obs_t,d1_t,d0_t,pi_weight_t,subject_size,fo,so,contrast,dict_group)
  print('OLS')  
  var_sim[3]= OLS_var_full_AS2(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint',fo,so,contrast,dict_group)
  print('WLS')
  var_sim[4]= OLS_var_full_AS2(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint',fo,so,contrast,dict_group)
  
  #WLS no harm
  var_sim[5]=  WLS_No_HARM[[2]]
  
  print('logit')
  var_sim[6]= Logit_var_full_AS2(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='unweighted_joint',fo,so,contrast,dict_group)
  
  #logit no harm
  var_sim[7]= Logit_No_HARM[[2]]  
  
  print('Optimal_Linear')  
  var_sim[8]= Optimal_Reg_var_full_AS2(y_obs_t,d1_t,d0_t,x,optimal_b)  
  print('Optimal_Logit') 
  var_sim[9]= Optimal_Logit_var_full_AS2(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_MS)  
  
  #generate regressors
  var_sim[10]=result_optimal_linear_genereated[[2]]
  var_sim[11]=result_optimal_logit_genereated[[2]]
  
  
  print('Variance Bound Estimation 2')
  #variance bound estimation
  var_sim2[1]= HT_var_full(y_obs_t,d1_t,d0_t,subject_size,fo,so,contrast,dict_group)
  var_sim2[2]= HA_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,subject_size,fo,so,contrast,dict_group)
  print('OLS')  
  var_sim2[3]= OLS_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='OLS_joint',fo,so,contrast,dict_group)
  print('WLS')
  var_sim2[4]= OLS_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='WLS_joint',fo,so,contrast,dict_group)
  
  print('NO_HARM_Linear')  
  WLS_No_HARM_AS=NO_HARM(y_obs_t,d1_t,d0_t,pi0,pi1,WLS_p1,WLS_p0)
  var_sim2[5]=  WLS_No_HARM_AS[[2]]
  
  print('logit')
  var_sim2[6]= Logit_var_full(y_obs_t,d1_t,d0_t,pi_weight_t,x,mode='unweighted_joint',fo,so,contrast,dict_group)
  print('Optimal_Linear')  

  
  print('NO_HARM_Logit')  
  Logit_No_HARM_AS=NO_HARM(y_obs_t,d1_t,d0_t,pi0,pi1,Logit_p1,Logit_p0)
  var_sim2[7]=  Logit_No_HARM_AS[[2]]  
  
  print('Optimal Reg')
  var_sim2[8]= Optimal_Reg_var_full(y_obs_t,d1_t,d0_t,x,optimal_b) 
  
  print('Optimal_Logit') 
  var_sim2[9]= Optimal_Logit_var_full(y_obs_t,d1_t,d0_t,x,optimal_b_Logit_MS)  

  
  #generate regressors
  var_sim2[10]= result_optimal_linear_genereated[[2]]
  var_sim2[11]= result_optimal_logit_genereated[[2]]
  
  
  return(list(est_sim,var_sim,var_sim2))
}