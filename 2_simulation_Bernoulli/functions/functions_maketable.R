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
  print('HT')
  z_HT = y * contrast_vec
  
  ####HA#####    
  print('HA')
  x_HA = data.frame(cbind(d0,d1))
  reg_HA =  lm(y~d1+d0-1)
  z_HA= y-predict(reg_HA,x_HA)
  z_HA = contrast_vec * z_HA
  
  print('OLS')
  ####OLS_joint####
  x_OLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_OLS_joint =  lm(y~d0+d1+x_predict-1,weights=fo_vec)
  OLS_pred=predict(reg_OLS_joint,x_OLS_joint)
  temp_OLS =  (y-predict(reg_OLS_joint,x_OLS_joint)) * fo_vec
  z_OLS_joint=  temp_OLS * as.matrix(x_OLS_joint) %*% solve( t(as.matrix(x_OLS_joint)) %*%diag(fo_vec)%*%as.matrix(x_OLS_joint)*(1/subject_size))
  z_OLS_joint = z_OLS_joint %*% contrast_OLS
  
  print('WLS')
  ####WLS_joint####
  x_WLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_WLS_joint =  lm(y~d0+d1+x_predict-1)
  temp_WLS = y-predict(reg_WLS_joint,x_WLS_joint)
  z_WLS_joint=  temp_WLS* as.matrix(x_WLS_joint) %*% solve( t(as.matrix(x_WLS_joint)) %*%as.matrix(x_WLS_joint)* (1/subject_size))
  z_WLS_joint = z_WLS_joint %*% contrast_OLS
  
  print('Noharm Linear')
  ###No_harm_WLS####
  
  WLS_pred=predict(reg_WLS_joint,x_WLS_joint)
  z_Noharm_WLS=Noharm_Reg(y1,y0,  WLS_pred)
  
  print('Logit')
  ###Logit_joint_Unweighted###
  x_Logit_joint_u = data.frame(cbind(d0,d1,x_predict))
  reg_Logit_joint_u =  glm(y~d0+d1+x_predict-1,family='quasibinomial',weights=fo_vec)
  z_Logit_joint_u = y-predict(reg_Logit_joint_u,x_Logit_joint_u,type='response')
  z_Logit_joint_u = z_Logit_joint_u * contrast_vec
  
  print('Noharm Logit')
  ##No_harm_Wlogit
  ULogit_pred=  predict(reg_Logit_joint_u,x_Logit_joint_u,type='response')
  z_Noharm_Logit=Noharm_Reg(y1,y0,ULogit_pred)  
  
  
  print('Optimal Linear')
  ####Optimal Coef####
  z_optimal = Optimal_Reg(y1,y0,x)
  
  
  print('Optimal Logit')
  ####Optimal Coef
  z_optimal_Logit_M = Optimal_Logit_M(y1,y0,x,Maxiter=1000,Optim=T)
  
  
  ####Optimal Coef
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
    
    #linear and no harm
    z_OLS_joint_t = z_OLS_joint[index_start:index_end]
    z_WLS_joint_t = z_WLS_joint[index_start:index_end]
    z_Noharm_WLS_t=z_Noharm_WLS[index_start:index_end]

    #logit and no harm
    z_Logit_joint_u_t =z_Logit_joint_u[index_start:index_end]
    z_Noharm_Logit_t=z_Noharm_Logit[index_start:index_end]  
    
    #optimal linear and logit
    z_optimal_t = z_optimal[index_start:index_end]
    z_optimal_Logit_M_t = z_optimal_Logit_M[index_start:index_end]
    
    #optimal imputed
    z_Noharm_WLS_Optimal_t=z_Noharm_WLS_Optimal[index_start:index_end]      
    z_Noharm_Logit_Optimal_t = z_Noharm_Logit_Optimal[index_start:index_end]   
    
    output[1]=output[1] +  t(z_HT_t) %*% normalized_cov %*% z_HT_t/(subject_size^2)
    output[2]=output[2] +  t(z_HA_t) %*% normalized_cov %*% z_HA_t/(subject_size^2)
    output[3]=output[3] +  t(z_OLS_joint_t) %*% normalized_cov %*% z_OLS_joint_t/(subject_size^2)
    output[4]=output[4] +  t(z_WLS_joint_t) %*% normalized_cov %*% z_WLS_joint_t/(subject_size^2)
    output[5]=output[5] +  t(z_Noharm_WLS_t) %*% normalized_cov %*%  z_Noharm_WLS_t/(subject_size^2)
    output[6]=output[6] +  t(z_Logit_joint_u_t) %*% normalized_cov %*% z_Logit_joint_u_t/(subject_size^2)
    output[7]=output[7] +  t(z_Noharm_Logit_t) %*% normalized_cov %*%  z_Noharm_Logit_t/(subject_size^2)    
    
    output[8]=output[8] +  t(z_optimal_t) %*% normalized_cov %*%  z_optimal_t/(subject_size^2)
    output[9]=output[9] +  t(z_optimal_Logit_M_t) %*% normalized_cov %*%  z_optimal_Logit_M_t/(subject_size^2)
    
    output[10]=output[10] +  t(z_Noharm_WLS_Optimal_t) %*% normalized_cov %*%  z_Noharm_WLS_Optimal_t/(subject_size^2)    
    output[11]=output[11] +  t(z_Noharm_Logit_Optimal_t) %*% normalized_cov %*%  z_Noharm_Logit_Optimal_t/(subject_size^2)
    
  }
  
  return(output)
}
theoretical_variance_bound_full2_AS2=function(y1,y0,x,fo,so,so_AS2,contrast,dict_group){
  
  
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
  print('HT')
  z_HT = y * contrast_vec
  
  ####HA#####    
  print('HA')
  x_HA = data.frame(cbind(d0,d1))
  reg_HA =  lm(y~d1+d0-1)
  z_HA= y-predict(reg_HA,x_HA)
  z_HA = contrast_vec * z_HA
  
  print('OLS')
  ####OLS_joint####
  x_OLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_OLS_joint =  lm(y~d0+d1+x_predict-1,weights=fo_vec)
  OLS_pred=predict(reg_OLS_joint,x_OLS_joint)
  temp_OLS =  (y-predict(reg_OLS_joint,x_OLS_joint)) * fo_vec
  z_OLS_joint=  temp_OLS * as.matrix(x_OLS_joint) %*% solve( t(as.matrix(x_OLS_joint)) %*%diag(fo_vec)%*%as.matrix(x_OLS_joint)*(1/subject_size))
  z_OLS_joint = z_OLS_joint %*% contrast_OLS
  
  print('WLS')
  ####WLS_joint####
  x_WLS_joint = data.frame(cbind(d0,d1,x_predict))
  reg_WLS_joint =  lm(y~d0+d1+x_predict-1)
  temp_WLS = y-predict(reg_WLS_joint,x_WLS_joint)
  z_WLS_joint=  temp_WLS* as.matrix(x_WLS_joint) %*% solve( t(as.matrix(x_WLS_joint)) %*%as.matrix(x_WLS_joint)* (1/subject_size))
  z_WLS_joint = z_WLS_joint %*% contrast_OLS
  
  print('Noharm Linear')
  ###No_harm_WLS####
  
  WLS_pred=predict(reg_WLS_joint,x_WLS_joint)
  z_Noharm_WLS=Noharm_Reg(y1,y0,  WLS_pred)
  
  print('Logit')
  ###Logit_joint_Unweighted###
  x_Logit_joint_u = data.frame(cbind(d0,d1,x_predict))
  reg_Logit_joint_u =  glm(y~d0+d1+x_predict-1,family='quasibinomial',weights=fo_vec)
  z_Logit_joint_u = y-predict(reg_Logit_joint_u,x_Logit_joint_u,type='response')
  z_Logit_joint_u = z_Logit_joint_u * contrast_vec
  
  print('Noharm Logit')
  ##No_harm_Wlogit
  ULogit_pred=  predict(reg_Logit_joint_u,x_Logit_joint_u,type='response')
  z_Noharm_Logit=Noharm_Reg(y1,y0,ULogit_pred)  
  
  
  print('Optimal Linear')
  ####Optimal Coef####
  z_optimal = Optimal_Reg(y1,y0,x)
  
  
  print('Optimal Logit')
  ####Optimal Coef
  z_optimal_Logit_M = Optimal_Logit_M(y1,y0,x,Maxiter=1000,Optim=T)
  
  
  ####Optimal Coef
  z_Noharm_WLS_Optimal=Noharm_Optimal_Reg(y1,y0,  OLS_pred)
  z_Noharm_Logit_Optimal=Noharm_Optimal_Reg(y1,y0,  ULogit_pred)
  
  output=rep(0,11)
  
  names(output)=c('var_HT','var_HA','var_OLS_joint','var_WLS_joint','var_Logit_joint_u','var_Logit_joint_w','var_OLS_DR','var_Linear_Optimal','var_Logit_Optimal_MS','no_harm_wls','no_harm_logit','var_Logit_Optimal_OS','no_harm_optimal_WLS')
  for (i in 1:nrow(dict_group)){
    
    c=dict_group[i,1]
    index_start=dict_group[i,4]
    index_end=dict_group[i,5]    
    
    cov_bound_t = so_AS2[[c]]
    #normalized_cov=so_t
    z_HT_t = z_HT[index_start:index_end]
    z_HA_t = z_HA[index_start:index_end]
    
    #linear and no harm
    z_OLS_joint_t = z_OLS_joint[index_start:index_end]
    z_WLS_joint_t = z_WLS_joint[index_start:index_end]
    z_Noharm_WLS_t=z_Noharm_WLS[index_start:index_end]
    
    #logit and no harm
    z_Logit_joint_u_t =z_Logit_joint_u[index_start:index_end]
    z_Noharm_Logit_t=z_Noharm_Logit[index_start:index_end]  
    
    #optimal linear and logit
    z_optimal_t = z_optimal[index_start:index_end]
    z_optimal_Logit_M_t = z_optimal_Logit_M[index_start:index_end]
    
    #optimal imputed
    z_Noharm_WLS_Optimal_t=z_Noharm_WLS_Optimal[index_start:index_end]      
    z_Noharm_Logit_Optimal_t = z_Noharm_Logit_Optimal[index_start:index_end]   
    
    output[1]=output[1] +  t(z_HT_t) %*% cov_bound_t %*% z_HT_t/(subject_size^2)
    output[2]=output[2] +  t(z_HA_t) %*% cov_bound_t %*% z_HA_t/(subject_size^2)
    output[3]=output[3] +  t(z_OLS_joint_t) %*% cov_bound_t %*% z_OLS_joint_t/(subject_size^2)
    output[4]=output[4] +  t(z_WLS_joint_t) %*% cov_bound_t %*% z_WLS_joint_t/(subject_size^2)
    output[5]=output[5] +  t(z_Noharm_WLS_t) %*% cov_bound_t %*%  z_Noharm_WLS_t/(subject_size^2)
    output[6]=output[6] +  t(z_Logit_joint_u_t) %*% cov_bound_t %*% z_Logit_joint_u_t/(subject_size^2)
    output[7]=output[7] +  t(z_Noharm_Logit_t) %*% cov_bound_t %*%  z_Noharm_Logit_t/(subject_size^2)    
    
    output[8]=output[8] +  t(z_optimal_t) %*% cov_bound_t %*%  z_optimal_t/(subject_size^2)
    output[9]=output[9] +  t(z_optimal_Logit_M_t) %*% cov_bound_t %*%  z_optimal_Logit_M_t/(subject_size^2)
    
    output[10]=output[10] +  t(z_Noharm_WLS_Optimal_t) %*% cov_bound_t %*%  z_Noharm_WLS_Optimal_t/(subject_size^2)    
    output[11]=output[11] +  t(z_Noharm_Logit_Optimal_t) %*% cov_bound_t %*%  z_Noharm_Logit_Optimal_t/(subject_size^2)
    
  }
  
  return(output)
}


HT_variance_full=function(y,sample_size,cov_bound,contrast){
  
  #make contrast estimator
  contrast_vec=rep(contrast,length(y)/2)
  
  #computing variance bound, accounting for contrasts
  var_est = 1/(sample_size)^2 * (contrast_vec*y) %*% cov_bound %*% (y*contrast_vec)
  return(var_est)
  
}

extract_result_imputed = function(expo1,expo2,y1,y0,fo,so,so_AS2,dict_group,subject_size){
  
  path=paste0('/home/hc654/palmer_scratch/final_analysis_Ds')
  setwd(path)
  
  print('Inside Function')
  setwd('/home/hc654/palmer_scratch/final_analysis_Ds/simulation_output')
  files=list.files()
  
  Sim_result_pattern=paste0(expo1,'_',expo2,'_Sim*')
  
  files_list1=files[grep(Sim_result_pattern,files)]
  #files_list2=files[grep(Sim_size_pattern,files)]
  
  print('loading data')
  
  est_sim_output=c()
  var_sim_output=c()
  var_sim2_output=c()
  sample_size_output=c()
  
  negative_var=0
  print('processing')
  start_time=proc.time()  
  
  for (i in 1:length(files_list1)){
    
    
    #load files
    print(i)
    file_name=files_list1[i]
    load(file_name)
    
    #store simulation information
    temp = output[[1]]
    temp2 = output[[2]]
    temp5=output[[3]]
    
    #tracking negative variance problem
    temp3=temp2<0
    temp4=any(temp3)
    negative_var=negative_var + sum(temp4)
      
      
    est_sim_output=rbind(est_sim_output,temp)
    var_sim_output=rbind(var_sim_output,temp2)
    var_sim2_output=rbind(var_sim2_output,temp5)
      
  }
  
  end_time=proc.time()  
  
  #load('result_Feb16.Rdata')
  #est_sim_output=est_sim_output[which(est_sim_output[,8]>-0.2),]
  #save(est_sim,var_sim,file='result_Feb16.Rdata')
  est_sim_output = est_sim_output[,1:11] - (mean(y1)-mean(y0))
  var_sim_output = var_sim_output[,1:11]
  var_sim2_output = var_sim2_output[,1:11]
  line1=apply(est_sim_output,2,mean)
  line2=apply(est_sim_output,2,var) 
  line3=sqrt(line1^2 + line2)
  line4=line2*subject_size
  line5=apply(var_sim_output,2,mean)  
  line5=subject_size*line5
  
  #theoretical variance and variance bound
  line6=theoretical_variance_full2(y1,y0,x,fo,so,contrast=c(-1,1),dict_group)
  line7=theoretical_variance_bound_full2_AS2(y1,y0,x,fo,so,so_AS2,contrast=c(-1,1),dict_group)
  
  #variance bound estimation quality
  line8= (apply(var_sim_output,2,mean)-line7)/line7 #bias 
  
  #coverage of the normal-based CI
  coverage = abs(est_sim_output/(sqrt(var_sim_output))) >= 1.96 
  
  #  bootstrap_coverage= ((est_sim_output/sqrt(var_sim_output)) >= quantile0025) * ((est_sim_output/sqrt(var_sim_output)) <= quantile0975)
  line9=1-apply(coverage,2,mean)
  
  #output the results
  result_EST=rbind(line1,line2,line3,line4,line5,line9,line6_temp, line7_temp,line8)  
  colnames(result_EST)=c('HT','HA','OLS','WLS','NO_HARM_WLS','LOGIT','NO_HARM_LOGIT','OPTIMAL_LINEAR','No_harm_Optimal_linear','No_harm_Optimal_logit')
  
  rownames(result_EST)=c('bias','true_var','MSE','normalized_var','normalized_var_b_est','normal coverage','linearized_var','variance bound','estimation error')  
  
  
  save_path='/home/hc654/palmer_scratch/final_results/'
  file_name = paste0(save_path,weight,'_Sim_',expo1,expo2,'_',date,'_Ds.Rdata')
  save(result_EST,negative_var,file=file_name)
  
}

load_data_imputed=function(expo1,expo2){
  
  nsim=1
  expo1=expo1
  expo2=expo2
  
  
  #initial membership informations
  

  prob=c(0.25,0.25,0.25,0.25)
  
  
  foso_file=paste0('/home/hc654/palmer_scratch/final_analysis_Bernoulli/','Bernoulli_',expo1,expo2,'_D.Rdata')
  load(foso_file)
  fo<<-mean_pre
  so<<-cov_pre
  
  so_bound_file=paste0('/home/hc654/palmer_scratch/final_analysis_Bernoulli/','Bernoulli_',expo1,expo2,'_Dbound.Rdata')
  load(so_bound_file)
  so_AS2<<-so_bound
  fo_vec<<-unlist(fo)
  
  ###basic simulation information
  #network size
  pol_size = network.size(net)
  
  #treatment status vectors/num of mappings
  status = 1:4
  num_mappings = 12
  
  #id information
  id_t = net %v% 'vertex.names'
  degree_t=degree(net,cmode='outdegree') #warning some units have only two friends so some exposure mappings are not defined
  
  
  #depending on exposure mapping some units must be dropped for identification reasons:
  if (expo1 == 7 | expo1 == 12 | expo2 ==7 | expo2==12){
    
    number_dropped_t =sum( (id_t[degree_t<3] %in% id_list))
    #delete units with less than three friends
    id_list=id_list[!id_list %in% id_t[degree_t<3]] 
  }else if (expo1 == 6 | expo1 == 11 | expo2 ==6 | expo2==11){
    
    number_dropped_t =sum( (id_t[degree_t<2] %in% id_list))
    #delete units with less than two friends
    id_list=id_list[!id_list %in% id_t[degree_t<2]]
    
  }else{
    
    number_dropped_t =sum( (id_t[degree_t<1] %in% id_list))
    #delete units with no friends
    id_list=id_list[!id_list %in% id_t[degree_t<1]] 
    
  }
  
  ###read from some input files
  pi0 = fo_vec[(1:length(fo_vec))%%2==1]
  pi1 = fo_vec[(1:length(fo_vec))%%2==0]
  
  ######################################################
  #EXTRACT SUBJECTS OF INTEREST#########################
  ######################################################
  ##divide into components####
  memberships_t=net %v% 'membership'
  
  #extract subjects of interest
  subjects_t = (1:pol_size)[id_t %in% id_list]
  #extract membership information
  subjects_memberships_t=memberships_t[id_t %in% id_list]
  #sort subjects according to component index
  subjects_t=subjects_t[order(subjects_memberships_t)]
  
  
  ######################################################
  #CREATE DICTIONARY FOR GROUP BELONGINGS###############
  ######################################################
  
  group_index = cumsum(table(subjects_memberships_t))
  dict_group=matrix(0,length(group_index),5)
  dict_group[,1]=names(group_index)
  dict_group[,2]= c(0,group_index[(1:(length(group_index))-1)])+1
  dict_group[,3]= group_index
  dict_group[,4]=c(0,2*group_index[(1:(length(group_index))-1)])+1
  dict_group[,5]= group_index*2
  dict_group<- matrix(as.numeric(dict_group),   ncol = ncol(dict_group))
  
  ###########################################################################################
  #SANITY CHECK: MAKE SURE THE ORDERING FOR FO AND SO MATRICES ARE CONSISTENT################
  ###########################################################################################
  
  print('Checking inconsisntey in id orderings of foso probabilities (no error msg is good)')
  
  
  compare1=compute_FOSO_all_components_sanity_check(net,expo1,expo2 ,id_list,option='Bernoulli')
  for (i in 1:nrow(dict_group)){
    
    #network index
    m = dict_group[i,1]
    
    #for those subjects who are in group m, what are their id?
    compare2_temp = id_t[subjects_t[dict_group[i,2]:dict_group[i,3]]]
    
    compare1_temp = compare1[[m]]  
    
    temp = (compare1_temp)==(compare2_temp)
    #if there is any inequality !temp will contain a true and this prompts an error msg
    if (any(!temp)==TRUE ){
      print(paste0('Inconsistency in network component m:', m))
    }
    
  }
  
  #################################################################
  ####SORT OUTCOMES################################################
  #################################################################
  #the outcomes are sorted according to the orderings of subjects_t
  
  Y=matrix(0,0,3)
  for (i in 1:nrow(dict_group)){
    
    #network index
    m = dict_group[i,1]
    
    #for those subjects who are in group m, what are their id?
    id_temp = id_t[subjects_t[dict_group[i,2]:dict_group[i,3]]]
    
    
    Y_temp = Y_impute[match(id_temp,Y_impute[,1]),c(1,expo1+1,expo2+1)]
    
    Y = rbind(Y,Y_temp)
    
  }
  
  
  #Y=Y[subjects_t,]
  #Y=Y_all[,c(1,2,3)]
  #weights_random = runif(subject_size)
  y0=Y[,2] #control
  y1=Y[,3] #treated
  #y0=y0_opt * (weights_random <= weight) + y0_imp * (weights_random > weight) 
  #y1=y1_opt * (weights_random <= weight) + y1_imp * (weights_random > weight) 
  
  ############################################################################
  ####SANITY CHECK AGAIN: MAKE THE OUTCOMES AND FOSO ARE ALIGNED##############
  ############################################################################
  
  print('Checking nnconsisntey in id orderings of outcomes (no error msg is good)')
  
  for (i in 1:nrow(dict_group)){
    
    
    #network index
    m = dict_group[i,1]
    
    #for those subjects who are in group m, what are their id?
    compare3_temp =Y[(dict_group[i,2]:dict_group[i,3]),1]
    
    
    compare1_temp = compare1[[m]]  
    temp = (compare1_temp)==(compare3_temp)
    
    #if there is any inequality !temp will contain a true and this prompts an error msg
    if (any(!temp)==TRUE ){
      
      print(paste0('Inconsistency in network component m:', m))
    }
    
  }
  
  
  
  #################################################################
  #######COVARIATES################################################
  #################################################################
  
  X=all_info2[match(id_t,all_info2$id),]
  X=X[,c('id','male','age','agpop','literacy','ricearea_2010','risk_averse','disaster_prob')] 
  x=c()
  for (i in 1:nrow(dict_group)){
    
    #network index
    m = dict_group[i,1]
    
    #for those subjects who are in group m, what are their id?
    id_temp = id_t[subjects_t[dict_group[i,2]:dict_group[i,3]]]
    
    
    X_temp = X[match(id_temp,X[,1]),]
    
    x = rbind(x,X_temp)
    
  }
  
  
  for (i in 2:ncol(x)){
    x[,i]=(x[,i]-mean(x[,i]))/sd(x[,i])
  }
  x=as.matrix(x)
  
  ############################################################################
  ####SANITY CHECK AGAIN: MAKE THE OUTCOMES AND COVARATES ALIGNED##############
  ############################################################################
  
  print('Checking nnconsisntey in id orderings of x (no error msg is good)')
  temp = (Y[,1]==x[,1])
  
  if (any(!temp)==TRUE ){
    
    print('Inconsistency in between outcomes and covariates')
  }
  
  x=x[,2:ncol(x)]
  
  #################################################################
  #######FRIENDSHIP################################################
  #################################################################
  
  #friendship information
  friends=matrix(FALSE,nrow=pol_size,ncol=5)
  for (i in 1:pol_size){
    friend_index_t=get.neighborhood(net,i,'out')
    if (length(friend_index_t)==0){
      next 
    }else
      friends[i,1:length(friend_index_t)] = friend_index_t
  }
  
  return(list(y1=y1,y0=y0,x=x,fo=fo,so=so,so_AS2=so_AS2,dict_group=dict_group))  
}
