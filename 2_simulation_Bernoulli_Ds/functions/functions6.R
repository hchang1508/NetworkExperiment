extract_result_imputed_Oct14 = function(expo1,expo2,date,weight){
  
  path=paste0('/home/hc654/palmer_scratch/final_analysis_Ds_',expo1,expo2,'_AS2/')
  #path=paste0('/home/hc654/palmer_scratch/final_analysis_simulation',expo1,expo2,'/')
  
  setwd(path)
  print('Inside Function')
  #setwd('/home/hc654/palmer_scratch/Unified_Complete_Simulation_Imputed/')
  files=list.files()
  
  Sim_result_pattern=paste0('^',weight,'_',expo1,expo2,'_Sim*')
  
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
    
    
    print(i)
    file_name=files_list1[i]
    #file_name2=files_list2[i]
    load(file_name)
    #load(file_name2)
    
    
    temp = output[[1]]
    if (length(temp)==11){
      temp2 = output[[2]]
      temp5=output[[3]]
      temp3=temp2<0
      
      #temp4=any(temp3)
      #negative_var=negative_var + sum(temp4)
      

        est_sim_output=rbind(est_sim_output,temp)
        var_sim_output=rbind(var_sim_output,temp2)
        var_sim2_output=rbind(var_sim2_output,temp5)
      
    }
  }  
  end_time=proc.time()  
  
  #load('result_Feb16.Rdata')
  #est_sim_output=est_sim_output[which(est_sim_output[,8]>-0.2),]
  #save(est_sim,var_sim,file='result_Feb16.Rdata')
  est_sim_output = est_sim_output[,1:9] - (mean(y1)-mean(y0))
  var_sim_output = var_sim_output[,1:9]
  var_sim2_output = var_sim2_output[,1:9]
  line1=apply(est_sim_output,2,mean)
  line2=apply(est_sim_output,2,var) 
  line3=sqrt(line1^2 + line2)
  line4=line2*subject_size
  line5=apply(var_sim_output,2,mean)  
  line5=subject_size*line5
  
  #need to load data
  
  line6=theoretical_variance_full(y1,y0,x,fo,so,contrast=c(-1,1),dict_group)
  line7=theoretical_variance_bound_full_AS2(y1,y0,x,fo,so,contrast=c(-1,1),dict_group)
  
  
  #variance_estimate
  line6_temp = line6[c(1:4,11,6,12,8,10)]
  line7_temp = line7[c(1:4,11,6,12,8,10)]
  line8= (apply(var_sim_output,2,mean)-line7_temp)/line7_temp
  line6_temp=line6_temp*subject_size
  line7_temp=line7_temp*subject_size
  
  
  coverage = abs(est_sim_output/(sqrt(var_sim_output))) >= 1.96 
  
  #  bootstrap_coverage= ((est_sim_output/sqrt(var_sim_output)) >= quantile0025) * ((est_sim_output/sqrt(var_sim_output)) <= quantile0975)
  line9=1-apply(coverage,2,mean)
  
  coverage2 = abs(est_sim_output/(sqrt(var_sim2_output))) >= 1.96
  line9_2=1-apply(coverage,2,mean)
  
  #  line10=apply(bootstrap_coverage,2,mean)
  #  line11= apply(sqrt(var_sim_output),2,mean) * qnorm(0.975)*2
  #  line12 = apply(sqrt(var_sim_output)*quantile0975,2,mean)  -apply(sqrt(var_sim_output)*quantile0025,2,mean)
  
  result_EST=rbind(line1,line2,line3,line4,line5,line9,line6_temp, line7_temp,line8)  
  colnames(result_EST)=c('HT','HA','OLS','WLS','NO_HARM_WLS','LOGIT','NO_HARM_LOGIT','OPTIMAL_LINEAR','OPTIMAL_LOGIT')
  
  rownames(result_EST)=c('bias','true_var','MSE','normalized_var','normalized_var_b_est','normal coverage','linearized_var','variance bound','estimation error')  

  
  save_path='/home/hc654/palmer_scratch/final_results/'
  file_name = paste0(save_path,weight,'_Sim_',expo1,expo2,'_',date,'_Ds.Rdata')
  save(result_EST,negative_var,file=file_name)
}
extract_result_imputed_Nov11 = function(expo1,expo2,date,weight){
  
  path=paste0('/home/hc654/palmer_scratch/final_analysis_Ds_',expo1,expo2,'_AS2/')
  #path=paste0('/home/hc654/palmer_scratch/final_analysis_simulation',expo1,expo2,'/')
  
  setwd(path)
  print('Inside Function')
  #setwd('/home/hc654/palmer_scratch/Unified_Complete_Simulation_Imputed/')
  files=list.files()
  
  Sim_result_pattern=paste0('^',weight,'_',expo1,expo2,'_Sim*')
  
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
    
    
    print(i)
    file_name=files_list1[i]
    #file_name2=files_list2[i]
    load(file_name)
    #load(file_name2)
    
    
    temp = output[[1]]
    if (length(temp)==11){
      temp2 = output[[2]]
      temp5=output[[3]]
      temp3=temp2<0
      
      #temp4=any(temp3)
      #negative_var=negative_var + sum(temp4)
      
      
      est_sim_output=rbind(est_sim_output,temp)
      var_sim_output=rbind(var_sim_output,temp2)
      var_sim2_output=rbind(var_sim2_output,temp5)
      
    }
  }  
  end_time=proc.time()  
  
  #load('result_Feb16.Rdata')
  #est_sim_output=est_sim_output[which(est_sim_output[,8]>-0.2),]
  #save(est_sim,var_sim,file='result_Feb16.Rdata')
  est_sim_output = est_sim_output[,1:10] - (mean(y1)-mean(y0))
  var_sim_output = var_sim_output[,1:10]
  var_sim2_output = var_sim2_output[,1:10]
  line1=apply(est_sim_output,2,mean)
  line2=apply(est_sim_output,2,var) 
  line3=sqrt(line1^2 + line2)
  line4=line2*subject_size
  line5=apply(var_sim_output,2,mean)  
  line5=subject_size*line5
  
  #need to load data
  
  line6=theoretical_variance_full2(y1,y0,x,fo,so,contrast=c(-1,1),dict_group)
  line7=theoretical_variance_bound_full2_AS2(y1,y0,x,fo,so,contrast=c(-1,1),dict_group)
  
  
 
  line6_temp = line6[c(1:4,9,5,10,8,11,12)]
  line7_temp = line7[c(1:4,9,5,10,8,11,12)]
  line8= (apply(var_sim_output,2,mean)-line7_temp)/line7_temp
  line6_temp=line6_temp*subject_size
  line7_temp=line7_temp*subject_size
  
  
  coverage = abs(est_sim_output/(sqrt(var_sim_output))) >= 1.96 
  
  #  bootstrap_coverage= ((est_sim_output/sqrt(var_sim_output)) >= quantile0025) * ((est_sim_output/sqrt(var_sim_output)) <= quantile0975)
  line9=1-apply(coverage,2,mean)
  
  coverage2 = abs(est_sim_output/(sqrt(var_sim2_output))) >= 1.96
  line9_2=1-apply(coverage,2,mean)
  
  #  line10=apply(bootstrap_coverage,2,mean)
  #  line11= apply(sqrt(var_sim_output),2,mean) * qnorm(0.975)*2
  #  line12 = apply(sqrt(var_sim_output)*quantile0975,2,mean)  -apply(sqrt(var_sim_output)*quantile0025,2,mean)
  
  result_EST=rbind(line1,line2,line3,line4,line5,line9,line6_temp, line7_temp,line8)  
  colnames(result_EST)=c('HT','HA','OLS','WLS','NO_HARM_WLS','LOGIT','NO_HARM_LOGIT','OPTIMAL_LINEAR','No_harm_Optimal_linear','No_harm_Optimal_logit')
  
  rownames(result_EST)=c('bias','true_var','MSE','normalized_var','normalized_var_b_est','normal coverage','linearized_var','variance bound','estimation error')  
  
  
  save_path='/home/hc654/palmer_scratch/final_results/'
  file_name = paste0(save_path,weight,'_Sim_',expo1,expo2,'_',date,'_Ds.Rdata')
  save(result_EST,negative_var,file=file_name)
}

load_data_imputed_Oct14=function(expo1,expo2,weight){
  
  nsim=1
  expo1=expo1
  expo2=expo2
  
  
  #initial membership informations
  

  prob=c(0.25,0.25,0.25,0.25)
  
  
  foso_file=paste0('/home/hc654/palmer_scratch/final_analysis_Ds/','Stratified_nat_',expo1,expo2,'_D.Rdata')
  load(foso_file)
  fo<<-mean_pre
  so<<-cov_pre
  
  so_bound_file=paste0('/home/hc654/palmer_scratch/final_analysis_Ds/','Stratified_nat_',expo1,expo2,'_Dbound.Rdata')
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
  
  
  memberships = net %v% 'membership'
  #extract subjects of interest
  subjects_t = (1:pol_size)[id_t %in% id_list]
  #extract membership information
  subjects_memberships_t=memberships[id_t %in% id_list]
  #sort subjects according to component index
  subjects_t=subjects_t[order(subjects_memberships_t)]
  
  #create dictionary for group belongings
  group_index = cumsum(table(subjects_memberships_t))
  dict_group=matrix(0,length(group_index),5)
  dict_group[,1]=names(group_index)
  dict_group[,2]= c(0,group_index[1:(length(group_index)-1)])+1
  dict_group[,3]= group_index
  dict_group[,4]=c(0,2*group_index[1:(length(group_index)-1)])+1
  dict_group[,5]= group_index*2
  dict_group<<- matrix(as.numeric(dict_group),   ncol = ncol(dict_group))
  contrast=c(-1,1)
  
  
  
  #prepare y and x informations, in the order of group beloings
  weights_random = runif(subject_size)
  y0_imp=Y_sim2[,expo1+1]
  y1_imp=Y_sim2[,expo2+1]
  y0<<-y0_opt * (weights_random <= weight) + y0_imp * (weights_random > weight) 
  y1<<-y1_opt * (weights_random <= weight) + y1_imp * (weights_random > weight) 
  Y_all=cbind(y0,y1,y1,y1,y1,y1,y1,y1,y1,y1,y1,y1)
  
  
  #covariates
  x=all_info2[match(id_t,all_info2$id),]
  x=x[subjects_t,]
  subject_size<<-length(subjects_t)
  #demean X
  #pre-treatment covarates and demeaning
  x=x[,c('male','age','agpop','literacy','ricearea_2010','risk_averse','disaster_prob')] 
  for (i in 1:ncol(x)){
    x[,i]=x[,i]-mean(x[,i])
  }
  
  x<<-as.matrix(x)
  
}
