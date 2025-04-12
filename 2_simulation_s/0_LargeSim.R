##header
source('0_header.R')

#index for parallelizing
input= commandArgs(trailingOnly=TRUE)
index=as.numeric(input[1])
expo1=as.numeric(input[2])
expo2=as.numeric(input[3])
weight=as.numeric(input[4])

#set random seed
set.seed(index)

#case_to_do=read.csv('/home/hc654/Unified/scripts_grace/simulation_network_optimal/case_to_do_14.csv')[,2]
#index=case_to_do[index]
#simulation parameters
nsim=1

time_start=Sys.time()


#initial membership informations

#drop cluster 7, 12 and 29, exposure mapping not defined
#comp_list=1:31
#comp_list=comp_list[-c(7,12,29)] #no exposure mapping for these two  
#comp_list=comp_list[-7]
#net_t = get.inducedSubgraph(net,which(memberships %in% comp_list))
prob=c(0.25,0.25,0.25,0.25)


#load foso informations

foso_file=paste0('/home/hc654/palmer_scratch/final_analysis_Ds/','Stratified_nat_',expo1,expo2,'_D.Rdata')
load(foso_file)
fo=mean_pre
so=cov_pre
# components_t=component.dist(net_t,connected='weak')
# memberships_t=components_t$membership
# csize_t = components_t$csize
# for (i in (unique(memberships_t))){
#  print(i)
#  net_tt=get.inducedSubgraph(net_t,which(memberships_t==i))
#  foso_t=FOSO(net_tt,individual_exposure,num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list,1,4)
#  fo[[i]]=foso_t[[1]]
#  so[[i]]=foso_t[[2]]
# }

#unfold fo for later calculations
fo_vec=unlist(fo)


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

##divide into components####
memberships_t=net %v% 'membership'



#extract subjects of interest
subjects_t = (1:pol_size)[id_t %in% id_list]
#extract membership information
subjects_memberships_t=memberships_t[id_t %in% id_list]
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
dict_group<- matrix(as.numeric(dict_group),   ncol = ncol(dict_group))
contrast=c(-1,1)



#method=c('HT','HA','OLS_joint','OLS_separate','WLS_joint','WLS_separate','unweighted_joint',
#         'unweighted_sepaprate','weighted_joint','weighted_separate')
method=c('HT','HA','OLS_joint','WLS_joint','No_harm_WS', 'unweighted_joint', 'No_harm_logit','Optimal_Linear','Optimal_Linear_fixed','Optimal_Logistic_MS','Optimal_Logit_mS')
est_sim=matrix(0,nrow=nsim,ncol=length(method)) 
var_sim=matrix(0,nrow=nsim,ncol=length(method)) 
sample_size = matrix(0,nrow=2,ncol = nsim)
optimization_ind=rep(nsim,0)
colnames(est_sim)=method
colnames(var_sim)=method

#prepare y and x informations, in the order of group beloings
#Y=Y_sim2[match(id_t,Y_sim2[,1]),c(1,expo1+1,expo2+1)]
#Y_all=Y_all[match(id_t,Y_all[,1]),]
#Y_all=Y_sim2[match(id_t,Y_sim2[,1]),]
#Y_all=Y_all[subjects_t,]
#Y=Y[subjects_t,]
#Y=Y_all[,c(1,2,3)]
weights_random = runif(subject_size)
y0=y0_opt * (weights_random <= weight) + y0_imp * (weights_random > weight) 
y1=y1_opt * (weights_random <= weight) + y1_imp * (weights_random > weight) 
Y_all=cbind(y0,y1,y1,y1,y1,y1,y1,y1,y1,y1,y1,y1)

#covariates
x=all_info2[match(id_t,all_info2$id),]
x=x[subjects_t,]
subject_size=length(subjects_t)
#demean X
#pre-treatment covarates and demeaning
x=x[,c('male','age','agpop','literacy','ricearea_2010','risk_averse','disaster_prob')] 
for (i in 1:ncol(x)){
  x[,i]=x[,i]-mean(x[,i])
}
x=as.matrix(x)


#friendship information
friends=matrix(FALSE,nrow=pol_size,ncol=5)
for (i in 1:pol_size){
  friend_index_t=get.neighborhood(net,i,'out')
  if (length(friend_index_t)==0){
    next 
  }else
    friends[i,1:length(friend_index_t)] = friend_index_t
}


for (i in 1:nsim){
  
  
  result=SIM_ONERUN_IMPUTED_AUG20(status,pol_size,prob,y1,y0,pi1,pi0,Y_all,index,expo1,expo2,option='Stratified')
  
  est_sim=result[[1]]
  var_sim=result[[2]]
    
}


time_end=Sys.time()

print(time_start-time_end)
#output=list(est_sim,var_sim, coef_sim_OLS,coef_sim_WLS,coef_sim_Logit,coef_sim_WLogit,coef_sim_OptimalOLS,coef_sim_OptimalLogit,coef_sim_OptimalLogitOS,coef_sim_OptimalLogitMS,coef_sim_OptimalLogitMSO)
output=list(est_sim,var_sim)



file_name=paste0('/home/hc654/palmer_scratch/final_analysis_simulation' ,expo1,expo2,'/',weight,'_',expo1,expo2,'_Sim_',index,'.Rdata')
#file_name2=paste0('/home/hc654/palmer_scratch/final_analysis_simulation/',weight,'_',expo1,expo2,'_SimSampleSize_',index,'.Rdata')
save(output,file=file_name)
#save(sample_size,file=file_name2)

#csv_name_est=paste0('/home/hc654/palmer_scratch/final_analysis_simulation/',expo1,expo2,'_Sim_est.csv')
#csv_name_var=paste0('/home/hc654/palmer_scratch/final_analysis_simulation/',expo1,expo2,'_Sim_var.csv')

#print(csv_name_est)
#print(csv_name_var)
#row1=matrix(c(index,est_sim),1,16)
#row2=matrix(c(index,var_sim),1,16)
#print(row1)
#print(row2)
#write.table(row1, file = csv_name_est, sep = ",", append = TRUE, quote = FALSE,
#            col.names = FALSE, row.names = FALSE)
#write.table(row2, file = csv_name_var, sep = ",", append = TRUE, quote = FALSE,
#            col.names = FALSE, row.names = FALSE)

