rm(list=ls())
source('0_header.R')
input= commandArgs(trailingOnly=TRUE)

expo1=as.numeric(input[1])
expo2=as.numeric(input[2])
prob_cases= as.integer(input[4])
case_to_do=read.csv(paste0('/home/hc654/NetworkExperiment/2_simulation_Bernoulli_Ds/case_to_do_b', prob_cases, '_',expo1,expo2,'.csv'))

print(input)
case_to_do=case_to_do[,2]

simulation_index=as.numeric(case_to_do[as.numeric(input[3])])

option='Bernoulli'

####################################################
#####Different assignment probabilities#############
####################################################


if (prob_cases==1){
  prob=rep(0.25,4)
}else if ( prob_cases ==2){
  prob=c(1/5,1/5,3/10,3/10)
}else if ( prob_cases ==3){
  prob=c(1/6,1/6,4/12,4/12)
}else if (prob_cases==4){
  prob=c(1/7,1/7,5/14,5/14)
}else if (prob_cases==5){
  prob=c(1/8,1/8,6/16,6/16)
}else if (prob_cases==6){
  prob=c(1/9,2/9,4/12,4/12)
}
save_path=paste0('/home/hc654/palmer_scratch/Bernoulli_compareDs/' ,prob_cases,'_Bernoulli_',expo1,'_',expo2,'_FOSO_',simulation_index,'.Rdata')

###########################################################################
#######Expo 1:  In First Round Simple (FRS)################################
#######Expo 2:  In First Round Intensive (FRI)#############################
#######Expo 3:  In Second Round Simple (SRS), no friend FR#################
#######Expo 4:  In SRS, no friends FRI#####################################
#######Expo 5:  In SRS, one friend in FI###################################
#######Expo 6:  In SRS, two friend in FRI##################################
#######Expo 7:  In SRS, >two friends in FRI################################
#######Expo 8:  In Second Round Intensive (SRI), no friend FR##############
#######Expo 9:  In SRI, no friend FRI######################################
#######Expo 10: In SRI, one friend FRI#####################################
#######Expo 11: In SRI, two friends FRI####################################
#######Expo 12: In SRI, >two friends FRI###################################
###########################################################################
print(expo1)
print(expo2)
print(simulation_index)
print(option)
#compoents for this iteration
id_t = net %v% 'vertex.names'
degree_t=degree(net,cmode='outdegree') #warning some units have only two friends so some exposure mappings are not defined

#depending on exposure mapping some units must be dropped for identification reasons:
if (expo1 == 7 | expo1 == 12 | expo2 ==7 | expo2==12){
  
  number_dropped_t =sum( (id_t[degree_t<3] %in% id_list))

  #delete units with less than three friends
  id_list=id_list[!(id_list %in% id_t[degree_t<3])] 

}else if (expo1 == 6 | expo1 == 11 | expo2 ==6 | expo2==11){
  
  number_dropped_t =sum( (id_t[degree_t<2] %in% id_list))

  #delete units with less than two friends
  id_list=id_list[!(id_list %in% id_t[degree_t<2])]
  
}else{
  
  number_dropped_t =sum( (id_t[degree_t<1] %in% id_list))
  #delete units with no friends
  id_list=id_list[!(id_list %in% id_t[degree_t<1])] 

  
}

#How many units dropped?
prompt_t = paste0(number_dropped_t,' units dropped.')
print(prompt_t)

#components=component.dist(net,connected='weak')
#memberships=components$membership

# comp_t = get.inducedSubgraph(net,which(memberships==network_index))
# id_t = comp_t %v% 'vertex.names'
# 
# ind = any(id_t %in% id_list)

#########################################################
#####compute FOSO for all components#####################
#########################################################
output=compute_FOSO_all_components(net,expo1,expo2,id_list,option,prob=prob)

#########################################################
#####save output#########################################
#########################################################
save(output,file=save_path)





