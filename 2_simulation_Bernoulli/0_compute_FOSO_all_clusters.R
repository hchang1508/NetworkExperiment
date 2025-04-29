rm(list=ls())
source('0_header.R')
input= commandArgs(trailingOnly=TRUE)

expo1=as.numeric(input[1])
expo2=as.numeric(input[2])
case_to_do=read.csv(paste0('/home/hc654/NetworkExperiment/2_simulation_Bernoulli/case_to_do_b',expo1,expo2,'.csv'))

print(input)
case_to_do=case_to_do[,2]

simulation_index=as.numeric(case_to_do[as.numeric(input[3])])

option='Bernoulli'

save_path=paste0('/home/hc654/palmer_scratch/final_analysis_Bernoulli/' ,option,'_nat_',expo1,'_',expo2,'_FOSO_',simulation_index,'.Rdata')

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
output=compute_FOSO_all_components(net,expo1,expo2,id_list,option)

#########################################################
#####save output#########################################
#########################################################
save(output,file=save_path)





