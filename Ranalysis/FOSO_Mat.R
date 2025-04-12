source('functions.R')
net = readRDS('net1.rds')


memberships = component.dist(net,connected=c('weak'))$membership
component_sizes = component.dist(net,connected=c('weak'))$csize
number_components = length(unique(memberships))

#drop small clusters: 10 is arbitrary here, 2 dropped
cluster_to_keep = (1:number_components)[component_sizes>10]
number_components-length(cluster_to_keep)

nrep=100

for (i in cluster_to_keep){
  
  prompt=paste0('calculating for cluster ', i)
  print(prompt)

  prob =  c(0.25,0.25,0.25,0.25)
  
  result = FOSO(net, FUN=individual_exposure,
                num_status=4,prob=prob,num_mapping=12,nrep=nrep,id_list=id_list,1,2)
  
  compare_ana = target_prob(prob,5)
  compare_com = result[[1]][1:12]
  
  diff = as.matrix(compare_ana-compare_com)
  norm(diff)
  }
  
  


