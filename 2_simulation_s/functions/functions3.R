
compute_FOSO_all_components=function(net,expo1,expo2,id_list,option){
  
  #compute first order and second order matrix of a network and save it as a list
  #Input includes
  # 1. net: a network object
  # 2. expo1: index for the first exposure mapping
  # 3. expo2: index for the second exposure mapping
  # 4. prob: probability of treatment arms (assuming Bernoulli)
  # 5. id_list: a list of units whose exposure values are of interest
  # 5a. The vertex name, an attribute of the network vertices, records id's.
  # 6. min_c_size: a threshold for minimum network size, used for dropping small networks
  
  #extract information for components
  memberships=net %v% "membership"
  num_components=length(unique(memberships))
  
  membership_list=unique(memberships)
  
  FO_all_components =list()
  SO_all_components=list()
  for( i in membership_list){
    prompt_t = paste0('Processing Component ',i)
    print(prompt_t)
    comp_t = get.inducedSubgraph(net,which(memberships==i))
    foso_t=compute_FOSO_per_cluster(comp_t,expo1,expo2,id_list,option)
    FO_all_components[[i]]=foso_t[[1]]
    SO_all_components[[i]]=foso_t[[2]]
    
  }
  
  output=list(FO_all_components,SO_all_components)
  
  return(output)
}

compute_FOSO_per_cluster=function(net,expo1,expo2,id_list,option){
  
  
  #First Order and Second Order Matrics
  expo1=expo1
  expo2=expo2
  
  foso=FOSO(net=net, num_status=4,prob=rep(0.25,4),num_mappings=12,nrep=10000,id_list,expo1,expo2,option=option)
  return(foso)
}