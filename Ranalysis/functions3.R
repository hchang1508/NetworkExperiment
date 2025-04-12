compute_FOSO_all_components=function(net,expo1,expo2,prob,id_list,min_c_size){
  
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
  components=component.dist(net,connected='weak')
  memberships=components$membership
  num_components=length(unique(memberships))
  csize =components$csize
  components_to_keep = (1:num_components)[csize>min_c_size]
  
  FO_all_components =list()
  SO_all_components=list()
  for( i in components_to_keep){
    prompt_t = paste0('Processing Component ',i)
    print(prompt_t)
    comp_t = get.inducedSubgraph(net,which(memberships==i))
    foso_t=compute_FOSO_per_cluster(comp_t,expo1,expo2,prob,id_list)
    FO_all_components[[i]]=foso_t[[1]]
    SO_all_components[[i]]=foso_t[[2]]
    
    }
  
  output=list(FO_all_components,SO_all_components)
  
  return(output)
}

compute_FOSO_per_cluster=function(net,expo1,expo2,prob,id_list){
  
  id_net = net %v% 'vertex.names'

  #First Order and Second Order Matrics
  prob=prob
  expo1=expo1
  expo2=expo2
  
  foso=FOSO(net=net, FUN=individual_exposure, num_status=4,prob=prob,num_mappings=12,nrep=10000,id_list,expo1,expo2)
  return(foso)
}