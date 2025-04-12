expo1=1
expo2=4

foso_file=paste0('/home/hc654/palmer_scratch/final_analysis_outputs/',expo1,expo2,'_D.Rdata')
load(foso_file)
fo=mean_pre
so=cov_pre
fo=unlist(fo)

problem_list=fo[order(fo)[1]]
problem_list

memberships = net %v% 'membership'
membership_list=unique(memberships)
membership_list=sort(memberships_list)
for( i in membership_list){
#  prompt_t = paste0('Processing Component ',i)
#  print(prompt_t)
  comp_t = get.inducedSubgraph(net,which(memberships==i))
  temp=unlist(mean_pre[i])
  network_subjects = comp_t %v% 'vertex.names'
  pol_size=network.size(comp_t)
  subjects = (1:pol_size)[network_subjects %in% id_list]
  subjects=rep(subjects,each=2)
  index=subjects[temp %in% problem_list]
  

  print(i)
  
  print(index)
  
}
#problem
#cluster 5
# index 190
#cluster 6
# index 26
#cluster 7
#index 91  518  573  594  611  851 1417
#cluster 8
#index 131
#cluster 30
#57

#cluster 5,6,7,30
net_5 = get.inducedSubgraph(net,which(memberships==5))
i=190
friend_index=get.neighborhood(net_5,i,'out')
group=net_5 %v% 'group'
(1:length(group))[group==group[i]]
friend_index
group[i]

net_6 = get.inducedSubgraph(net,which(memberships==6))
i=26
friend_index=get.neighborhood(net_6,i,'out')
group=net_6 %v% 'group'
(1:length(group))[group==group[i]]
friend_index
group[i]

net_7 = get.inducedSubgraph(net,which(memberships==7))
i=1377
friend_index=get.neighborhood(net_7,i,'out')
group=net_7 %v% 'group'
assignment=net_7 %v% 'assignment'
(assignment)[group==group[i]]
friend_index
group[i]

net_7 = get.inducedSubgraph(net,which(memberships==7))
i=518
friend_index=get.neighborhood(net_7,i,'out')
group=net_7 %v% 'group'
(1:length(group))[group==group[i]]
friend_index
group[i]

net_7 = get.inducedSubgraph(net,which(memberships==7))
i=573
friend_index=get.neighborhood(net_7,i,'out')
group=net_7 %v% 'group'
(1:length(group))[group==group[i]]
friend_index
group[i]

net_7 = get.inducedSubgraph(net,which(memberships==7))
i=852
friend_index=get.neighborhood(net_7,i,'out')
group=net_7 %v% 'group'
(1:length(group))[group==group[i]]
friend_index
group[i]

net_7 = get.inducedSubgraph(net,which(memberships==7))
i=1481
friend_index=get.neighborhood(net_7,i,'out')
group=net_7 %v% 'group'
(1:length(group))[group==group[i]]
friend_index
group[i]

net_8 = get.inducedSubgraph(net,which(memberships==8))
i=131
friend_index=get.neighborhood(net_8,i,'out')
group=net_8 %v% 'group'
(1:length(group))[group==group[i]]
friend_index
group[i]

net_9 = get.inducedSubgraph(net,which(memberships==30))
i=57
friend_index=get.neighborhood(net_9,i,'out')
group=net_9 %v% 'group'
(1:length(group))[group==group[i]]
friend_index
group[i]

# [1] 6
# [1] 26 27 30 32 33
# [1] 7
# [1] 567 734 851
# [1] 8
# [1] 136 154
# [1] 10
net_6 = get.inducedSubgraph(net,which(memberships==6))
i=27
friend_index=get.neighborhood(net_6,i,'out')
group=net_6 %v% 'group'
(1:length(group))[group==group[i]]
friend_index
group[i]

net_7 = get.inducedSubgraph(net,which(memberships==7))
i=567
friend_index=get.neighborhood(net_7,i,'out')
group=net_7 %v% 'group'
(1:length(group))[group==group[i]]
friend_index
group[i]

net_8 = get.inducedSubgraph(net,which(memberships==8))
i=154
friend_index=get.neighborhood(net_8,i,'out')
group=net_8 %v% 'group'
(1:length(group))[group==group[i]]
friend_index
group[i]

net_7 = get.inducedSubgraph(net,which(memberships==7))
i=39
friend_index=get.neighborhood(net_7,i,'out')
group=net_7 %v% 'group'
(1:length(group))[group==group[i]]
friend_index
group[i]



net_7 = get.inducedSubgraph(net,which(memberships==7))
i=509
friend_index=get.neighborhood(net_7,i,'out')
group=net_7 %v% 'group'
(1:length(group))[group==group[i]]
friend_index
group[i]


net_7 = get.inducedSubgraph(net,which(memberships==7))
i=1652
friend_index=get.neighborhood(net_7,i,'out')
group=net_7 %v% 'group'
name=net_7 %v% 'vertex.names'
(1:length(group))[group==group[i]]
name[group==group[i]]
friend_index
group[i]