rm(list=ls())
library(statnet)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

#setwd("/Users/haoge/Library/CloudStorage/GoogleDrive-hgchang1508@gmail.com/My Drive/UnifiedDesignBasedInference/Cai(2015)Insurance/FinalData/")
setwd("/Volumes/GoogleDrive/My Drive/UnifiedDesignBasedInference/Cai(2015)Insurance/FinalData/")

#load data
network_info=read.csv("network_info.csv")
length(unique(network_info$id))

all_info=read.csv("all_info.csv")

#create a net object
relation_data <- network_info[,c('id','network_id')]
net <- network(relation_data,matrix.type="edgelist")

#administrative villages
id_list=net %v% "vertex.names"
villages= all_info[match(id_list,all_info$id),"village"] 
sum(is.na(villages)) 
net %v% "village" = villages

#natural village
id_list=net %v% "vertex.names"
nat_villages= all_info[match(id_list,all_info$id),"address"] 
sum(is.na(nat_villages)) 
net %v% "nat_village" = nat_villages


#47 villages
length(unique(net%v%"village"))
#number of households
network.size(net)


#max in and out degree degree
max(degree(net,cmode="outdegree"))
max(degree(net,cmode="indegree"))
max(degree(net,cmode="freeman"))

#inspect graph
components(net,'weak')

#number of households
network.size(net)

#max degree
max(degree(net,gmode="graph"))

#plotting the graph
par(mfrow=c(1,1))

#for coloring
village_color = as.factor(get.vertex.attribute(net,'village'))
coul <- brewer.pal(11,'Spectral')
coul <- colorRampPalette(coul)(length(levels(village_color)))
legend_village = as.numeric(as.character(factor(village_color)))
temp=cbind(coul[village_color],legend_village )
temp=data.frame(temp)
color_list=distinct(temp)
color_list=color_list[order(as.numeric(color_list[,2])),]

#ploting
pdf(file="Friendships_Complete_Village_Info.pdf",
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")

gplot(net,gmode='graph',mode="kamadakawai",vertex.col = coul[village_color], vertex.cex = 1.5,main="Reported Friendships, Complete Village Info")
legend("bottomleft",legend=color_list[order(as.numeric(color_list[,2]))[25:47],2],col=color_list[order(as.numeric(color_list[,2]))[25:47],1],pch=19,pt.cex=1.5,bty="n")
legend("bottomright",legend=color_list[order(as.numeric(color_list[,2]))[1:24],2],col=color_list[order(as.numeric(color_list[,2]))[1:24],1],pch=19,pt.cex=1.5,bty="n",
       title="Admin Villages")
dev.off()


#plotting a component
membership=component.dist(net,connected='weak')$membership
membership_info= cbind(1:length(membership),membership)#info for graphing

compare=cbind(get.vertex.attribute(net,"membership"),get.vertex.attribute(net,"nat_village"))
colnames(compare)=c('membership','nat_village')

##check village belongings
nat_villages= unique(compare[,2])
result=matrix(0,length(nat_villages),2)
result[,1]=nat_villages
for (i in 1:length(nat_villages)){

  temp_village=nat_villages[i]
  membership_temp=compare[which(compare[,2]==temp_village),1]
  result[i,2]=length(unique(membership_temp))
}

##check membership sizes
nat_cluster= unique(compare[,1])
result_c=matrix(0,length(nat_cluster),2)
result_c[,1]=nat_cluster
for (i in 1:length(nat_cluster)){
  
  temp_cluster=nat_cluster[i]
  membership_temp=compare[which(compare[,1]==temp_cluster),2]
  result_c[i,2]=length(unique(membership_temp))
}

#No.7 component has the largest size 1914
#No.8 component has the second largest size 577 
table(membership)
####Visulize the 8th component
pdf("directed_graphs_per_component.pdf" ,
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")
progress=paste0('ploting ',7, 'th component....')
print(progress)
#extract ith subgraph
subgraph_info = membership_info[which(membership_info[,2]==8),]
subgraph_vertex = as.numeric(subgraph_info[,1])
subgraph=get.inducedSubgraph(net, subgraph_vertex, alters = NULL, eid = NULL)

#for coloring
village_color = as.factor(get.vertex.attribute(subgraph,'nat_village'))
coul <- brewer.pal(11,'Spectral')
coul <- colorRampPalette(coul)(length(levels(village_color)))
coul = sample(coul,length(coul)) #break color orders randomly
legend_village = as.numeric(as.character(factor(village_color)))
temp=cbind(coul[village_color],legend_village )
temp=data.frame(temp)
color_list=distinct(temp)
color_list=color_list[order(as.numeric(color_list[,2])),]

#title=paste0("Reported Friendships",", 8th Component")
title=paste0("Reported Friendships",", Second Largest Component")
gplot(subgraph,gmode='graph',mode="kamadakawai", vertex.col=coul[village_color],
      vertex.cex = 1.5,main=title)
legend("bottomleft",legend=color_list[,2],col=color_list[,1],pch=19,pt.cex=1.5,bty="n",title='Natural Village')
dev.off()

####Check correct match village assignment

index=sample(1:4000,20)
names = net %v% "vertex.names"
names = names[index]
village = net %v% "village"
village = village[index]
expect=cbind(names,village)
correct=all_info[which(all_info$id %in% names),c('id',"village")]
cbind(correct[match(expect[,'names'],correct$id),],expect)



