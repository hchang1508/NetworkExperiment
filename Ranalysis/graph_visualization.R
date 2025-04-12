library(statnet)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
rm(list=ls())
load('data/truncated_graph/net_complete.RData')
#plotting
id_list=net %v% "vertex.names"
villages= X[match(id_list,X$id),"village"] 
villages[is.na(villages)] = 52 
net %v% "village" = villages
  
#47 villages
length(unique(net%v%"village"))
#number of households
network.size(net)

#max in and out degree degree
max(degree(net,cmode="outdegree"))
max(degree(net,cmode="indegree"))
max(degree(net,cmode="freeman"))
#We first analyze directed graph
components(net,'weak')

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
color_list[,2]=as.character(color_list[,2])
color_list[48,2]='Missing'
#for ploting
pdf(file="Friendships_All_Units.pdf",
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")

gplot(net,gmode='graph',mode="kamadakawai",vertex.col = coul[village_color], vertex.cex = 1.5,main="Reported Friendships, All Units")
legend("bottomleft",legend=color_list[order(as.numeric(color_list[,2]))[25:48],2],col=color_list[order(as.numeric(color_list[,2]))[25:48],1],pch=19,pt.cex=1.5,bty="n")
legend("bottomright",legend=color_list[order(as.numeric(color_list[,2]))[1:24],2],col=color_list[order(as.numeric(color_list[,2]))[1:24],1],pch=19,pt.cex=1.5,bty="n",
       title="Villages")
dev.off()


#Analyze undirected graph

net2=get.inducedSubgraph(net,which(net %v% 'village'!=52))
#number of households
network.size(net2)

#max degree
max(degree(net2,gmode="graph"))

#We analyze undirected graph
components(net2,'weak')

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

gplot(net2,gmode='graph',mode="kamadakawai",vertex.col = coul[village_color], vertex.cex = 1.5,main="Reported Friendships, Complete Village Info")
legend("bottomleft",legend=color_list[order(as.numeric(color_list[,2]))[25:47],2],col=color_list[order(as.numeric(color_list[,2]))[25:47],1],pch=19,pt.cex=1.5,bty="n")
legend("bottomright",legend=color_list[order(as.numeric(color_list[,2]))[1:24],2],col=color_list[order(as.numeric(color_list[,2]))[1:24],1],pch=19,pt.cex=1.5,bty="n",
       title="Villages")
dev.off()


####Visulize the 8th component
pdf("directed_graphs_per_component.pdf" ,
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")
progress=paste0('ploting ',8, 'th component....')
print(progress)
#extract ith subgraph
subgraph_info = membership_info[which(membership_info[,2]==8),]
subgraph_vertex = as.numeric(subgraph_info[,1])
subgraph=get.inducedSubgraph(net2, subgraph_vertex, alters = NULL, eid = NULL)

#for coloring
village_color = as.factor(get.vertex.attribute(subgraph,'village'))
coul <- brewer.pal(11,'Spectral')
coul <- colorRampPalette(coul)(length(levels(village_color)))
coul = sample(coul,length(coul)) #break color orders randomly
legend_village = as.numeric(as.character(factor(village_color)))
temp=cbind(coul[village_color],legend_village )
temp=data.frame(temp)
color_list=distinct(temp)
color_list=color_list[order(as.numeric(color_list[,2])),]

title=paste0("Reported Friendships",", 8th Component")
gplot(subgraph,gmode='graph',mode="kamadakawai", vertex.col=coul[village_color],
      vertex.cex = 1.5,main=title)
legend("bottomleft",legend=color_list[,2],col=color_list[,1],pch=19,pt.cex=1.5,bty="n",title='village index')
dev.off()

####Check correct match village assignment
index=sample(1:4000,20)
names = net %v% "vertex.names"
names = names[index]
village = net %v% "village"
village = village[index]

expect=cbind(names,village)

correct=covar[which(covar$id %in% names),c('id',"village")]

cbind(correct[match(expect[,'names'],correct$id),],expect)

