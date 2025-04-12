library(statnet)
library(UserNetR)
library(RColorBrewer)

#this requires data from "graph_visualization.R"

#extract information about memberships
membership=component.dist(net2,connected='weak')$membership
membership_info= cbind(1:length(membership),membership)#info for graphing

#No.1 compoent has the largest size 1914
#No.8 component has the second largest size 577 
table(membership)



pdf("directed_graphs_per_component.pdf" ,
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")

for(i in 1:length(unique(membership))) {
  
  progress=paste0('ploting ',i, 'th component....')
  print(progress)
  #extract ith subgraph
  subgraph_info = membership_info[which(membership_info[,2]==i),]
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
  
  title=paste0("Component " ,i, ", (Direct)")
  gplot(subgraph,gmode='graph',mode="kamadakawai", vertex.col=coul[village_color],
        vertex.cex = 1.5,main=title)
  legend("bottomleft",legend=color_list[,2],col=color_list[,1],pch=19,pt.cex=1.5,bty="n",title='village index')
  
  
  
}

dev.off()


