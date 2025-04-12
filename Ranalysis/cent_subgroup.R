detach(package:statnet)
library(igraph)
library(intergraph)

id_list=net %v% "vertex.names"
net %v% "village" = covar[match(id_list,covar$id),"village"]
net1<- asIgraph(net)
V(net1)$village=net %v% "village" #include attributes
net1=as.undirected(net1,mode=c('collapse'))
membership=components(net1)$membership


#####Trying k-core algorithms#########
component14 = induced_subgraph(net1,which(membership==14))
ecount(component14)
component14=as.undirected(component14,mode=c('collapse'))
ecount(component14)
plot(component14)


coreness <- graph.coreness(component7)
table(coreness)
maxCoreness=max(coreness)
V(component14)$color <- coreness + 1
op <- par(mar = rep(0, 4))
plot(component14,vertex.label.cex=0.5)

colors <- rainbow(maxCoreness)
op <- par(mar = rep(0, 4))
plot(component14,vertex.label=coreness,
     vertex.color=colors[coreness])

V(component14)$color <- colors[coreness]
component141_5 <- component14
component142_5 <- induced.subgraph(component14,
                             vids=which(coreness > 1))
component143_5 <- induced.subgraph(component14,
                             vids=which(coreness > 2))
component144_5 <- induced.subgraph(component14,
                             vids=which(coreness > 3))
component145_5 <- induced.subgraph(component14,
                             vids=which(coreness > 4))

lay <- layout.fruchterman.reingold(component14)
op <- par(mfrow=c(3,2),mar = c(3,0,2,0))
plot(component141_5 ,layout=lay,main="All k-cores")
plot(component142_5 ,layout=lay[which(coreness > 1),],
     main="k-cores 2-5")
plot(component143_5 ,layout=lay[which(coreness > 2),],
     main="k-cores 3-5")
plot(component144_5 ,layout=lay[which(coreness > 3),],
     main="k-cores 4-5")
plot(component145_5 ,layout=lay[which(coreness > 4),],
     main="k-cores 5-5")


#####Centrality#######
detach(package:igraph)
library(statnet)
max(degree(net,gmode='digraph',cmode='indegree'))
max(degree(net,gmode='digraph',cmode='outdegree'))

#######Trying Modulairty##########
modularity(net1,V(net1)$village)



######Trying differecent centrality on Component14#####
membership=component.dist(net,connected='weak')$membership
c14 = get.inducedSubgraph(net,which(membership==14))
coul <- brewer.pal(5,'Spectral')
degree=degree(c14,gmode='graph')
degree_col=cut_number(degree,5)
btw = betweenness(c14,gmode='graph')
btw_col=cut_number(btw,5)
close = closeness(c14,gmode='graph')
close_col=cut_number(close, 5)
eigen=evcent(c14)
eigen_col=cut_number(eigen, 5)

legend_text=c('1st','2nd','3rd','4th','5th')
coords <- gplot(c14,mode='kamadakawai',vertex.cex = 1.0,gmode='graph')
pdf("compare_cent.pdf" ,
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")
par(mfrow=c(2,2))
gplot(c14,vertex.col = coul[degree_col],coord=coords,main='degree centrality',gmode='graph',label = 1:56,label.cex = 0.5)
legend('bottomleft',legend=legend_text,col=coul,pch=19,pt.cex=1.5,bty="n",title='quntile')
gplot(c14,vertex.col = coul[btw_col],coord=coords,main='between centrality',gmode='graph',label = 1:56,label.cex = 0.5)
legend('bottomleft',legend=legend_text,col=coul,pch=19,pt.cex=1.5,bty="n",title='quntile')
gplot(c14,vertex.col = coul[close_col],coord=coords,main='close centrality',gmode='graph',label = 1:56,label.cex = 0.5)
legend('bottomleft',legend=legend_text,col=coul,pch=19,pt.cex=1.5,bty="n",title='quntile')
gplot(c14,vertex.col = coul[eigen_col],coord=coords,main='eigen centrality',gmode='graph',label = 1:56,label.cex = 0.5)
legend('bottomleft',legend=legend_text,col=coul,pch=19,pt.cex=1.5,bty="n",title='quntile')
dev.off()