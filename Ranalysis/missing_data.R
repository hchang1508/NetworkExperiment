covars = c('id','male','age', 'agpop',
           'educ', 'ricearea_2010',
           'rice_inc', 'disaster_yes', 'disaster_loss' ,
           'risk_averse', 'disaster_prob', 'literacy', 'understanding')
X=covar[,covars]
percent = rep(0,ncol(X))
for (i in 1:ncol(X)){
  percent[i]=sum(is.na(X[,i]))/nrow(X)
  
}
names(percent)=covars
