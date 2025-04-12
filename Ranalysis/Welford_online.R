Welford_online=function(mean_pre,cov_pre,n,new_data){
 
  mean_new = (n-1)*mean_pre + new_data
  mean_new = mean_new / n
  diff= mean_pre-mean_new
  cov_new = cov_pre + (n-1) * (diff %*% t(diff)) + (new_data-mean_new) %*% t(new_data-mean_new)
  
  return(list(mean_new,cov_new))
}
  

Welford_online2=function(mean_pre,cov_pre,n,new_data){
  
  mean_new = (n-1)*mean_pre + new_data
  mean_new = mean_new / n
  diff= mean_pre-mean_new
  cov_new = cov_pre + (new_data-mean_pre) %*% t(new_data-mean_new)
  
  return(list(mean_new,cov_new))
}
data=cbind(rnorm(100,0,1),rnorm(100,0,1))
data[,1] = data[,1]+ data[,2]

mean_pre=0
cov_pre=matrix(0,2,2)
for (i in 1:100){
  
  new_data=data[i,]
  
  result = Welford_online(mean_pre,cov_pre,i,new_data)
  result2= Welford_online2(mean_pre,cov_pre,i,new_data)
  
  mean_pre=result[[1]]
  cov_pre=result[[2]]
  
  mean_pre2=result2[[1]]
  cov_pre2=result2[[2]]
  
}

print((cov_pre/99)-cov(data))
print((cov_pre2/99)-cov(data))

mean_pre-apply(data,2,mean)
