
temp=colnames(all_info2)

temp2=all_info2[,c('id','group')]
temp2['assignment']=99

group_index=unique(temp$group)

dictionary=c()
for (i in group_index){
  
  ind=temp2$group==i
  nobs=sum(ind)
  
  integer_part= nobs %/% 4
  
  reminder_part = nobs - 4 * integer_part
  
  temp3=4:1
  temp=rep(c(4,3,2,1),integer_part+1)
  temp=temp[1:nobs]
  temp=sample(temp,length(temp),replace=FALSE)
  
  temp2[ind,'assignment']=temp
  
}