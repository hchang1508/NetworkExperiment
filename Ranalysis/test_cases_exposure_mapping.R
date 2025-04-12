#test cases 1
exp_t1=c(1,0,0,0,0,0)
output_t1=individual_exposure(exp_t1)
anw_t1=c(1,rep(0,11))
print(all(anw_t1==output_t1))

#test cases 2
exp_t2=c(1,0,1,0,0,0)
output_t2=individual_exposure(exp_t2)
anw_t2=c(1,rep(0,11))
print(all(anw_t2==output_t2))

#test cases 3
exp_t3=c(1,0,0,0,1,1)
output_t3=individual_exposure(exp_t3)
anw_t3=c(1,rep(0,11))
print(all(anw_t3==output_t3))

#test cases 4
exp_t4=c(0,1,0,0,0,0)
output_t4=individual_exposure(exp_t4)
anw_t4=c(0,1,rep(0,10))
print(all(anw_t4==output_t4))

#test cases 5
exp_t5=c(0,1,0,0,1,1)
output_t5=individual_exposure(exp_t5)
anw_t5=c(0,1,rep(0,10))
print(all(anw_t5==output_t5))

#test cases 6
exp_t6=c(0,0,1,0,0,0)
output_t6=individual_exposure(exp_t6)
anw_t6=c(0,0,1,rep(0,9))
print(all(anw_t6==output_t6))

#test cases 7
exp_t7=c(0,0,1,0,1,0)
output_t7=individual_exposure(exp_t7)
anw_t7=c(0,0,0,1,rep(0,8))
print(all(anw_t7==output_t7))

#test cases 8
exp_t8=c(0,0,1,0,1,1)
output_t8=individual_exposure(exp_t8)
anw_t8=c(0,0,0,0,1,rep(0,7))
print(all(anw_t8==output_t8))

#test cases 9
exp_t9=c(0,0,1,0,2,1)
output_t9=individual_exposure(exp_t9)
anw_t9=c(0,0,0,0,1,rep(0,7))
print(all(anw_t9==output_t9))

#test cases 10
exp_t10=c(0,0,1,0,2,2)
output_t10=individual_exposure(exp_t10)
anw_t10=c(0,0,0,0,0,1,rep(0,6))
print(all(anw_t10==output_t10))

#test cases 11
exp_t11=c(0,0,1,0,2,4)
output_t11=individual_exposure(exp_t11)
anw_t11=c(0,0,0,0,0,0,1,rep(0,5))
print(all(anw_t11==output_t11))

#test cases 12
exp_t12=c(0,0,0,1,0,0)
output_t12=individual_exposure(exp_t12)
anw_t12=c(0,0,0,0,0,0,0,1,rep(0,4))
print(all(anw_t12==output_t12))

#test cases 13
exp_t13=c(0,0,0,1,1,0)
output_t13=individual_exposure(exp_t13)
anw_t13=c(0,0,0,0,0,0,0,0,1,rep(0,3))
print(all(anw_t13==output_t13))

#test cases 14
exp_t14=c(0,0,0,1,1,1)
output_t14=individual_exposure(exp_t14)
anw_t14=c(0,0,0,0,0,0,0,0,0,1,rep(0,2))
print(all(anw_t14==output_t14))

#test cases 15
exp_t15=c(0,0,0,1,2,1)
output_t15=individual_exposure(exp_t15)
anw_t15=c(0,0,0,0,0,0,0,0,0,1,rep(0,2))
print(all(anw_t15==output_t15))

#test cases 16
exp_t16=c(0,0,0,1,2,2)
output_t16=individual_exposure(exp_t16)
anw_t16=c(0,0,0,0,0,0,0,0,0,0,1,rep(0,1))
print(all(anw_t16==output_t16))

#test cases 17
exp_t17=c(0,0,0,1,2,3)
output_t17=individual_exposure(exp_t17)
anw_t17=c(0,0,0,0,0,0,0,0,0,0,0,1)
print(all(anw_t17==output_t17))