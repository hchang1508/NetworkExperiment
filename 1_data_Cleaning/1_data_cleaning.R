library(dplyr)

################################################
####Set working directory########################
################################################

setwd("/Users/haoge/Dropbox/Research/AnalysisExperiment/Cai(2015)/FinalData")

################################################
####To save processed data######################
################################################
target_people_info_path = '/Users/haoge/Dropbox/Research/AnalysisExperiment/Cai(2015)/FinalData/target_info.csv'
all_people_info_path = '/Users/haoge/Dropbox/Research/AnalysisExperiment/Cai(2015)/FinalData/all_info.csv'
network_info_path = '/Users/haoge/Dropbox/Research/AnalysisExperiment/Cai(2015)/FinalData/network_info.csv'
data_table2_path = '/Users/haoge/Dropbox/Research/AnalysisExperiment/Cai(2015)/FinalData/data_table2.csv'

################################################
####Import Dataset##############################
################################################
basic_info=read.csv("0422survey.csv")
network_info=read.csv("0422allinforawnet.csv")


################################################
####4902 households in the dataset##############
################################################
basic_info_id=unique(basic_info$id)
length(unique(basic_info$id))

# Only individuals in the second round + NoInfo + have 5 friends
data_table2 = basic_info[basic_info$delay == 1 & basic_info$info_none == 1, ] #data for households whose results are reported in Table 2 of Cat et al paper
########################################################
#Step 1: Dropping missing agpop and ricearea values#####
########################################################
basic_info=basic_info[ !(is.na(basic_info$agpop) | is.na(basic_info$ricearea_2010)) ,] #basic_info
basic_info_id=unique(basic_info$id)
length(basic_info_id)

########################################################
#Step 2: Dropping rows with missing network infroamtion#
########################################################
#remove NA entries for nominated friends
network_info2=network_info[!is.na(network_info$network_id),]
length(unique(network_info2$id))
nrow(network_info2)

########################################################
#Step 3: Remove self-nominations########################
########################################################
network_info3=network_info2[which(network_info2$id != network_info2$network_id ),]
length(unique(network_info3$id))
length(unique(network_info3$network_id))
nrow(network_info3)

########################################################
#Step 4: Remove repeated nominatos######################
########################################################
network_info4=distinct(network_info3,network_info3$id,network_info3$network_id,.keep_all=TRUE)
length(unique(network_info4$network_id))
length(unique(network_info4$id))

#remaining people
nominator=unique(network_info4$id)
nominated=unique(network_info4$network_id)


length(setdiff(nominator,basic_info_id))
length(setdiff(nominated,basic_info_id))

########################################################
#Step 5: Remove people with no information##############
########################################################
network_info5=network_info4[ (network_info4$id %in% basic_info_id),]
network_info6=network_info5[ (network_info5$network_id %in% basic_info_id),]  #notice missing value (99) is automatically excluded
(nrow(network_info3)-nrow(network_info5))/nrow(network_info3)
(nrow(network_info3)-nrow(network_info6))/nrow(network_info3)
length(unique(network_info5$id))
length(unique(network_info6$id))


#########################################################################
######Some Investigation#################################################
#How many nomiated households with treatment information are dropped#####
#########################################################################
network_info_check=network_info5[ !(network_info5$network_id %in% basic_info_id),] #nominated friends have no record in the basic_info
before=length(unique(network_info_check$network_id)) #2757 nominations havve no information
network_info_check=distinct(network_info_check,network_info_check$network_id,network_info_check$intensive,network_info_check$delay,.keep_all=FALSE) #remove repeated records

network_info_check1=network_info_check[(is.na(network_info_check$`network_info_check$intensive`) | is.na(network_info_check$`network_info_check$delay`)) ,] #some people has treatment record in the system
after1=length(unique(network_info_check1$`network_info_check$network_id`))
before-after1

network_info_check2=network_info_check[(!is.na(network_info_check$`network_info_check$intensive`) & !is.na(network_info_check$`network_info_check$delay`)) ,] #some people has treatment record in the system
after2=length(unique(network_info_check2$`network_info_check$network_id`))
after2 #two households with both missing and complete treatment records

#########################################################################
######Final Data#########################################################
#########################################################################

friend=network_info6[,c('id','network_id')]
friend %>% 
  group_by(id) %>%
  summarise(no_rows = length(network_id))

#final list of candidates
all_people=union(network_info6$id,network_info6$network_id) #we have covariates information for all people here
target_people=unique(network_info6$id) #we have outcome informaton, covariate information and frendship information for all these people; this is the population of interest
basic_info_target = basic_info[basic_info$id %in% target_people,]

#merging covariate data with all people and target people
all_people=as.data.frame(all_people)
colnames(all_people)=c('id')
all_people_info=merge(all_people,basic_info,by='id')

target_people=as.data.frame(target_people)
colnames(target_people)=c('id')
target_people_info=merge(target_people,basic_info,by='id')


#filtering data_table2 based on the availablity of network information
data_table2_filtered=data_table2[ which(data_table2$id %in% as.matrix(target_people)), ]
nrow(data_table2_filtered)
nrow(data_table2) - nrow(data_table2_filtered)

#alternative way to extract the same data (for sanity check)
#data_table3 = target_people_info[target_people_info$delay == 1 & target_people_info$info_none == 1, ] #data for households whose results are reported in Table 3 of Cat et al paper


#########################################################################
######Some Investigation#################################################
#Who are the subjects dropped with the filtered data#####################
#########################################################################
#data_table2_filtered=data_table2[ which( !(data_table2$id %in% target_people)), ]
#View(data_table2_filtered)
#id=data_table2_filtered[,'id']

#many dropped households did not nominate any friends
#for (i in 1:length(id)){
#  print(i)
#  print(network_info[which(network_info$id==id[i]),])
#}

#id_with_record=c()
#for (i in 1:length(id)){
#  if (any(!is.na(network_info[which(network_info$id==id[i]),'network_id']))==TRUE){
#    id_with_record=c(id_with_record,id[i])
#  }
#}
#network_info[which(network_info$id==id_with_record[i]),]


######################################################################################################################################################
####Households may be dropped if they or their nominated friends don't have the relevant information for stratification###############################
######################################################################################################################################################
#nomination=network_info[which(network_info$id %in% id_with_record),c('id','network_id')]
#nominated_friend=network_info[which(network_info$id %in% id_with_record),c('network_id')]

#####Some households are dropped because their nominated friends don't have the relevant information####################
#####But three ohuseholds are dropped because they don't have the relevant information themselves (rice area 2010) #####
######################################################################################################################################################
#nominated_friend_with_record=nominated_friend[nominated_friend %in% basic_info_id]
#nomination_with_record=nomination[which(nomination$network_id %in% nominated_friend_with_record),]
#basic_info[which(basic_info$id%in% nomination_with_record$id),]


#Check if have missing entries
sum(is.na(all_people_info$ricearea_2010))
sum(is.na(all_people_info$agpop))


write.csv(target_people_info,target_people_info_path)
write.csv(all_people_info,all_people_info_path)
write.csv(network_info6, network_info_path)
write.csv(data_table2_filtered, data_table2_path)

