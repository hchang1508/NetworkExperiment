#packages 
library(statnet)

#set working path
working_path='/Volumes/GoogleDrive/My Drive/UnifiedDesignBasedInference/Cai(2015)Insurance/Ranalysis'
setwd(working_path)

#load data
data_path='/Volumes/GoogleDrive/My Drive/UnifiedDesignBasedInference/Cai(2015)Insurance/Ranalysis/data/truncated_graph'
data_path=paste0(data_path,'/net.RData')
load(data_path)

#functions
source('functions.R')
source('functions2.R')
source('functions3.R')


