#preprocessing behavioural data from ambiguous voices fMRI study

#library
library(dplyr)
library(tidyr)
library(stringr)
#directory
setwd("~/projects/ambiguous_voices/FMRI_study/data/behavioral")

path<-paste(getwd())


#read version A files

sub.folders <- list.files(path=getwd(), pattern="av[0-9]{2}[a-z]{1}",include.dirs = TRUE)
for (j in sub.folders) {
  path<-paste(getwd(),"/",j,sep="")
  filenames = list.files(path=path,pattern="AV_[A-Za-z]{3}_av[0-9]{2}[a-z]{1}.csv")
  for (i in filenames) {
    name <- gsub(".csv","",i)
    name <-gsub("/","_",name)
    assign(name,read.csv(paste(path,"/",i,sep=""))) #read in the table and name as "name"
  }
}

#recode button presses for baseline condition
#version a
file_list<-ls(pattern="AV_Bas_av[0-9]{2}a")
for (i in 1:length(file_list)) {
  name=file_list[i]
  thistable<-get(file_list[i])
  thistable$key_num<-dplyr::recode(thistable$key, b = 1, z = 2, g = 3, r = 4, .default=99)
  thistable$key_num[thistable$key_num==99]<-NA
  thistable$X<-NULL
  assign(name,thistable)
  i=i+1
}

#version b
file_list<-ls(pattern="AV_Bas_av[0-9]{2}b")
for (i in 1:length(file_list)) {
  name=file_list[i]
  print(i)
  thistable<-get(file_list[i])
  thistable$key_num<-dplyr::recode(thistable$key, b = 4, z = 3, g = 2, r = 1, .default=99)
  thistable$key_num[thistable$key_num==99]<-NA
  thistable$X<-NULL
  assign(name,thistable)
  i=i+1
}


#recode button presses to numbers 0 and 1 (opposite assignment for versions a and b)
#version a
file_list<-ls(pattern="AV_[A-Za-z]{3}_av[0-9]{2}a")
file_list<-file_list[str_detect(file_list, "Bas")==FALSE]
for (i in 1:length(file_list)) {
  name=file_list[i]
  thistable<-get(file_list[i])
  thistable$key_num<-dplyr::recode(thistable$key, z = 1, b = 0, .default=99)
  thistable$key_num[thistable$key_num==99]<-NA
  thistable$X<-NULL
  assign(name,thistable)
  i=i+1
}

#version b
file_list<-ls(pattern="AV_[A-Za-z]{3}_av[0-9]{2}b")
file_list<-file_list[str_detect(file_list, "Bas")==FALSE]
for (i in 1:length(file_list)) {
  name=file_list[i]
  print(i)
  thistable<-get(file_list[i])
  thistable$key_num<-dplyr::recode(thistable$key, z = 0, b = 1, .default=99)
  thistable$key_num[thistable$key_num==99]<-NA
  thistable$X<-NULL
  assign(name,thistable)
  i=i+1
}


