library(dplyr)
library(tidyr)
library(readr)

# in the following line, change 200713 to name of data from sosci survey
data <- read_delim("~/projects/ambiguous_voices/FMRI_study/learning/data/200902.csv", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)

#only take participants from the interview version (no pretest/admin)
data<-filter(data, data$MODE=="interview")



#---------------------------------------------
#####inclusion criteria#####
#---------------------------------------------
demographics<-select(data, starts_with("DE")) 
health<-select(data, starts_with("DE08")) 
music<-select(data, starts_with("DE09")) 
inclusion<-cbind(health, music)

for (i in 1:nrow(inclusion)){
  name<-data$DE01_01
  print(name[i])
  if (is.na(all(inclusion[i,])))
    print("do not inlcude/inspect")
  else if (all(inclusion[i,]==2 |3))
    print("include")
  else
    print("do not include/inspect")
  
}


#---------------------------------------------
#####learning task#####
#---------------------------------------------

#timbre (woodwind vs. string instrument)
woodwind<-select(data, num_range("IG", 31:38)) # answer woodwind ==1
string<-select(data, num_range("IG", 24:30)) #answer string==2
string<-cbind(string, data$IG01)

results<-data.frame()
results[1:nrow(data),1]<-data$DE01_01 

for (i in 1:nrow(string)){
corr_s<-sum(string[i,]==2)
corr_w<-sum(woodwind[i,]==1)
perc<-(corr_s+corr_w)/16
results[i,2]<-corr_s
results[i,3]<-corr_w
results[i,4]<-perc
}


#instrument

cello<-cbind(data$IG02,data$IG08, data$IG09, data$IG10)
viola<-select(data, num_range("IG", 11:14))
clarinet<-select(data, num_range("IG", 15:18))
bclarinet<-select(data, num_range("IG", 19:22))

for (i in 1:nrow(string)){
  corr_vc<-sum(cello[i,]==2)
  corr_va<-sum(viola[i,]==1)
  corr_cl<-sum(clarinet[i,]==3)
  corr_bcl<-sum(bclarinet[1,]==4)
  perc<-(corr_vc+corr_va+corr_cl+corr_bcl)/16
  results[i,5]<-corr_vc
  results[i,6]<-corr_va
  results[i,7]<-corr_cl
  results[i,8]<-corr_bcl
  results[i,9]<-perc
}


#baseline
str<-cbind(data$VG03,data$VG04, data$VG05, data$VG06)
ww<-cbind(data$VG07, data$VG08, data$VG09, data$VG10)
male<-select(data, num_range("VG", 11:14))
female<-select(data, num_range("VG", 15:18))

for (i in 1:nrow(string)){
  corr_m<-sum(male[i,]==1)
  corr_f<-sum(female[i,]==2)
  corr_str<-sum(str[i,]==3)
  corr_ww<-sum(ww[1,]==4)
  perc<-(corr_m+corr_f+corr_str+corr_ww)/16
  results[i, 10:14]<- c(corr_m, corr_f, corr_str, corr_ww, perc)
}

#format results table

colnames(results)<-c("code","corr_s","corr_w", "perc_T", 
                     "corr_vc","corr_va","corr_cl", "corr_bcl", "perc_I",
                     "corr_m", "corr_f", "corr_str", "corr_ww","perc_B")

#------------------------------
#function that outputs results per participant
#------------------------------
#example use: incltest("av32",data, results)
incltest<-function(vpcode, ds, perc){
  vp_data<- filter(ds, ds$DE01_01 == vpcode)
  vp_data<-as.data.frame((vp_data))
  vp_health<-select(vp_data, starts_with("DE08"))
  print(vp_health)
  vp_music<-select(vp_data, starts_with("DE09"))
  if (is.na(all(vp_health[i,])))
    print("health: na responses")
  else if (all(vp_health[i,]==2 |3))
    print("health criteria fulfilled")
  else
    print("health criteria flagged")
  
  if (is.na(all(vp_music[i,])))
    print("music: na responses")
  else if (all(vp_music[i,]==2 |3))
    print("music criteria fulfilled")
  else
    print("music criteria flagged")
  
  cat('age', vp_data$DE03_01, sep= ": ")
  cat(sep="\n")
  cat('sex (1=f, 2=m)', vp_data$DE04, sep = ": ")
  cat(sep="\n")
  cat('handedness (1=r, 2=l)', vp_data$DE05, sep = ": ")
  cat(sep="\n")
  cat('native language', vp_data$DE06_01, sep = ": ")
  cat(sep="\n")
  cat('englisch proficiency (native 1-7 bad)', vp_data$DE07, sep = ": ")
  cat(sep="\n")
  cat('vision', vp_data$DE15, 'drugs',vp_data$DE17_01, sep = "  ")
  cat(sep="\n")
  cat('nr instruments', vp_data$DE10, sep=": " )
  cat(sep="\n")
  cat('age instruments', vp_data$DE11x01, vp_data$DE11x02, vp_data$DE11x03,vp_data$DE11x04, sep = "  ")
  cat(sep="\n")
  cat('instr3/voc3/listening3', vp_data$DE12_01, vp_data$DE13_01, vp_data$DE14_01, sep = "  ")
  cat(sep="\n")
  
  vp_perc<- filter(perc, perc$code == vpcode)
  return(vp_perc)
}

