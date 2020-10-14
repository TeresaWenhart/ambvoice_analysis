library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(reshape2)

#first, create tables for conditions and subjects with preprocessing

datalist<-ls(pattern="AV_Bas_av[0-9]{2}[a-z]{1}")
mylist<-list()
for (i in 1:length(datalist)){
  thisdata<-get(datalist[[i]])
  thisdata$morphrate1<-as.numeric(as.character(thisdata$morphrate1))
  thisdata$morphrate2<-as.numeric(as.character(thisdata$morphrate2))
  thisdata$Category1<-as.factor(thisdata$Category1)
  thisdata$Category2<-as.factor(thisdata$Category2)
  thisdata$Prototype1<-as.factor(thisdata$Prototype1)
  thisdata$Prototype2<-as.factor(thisdata$Prototype2)
  thisdata$ISI_jitter<-as.factor(thisdata$ISI_jitter)
  thisdata$rt<-as.factor(thisdata$rt)
  thisdata$trial<-as.factor(thisdata$trial)
  thisdata$onset<-as.factor(thisdata$onset)
  thisdata$keycorr<-dplyr::recode(thisdata$Category1, 'Fem' = 1, 'Mal' = 2, 'Vla' = 3, 'Vlc' = 3, 'Cla'=4, 'Bcl'=4, .default=99)
  thisdata$correct<-thisdata$key_num==thisdata$keycorr
  thisdata$correct<-as.character(thisdata$correct)
  thisdata$correct<-dplyr::recode(thisdata$correct, 'TRUE' = 1, 'FALSE' = 0)
  thisdata$use.trigger<-as.logical(thisdata$use.trigger)
  thisdata$frameRate<-as.factor(thisdata$frameRate)
  mylist[[i]]<-thisdata
}
allsubjects_Bas<-dplyr::bind_rows(mylist)
allsubjects_Bas$morphrate1<-as.numeric(as.character(allsubjects_Bas$morphrate1))

#performance on each category
df<-data.frame()
rm(name)
name<-c()
for (i in 1:length(mylist)){
  thisdata<-mylist[[i]]
  name[i]<-thisdata$participant[1]
  FEM<-sum(thisdata$correct[thisdata$Category1=="Fem"], na.rm= TRUE)/sum(thisdata$Category1=="Fem", na.rm= TRUE)
  MAL<-sum(thisdata$correct[thisdata$Category1=="Mal"], na.rm= TRUE)/sum(thisdata$Category1=="Mal", na.rm= TRUE)
  STRING<-sum(thisdata$correct[thisdata$keycorr==3], na.rm= TRUE)/sum(thisdata$keycorr==3, na.rm= TRUE)
  WOOD<-sum(thisdata$correct[thisdata$keycorr==4], na.rm= TRUE)/sum(thisdata$keycorr==4, na.rm= TRUE)
  perc<-c(FEM,MAL, STRING,WOOD)
  df[i, 1:4]<-perc
  }
row.names(df)<-name
colnames(df)<-c("fem","mal", "string", "wood")
df$subj<-rownames(df)

# bring your data to long format as needed by ggplot
df.molten <- melt(df, value.name="Count", variable.name="Variable", na.rm=TRUE)
#plot
ggplot(df.molten, aes(y=Count, x=Variable, fill=subj)) +  
  geom_bar(stat="identity", position="dodge")+ theme_bw()+
  ggtitle("Baseline performance")+xlab("sound category")+ylab("percent correct")+theme(text = element_text(size=20))

