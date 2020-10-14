#prepare analysis averaging over participants

library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)

#first, create tables for conditions and subjects with preprocessing

datalist<-ls(pattern="AV_Gen_av[0-9]{2}[a-z]{1}")
mylist<-list()
for (i in 1:length(datalist)){
  thisdata<-get(datalist[[i]])
  thisdata$pairs<-paste(thisdata$Category1, thisdata$Prototype1,# new conditions
                        thisdata$Category2,thisdata$Prototype2)
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
  thisdata$pairs<-as.factor(thisdata$pairs)
  thisdata$use.trigger<-as.logical(thisdata$use.trigger)
  thisdata$frameRate<-as.factor(thisdata$frameRate)
  mylist[[i]]<-thisdata
}
allsubjects_Gender<-dplyr::bind_rows(mylist)
allsubjects_Gender$morphrate1<-as.numeric(as.character(allsubjects_Gender$morphrate1))


############################
#fit and plot psychometric function
############################
library(quickpsy)
library(DEoptim)
library(parallel)
library(gridExtra)

allsubjects_Gender <- na.omit(allsubjects_Gender) 
fit <- quickpsy(d=allsubjects_Gender, x=morphrate2, k=key_num, grouping = .(pairs), B = 100)
pcurve<-plotcurves(fit)
pcurve<-pcurve+scale_x_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
labels=c("0% male, 100% female","20%","40%","60%", "80%", "100% male, 0% female"),limits=c(0,1)) +
  scale_y_continuous(name ="probability for 'male'")+
  ggtitle("Gender-task (all subjects)")+theme_bw()  +theme(text = element_text(size=20))


parplot<-plotpar(fit)
parplot<-parplot+ggtitle("Gender (all subjects) - parameters")+theme_bw()+scale_y_continuous(name ="morphing steps")

thresplot<-plotthresholds(fit, geom = 'point')
thresplot<-thresplot+ggtitle("Gender-task (all subjects)")+theme_bw()+
  scale_y_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                     labels=c("0% male, 100% female","20%","40%","60%", "80%", "100% male, 0% female"),limits=c(0,1))+
  theme(text = element_text(size=20))
thresplot
pcurve
parplot

grid.arrange(parplot, pcurve, thresplot, ncol = 2)

#one curve per subject, all data
allsubjects_Gender <- na.omit(allsubjects_Gender) 
fit <- quickpsy(d=allsubjects_Gender, x=morphrate2, k=key_num, grouping = .(participant), B = 100)
pcurve<-plotcurves(fit)
pcurve<-pcurve+scale_x_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                                  labels=c("0% male, 100% female","20%","40%","60%", "80%", "100% male, 0% female"),limits=c(0,1)) +
  scale_y_continuous(name ="probability for 'male'")+
  ggtitle("Gender-task (all subjects)")+theme_bw()  +theme(text = element_text(size=20))


parplot<-plotpar(fit)
parplot<-parplot+ggtitle("Gender (all subjects) - parameters")+theme_bw()+scale_y_continuous(name ="morphing steps")

thresplot<-plotthresholds(fit, geom = 'point')
thresplot<-thresplot+ggtitle("Gender-task (all subjects)")+theme_bw()+
  scale_y_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                     labels=c("0% male, 100% female","20%","40%","60%", "80%", "100% male, 0% female"),limits=c(0,1))+
  theme(text = element_text(size=20))
thresplot
pcurve
parplot

grid.arrange(parplot, pcurve, thresplot, ncol = 2)



###########################
#fit and plot psychometric function per subject
###########################
for (i in 1:length(mylist)){
  thisdata<-mylist[i]
  thisdata<-as.data.frame(thisdata)
  name<-datalist[[i]]
  thisdata <- na.omit(thisdata) 
  fit <- quickpsy(d=thisdata, x=morphrate2, k=key_num, grouping = .(pairs), B = 100)
  pcurve<-plotcurves(fit)
  pcurve<-pcurve+scale_x_continuous(limits=c(0.0,1))+ggtitle(name)+ggtitle(name)
  print(pcurve)
  #parplot<-plotpar(fit)
  #parplot<-parplot+scale_y_continuous(limits=c(0.0,1))+ggtitle(name)
  #print(parplot)
  #thresplot<-plotthresholds(fit, geom = 'point')
  #thresplot<-thresplot+scale_y_continuous(limits=c(0.0,1))+ggtitle(name)
  #print(thresplot)
}


