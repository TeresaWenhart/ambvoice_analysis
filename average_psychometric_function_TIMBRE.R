#prepare analysis averaging over participants

library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)

#first, create tables for conditions and subjects with preprocessing

datalist<-ls(pattern="AV_Tim_av[0-9]{2}[a-z]{1}")
mylist<-list()
for (i in 1:length(datalist)){
  thisdata<-get(datalist[[i]])
  thisdata$pairs<-paste(thisdata$Category1, thisdata$Prototype1,# new conditions
                        thisdata$Category2,thisdata$Prototype2)
  thisdata$pitch<-paste(thisdata$Prototype1,# new conditions
                        thisdata$Prototype2)
  thisdata$instrument<-paste(thisdata$Category1, # new conditions
                        thisdata$Category2)
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
allsubjects_TIMBRE<-dplyr::bind_rows(mylist)
allsubjects_TIMBRE$morphrate1<-as.numeric(as.character(allsubjects_TIMBRE$morphrate1))


############################
#fit and plot psychometric function
############################
library(quickpsy)
library(DEoptim)
library(parallel)
library(gridExtra)
#all pairs
allsubjects_TIMBRE <- na.omit(allsubjects_TIMBRE) 
fit <- quickpsy(d=allsubjects_TIMBRE, x=morphrate1, k=key_num, grouping = .(pairs), B = 100)

pcurve<-plotcurves(fit)
pcurve<-pcurve+scale_x_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                                  labels=c("0% woodwind, 100% string","20%","40%","60%", "80%", "100% woodwind, 0% string"),limits=c(0,1)) +
  scale_y_continuous(name ="probability for 'woodwind'")+
  ggtitle("Timbre-task (all subjects)")+theme_bw()  +theme(text = element_text(size=20))

parplot<-plotpar(fit)
parplot<-parplot+scale_y_continuous(name ="morphing steps")+ggtitle("Timbre-task (all subjects)")+theme_bw()

thresplot<-plotthresholds(fit, geom = 'point')
thresplot<-thresplot+ggtitle("Timbre-task (all subjects)")+theme_bw()+
  scale_y_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                     labels=c("0% woodwind, 100% string","20%","40%","60%", "80%", "100% woodwind, 0% string"),limits=c(0,1))  +theme(text = element_text(size=20))


pcurve
parplot
thresplot

grid.arrange(parplot, pcurve, thresplot, ncol = 2)

#one curve per subject, all trials

allsubjects_TIMBRE <- na.omit(allsubjects_TIMBRE) 
fit <- quickpsy(d=allsubjects_TIMBRE, x=morphrate1, k=key_num, grouping = .(participant), B = 100)

pcurve<-plotcurves(fit)
pcurve<-pcurve+scale_x_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                                  labels=c("0% woodwind, 100% string","20%","40%","60%", "80%", "100% woodwind, 0% string"),limits=c(0,1)) +
  scale_y_continuous(name ="probability for 'woodwind'")+
  ggtitle("Timbre-task (all subjects)")+theme_bw()  +theme(text = element_text(size=20))

parplot<-plotpar(fit)
parplot<-parplot+scale_y_continuous(name ="morphing steps")+ggtitle("Timbre-task (all subjects)")+theme_bw()

thresplot<-plotthresholds(fit, geom = 'point')
thresplot<-thresplot+ggtitle("Timbre-task (all subjects)")+theme_bw()+
  scale_y_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                     labels=c("0% woodwind, 100% string","20%","40%","60%", "80%", "100% woodwind, 0% string"),limits=c(0,1))  +theme(text = element_text(size=20))


pcurve
parplot
thresplot

grid.arrange(parplot, pcurve, thresplot, ncol = 2)



# pitch pairs
allsubjects_TIMBRE <- na.omit(allsubjects_TIMBRE) 
fit <- quickpsy(d=allsubjects_TIMBRE, x=morphrate1, k=key_num, grouping = .(pitch), B = 100)


pcurve<-plotcurves(fit)
pcurve<-pcurve+scale_x_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                                  labels=c("0% woodwind, 100% string","20%","40%","60%", "80%", "100% woodwind, 0% string"),limits=c(0,1)) +
  scale_y_continuous(name ="probability for 'woodwind'")+
  ggtitle("Timbre-task (all subjects)")+theme_bw()  +theme(text = element_text(size=20))

parplot<-plotpar(fit)
parplot<-parplot+scale_y_continuous(name ="morphing steps")+ggtitle("Timbre-task (all subjects)")+theme_bw()

thresplot<-plotthresholds(fit, geom = 'point')
thresplot<-thresplot+ggtitle("Timbre-task (all subjects)")+theme_bw()+
  scale_y_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                     labels=c("0% woodwind, 100% string","20%","40%","60%", "80%", "100% woodwind, 0% string"),limits=c(0,1))  +theme(text = element_text(size=20))


pcurve
parplot
thresplot

grid.arrange(parplot, pcurve, thresplot, ncol = 2)

#only instrument categories
allsubjects_TIMBRE <- na.omit(allsubjects_TIMBRE) 
fit <- quickpsy(d=allsubjects_TIMBRE, x=morphrate1, k=key_num, grouping = .(instrument), B = 100)


pcurve<-plotcurves(fit)
pcurve<-pcurve+scale_x_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                                  labels=c("0% woodwind, 100% string","20%","40%","60%", "80%", "100% woodwind, 0% string"),limits=c(0,1)) +
  scale_y_continuous(name ="probability for 'woodwind'")+
  ggtitle("Timbre-task (all subjects)")+theme_bw()

parplot<-plotpar(fit)
parplot<-parplot+scale_y_continuous(name ="morphing steps")+ggtitle("Timbre-task (all subjects)")+theme_bw()

thresplot<-plotthresholds(fit, geom = 'point')
thresplot<-thresplot+ggtitle("Timbre-task (all subjects)")+theme_bw()+
  scale_y_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                     labels=c("0% woodwind, 100% string","20%","40%","60%", "80%", "100% woodwind, 0% string"),limits=c(0,1))


pcurve
parplot
thresplot

grid.arrange(parplot, pcurve, thresplot, ncol = 2)

###########################
#fit and plot psychometric function per subject
###########################
for (i in 1:length(mylist)){
  thisdata<-mylist[i]
  thisdata<-as.data.frame(thisdata)
  name<-datalist[[i]]
  thisdata <- na.omit(thisdata) 
  fit <- quickpsy(d=thisdata, x=morphrate1, k=key_num, grouping = .(pairs), B = 10)
  pcurve<-plotcurves(fit)
  pcurve<-pcurve+scale_x_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                                    labels=c("0% woodwind, 100% string","20%","40%","60%", "80%", "100% woodwind, 0% string"),limits=c(0,1)) +
    scale_y_continuous(name ="probability for 'woodwind'")+
    ggtitle(name)+theme_bw() +theme(text = element_text(size=20))
  print(pcurve)
  parplot<-plotpar(fit)
  parplot<-parplot+scale_y_continuous(limits=c(0.0,1))+ggtitle(name)
  print(parplot)
  thresplot<-plotthresholds(fit, geom = 'point')
  thresplot<-thresplot+scale_y_continuous(limits=c(0.0,1))+ggtitle(name)
  print(thresplot)
}


###########################
#all subjects but only morphs with SAME PITCH
###########################
T1<-filter(allsubjects_TIMBRE,allsubjects_TIMBRE$Category1=="Cla"&allsubjects_TIMBRE$Category2=="Vla")
T2<-filter(allsubjects_TIMBRE,allsubjects_TIMBRE$Category1=="Bcl"&allsubjects_TIMBRE$Category2=="Vlc")
allsubjects_TIMBRE_SAME<-rbind(T1,T2)

fit <- quickpsy(d=allsubjects_TIMBRE_SAME, x=morphrate1, k=key_num, grouping = .(pairs), B = 100)
pcurve<-plotcurves(fit)
pcurve+scale_x_continuous(limits=c(0.0,1))+ggtitle("all subjects")
parplot<-plotpar(fit)
parplot+scale_y_continuous(limits=c(0.0,1))+ggtitle("all subjects")
thresplot<-plotthresholds(fit, geom = 'point')
thresplot+scale_y_continuous(limits=c(0.0,1))+ggtitle("all subjects")


###########################
#all subjects but only morphs with Different PITCH
###########################
T3<-filter(allsubjects_TIMBRE,allsubjects_TIMBRE$Category1=="Cla"&allsubjects_TIMBRE$Category2=="Vlc")
T4<-filter(allsubjects_TIMBRE,allsubjects_TIMBRE$Category1=="Bcl"&allsubjects_TIMBRE$Category2=="Vla")
allsubjects_TIMBRE_DIFF<-rbind(T3,T4)

fit <- quickpsy(d=allsubjects_TIMBRE_DIFF, x=morphrate1, k=key_num, grouping = .(pairs), B = 100)
pcurve<-plotcurves(fit)
pcurve+scale_x_continuous(limits=c(0.0,1))+ggtitle("all subjects")
parplot<-plotpar(fit)
parplot+scale_y_continuous(limits=c(0.0,1))+ggtitle("all subjects")
thresplot<-plotthresholds(fit, geom = 'point')
thresplot+scale_y_continuous(limits=c(0.0,1))+ggtitle("all subjects")



###########################
#all subjects but only morphs with Different PITCH (according to gender task) and only Vlc and Cl
###########################
T3<-filter(allsubjects_TIMBRE,allsubjects_TIMBRE$Category1=="Cla"&allsubjects_TIMBRE$Category2=="Vlc")
allsubjects_TIMBRE_DIFF<-T3

fit <- quickpsy(d=allsubjects_TIMBRE_DIFF, x=morphrate1, k=key_num, grouping = .(pairs), B = 100)
pcurve<-plotcurves(fit)
pcurve<-pcurve+scale_x_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                                  labels=c("0% woodwind, 100% string","20%","40%","60%", "80%", "100% woodwind, 0% string"),limits=c(0,1)) +
  scale_y_continuous(name ="probability for 'woodwind'")+
  ggtitle("Timbre-task (all subjects)")+theme_bw()

parplot<-plotpar(fit)
parplot<-parplot+scale_y_continuous(name ="morphing steps")+ggtitle("Timbre-task (all subjects)")+theme_bw()

thresplot<-plotthresholds(fit, geom = 'point')
thresplot<-thresplot+ggtitle("Timbre-task (all subjects)")+theme_bw()+
  scale_y_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                     labels=c("0% woodwind, 100% string","20%","40%","60%", "80%", "100% woodwind, 0% string"),limits=c(0,1))


pcurve
parplot
thresplot

grid.arrange(parplot, pcurve, thresplot, ncol = 2)

