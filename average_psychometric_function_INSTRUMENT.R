#prepare analysis averaging over participants

library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
#first, create tables for conditions and subjects with preprocessing

datalist<-ls(pattern="AV_Ins_av[0-9]{2}[a-z]{1}")
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
  thisdata$ISI_jitter<-as.factor(thisdata$rt)
  thisdata$rt<-as.factor(thisdata$ISI_jitter)
  thisdata$trial<-as.factor(thisdata$trial)
  thisdata$onset<-as.factor(thisdata$onset)
  thisdata$pairs<-as.factor(thisdata$pairs)
  thisdata$key_num<-dplyr::recode(thisdata$key_num, '0' = '1', '1' = '0')
  thisdata$key_num<-as.numeric(thisdata$key_num)
  thisdata$use.trigger<-as.logical(thisdata$use.trigger)
  thisdata$frameRate<-as.factor(thisdata$frameRate)
  mylist[[i]]<-thisdata
}
allsubjects_INSTRUMENT<-dplyr::bind_rows(mylist)
allsubjects_INSTRUMENT$morphrate1<-as.numeric(as.character(allsubjects_INSTRUMENT$morphrate1))


############################
#fit and plot psychometric function
############################
library(quickpsy)
library(DEoptim)
library(parallel)
library(gridExtra)
#all pairs
allsubjects_INSTRUMENT <- na.omit(allsubjects_INSTRUMENT) 
fit <- quickpsy(d=allsubjects_INSTRUMENT, x=morphrate1, k=key_num, grouping = .(pairs), B = 100)
pcurve<-plotcurves(fit)
pcurve<-pcurve+scale_x_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                                  labels=c("0% instrument, 100% voice","20%","40%","60%", "80%", "100% instrument, 0% voice"),limits=c(0,1)) +
  scale_y_continuous(name ="probability for 'instrument'")+
  ggtitle("Instrument-task (all subjects)")+theme_bw()+  theme(text = element_text(size=20))

parplot<-plotpar(fit)
parplot<-parplot+scale_y_continuous(name ="morphing steps")+ggtitle("Instrument-task (all subjects)")+theme_bw()

thresplot<-plotthresholds(fit, geom = 'point')
thresplot<-thresplot+ggtitle("Instrument-task (all subjects)")+theme_bw()+
  scale_y_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                     labels=c("0% instrument, 100% voice","20%","40%","60%", "80%", "100% instrument, 0% voice"),limits=c(0,1))+  theme(text = element_text(size=20))

pcurve
parplot
thresplot

grid.arrange(parplot, pcurve, thresplot, ncol = 2)

#one curve per subject, all trial

allsubjects_INSTRUMENT <- na.omit(allsubjects_INSTRUMENT) 
fit <- quickpsy(d=allsubjects_INSTRUMENT, x=morphrate1, k=key_num, grouping = .(participant), B = 100)
pcurve<-plotcurves(fit)
pcurve<-pcurve+scale_x_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                                  labels=c("0% instrument, 100% voice","20%","40%","60%", "80%", "100% instrument, 0% voice"),limits=c(0,1)) +
  scale_y_continuous(name ="probability for 'instrument'")+
  ggtitle("Instrument-task (all subjects)")+theme_bw()+  theme(text = element_text(size=20))

parplot<-plotpar(fit)
parplot<-parplot+scale_y_continuous(name ="morphing steps")+ggtitle("Instrument-task (all subjects)")+theme_bw()

thresplot<-plotthresholds(fit, geom = 'point')
thresplot<-thresplot+ggtitle("Instrument-task (all subjects)")+theme_bw()+
  scale_y_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                     labels=c("0% instrument, 100% voice","20%","40%","60%", "80%", "100% instrument, 0% voice"),limits=c(0,1))+  theme(text = element_text(size=20))

pcurve
parplot
thresplot

grid.arrange(parplot, pcurve, thresplot, ncol = 2)


#only instrument categories
allsubjects_INSTRUMENT <- na.omit(allsubjects_INSTRUMENT) 
fit <- quickpsy(d=allsubjects_INSTRUMENT, x=morphrate1, k=key_num, grouping = .(Category1), B = 100)

pcurve<-plotcurves(fit)
pcurve<-pcurve+scale_x_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                                  labels=c("0% instrument, 100% voice","20%","40%","60%", "80%", "100% instrument, 0% voice"),limits=c(0,1)) +
  scale_y_continuous(name ="probability for 'instrument'")+
  ggtitle("Instrument-task (all subjects)")+theme_bw()+  theme(text = element_text(size=20))


parplot<-plotpar(fit)
parplot<-parplot+scale_y_continuous(name ="morphing steps")+ggtitle("Instrument-task (all subjects)")+theme_bw()

thresplot<-plotthresholds(fit, geom = 'point')
thresplot<-thresplot+ggtitle("Instrument-task (all subjects)")+theme_bw()+
  scale_y_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                     labels=c("0% instrument, 100% voice","20%","40%","60%", "80%", "100% instrument, 0% voice"),limits=c(0,1))+  theme(text = element_text(size=20))

pcurve
parplot
thresplot

grid.arrange(parplot, pcurve, thresplot, ncol = 2)

#only gender
allsubjects_INSTRUMENT <- na.omit(allsubjects_INSTRUMENT) 
fit <- quickpsy(d=allsubjects_INSTRUMENT, x=morphrate1, k=key_num, grouping = .(Category2), B = 100)

pcurve<-plotcurves(fit)
pcurve<-pcurve+scale_x_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                                  labels=c("0% instrument, 100% voice","20%","40%","60%", "80%", "100% instrument, 0% voice"),limits=c(0,1)) +
  scale_y_continuous(name ="probability for 'instrument'")+
  ggtitle("Instrument-task (all subjects)")+theme_bw()

parplot<-plotpar(fit)
parplot<-parplot+scale_y_continuous(name ="morphing steps")+ggtitle("Instrument-task (all subjects)")+theme_bw()

thresplot<-plotthresholds(fit, geom = 'point')
thresplot<-thresplot+ggtitle("Instrument-task (all subjects)")+theme_bw()+
  scale_y_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                     labels=c("0% instrument, 100% voice","20%","40%","60%", "80%", "100% instrument, 0% voice"),limits=c(0,1))

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
  pcurve<-pcurve+scale_x_continuous(limits=c(0.0,1))+ggtitle(name)+ggtitle(name)
  print(pcurve)
  parplot<-plotpar(fit)
  parplot<-parplot+scale_y_continuous(limits=c(0.0,1))+ggtitle(name)
  print(parplot)
  thresplot<-plotthresholds(fit, geom = 'point')
  thresplot<-thresplot+scale_y_continuous(limits=c(0.0,1))+ggtitle(name)
  print(thresplot)
}

###########################
#all subjects but only morphs with Different PITCH (according to gender task) and only Vlc and Cl
###########################
T3<-filter(allsubjects_INSTRUMENT,allsubjects_INSTRUMENT$Category1=="Cla" | allsubjects_INSTRUMENT$Category1=="Vlc")
allsubjects_Instr_DIFF<-T3

fit <- quickpsy(d=allsubjects_Instr_DIFF, x=morphrate1, k=key_num, grouping = .(pairs), B = 100)

pcurve<-plotcurves(fit)
pcurve<-pcurve+scale_x_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                                  labels=c("0% instrument, 100% voice","20%","40%","60%", "80%", "100% instrument, 0% voice"),limits=c(0,1)) +
  scale_y_continuous(name ="probability for 'instrument'")+
  ggtitle("Instrument-task (all subjects)")+theme_bw()

parplot<-plotpar(fit)
parplot<-parplot+scale_y_continuous(name ="morphing steps")+ggtitle("Instrument-task (all subjects)")+theme_bw()

thresplot<-plotthresholds(fit, geom = 'point')
thresplot<-thresplot+ggtitle("Instrument-task (all subjects)")+theme_bw()+
  scale_y_continuous(name ="morphing steps", breaks =c(0,0.2,0.4,0.6,0.8,1), 
                     labels=c("0% instrument, 100% voice","20%","40%","60%", "80%", "100% instrument, 0% voice"),limits=c(0,1))
pcurve
parplot
thresplot

grid.arrange(parplot, pcurve, thresplot, ncol = 2)
