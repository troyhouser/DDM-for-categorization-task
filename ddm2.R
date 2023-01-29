tpc = read.csv("~/Dropbox (University of Oregon)/TPC_DDM/TPC_ddm.csv")
tpc = na.omit(tpc)
library(RWiener)
library(tidyverse)
library(tidyr)
library(reshape2)
library(dplyr)
many_drifts <- function (x, datlist) {
  loss=0
  for (myVar in 1:length(datlist))
    loss = loss - logLik(wdm(datlist[[myVar]], 
                             alpha=x[myVar], tau=x[myVar+5], beta=0.5, delta=x[myVar+10]))
  return(loss)
}
p_ddm5 = rep(c(1,0.1,0),each=5)

data1.ddm5 = list()
subs = unique(tpc$subj_idx)
chi_res = list()
for(sub in 1:length(subs)){
  df = tpc[tpc$subj_idx==subs[sub],]

  q = df$rt
  resp = ifelse(df$response==1,"lower","upper")
  cond = df$dist_newold
  
  data1.w = data.frame(q,resp,cond)
  table(data1.w$cond,data1.w$resp)
  data1.w$resp = factor(data1.w$resp)
  data1.w$resp = relevel(data1.w$resp,ref = "upper")
  data1.w = as.wiener(data1.w)
  data1.s1 = subset(data1.w,cond=="new0")
  data1.s2 = subset(data1.w,cond=="new1")
  data1.s3 = subset(data1.w,cond=="new2")
  data1.s4 = subset(data1.w,cond=="new3")
  data1.s5 = subset(data1.w,cond=="old2")

  data1.both = list(data1.s1,data1.s2,data1.s3,data1.s4,data1.s5)
  data1.ddm5[[sub]] = optim(par=p_ddm5,f=many_drifts,dat=data1.both)
  
  parm_recov_s1 <- cbind(rwiener(n=sum(cond=="new0"), alpha=data1.ddm5[[sub]]$par[1],
                                 tau=data1.ddm5[[sub]]$par[6],
                                 beta=.5,
                                 delta=data1.ddm5[[sub]]$par[11]), # v1 for s1
                         cond="new0")
  parm_recov_s2 <- cbind(rwiener(n=sum(cond=="new1"), alpha=data1.ddm5[[sub]]$par[2],
                                 tau=data1.ddm5[[sub]]$par[7],
                                 beta=.5,
                                 delta=data1.ddm5[[sub]]$par[12]), # v2 for s1
                         cond="new1")
  parm_recov_s3 <- cbind(rwiener(n=sum(cond=="new2"), alpha=data1.ddm5[[sub]]$par[3],
                                 tau=data1.ddm5[[sub]]$par[8],
                                 beta=.5,
                                 delta=data1.ddm5[[sub]]$par[13]), # v3 for s1
                         cond="new2")
  parm_recov_s4 <- cbind(rwiener(n=sum(cond=="new3"), alpha=data1.ddm5[[sub]]$par[4],
                                 tau=data1.ddm5[[sub]]$par[9],
                                 beta=.5,
                                 delta=data1.ddm5[[sub]]$par[14]), # v4 for s1
                         cond="new3")
  parm_recov_s5 <- cbind(rwiener(n=sum(cond=="old2"), alpha=data1.ddm5[[sub]]$par[5],
                                 tau=data1.ddm5[[sub]]$par[10],
                                 beta=.5,
                                 delta=data1.ddm5[[sub]]$par[15]), # v5 for s1
                         cond="old2")
  
  parm_recov.w <- as.wiener(rbind(parm_recov_s1, parm_recov_s2, parm_recov_s3, parm_recov_s4, parm_recov_s5))
  parm_recov.both <- list(subset(parm_recov.w, cond=="new0"), subset(parm_recov.w, cond=="new1"),
                          subset(parm_recov.w, cond=="new2"), subset(parm_recov.w, cond=="new3"),
                          subset(parm_recov.w, cond=="old2"))
  (x <- table(parm_recov.w$resp, parm_recov.w$cond))
  (x2 = table(df$response,df$dist_newold))
  chi_res[[sub]] = chisq.test(as.vector(x),as.vector(x2))
  
}
data1.ddm5

pars = matrix(0,length(data1.ddm5),15)
for(i in 1:length(data1.ddm5)){
  pars[i,] = data1.ddm5[[i]]$par
}
colnames(pars) = c("thr_0","thr_1","thr_2","thr_3","thr_old2",
                   "non_0","non_1","non_2","non_3","non_old2",
                   "drift_0","drift_1","drift_2","drift_3","drift_old2")
pars = as.data.frame(pars)
barplot(c(mean(pars$drift_0),mean(pars$drift_1),mean(pars$drift_2),mean(pars$drift_3)))


my.quantiles <- c(0.1, 0.3, 0.5, 0.7, 0.9) # copied from above
quantiles.dataR <- parm_recov.w %>% group_by(cond, resp) %>% nest() %>%
  mutate(quantiles=map(data, 
                       ~ as.data.frame(t(quantile(.x$q, probs=my.quantiles))))) %>%
  unnest(quantiles) %>%
  gather("quantile", "RT", '10%':'90%') %>%
  arrange(cond, resp)
par(mfrow=c(1,1))
tempX = list((100*x[2,1])/(x[1,1]+x[2,1]),(100*x[2,2])/(x[1,2]+x[2,2]),(100*x[2,3])/(x[1,3]+x[2,3]),(100*x[2,4])/(x[1,4]+x[2,4]),(100*x[2,5])/(x[1,5]+x[2,5]))
names(tempX)=c("s1","s2","s3","s4","s5")
(b=barplot(unlist(tempX),ylim=c(0,100),ylab="percent correct",xlab="stimulus",main="predicted accuracy"))
points(x=b,unlist(tempX),col="red")

library(tidyr)

long = gather(pars,condition,value,factor_key = T)[1:1080,]
long$sub = rep(1:72,15)
long$dist = rep(rep(c(0:3,2),each=72),3)
long$newold = rep(c(rep("new",288),rep("old",72)),3)
long$par = rep(c("thr","non","drift"),each=360)
library(lme4)
mThr = bayesglm(value~dist+newold+(1|sub),family="gaussian",long[long$par=="thr",])
summary(mThr)
mNon = bayesglm(value~dist+newold+(1|sub),family="gaussian",long[long$par=="non",])
summary(mNon)
mDrift = bayesglm(abs(value)~dist+newold+(1|sub),family="gaussian",long[long$par=="drift",])
summary(mDrift)
##############
library(DMCfun)
setwd("~/Downloads/DMC_190819/")
source("dmc/dmc.R")
load_model("DDM","ddm.R")
factors = list(S=c("s1","s2","s3","s4","s5"))
responses = c("upper","lower")
match.map = list(M=list(s1="upper",s2="upper",s3="upper",s4="upper",s5="lower"))
p.map.ddm15 = list(a="S",v="S",z="1",d="S",sz="1",sv="1",t0="1",st0="1")                 
constants.ddm15 = c(sv=1,sz=0.2,st0=0,z=0)
df$resp = resp
df[df$dist_newold=="new0","dist_newold"]="s1"
df[df$dist_newold=="new1","dist_newold"]="s2"
df[df$dist_newold=="new2","dist_newold"]="s3"
df[df$dist_newold=="new3","dist_newold"]="s4"
df[df$dist_newold=="old2","dist_newold"]="s5"
colnames(df) = c("id","R","RT","dist","S","resp")
B.model.ddm15 = model.dmc(p.map.ddm15,match.map=match.map,constants=constants.ddm15,
                          factors=factors,responses=responses,type="rd")
data1.Bmodel.ddm15 <- data.model.dmc(df, B.model.ddm15)
######################################
ddm3 <- function (x, datlist) {
  loss=0
  for (myVar in 1:length(datlist))
    loss = loss - logLik(wdm(datlist[[myVar]], 
                             alpha=x[myVar], tau=x[myVar+5], beta=0.5, delta=x[myVar+10]))
  return(loss)
}
ddm_thr.dri =function (x, datlist) {
  loss=0
  for (myVar in 1:length(datlist))
    loss = loss - logLik(wdm(datlist[[myVar]], 
                             alpha=x[myVar], tau=x[11], beta=0.5, delta=x[myVar+5]))
  return(loss)
}
ddm_thr.non =function (x, datlist) {
  loss=0
  for (myVar in 1:length(datlist))
    loss = loss - logLik(wdm(datlist[[myVar]], 
                             alpha=x[myVar], tau=x[myVar+5], beta=0.5, delta=x[11]))
  return(loss)
}
ddm_dri.non =function (x, datlist) {
  loss=0
  for (myVar in 1:length(datlist))
    loss = loss - logLik(wdm(datlist[[myVar]], 
                             alpha=x[11], tau=x[myVar+5], beta=0.5, delta=x[myVar]))
  return(loss)
}
ddm_dri =function (x, datlist) {
  loss=0
  for (myVar in 1:length(datlist))
    loss = loss - logLik(wdm(datlist[[myVar]], 
                             alpha=x[6], tau=x[7], beta=0.5, delta=x[myVar]))
  return(loss)
}
ddm_thr =function (x, datlist) {
  loss=0
  for (myVar in 1:length(datlist))
    loss = loss - logLik(wdm(datlist[[myVar]], 
                             alpha=x[myVar], tau=x[7], beta=0.5, delta=x[6]))
  return(loss)
}
ddm_non =function (x, datlist) {
  loss=0
  for (myVar in 1:length(datlist))
    loss = loss - logLik(wdm(datlist[[myVar]], 
                             alpha=x[7], tau=x[myVar], beta=0.5, delta=x[6]))
  return(loss)
}
p_ddm3 = rep(c(1,0.1,0),each=5)
p_ddm_thr.dri = c(rep(c(1,0),each=5),0.1)
p_ddm_thr.non = c(rep(c(1,0.1),each=5),0)
p_ddm_dri.non = c(rep(c(0,0.1),each=5),1)
p_ddm_dri = c(rep(0,each=5),1,0.1)
p_ddm_thr = c(rep(1,each=5),0,0.1)
p_ddm_non = c(rep(0.1,each=5),0,1)

data1_ddm3 = list()
data1_thr.dri = list()
data1_thr.non = list()
data1_dri.non = list()
data1_dri = list()
data1_thr = list()
data1_non = list()

subs = unique(tpc$subj_idx)
for(sub in 1:length(subs)){
  df = tpc[tpc$subj_idx==subs[sub],]
  q = df$rt
  resp = ifelse(df$response==1,"lower","upper")
  cond = df$dist_newold
  data1.w = data.frame(q,resp,cond)
  table(data1.w$cond,data1.w$resp)
  data1.w$resp = factor(data1.w$resp)
  data1.w$resp = relevel(data1.w$resp,ref = "upper")
  data1.w = as.wiener(data1.w)
  data1.s1 = subset(data1.w,cond=="new0")
  data1.s2 = subset(data1.w,cond=="new1")
  data1.s3 = subset(data1.w,cond=="new2")
  data1.s4 = subset(data1.w,cond=="new3")
  data1.s5 = subset(data1.w,cond=="old2")
  data1.both = list(data1.s1,data1.s2,data1.s3,data1.s4,data1.s5)
  
  data1_ddm3[[sub]] = optim(par=p_ddm3,f=ddm3,dat=data1.both)
  data1_thr.dri[[sub]] = optim(par=p_ddm_thr.dri,f=ddm_thr.dri,dat=data1.both)
  data1_thr.non[[sub]] = optim(par=p_ddm_thr.non,f=ddm_thr.non,dat=data1.both)
  data1_dri.non[[sub]] = optim(par=p_ddm_dri.non,f=ddm_dri.non,dat=data1.both)
  data1_dri[[sub]] = optim(par=p_ddm_dri,f=ddm_dri,dat=data1.both)
  data1_thr[[sub]] = optim(par=p_ddm_thr,f=ddm_thr,dat=data1.both)
  data1_non[[sub]] = optim(par=p_ddm_non,f=ddm_non,dat=data1.both)
}
ntrials = c()
for(sub in 1:length(subs)){
  df = tpc[tpc$subj_idx==subs[sub],]
  ntrials[sub] = nrow(df)
}
bicmat = aicmat = matrix(NA,72,7,byrow=F)
for(i in 1:72){
  bicmat[i,1] = 2*data1_ddm3[[i]]$value+log(ntrials[i])*15
  bicmat[i,2] = 2*data1_thr.dri[[i]]$value+log(ntrials[i])*11
  bicmat[i,3] = 2*data1_thr.non[[i]]$value+log(ntrials[i])*11
  bicmat[i,4] = 2*data1_dri.non[[i]]$value+log(ntrials[i])*11
  bicmat[i,5] = 2*data1_dri[[i]]$value+log(ntrials[i])*7
  bicmat[i,6] = 2*data1_thr[[i]]$value+log(ntrials[i])*7
  bicmat[i,7] = 2*data1_non[[i]]$value+log(ntrials[i])*7
  
  aicmat[i,1] = 2*data1_ddm3[[i]]$value+30
  aicmat[i,2] = 2*data1_thr.dri[[i]]$value+22
  aicmat[i,3] = 2*data1_thr.non[[i]]$value+22
  aicmat[i,4] = 2*data1_dri.non[[i]]$value+22
  aicmat[i,5] = 2*data1_dri[[i]]$value+14
  aicmat[i,6] = 2*data1_thr[[i]]$value+14
  aicmat[i,7] = 2*data1_non[[i]]$value+14
}
library(bmsR)
library(qpcR)
library(viridis)
library(ggbeeswarm)
m = akaike.weights(bicmat[cat$strategy==2,])$weights
m = matrix(m,sum(cat$strategy==2),7,byrow=F)
pxp = VB_bms(log(m))
barplot(pxp$pxp)
aic = data.frame(AIC = c(aicmat[,1],aicmat[,2],aicmat[,3],aicmat[,4],aicmat[,5],aicmat[,6],aicmat[,7]),
                 subs = rep(1:72,7),
                 model = rep(c("DTN","DT","TN","DN","D","T","N"),each=72))
bic = data.frame(BIC = c(bicmat[,1],bicmat[,2],bicmat[,3],bicmat[,4],bicmat[,5],bicmat[,6],bicmat[,7]),
                 subs = rep(1:72,7),
                 model = rep(c("DTN","DT","TN","DN","D","T","N"),each=72))

ggplot(aic,aes(x=model,y=AIC,fill=model))+
  geom_boxplot()+
  geom_beeswarm(aes(fill=model,colour=model),color="black",shape=21,dodge.width = 1.1,
                priority = "density")+
  scale_fill_viridis(discrete=T)+
  scale_y_continuous(limits = c(),expand= c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(bic,aes(x=model,y=BIC,fill=model))+
  geom_boxplot()+
  geom_beeswarm(aes(fill=model,colour=model),color="black",shape=21,dodge.width = 1.1,
                priority = "density")+
  scale_fill_viridis(discrete=T)+
  scale_y_continuous(limits = c(),expand= c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


cat = read.csv("~/Dropbox (University of Oregon)/Troy&Dasa/TPC/version2022/cat_models_20221204.csv")
drift_mod_thr = c()
drift_mod_non = c()
for(i in 1:72){
  drift_mod_thr[i] = data1_dri[[i]]$par[6]
  drift_mod_non[i] = data1_dri[[i]]$par[7]
}
t.test(drift_mod_thr[cat$strategy==1],drift_mod_thr[cat$strategy==2])
t.test(drift_mod_non[cat$strategy==1],drift_mod_non[cat$strategy==2])
