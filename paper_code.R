library(here)
library(pROC)
library(ggplot2)
library(pscl)
library(lmtest)
library(MASS)

source(here::here("functions.R"))
dfall <- readRDS(here::here("mosdf.rds"))

##### Split the datset into train and test datasets  #####
exactdate <- as.Date(dfall[nrow(dfall)*0.75,"date"])
weekstart <- seq(as.Date(dfall[1,"date"]), as.Date(dfall[nrow(dfall),"date"]), 7)
splitdate <- weekstart[which.min(abs(weekstart-exactdate))]
train <- dfall[as.Date(dfall[,"date"])<as.Date(splitdate),]
test <- dfall[as.Date(dfall[,"date"])>=as.Date(splitdate),]

#############################################
#######   S T A G E     O N E   #############
#############################################

##### Fit a regression model for the train data set  #####
dfa <- data.frame("larva"=train[,"larva"], 
                  MaxTemp=train[,"MaxTemp"], 
                  MaxPrecip=train[,"MaxPrecip"],
                  AvgRHumi=train[,"AvgRHumi"],
                  AvgHGT=train[,"AvgHGT"],
                  Con=train[,"Con"], 
                  county=train[,"county"],
                  date=train[,"date"],
                  lat=train[,"lat"],
                  lon=train[,"lon"])
dfa <- dfa[!is.na(dfa[,"Con"])&dfa[,"Con"]>-Inf,]
df2 <- cbind.data.frame(larvaB=ifelse(dfa[,1]<1,0,1),dfa[,2:(ncol(dfa)-4)])
reg3 = step(glm(larvaB~MaxTemp+MaxPrecip+AvgRHumi+AvgHGT+offset(Con),family=binomial(link='logit'),data=df2), 
            scope=list(upper = ~ .+.^2,lower=~.), trace=FALSE)
summary(reg3)
t2 = threshold1(reg3$fitted.values, df2[,1])
print(summary(reg3))
print(paste('threshold',round(t2,3)))


##########################################################
######## accuracy measures for the train set   ###########
B <- cbind.data.frame(larvaB=df2[,"larvaB"],
                      pred = ifelse(reg3$fitted.values>t2,"p_Yes","p_No"),
                      larva=dfa[,"larva"],
                      county=dfa[,"county"],
                      predp = reg3$fitted.values,
                      lat=dfa[,"lat"],lon=dfa[,"lon"])
print(paste('Accuracy:',round(1-mean(ifelse(reg3$fitted.values>t2,1,0) != B[,"larvaB"]),2)))
print(paste('Sensitivity(=tpr):',round(sum(B[,"larvaB"]==1& B[,"pred"]=="p_Yes")/sum(B[,"larvaB"]==1),2)))
print(paste('specificity(=1-fpr):',1-round(sum(B[,"larvaB"]==0& B[,"pred"]=="p_Yes")/sum(B[,"larvaB"]==0),2)))
print(paste('AUC:',round(pROC::auc(B$larvaB, ifelse(B$pred=="p_Yes",1,0)),3)))


## predicted mosquito percentages for the train period by county ##
bb <- as.matrix(table(B[,"county"],B[,"pred"]))
bb2 <- cbind.data.frame(rownames(bb),do.call(rbind,lapply(1:nrow(bb),function(i) round(bb[i,]/sum(bb[i,]),3) )))
bb2[order(bb2[,"p_Yes"],decreasing=TRUE),]
#replace the county names required for plotting
dfpos <- cbind.data.frame(name_replace(bb2[,1]),bb2[,2:3]) 

##########################################################
######## accuracy measures for the test set   ###########
dfa_test <- data.frame("larva"=test[,"larva"], 
                       MaxTemp=test[,"MaxTemp"], 
                       MaxPrecip=test[,"MaxPrecip"],
                       AvgRHumi=test[,"AvgRHumi"],AvgHGT=test[,"AvgHGT"],
                       Con=test[,"Con"], 
                       county=test[,"county"],
                       date=test[,"date"])
dfa_test <- dfa_test[!is.na(dfa_test[,"Con"])&dfa_test[,"Con"]>-Inf,]
df2_test <- cbind.data.frame(larvaB=ifelse(dfa_test[,1]<1,0,1),dfa_test[,2:(ncol(dfa_test)-2)])
predp = predict(reg3, df2_test[,-1], type="response")
B_test <- cbind.data.frame(larvaB=df2_test[,"larvaB"],pred = ifelse(predp>t2,"p_Yes","p_No"),
                           larva=dfa_test[,"larva"],county=dfa_test[,"county"], predp = predp)
print(paste('Accuracy:',round(1-mean(ifelse(B_test[,"predp"]>t2,1,0) != B_test[,"larvaB"]),3)))
print(paste('sensitivity(=tpr):',round(sum(B_test[,"larvaB"]==1& B_test[,"pred"]=="p_Yes")/sum(B_test[,"larvaB"]==1),3)))
print(paste('specificity(=1-fpr):',1-round(sum(B_test[,"larvaB"]==0& B_test[,"pred"]=="p_Yes")/sum(B_test[,"larvaB"]==0),3)))
print(paste('AUC:',round(pROC::auc(B_test$larvaB, ifelse(B_test$pred=="p_Yes",1,0)),3)))


## predicted mosquito percentages for the train period by county ##
bb_test <- as.matrix(table(B_test[,"county"],B_test[,"pred"]))
bb2_test <- cbind.data.frame(rownames(bb_test),do.call(rbind,lapply(1:nrow(bb_test),function(i) round(bb_test[i,]/sum(bb_test[i,]),3))))
#replace the county names required for plotting
dfpos_test <- cbind.data.frame(name_replace(bb2_test[,1]),bb2_test[,2:3])

##########
### Wilcoxon test to compare the rank order of predictions in train and test sets
wilcox.test(bb2$p_Yes, bb2_test$p_Yes, paired=FALSE) 

## observed mosquito percentages by county ##
bb_obs <- as.matrix(table(dfall[,"county"],ifelse(dfall[,"larva"]>0,"p_Yes","p_No")))
bb2_obs <- cbind.data.frame(rownames(bb_obs),do.call(rbind,lapply(1:nrow(bb_obs),function(i) round(bb_obs[i,]/sum(bb_obs[i,]),3) )))
#replace the county names required for plotting
dfpos_obs <- cbind.data.frame(name_replace(bb2_obs[,1]),bb2_obs[,2:3])



############################################################
########  THRESHOLD FOR CLASSITYING COUNTIES ###############
### first-time only
## seed_seq <- 1:1000
## out <- do.call(rbind,lapply(1:length(seed_seq), pos_boot,df=train,sampling_r=boot_samples))
## saveRDS(out,file="pos_percentage.rds")

pos_train <- readRDS(here::here("pos_percentage.rds"))
county <-  unique(pos_train[,"county"])

pos_coun <- bb2_obs[bb2_obs[,3]>0.01,1]
pos_train$response <- ifelse(pos_train$county%in%pos_coun==TRUE,1,0)

t3 <- threshold1(pos_train$p_Y,pos_train$response)
print(paste('threshold to classify counties: ',round(t3,3)))


######## accuracy measures for classifying the counties   ###########
cc <- cbind.data.frame(res=pos_train$response,pred = ifelse(pos_train$p_Y>t3,"p_Yes","p_No"))
print(paste('Threshold level:', round(t3,3)))
print(paste('Accuracy:', round(1-mean(ifelse(pos_train$p_Y>t3,1,0) != cc[,1]) ,3)))
print(paste('sensitivity(=tpr):', round(sum(cc[,1]==1& cc[,2]=="p_Yes")/sum(cc[,1]==1),3)))
print(paste('specificity(=1-fpr):',1-round(sum(cc[,1]==0& cc[,2]=="p_Yes")/sum(cc[,1]==0),3)))

### boxplot of larvae positive locations  ###
fig4 <- ggplot(pos_train, aes(reorder(county, p_Y, FUN=median),p_Y*100)) + geom_boxplot() + coord_flip() +
  xlab("Administrative division")  +ylab("Percentage of larvae positive locations") +  
  geom_hline(aes(yintercept=t3*100), colour="#990000",size=1) + 
  annotate(geom="text", x=6, y=37,label=paste0("threshold level = ",t3*100,"%"),size=4)+
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
fig4
#ggsave("../Figs/Fig4.tiff",fig4)


##################################################################
####  PLOTTING THE POSITIVE PERCENTAGES OF COUNTIES (FIG 3)  #####

# detrmine the breakpoints
if(abs(t3*100-round(t3*100,-1))<5)
  {brksq <- c(0,round(t3*100,-1)/2, t3*100,seq(round(t3*100,-1)+10,100,by=10))}else
  {brksq <- c(0,round(t3*100,-1)/2, t3*100,seq(round(t3*100,-1)+10,100,by=10))}

fig3a <- map_pos(dfpos,brk=brksq,plt='OrRd')
fig3a 
#ggsave("../Figs/Fig3a.tiff",fig3a)

fig3b <- map_pos(dfpos_test,brk=brksq,plt='OrRd')
fig3b
#ggsave("../Figs/Fig3b.tiff",fig3b)

fig3c <- map_pos(dfpos_obs,brk=c(0,10,20,30,40,50),plt='YlGn')
fig3c
#ggsave("../Figs/Fig3c.tiff",fig3c)





#############################################
#######   S T A G E     T W O   #############
#############################################

### ### ### 
### T R A I N    S E T ###

## filter the datset to contain only larvae positive counties - train set
ispos <- bb2[bb2[,"p_Yes"]>t3,1]
cou <- function(i) grep(paste("^",ispos[i],"$", sep=""),dfa[,"county"])
index <- do.call(c,apply(matrix(1:length(ispos),length(ispos),1) ,1 , cou))
dfb <- dfa[index,]
df <- dfb[,-which(colnames(dfb) %in% c("county","date"))]

### Fit a zero-inflated regression model for the train data set 
m.znb <- step(zeroinfl(data=df,larva~MaxTemp+MaxPrecip+AvgRHumi+AvgHGT+offset(Con),
                       dist = "negbin", EM = TRUE) , scope=list(upper = ~ .,lower=~MaxTemp+AvgRHumi), trace=FALSE)
print(summary(m.znb))


### Prediction intervals for the train set 
znbint <- pz(m.znb,se=TRUE,type="response")
indx2 <- order(znbint$out$upper.y)   #order by upper value
n <- length(indx2)
zpltnb <- cbind.data.frame(x=1:n,obs= df[,1][indx2], fits=znbint$fitted[indx2],
                           mu0=znbint$out$lower.yhat[indx2],mu1=znbint$out$upper.yhat[indx2],
                           y0_3=znbint$out$lower.y[indx2],y1_3=znbint$out$upper.y[indx2])
A <- cbind.data.frame(x=zpltnb[,"x"], obs=zpltnb[,"obs"], znb_fit=round(zpltnb[,"fits"],3), 
                      znb0=zpltnb[,"y0_3"],znb1=zpltnb[,"y1_3"])
fig7a <- ggplot(data=A, ylim=range(cbind("znb0","znb1")), aes(x=A[,"x"], y=A[,"obs"]))+xlab("index")+
  ylab(expression(paste("Number of ", italic("Ae. aegypti"), " larvae")))+
  geom_linerange(aes(ymin =A[,"znb0"], ymax = A[,"znb1"]),col="green")+
  geom_point(cex=0.2,col="blue",alpha=0.2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
fig7a
#ggsave("../Figs/Fig7a.tiff",fig7a)

### coverage train set
b=0
for(i in 1:nrow(A))
{
  y=A[i,"obs"]
  if((y >= as.numeric(A[i,"znb0"])) & (y <= as.numeric(A[i,"znb1"]))) b=b+1;  
}
print(paste("coverage - zero inflated neg bin model (train set): ", round(b/nrow(zpltnb)*100,2))) 

### MSE train set
print(paste("MSE (train set): ",round(sum((zpltnb[,"obs"]-zpltnb[,"fits"])^2)/nrow(dfa),1)))

### ### ### 
### T E S T   S E T ###
### filter the datset to contain only larvae positive counties - test set
cou_test <- function(i) grep(paste("^",ispos[i],"$", sep=""),dfa_test[,"county"])
index_test <- do.call(c,apply(matrix(1:length(ispos),length(ispos),1) ,1 , cou_test))
dfb_test <- dfa_test[index_test,]
df_test <- dfb_test[,-which(colnames(dfb_test) %in% c("county","date"))]

###  Prediction intervals for the test set
znbint_test <- pz(m.znb,newdata=df_test,se=TRUE,type="response")
indx2_test <- order(znbint_test$out$upper.y)  #order by upper value
n <- length(indx2_test)
zpltnb_test <- cbind.data.frame(x=1:n,obs= df_test[,1][indx2_test], fits=znbint_test$fitted[indx2_test],
                                mu0=znbint_test$out$lower.yhat[indx2_test],mu1=znbint_test$out$upper.yhat[indx2_test],
                                y0_3=znbint_test$out$lower.y[indx2_test],y1_3=znbint_test$out$upper.y[indx2_test])
A_test <- cbind.data.frame(x=zpltnb_test[,"x"], obs=zpltnb_test[,"obs"], znb_fit=round(zpltnb_test[,"fits"],3), 
                           znb0=zpltnb_test[,"y0_3"],znb1=zpltnb_test[,"y1_3"])
fig7b <- ggplot(data=A_test, ylim=range(cbind("znb0","znb1")), aes(x=A_test[,"x"], y=A_test[,"obs"]))+xlab("index")+
  ylab(expression(paste("Number of ", italic("Ae. aegypti"), " larvae")))+
  geom_linerange(aes(ymin =A_test[,"znb0"], ymax = A_test[,"znb1"]),col="green")+
  geom_point(cex=0.2,col="blue",alpha=0.2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
fig7b
#ggsave("../Figs/Fig7b.tiff",fig7b)

### coverage test set
b=0
for(i in 1:nrow(A_test))
{
  y=A_test[i,"obs"]
  #if((y >= as.numeric(A_test[i,"nb0"])) & (y <= as.numeric(A_test[i,"nb1"]))) a=a+1;  
  if((y >= as.numeric(A_test[i,"znb0"])) & (y <= as.numeric(A_test[i,"znb1"]))) b=b+1;  
}
#print(paste("coverage - neg bin model : ", round(a/nrow(pltnb)*100,2))) 
print(paste("coverage - zero inflated neg bin model (test set): ", round(b/nrow(zpltnb_test)*100,2))) 

### MSE test set
print(paste("MSE (train set): ",round(sum((zpltnb_test[,"obs"]-zpltnb_test[,"fits"])^2)/nrow(dfa_test),1)))

###############################
#### weekly predictions 
# train set
df_c <- cbind.data.frame(dfb, larva_st=dfb[,"larva"], fit_znb=m.znb$fitted.values)
c_name="all"
firstday <- seq(as.Date(df_c[1,"date"]), as.Date(df_c[nrow(df_c),"date"]), 7)
firstday <- firstday[-length(firstday)]
division_weekly  <- divisional(df_c,c_name,firstday,c("larva_st","fit_znb"))

# test set
df_c_test <- cbind.data.frame(dfb_test, larva_st = dfb_test[,"larva"], fit_znb = znbint_test$fitted)
firstday_test <- seq(as.Date(df_c_test[1,"date"]), as.Date(df_c_test[nrow(df_c_test),"date"]), 7)
firstday_test <- firstday_test[-length(firstday_test)]
division_weekly_test  <- divisional(df_c_test,c_name,firstday_test,c("larva_st","fit_znb"))
division_weekly_test$`Week number` <- (nrow(division_weekly)+1):(nrow(division_weekly)+nrow(division_weekly_test))

df_tt <- rbind.data.frame(cbind.data.frame(week=division_weekly$`Week number`,num=division_weekly$sum_larva_st,set="Observed for train set"),
                          cbind.data.frame(week=division_weekly$`Week number`,num=division_weekly$sum_fit_znb,set="Predicted for train set"),
                          cbind.data.frame(week=division_weekly_test$`Week number`,num=division_weekly_test$sum_larva_st,set="Observed for test set"),
                          cbind.data.frame(week=division_weekly_test$`Week number`,num=division_weekly_test$sum_fit_znb,set="Predicted for test set"),
                          cbind.data.frame(week=division_weekly$`Week number`[nrow(division_weekly)],num=division_weekly$sum_larva_st[nrow(division_weekly)],set="Observed for test set"),
                          cbind.data.frame(week=division_weekly$`Week number`[nrow(division_weekly)],num=division_weekly$sum_fit_znb[nrow(division_weekly)],set="Predicted for test set"))

fig5 <- ggplot(df_tt, aes(week, num, colour = set, linetype=set)) +
  geom_line() + xlab("Week") + 
  ylab(expression(paste("Weekly number of ", italic("Ae. aegypti"), " larvae"))) +
  scale_colour_manual("",values=c("blue","green","blue","green")) +
  scale_linetype_manual("", values=c(1,1,2,2))+
  theme(legend.text=element_text(size=8),legend.title=element_blank(),
        legend.key.width = unit(0.5, "cm"),legend.key.height = unit(0.8, "cm"),
        legend.position = "bottom", 
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
fig5
#ggsave("../Figs/Fig5.tiff", fig5)



##############################
########   DICUSSION  ########

## COMPARISON WITH ONE-STAGE MODEL
df_z <- dfa[,-which(colnames(dfa) %in% c("county","date"))]

m.znb_z <- step(zeroinfl(data=df_z,larva~MaxTemp+MaxPrecip+AvgRHumi+AvgHGT+offset(Con),
                         dist = "negbin", EM = TRUE) , scope=list(upper = ~ .,lower=~MaxTemp+AvgRHumi), trace=FALSE)

##### MSE#
############################
sum((zpltnb_z[,"obs"]-zpltnb_z[,"fits"])^2)/nrow(zpltnb_z)
mean((df_z$larva - m.znb_z$fitted.values)^2)

##### The sn/sp for the one stage model ###
##############
D <- cbind.data.frame(dfa,fits=round(m.znb_z$fitted.values,3),
                      larvaB=ifelse(dfa[,"larva"]<1,0,1),
                      pred = ifelse(round(m.znb_z$fitted.values,1)<1,"p_No","p_Yes"))
#In determining yes or no the predictions are rounded to closest integer
#b <- bb2[bb2[,"p_Yes"]<0.1,1]
#couN <- function(i) grep(paste("^",b[i],"$", sep=""),D[,"county"])
#index <- do.call(c,apply(matrix(1:length(b),length(b),1) ,1 , couN))
#Db <- D[index,]  #D for counties we label NO
#(sum((zpltnb[,"obs"]-zpltnb[,"fits"])^2)+sum((aaa[,"larva"])^2))/nrow(dfa)
#print(table(D[,"larvaB"],D[,"pred"]))

print(paste('Accuracy',round(1-mean(ifelse(round(m.znb_z$fitted.values,1)>1,1,0) != D[,"larvaB"]),3)))
print(paste('tpr(=sensitivity)',round(sum(D[,"larvaB"]==1& D[,"pred"]=="p_Yes")/sum(D[,"larvaB"]==1),3)))
print(paste('fpr(=1-specificity)',round(sum(D[,"larvaB"]==0& D[,"pred"]=="p_Yes")/sum(D[,"larvaB"]==0),3)))



############################################
#######   OTHER SUPPLEMENTARY  #############

### dotplot of larve  ###
fig6 <- ggplot(dfall, aes(larva,col=larva)) +stat_ecdf(geom = "step")+
  xlab("Larvae count")+ylab("Cumulative density")+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
fig6
#ggsave("../Figs/fig6.tiff", fig6)

### model comparison ###
m.lm <- step(lm(data=df,larva~MaxTemp+MaxPrecip+AvgRHumi+AvgHGT+offset(Con)) , scope=list(upper = ~ .,lower=~MaxTemp+AvgRHumi), trace=FALSE)
m.p <- step(glm(data=df,larva~MaxTemp+MaxPrecip+AvgRHumi+AvgHGT+offset(Con),family=poisson), scope=list(upper = ~ .,lower=~MaxTemp+AvgRHumi), trace=FALSE)
m.nb <- step(glm.nb(data=df,larva~MaxTemp+MaxPrecip+AvgRHumi+AvgHGT+offset(Con)), scope=list(upper = ~ .,lower=~MaxTemp+AvgRHumi), trace=FALSE)
m.zp <- step(zeroinfl(data=df,larva~MaxTemp+MaxPrecip+AvgRHumi+AvgHGT+offset(Con), dist = "poisson", EM = TRUE), scope=list(upper = ~ .,lower=~MaxTemp+AvgRHumi), trace=FALSE)

print(lrtest(m.p, m.lm))
print(vuong(m.nb, m.p))
print(vuong(m.zp, m.p))
print(vuong(m.znb, m.nb))
print(vuong(m.znb, m.zp))