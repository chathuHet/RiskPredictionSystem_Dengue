threshold1 <- function(predict, response) {
  perf <- ROCR::performance(ROCR::prediction(predict, response), "sens", "spec")
  df <- data.frame(cut = perf@alpha.values[[1]], sens = perf@x.values[[1]], spec = perf@y.values[[1]])
  df[which.max(df$sens + df$spec), "cut"]
}



name_replace <- function(A)
{
  A <- gsub("Changhua County", "Changhua", A)
  A <- gsub("Hualien County", "Hualien", A)
  A <- gsub("Kaohsiung City","Kaohsiung", A)
  A <- gsub("Keelung City","Keelung", A)
  A <- gsub("Miaoli County","Miaoli", A)
  A <- gsub("Nantou County","Nantou", A)
  A <- gsub("Penghu County","Penghu", A)
  A <- gsub("Pingtung County", "Pingtung", A)
  A <- gsub("Taichung City", "Taichung", A)
  A <- gsub("Tainan City", "Tainan" , A)
  A <- gsub("Taipei City", "Taipei", A)
  A <- gsub("New Taipei","New Taipei City",  A)
  A <- gsub("Taitung County" ,"Taitung", A)
  A <- gsub("Taoyuan City", "Taoyuan", A)
  A <- gsub("Yilan County","Yilan" , A)
  A <- gsub("Yunlin County","Yulin", A)
  A
}



boot_samples <- function(i,df,a,prop=0.8,seed1=2)
{
  set.seed(seed1+21)
  df_c <- subset(df,county==a[i])
  trainRows = sample(1:nrow(df_c),prop*nrow(df_c),replace = TRUE)
  train = df_c[trainRows,]
  train
}


pos_boot <- function(seed2,df,sampling_r)
{
  a <- unique(df[,"county"])
  train_tr <- do.call(rbind,lapply(1:length(a), sampling_r,df=train,a=a,prop=1,seed1=seed2)) 
  train_tr <- train_tr[order(train_tr[,"date"]),]
  
  dfa <- data.frame("larva"=train_tr[,"larva"], 
                    MaxTemp=train_tr[,"MaxTemp"], 
                    MaxPrecip=train_tr[,"MaxPrecip"],
                    AvgRHumi=train_tr[,"AvgRHumi"],
                    AvgHGT=train_tr[,"AvgHGT"],
                    Con=train_tr[,"Con"], 
                    county=train_tr[,"county"],
                    date=train_tr[,"date"])
  dfa <- dfa[!is.na(dfa[,"Con"])&dfa[,"Con"]>-Inf,]
  df2 <- cbind.data.frame(larvaB=ifelse(dfa[,1]<1,0,1),dfa[,2:(ncol(dfa)-2)])
  reg3 = step(glm(larvaB~MaxTemp+MaxPrecip+AvgRHumi+AvgHGT+offset(Con),family=binomial(link='logit'),data=df2), scope=list(upper = ~ .+.^2,lower=~.), trace=FALSE)
  summary(reg3)
  t2 = threshold1(reg3$fitted.values, df2[,1])
  threshold = round(t2,3)  
  
  B <- cbind.data.frame(larvaB=df2[,"larvaB"],
                        pred = ifelse(reg3$fitted.values>t2,"p_Yes","p_No"),
                        larva=dfa[,"larva"],
                        county=dfa[,"county"],
                        dfa[,2:6],
                        predp = reg3$fitted.values)
  misClasificError_train <- mean(ifelse(reg3$fitted.values>t2,1,0) != B[,"larvaB"])
  acc_train = round(1-misClasificError_train,3)
  sn_train = round(sum(B[,"larvaB"]==1& B[,"pred"]=="p_Yes")/sum(B[,"larvaB"]==1),3)
  sp_train = round(1- sum(B[,"larvaB"]==0& B[,"pred"]=="p_Yes")/sum(B[,"larvaB"]==0),3)
  bb_train <- as.matrix(table(B[,"county"],B[,"pred"]))
  bb2_train <- cbind.data.frame(rownames(bb_train),do.call(rbind,lapply(1:nrow(bb_train),function(i) round(bb_train[i,]/sum(bb_train[i,]),3) )))
  
  out <- cbind.data.frame(threshold,acc=acc_train,sn=sn_train,sp=sp_train,p_Y= bb2_train[,"p_Yes"],county=bb2_train[,1])
  out
}


####  FUNCTION FOR PLOTTING THE POSITIVE PERCENTAGES OF COUNTIES (FIG 3)  #####
map_pos <- function(df,brk,plt='PuRd')
{
  gadm <- readRDS("GADM_2.8_TWN_adm2.rds")
  taiwan.adm2.spdf <- get("gadm")
  taiwan.adm2.df <- fortify(taiwan.adm2.spdf, region = "NAME_2")
  
  pos.df  <- cbind.data.frame(id=unique(taiwan.adm2.df[,'id']),pos=NA)
  for(k in 1:nrow(pos.df ))
  {if(pos.df [k,'id']%in%df[,1]) pos.df [k,"pos"]= round(df[which(as.character(pos.df [k,'id'])==as.character(df[,1])),"p_Yes"]*100,0)}
  
  taiwan.adm2.df <- merge(taiwan.adm2.df, pos.df, by.y = 'id', all.x = TRUE)
  taiwan.adm2.df <- taiwan.adm2.df[findInterval(taiwan.adm2.df$long,c(119, 123))==1 & findInterval(taiwan.adm2.df$lat,c(21.7, 25.5))==1,]
  taiwan.adm2.df <- taiwan.adm2.df[!is.na(taiwan.adm2.df[,"pos"]),]
  taiwan.adm2.centroids.df <- data.frame(long = coordinates(taiwan.adm2.spdf)[, 1], 
                                         lat = coordinates(taiwan.adm2.spdf)[, 2]) 
  
  # Get names and id numbers corresponding to administrative areas
  taiwan.adm2.centroids.df[, 'ID_2'] <- taiwan.adm2.spdf@data[,'ID_2']          
  taiwan.adm2.centroids.df[, 'NAME_2'] <- taiwan.adm2.spdf@data[,'NAME_2']
  labtext <-  taiwan.adm2.centroids.df[-which(taiwan.adm2.centroids.df$NAME_2%in%c("Kinmen", "Lienkiang (Matsu Islands)", "Penghu")),]
  
  ggplot(taiwan.adm2.df, aes(x = long, y = lat, group = group)) + 
    geom_polygon(aes(fill = cut(pos,breaks = brk, include.lowest = TRUE))) + 
    geom_path(colour = 'gray', linetype = "solid")+ #"dashed") +
    #geom_text(data = taiwan.adm2.centroids.df, aes(label = NAME_2, x = long, y = lat, group = NAME_2), size = 2) + 
    geom_text(data = labtext, aes(label = NAME_2, x = long, y = lat, group = NAME_2), size = 3) + 
    labs(x=" ", y=" ") + 
    theme_bw() + 
    scale_fill_brewer('% of postive', palette  = plt) + 
    coord_map() + 
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
          panel.border = element_blank(),#legend.position = "none")
          legend.text=element_text(size=12),legend.title=element_text(size=12),
          legend.key.width = unit(0.5, "cm"),legend.key.height = unit(1, "cm"),
          legend.position = c(1.3,0.4))
}

##  function for calculating prediction interavals  ##
pz <- function(object, newdata, type = c("response", "prob"),
               se=FALSE,MC=1000,level=.95,
               na.action = na.pass, ...)
{
  type <- match.arg(type)
  
  ## if no new data supplied
  if(missing(newdata))
  {
    rval <- object$fitted.values
    if(!is.null(object$x)) {
      X <- object$x$count
      Z <- object$x$zero
    }
    else if(!is.null(object$model)) {
      X <- as.matrix(cbind(model.matrix(object$terms$count, object$model, contrasts= object$contrasts$count),Con=df[,"Con"]))
      Z <- as.matrix(cbind(model.matrix(object$terms$zero,  object$model, contrasts = object$contrasts$zero),Con=df[,"Con"]))
      mu <- exp(X %*% c(object$coefficients$count,1))[,1]
      phi <- object$linkinv(Z %*% c(object$coefficients$zero,1))[,1]
    }
    else {
      stop("no X and/or Z matrices can be extracted from fitted model")
    }
    if(type == "prob"){
      mu <- exp(X %*% object$coefficients$count)[,1]
      phi <- object$linkinv(Z %*% object$coefficients$zero)[,1]
    }
  }
  else {
    mf <- model.frame(delete.response(object$terms$full), newdata, xlev = object$levels)
    X <- cbind(model.matrix(delete.response(object$terms$count), mf, contrasts = object$contrasts$count),Con=newdata[,"Con"])
    Z <- cbind(model.matrix(delete.response(object$terms$zero),  mf, contrasts = object$contrasts$zero),Con=newdata[,"Con"])
    mu <- exp(X %*% c(object$coefficients$count,1))[,1]
    phi <- object$linkinv(Z %*% c(object$coefficients$zero,1))[,1]
    rval <- (1-phi) * mu
  }   
  
  if(se & !is.null(X) & !is.null(Z))
  {
    require(mvtnorm)
    require(ZIM)
    vc <- -solve(object$optim$hessian)
    kx <- length(object$coefficients$count)
    kz <- length(object$coefficients$zero)
    parms <- object$optim$par
    
    if(type!="prob")
    {
      yhat.sim <- matrix(NA,MC,dim(X)[1])
      zpr.sim <- matrix(NA,MC,dim(X)[1])
      nonz.sim <- matrix(NA,MC,dim(X)[1])
      y.sim <- matrix(NA,MC,dim(X)[1])
      
      for(i in 1:MC)
      {
        set.seed(22+i)
        parms.sim <- rmvnorm(n=1,mean=parms,sigma=vc)
        beta <- c(parms.sim[1:kx],1)
        gamma <- c(parms.sim[(kx+1):(kx+kz)],1)
        mu.sim <- exp(X%*%beta)[,1]
        phi.sim <- object$linkinv(Z%*%gamma)[,1]
        yhat.sim[i,] <- (1-phi.sim)*mu.sim
        if(object$dist== "negbin") y.sim[i,] <- rzinb(length(mu),k=object$theta,lambda=mu,omega=phi)
        if(object$dist=="poisson") y.sim[i,] <- rzip(length(mu),lambda=mu,omega=phi)
      }
    }
    out <- list()
    out$lower.yhat <- apply(yhat.sim,2,quantile,(1-level)/2)
    out$upper.yhat  <- apply(yhat.sim,2,quantile,1-((1-level)/2))
    out$se <- apply(yhat.sim,2,sd)
    out$lower.y <- apply(y.sim,2,quantile,(1-level)/2)
    out$upper.y  <- apply(y.sim,2,quantile,1-((1-level)/2))
  }
  if(se) rval <- list(fitted=rval,out=out,phi=phi,mu=mu)
  rval
}


### function for accumulating the predcited values by week  ###]
divisional <- function(A_division,division,firstday,varnames)
{
  weekly <- function(j,division){
    c(j,sum(as.Date(A_division[,"date"])%in%seq(firstday[j],firstday[j]+6,1)),
      as.numeric(colSums(A_division[as.Date(A_division[,"date"])%in%seq(firstday[j],firstday[j]+6,1),varnames])),
      as.character(division))
  }
  #a function to compoute the weekly sums for the interested fields
  jj=matrix(1:length(firstday),length(firstday),1)
  weeklysums <- t(apply(jj,1,weekly,division=division))  ## calls the function weekly for all of the weeks
  division_weekly <- data.frame(
    as.numeric(weeklysums[,1]),as.Date(firstday[as.numeric(weeklysums[,1])]),
    apply(weeklysums[,2:(ncol(weeklysums)-1)], 2, function(x) as.numeric(x)),
    weeklysums[,ncol(weeklysums)])
  colnames(division_weekly)<-c("Week number","First day of week","Number of samples",paste0("sum_",varnames),"county")
  I_missing=max(as.Date(A_division[,"date"])) < division_weekly[,"First day of week"]  #the weeks where the data is not available for
  division_weekly[I_missing,3:ncol(division_weekly)] <- NA
  division_weekly
}

