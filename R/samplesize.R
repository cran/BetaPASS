betapwr2.base <- function(mu0,sd0,mu1,sampsize,trials,seed,link.type,equal.precision,sd1,sig.level){
  #Set seed
  set.seed(seed)
  
  #Set parameters
  if(equal.precision == TRUE){
    phi<- ((mu0*(1-mu0))/(sd0*sd0))-1
    a0<- mu0*phi
    b0<- (1-mu0)*phi
    a1<- mu1*phi
    b1<- (1-mu1)*phi
  }
  else{
    if(is.null(sd1)==TRUE){
      stop("miss sd1 with equal dispersion parameter assumption")
    }
    else{
      phi <- ((mu0*(1-mu0))/(sd0*sd0))-1
      phi1 <- ((mu1*(1-mu1))/(sd1*sd1))-1
      if(phi < 0){
        stop("phi must be greater than 0")
      }
      if(phi1 < 0){
        stop("phi1 must be greater than 0")
      }
      a0<- mu0*phi
      b0<- (1-mu0)*phi
      a1<- mu1*phi1
      b1<- (1-mu1)*phi1
    }
  }
  
  #n is trial, trials is number of trials
  #Y.H0.1 is simulated Y under H0, version 1
  #Y.Ha.1 is simulated Y under Ha, version 1
  Y.H0.1 <- Y.Ha.1 <- matrix(NA,sampsize,trials)
  for (n in 1:trials) {
    for (s in 1:sampsize) {
      Y.H0.1[s,n] <- rbeta(1,a0,b0)
      Y.Ha.1[s,n] <- rbeta(1,a1,b1)
    }
  }
  
  #Reshape simulated Y
  #Switch y value in matrix, into a 'y' column in data.frame
  Y.H0.2 <- reshape::melt(Y.H0.1)  
  Y.Ha.2 <- reshape::melt(Y.Ha.1)
  
  #Combine Y.H0 and Y.H1
  Y.mat <- rbind(Y.H0.2,Y.Ha.2)
  names(Y.mat) <- c( "sample","trials","y")
  tmt <-c(rep(0,(trials*sampsize)),rep(1,(trials*sampsize)))
  #Combine "sample trial y" with "tmt"(0,1)
  #Set simulation matrix as sim, ordered by trials
  sim <- data.frame(Y.mat,tmt)  
  sim <- sim[order(sim$trials),]
  
  if(max(sim[,3]) > (1-1e-16) | min(sim[,3]) < 1e-16){
    sim[,3] <- (sim[,3] * (sampsize - 1) + 0.5) / sampsize
  }
  
  if(link.type=="wilcoxon"){
    outtest <- matrix(NA,nrow=trials,ncol=1)
    for(i in 1:trials){
      sub.sim <-  data.frame(sim[which(sim$trial==i),])
      out.wil <- wilcox.test(sub.sim[which(sub.sim[,4]==0),3],sub.sim[which(sub.sim[,4]==1),3])
      outtest[i,1] <- out.wil$p.value
    }
    Power = mean(as.numeric(outtest<sig.level))
  }
  else{
    outtest <- matrix(NA,nrow=trials,ncol=1)
    for(i in 1:trials){
      sub.sim <-  data.frame(sim[which(sim$trial==i),])
      fit1 <- suppressWarnings(do.call(betareg::betareg,list(formula=(sub.sim[,3] ~ sub.sim[,4]), link = link.type,type ="ML", model = T )))
      # link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog")
      out <- suppressWarnings(lmtest::waldtest(fit1,test = "F"))
      outtest[i,1] <- out$`Pr(>F)`[2]
    }
    Power = mean(as.numeric(outtest<sig.level))
  }
  return(Power)
}

betapwr2 <- function(mu0,sd0,mu1,sampsize,trials,seed,link.type,equal.precision,sd1,sig.level){
  seed.new <- seed
  Power <- tryCatch(betapwr2.base(mu0,sd0,mu1,sampsize,trials,seed.new,link.type,equal.precision,sd1,sig.level),error=function(e){return(NA)})
  while(is.na(Power[1])){
    seed.new <- seed.new + 1
    Power <- tryCatch(betapwr2.base(mu0,sd0,mu1,sampsize,trials,seed.new,link.type,equal.precision,sd1,sig.level),error=function(e){return(NA)})
  }
  return(Power)
}

sample.size.mid <- function(mu0,sd0,mu1,power.min,sig.level,trials,delta,seed,link.type,equal.precision,sd1){
  sample.size.starting <- round(do.call("power.t.test",list(delta = (mu1-mu0), sd = sd0, sig.level = sig.level,power = power.min))$n,0)
  if(is.null(delta)){
    delta <- 1
  }
  delta.current <- max(floor(sample.size.starting/2),delta)
  power.starting <- betapwr2(mu0,sd0,mu1,sample.size.starting,trials,seed,link.type,equal.precision,sd1,sig.level)
  if(power.starting < power.min){
    sample.size.ending <- sample.size.starting + delta.current
    power.ending <- betapwr2(mu0,sd0,mu1,sample.size.ending,trials,seed,link.type,equal.precision,sd1,sig.level)
    while(power.ending < power.min){
      power.starting <- power.ending
      sample.size.starting <- sample.size.ending
      sample.size.ending <- sample.size.starting + delta.current
      power.ending <- betapwr2(mu0,sd0,mu1,sample.size.ending,trials,seed,link.type,equal.precision,sd1,sig.level)
    }
    sample.size.lower <- sample.size.starting
    power.lower <- power.starting
    sample.size.upper <- sample.size.ending
    power.upper <- power.ending
  }
  if(power.starting >= power.min){
    sample.size.ending <- max((sample.size.starting - delta.current),3)
    power.ending <- betapwr2(mu0,sd0,mu1,sample.size.ending,trials,seed,link.type,equal.precision,sd1,sig.level)
    while((power.ending >= power.min)&(sample.size.ending > 3)){
      power.starting <- power.ending
      sample.size.starting <- sample.size.ending
      sample.size.ending <- max((sample.size.starting - delta.current),3)
      power.ending <- betapwr2(mu0,sd0,mu1,sample.size.ending,trials,seed,link.type,equal.precision,sd1,sig.level)
    }
    sample.size.lower <- sample.size.ending
    power.lower <- power.ending
    sample.size.upper <- sample.size.starting
    power.upper <- power.starting
  }
  reach.flag <- 0
  while(reach.flag == 0){
    if((sample.size.upper-sample.size.lower)<=delta){
      reach.flag <- 1
      sample.size.min <- sample.size.upper
      power.output <- power.upper
    }
    else{
      sample.size.midpt <- (sample.size.upper+sample.size.lower)%/%2
      power.midpt <- betapwr2(mu0,sd0,mu1,sample.size.midpt,trials,seed,link.type,equal.precision,sd1,sig.level)
      if(power.midpt < power.min){
        sample.size.lower <- sample.size.midpt
        power.lower <- power.midpt
      }
      else{
        sample.size.upper <- sample.size.midpt
        power.upper <- power.midpt
      }
    }
  }
  output <- matrix(c(sample.size.min,power.output),nrow = 1)
  colnames(output) <- c("minimum sample size","minimum power")
  return(output)
}

doit2 <- function(mu0,sd0,mu1.start, mu1.end, mu1.by, power.start, power.end, power.by, sig.level,trials,delta,seed,link.type,sample.size.matrix,equal.precision,sd1){
  loops <- length(seq(mu1.start,mu1.end,mu1.by))
  trials.parts <- 5*((nrow(sample.size.matrix)>10)+1)
  trials.total <- c(as.numeric(quantile(1:nrow(sample.size.matrix),(1:trials.parts)/(5*((nrow(sample.size.matrix)>10)+1)))),nrow(sample.size.matrix)+1)
  
  if(nrow(sample.size.matrix)==1){
    trials.current <- length(trials.total)
  }
  else{
    trials.current <- 2
  }
  
  Tmp.loc <- 1
  seed.start <- seed
  for (power.target in seq(power.start,power.end,power.by)) {
    mu1 <- mu1.start
    sample.size.unit <- matrix(NA, nrow = 2, ncol = length(link.type))
    for(i in 1:loops){
      while(Tmp.loc>=trials.total[trials.current]){
        print(noquote(paste0((trials.current-1)*(100/trials.parts),"% completed")))
        trials.current <- trials.current+1
      }
      for(j in 1:length(link.type)){
        sample.size.unit[,j] <- as.numeric(do.call("sample.size.mid",list(mu0 = mu0, sd0 = sd0, mu1 = mu1, 
                                                                          power.min = power.target,sig.level = sig.level,
                                                                          trials = trials,delta = delta,seed = seed.start,
                                                                          link.type = link.type[j],equal.precision=equal.precision,sd1=sd1)))
      }
      sample.size.matrix[Tmp.loc,] <- c(as.numeric(sample.size.unit),power.target,mu1)
      mu1 <- mu1+ mu1.by
      Tmp.loc <- Tmp.loc+1
      seed.start <- seed.start + 1
    }
  }
  return(sample.size.matrix)
}

#' @title Find minimum sample size with Beta distribution
#' @description  Find minimum sample sizes with Beta distribution and given mu0,sd0,mu1 and target powers.
#' @details The samplesize function allows you to control the number of trials in the simulation, 
#' the target power, delta, and the alternative means.
#' You can fix the alternative and vary power to match a desired sample size; 
#' Use default values for the number of trials for a quick view;
#' Use a larger number of trials (say 1000) and a smaller delta (say 1) to get better estimates.
#' @usage samplesize(mu0, sd0, mu1.start, mu1.end = NULL, mu1.by = NULL, 
#' power.start, power.end = NULL, power.by = NULL, sig.level = 0.05, 
#' trials = 100, delta = 1, seed = 1, link.type = "logit", 
#' equal.precision = TRUE, sd1 = NULL)
#' @param mu0 the mean for the control group
#' @param sd0 the standard deviation for the control group
#' @param mu1.start the starting value of mean for the treatment group under the alternative mu1
#' @param mu1.end the ending value of mean for the treatment group under the alternative mu1
#' @param mu1.by the step length of mean for the treatment group under the alternative mu1
#' @param power.start the starting value of target power
#' @param power.end the ending value of target power
#' @param power.by the step length of target power
#' @param sig.level significant level; default value is 0.05
#' @param trials the number of trials; default value is 100
#' @param delta the accuracy of the result; must be integer
#' @param seed the seed used in the simulation
#' @param link.type default link is "logit". Other link options include: "logit", "probit", "cloglog", "log", "loglog", "wilcoxon", or you can use "all" for all types of link
#' @param equal.precision equal dispersion parameter assumption in simulation
#' @param sd1 the standard deviation for the treatment group. Only applicable when equal.precision = FALSE
#' @return Return a table including minimum sample size and power, as well as the target power and mu1:
#' \item{minimum sample size: link type:}{minimum sample size for given given mu0, sd0, mu1, target power and type of link.}
#' \item{minimum power: link type:}{the minimum power greater than or equal to target power.}
#' \item{target power:}{the target power.}
#' \item{mu1:}{mean for the treatment group under the alternative.}
#' \item{mu0:}{the mean for the control group.}
#' \item{sd0:}{the standard deviation for the control group.}
#' @examples 
#' samplesize(mu0=0.56, sd0=0.255, mu1.start = 0.8, power.start =  0.9, trials = 25,
#' link.type = c("logit","wilcoxon"))
#' @importFrom stats rbeta wilcox.test quantile
#' @export

samplesize <- function(mu0, sd0, mu1.start, mu1.end = NULL, mu1.by = NULL,
                        power.start, power.end = NULL, power.by = NULL, 
                        sig.level=0.05, trials=100, delta=1, seed=1, 
                        link.type="logit", equal.precision=TRUE, sd1=NULL){
  if(link.type[1]=="all"){
    link.type <- c("logit", "probit", "cloglog", "log", "loglog","wilcoxon")
  }
  if(is.null(mu1.end) & is.null(mu1.by)){
    mu1.end <- mu1.start
    mu1.by <- 0
  }
  if(is.null(power.end) & is.null(power.by)){
    power.end <- power.start
    power.by <- 0
  }
  Power.matrix <- matrix(nrow=(length(seq(power.start,power.end,power.by))*length(seq(mu1.start,mu1.end,mu1.by))),ncol=(2*length(link.type)+2),NA)
  Power.matrix <- data.frame(Power.matrix)
  Power.matrix <- do.call("doit2",list(mu0 = mu0, sd0 = sd0, 
                                       mu1.start = mu1.start, mu1.end = mu1.end, mu1.by = mu1.by,
                                       power.start = power.start, power.end = power.end, power.by = power.by,
                                       sig.level = sig.level, trials = trials, delta = delta, seed = seed,
                                       link.type = link.type, sample.size.matrix = Power.matrix,equal.precision=equal.precision,sd1=sd1))
  Power.matrix <- cbind(Power.matrix,rep(mu0,nrow(Power.matrix)),rep(sd0,nrow(Power.matrix)))
  Power.names <- rep(NA,2*length(link.type))
  for(i in 1:length(link.type)){
    Power.names[2*i-1] <- paste("minimum sample size:",link.type[i])
    Power.names[2*i] <- paste("minimum power:",link.type[i])
  }
  colnames(Power.matrix) <- c(Power.names, "target power","mu1","mu0","sd0")
  Power.matrix <- Power.matrix[order(Power.matrix[,"target power"]),]
  return(Power.matrix)
}