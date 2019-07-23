betapwr.PAR.base <- function(mu0, sd0, mu1, sampsize, N, seed, link.type, equal.precision, sd1, sig.level){
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
  
  #n is trial, N is number of trials
  #Y.H0.1 is simulated Y under H0, version 1
  #Y.Ha.1 is simulated Y under Ha, version 1
  Y.H0.1 <- Y.Ha.1 <- matrix(NA,sampsize,N)
  for (n in 1:N) {
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
  tmt <-c(rep(0,(N*sampsize)),rep(1,(N*sampsize)))
  #Combine "sample trial y" with "tmt"(0,1)
  #Set simulation matrix as sim, ordered by trials
  sim <- data.frame(Y.mat,tmt)  
  sim <- sim[order(sim$trials),]
  
  #Define output matrix outtest
  outtest <- rep(NA,N)
  if(max(sim[,3]) > (1-1e-16) | min(sim[,3]) < 1e-16){
    sim[,3] <- (sim[,3] * (sampsize - 1) + 0.5) / sampsize
  }
  for(i in 1:N){
    sub.sim <-  data.frame(sim[which(sim$trial==i),])
    fit1 <- suppressWarnings(do.call(betareg::betareg,list(formula=(sub.sim[,3] ~ sub.sim[,4]), link = link.type,type ="ML", model = T , subset = sub.sim[,4] != "c1")))
    out <- suppressWarnings(lmtest::waldtest(fit1,test = "F"))
    outtest[i] <- as.numeric(out$`Pr(>F)`[2])
  }
  power <- mean(outtest<sig.level)
  
  return(power)
}

betapwr.PAR <- function(mu0, sd0, mu1, sampsize, N, seed, link.type, equal.precision, sd1, sig.level){
  seed.new <- seed
  Power <- tryCatch(do.call("betapwr.PAR.base", list(mu0, sd0, mu1, sampsize, N, seed.new, link.type, equal.precision, sd1, sig.level)),error=function(e){return(NA)})
  while(is.na(Power[1])){
    seed.new <- seed.new + 1
    Power <- tryCatch(do.call("betapwr.PAR.base", list(mu0, sd0, mu1, sampsize, N, seed.new, link.type, equal.precision, sd1, sig.level)),error=function(e){return(NA)})
  }
  return(Power)
}

betapwr.NPAR.base <- function(mu0, sd0, mu1, sampsize, N, seed, equal.precision, sd1, sig.level){
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
  
  #n is trial, N is number of trials
  #Y.H0.1 is simulated Y under H0, version 1
  #Y.Ha.1 is simulated Y under Ha, version 1
  Y.H0.1 <- Y.Ha.1 <- matrix(NA,sampsize,N)
  for (n in 1:N) {
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
  tmt <-c(rep(0,(N*sampsize)),rep(1,(N*sampsize)))
  #Combine "sample trial y" with "tmt"(0,1)
  #Set simulation matrix as sim, ordered by trials
  sim <- data.frame(Y.mat,tmt)  
  sim <- sim[order(sim$trials),]
  
  #Define output matrix outtest
  outtest <- rep(NA,N)
  if(max(sim[,3]) > (1-1e-16) | min(sim[,3]) < 1e-16){
    sim[,3] <- (sim[,3] * (sampsize - 1) + 0.5) / sampsize
  }
  for(i in 1:N){
    sub.sim <-  data.frame(sim[which(sim$trial==i),])
    out.wil <- wilcox.test(sub.sim[which(sub.sim[,4]==0),3],sub.sim[which(sub.sim[,4]==1),3])
    outtest[i] <- as.numeric(out.wil$p.value)
  }
  power <- mean(outtest<sig.level)
  
  return(power)
}

betapwr.NPAR <- function(mu0, sd0, mu1, sampsize, N, seed, link.type, equal.precision, sd1, sig.level){
  seed.new <- seed
  Power <- tryCatch(do.call("betapwr.NPAR.base",list(mu0, sd0, mu1, sampsize, N, seed.new, equal.precision, sd1, sig.level)),error=function(e){return(NA)})
  while(is.na(Power[1])){
    seed.new <- seed.new + 1
    Power <- tryCatch(do.call("betapwr.NPAR.base",list(mu0, sd0, mu1, sampsize, N, seed.new, equal.precision, sd1, sig.level)),error=function(e){return(NA)})
  }
  return(Power)
}

doit <- function(mu0, sd0, mu1.start, mu1.end, mu1.by, sd1, ss.start, ss.end, ss.by, trials, seed, Power.matrix, link.type, equal.precision, sig.level){
  # prepare the processing bar
  loops <- length(seq(mu1.start,mu1.end,mu1.by))
  N.parts <- 5*((nrow(Power.matrix)>10)+1)
  N.total <- c(as.numeric(quantile(1:nrow(Power.matrix),(1:N.parts)/(5*((nrow(Power.matrix)>10)+1)))),nrow(Power.matrix)+1)
  if(nrow(Power.matrix)==1){
    N.current <- length(N.total)
  }
  else{
    N.current <- 2
  }
  
  # define total # of loops
  Tmp.loc <- 1
  for (sampsize in seq(ss.start,ss.end,ss.by)) {
    mu1 <- mu1.start
    seed.start <- seed
    for(i in 1:loops){
      while(Tmp.loc>=N.total[N.current]){
        print(noquote(paste0((N.current-1)*(100/N.parts),"% completed")))
        N.current <- N.current+1
      }
      seed.start <- seed.start +1
      #Plug in initials into defined function's variables
      Power.PAR <- rep(NA,length(link.type))
      for(j in 1:length(link.type)){
        Power.PAR[j] <- do.call("betapwr.PAR",list(mu0 = mu0, sd0 = sd0, mu1 = mu1,sampsize = sampsize, 
                                                   N = trials, seed = seed.start, link.type = link.type[j],
                                                   equal.precision = equal.precision, sd1 = sd1, sig.level = sig.level))
      }
      Power.NPAR <- do.call("betapwr.NPAR",list(mu0 = mu0, sd0 = sd0, mu1 = mu1,sampsize = sampsize, 
                                                N = trials, seed = seed.start, equal.precision = equal.precision, 
                                                sd1 = sd1, sig.level = sig.level))
      Power.matrix[Tmp.loc,] <- c(Power.PAR,Power.NPAR,sampsize,mu1)
      mu1 <- mu1+ mu1.by
      Tmp.loc <- Tmp.loc+1
    }
  }
  return(Power.matrix)
}



#' @title Find Power with Beta distribution
#' @description  Find the power for a given sample size when testing the null hypothesis that the means for the control and treatment groups are equal against a two-sided alternative.
#' @details betapower function allows you to control the number of trials in the simulation, 
#' the sample sizes used, and the alternative means. 
#' You can fix the alternative and vary sample size to match a desired power;
#' You can fix the sample size and vary the alternative to see which will match a desired power;
#' You can vary both;
#' Start with a small number of trials (say 100) to determine the rough range of sample sizes or alternatives;
#' Use a larger number of trials (say 1000) to get better estimates.
#' @usage betapower(mu0, sd0, mu1.start, mu1.end = NULL, mu1.by = NULL, 
#' ss.start, ss.end = NULL, ss.by = NULL, sig.level = 0.05,
#' trials = 100, seed = 1, link.type="logit",
#' equal.precision=TRUE, sd1 = NULL)
#' @param mu0 the mean for the control group
#' @param sd0 the standard deviation for the control group
#' @param mu1.start the starting value of mean for the treatment group under the alternative mu1
#' @param mu1.end the ending value of mean for the treatment group under the alternative mu1
#' @param mu1.by the step length of mean for the treatment group under the alternative mu1
#' @param ss.start the starting value of sample size
#' @param ss.end the ending value of sample size
#' @param ss.by the step length of sample size
#' @param sig.level significant level; default value is 0.05
#' @param trials the number of trials
#' @param seed the seed used in the simulation
#' @param link.type the type of link used in the beta regression. Default value is "logit", or you can use "all" or choose one or more of the following: "logit", "probit", "cloglog", "cauchit", "log", "loglog"
#' @param equal.precision equal dispersion parameter assumption in simulation
#' @param sd1 the standard deviation for the treatment group. Only applicable when equal.precision = FALSE
#' @return Return a matrix with 7 to 12 columns: 
#' \item{power.of.GLM: link name}{the power using regression method; it will return the power with every links if you use link.type = "all" statement.}
#' \item{power.of.Wilcoxon.test}{the power from Wilcoxon Rank sum test.}
#' \item{sample size}{sample size.} 
#' \item{mu1}{the mean for the treatment group under the alternative.}
#' \item{mu0}{the mean for the control group.}
#' \item{sd0}{the standard deviation for the control group.}
#' \item{trials}{the number of trials.}
#' @examples 
#' betapower(mu0 = 0.56, sd0 = 0.255, mu1.start = .70, mu1.end = .75, mu1.by = .05, 
#' ss.start = 30, ss.end = 50, ss.by = 20, trials = 20)
#' @importFrom stats rbeta wilcox.test quantile
#' @export

betapower <-function(mu0, sd0, mu1.start, mu1.end = NULL, mu1.by = NULL, 
                     ss.start, ss.end = NULL, ss.by = NULL, sig.level = 0.05,
                     trials = 100, seed = 1, link.type="logit", 
                     equal.precision=TRUE, sd1 = NULL){
  # define link.type = "all"
  if(link.type[1]=="all"){
    link.type <- c("logit", "probit", "cloglog", "log", "loglog")
  }
  # if mu1.end & mu1.by = NULL, set mu1.end as mu1.start
  if(is.null(mu1.end) & is.null(mu1.by)){
    mu1.end <- mu1.start
    mu1.by <- 0
  }
  # if ss.end & ss.by = NULL, set ss.end as ss.start
  if(is.null(ss.end) & is.null(ss.by)){
    ss.end <- ss.start
    ss.by <- 0
  }
  # initialize power table
  Power.matrix <- matrix(nrow=(length(seq(ss.start,ss.end,ss.by))*length(seq(mu1.start,mu1.end,mu1.by))),ncol=(length(link.type)+3),NA)
  Power.matrix <- data.frame(Power.matrix)
  # invoke doit function
  Power.matrix <- do.call("doit",list(mu0 = mu0, sd0 = sd0, mu1.start = mu1.start, mu1.end = mu1.end, 
                                      mu1.by = mu1.by, sd1 = sd1, ss.start = ss.start, ss.end = ss.end, 
                                      ss.by = ss.by, trials = trials, seed = seed, Power.matrix = Power.matrix, 
                                      link.type = link.type, equal.precision = equal.precision, sig.level = sig.level))
  # add in parameter values
  Power.matrix <- cbind(Power.matrix,rep(mu0,nrow(Power.matrix)),rep(sd0,nrow(Power.matrix)),rep(trials,nrow(Power.matrix)))
  # rename the columns of table
  Power.names <- rep(NA,length(link.type))
  for(i in 1:length(link.type)){
    Power.names[i] <- paste("power of GLM:",link.type[i])
  }
  colnames(Power.matrix) <- c(Power.names, "power of Wilcoxon test","sample size","mu1","mu0","sd0","trials")
  Power.matrix <- Power.matrix[order(Power.matrix[,"sample size"]),]
  # output power table
  return(Power.matrix)
}