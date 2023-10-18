## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  if (!require(BetaPASS)){
#    devtools::install_github("CastleLi/draft/BetaPASS")
#    Needed_packages <- c("Rcpp","betareg","ggplot")
#    install.packages(Needed_packages)
#  }

## ---- results='hide'----------------------------------------------------------
library(BetaPASS)
Power.mat <- betapower(mu0 = 0.56, sd0 = 0.255,
                       mu1.start = .70, mu1.end = .75, mu1.by = .05,
                       ss.start = 30, ss.end = 50, ss.by = 5,
                       trials = 40, seed = 1, 
                       link.type = c("logit"))

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(Power.mat$Power.matrix)

## ---- fig.show='hold', fig.width = 9, fig.height =6---------------------------
plot(Power.mat, link.type = "logit", by = "mu1")

## -----------------------------------------------------------------------------
samplesize(mu0=0.56, sd0=0.255, mu1.start = 0.75, power.start = 0.8, trials = 40,
           link.type = c("logit","wilcoxon"))

## -----------------------------------------------------------------------------
SS.mat <- samplesize(mu0=0.56, sd0=0.255, 
                       mu1.start = 0.70, mu1.end = 0.75, mu1.by = 0.05, 
                       power.start = 0.8, power.end = 0.9, power.by = 0.1, 
                       trials = 40, link.type = c("logit","wilcoxon"))

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(SS.mat$Power.matrix)

## ---- fig.show='hold', fig.width = 9, fig.height =6---------------------------
plot(SS.mat, link.type = c("logit","wilcoxon"))

