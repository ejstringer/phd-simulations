                          # Simulation Models 
# libraries --------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(dartR)
library(emphd)
library(lme4)
library(scales)
        
# load -------------------------------------------------------------------------
simhoInd <- read.csv("./output/sim_ind_heterozygosity.csv") 
simho <- read.csv("./output/sim_heterozygosity.csv")
simhe <- simho
simhe$ho <- simhe$he 

simhoInd <- read.csv("./output/sim_ind_heterozygosity_ddd_n.csv")
simho <- read.csv("./output/sim_heterozygosity_ddd_n.csv")
simhe <- simho
simhe$ho <- simhe$he 

simhoInd <- read.csv("./output/sim_fst_ddd_n.csv")
simho <- read.csv("./output/sim_fst_ddd_p.csv")
simhe <- simho
simhe$ho <- simhe$he 

simhoInd$ho <- simhoInd$Fst
simho$ho <- simho$Fst
simhe$ho <- simhe$Fst
simhe$he <- simhe$Fst

# models -----------------------------------------------------------------------
simmodel <- function(simhoInd){
  
  dataName <- (deparse(substitute(simhoInd)))
  simhoInd$gen <- as.numeric(as.character(sub("gen", "", simhoInd$trip)))
  simhoInd$sinceEvent <- "event 2"
  simhoInd$sinceEvent[simhoInd$gen < 5] <- "event 1"
  simhoInd$sinceEvent[simhoInd$gen > 14]    <- "event 3"
  
  event <- levels(factor(simhoInd$sinceEvent))
  
  mlinear <- lapply(event,
                    function(x) lm(log(ho) ~ gen, 
                                   data = filter(simhoInd, 
                                                 sinceEvent == x)))
  names(mlinear) <- paste(dataName, event) 
  sjPlot::tab_model(mlinear, title = dataName,
                    dv.labels = event,
                    collapse.ci = T,
                    digits = 4)
  #plot(simhoInd$gen, simhoInd$ho, col = factor(simhoInd$sinceEvent))
  return(mlinear)
  
}
mind <- simmodel(simhoInd)
mho <- simmodel(simho)
mhe <- simmodel(simhe)

# sim table --------------------------------------------------------------------
m.fun <- function(x, roundTo = 4, coeftype = "slope"){
  
  s <- summary(x)
  int <- s$coefficients[1,1]
  slop <- s$coefficients[2,1]
  intse <- s$coefficients[1,2]
  slopse <- s$coefficients[2,2]
  
  intci <- c(int - intse*2, int + intse*2)
  slopci <- c(slop - slopse*2, slop + slopse*2)
  
  if(coeftype == "slope") dd <- slopci
  if(coeftype == "intercept") dd <- intci
  return(dd)
} 
cIntervals <- lapply(c(mind, mho, mhe), m.fun) %>% 
  do.call("rbind", .)
cIntervals

cIntervalsInt <- lapply(c(mind, mho, mhe), function(x) m.fun(x, 4, "intercept")) %>% 
  do.call("rbind", .)
cIntervalsInt

# plot coefficients ------------------------------------------------------------
mcoef <- sapply(c(mind, mho, mhe), coef)[2,]
mcoefint <- sapply(c(mind, mho, mhe), coef)[1,]

hetcoef <- data.frame(model = factor(gsub(" .*$", "", names(mcoef))),
           event = as.numeric(gsub("[^0-9.-]", "", names(mcoef))),
           het = mcoef,
           row.names = NULL)
hetcoef$lower <- cIntervals[,1]
hetcoef$upper <- cIntervals[,2]
hetcoef$R2 <- sapply(c(mind, mho, mhe), function(x) summary(x)$r.squared)
hetcoef$gen <- rep(c(5,10,5), 3)
hetcoef$eventgen <- paste0("event ", hetcoef$event, " (", hetcoef$gen, 
                           " generations)")
hetcoef$int <- mcoefint
hetcoef$intlower <- cIntervalsInt[,1]
hetcoef$intupper <- cIntervalsInt[,2]
hetcoef <- filter(hetcoef, model == 'simhoInd')
pp1 <- ggplot(hetcoef, aes(model, het, colour = model))+
  geom_point(size = 3) +
  facet_grid(~eventgen) +
  theme_bw()+
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  ylab("slope coefficient")+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.0);pp1

pp2 <- ggplot(hetcoef, aes(model, int, colour = model))+
  geom_point(size = 3) +
  facet_grid(~eventgen) +
  theme_bw() +
  ylab("intercept coefficient")+
  theme(legend.position = "none",
        axis.title.x = element_blank())+
  geom_errorbar(aes(ymin = intlower, ymax = intupper), width = 0.0);pp2

pp3 <- ggplot(hetcoef, aes(model, R2, colour = model))+
  geom_point(size = 3) +
  facet_grid(~eventgen) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.title.x = element_blank()) +
  ylab(expression(italic(R^2))); pp3


gpp <- gridExtra::grid.arrange(pp1, pp2, pp3, nrow = 3)


ggsave("./figures/heterozygosity_simulations_model.png", gpp)
 
# plot data --------------------------------------------------------------------
mancol <- c(he = "#72D075", ho = "#2553DA", indHo = "#E05BA4")

simhoInd$gen <- as.numeric(as.character(sub("gen", "", simhoInd$trip)))
simho$gen <- as.numeric(as.character(sub("gen", "", simho$trip)))
simhe$gen <- as.numeric(as.character(sub("gen", "", simhe$trip)))

df <- simhoInd %>% 
  group_by(gen) %>% 
  summarise(ho = mean(ho))

pp <- ggplot(simhoInd, aes(x = gen, y = ho))+
  geom_jitter(width = 0.1, alpha = 0.2) +
  geom_line(data = df, colour = 'red', size = 3)+
#  geom_point(data = simhe, aes(x = gen, y = he), size = 3, fill = 3, pch = 21)+
#  geom_point(data = simho, aes(x = gen, y = ho), size = 3, fill = 2, pch = 21)+
  theme_bw()+
  ylab("heterozygosity")+
  xlab("generation") +
  scale_x_continuous(breaks = 0:20);
pp  
ggsave("./figures/heterozygosity_simulations_dp.png", pp)

# simulate fis vs fst ----------------------------------------------------------
cal.heterozygosity <- function(pop, metric = "he"){
  N <- ncol(pop)
  #print(N)
  p <- sum(pop)/(N*2)
  q <- 1-p
  he <- p*q*2
  ho <- sum(colSums(pop) == 1)/(N)
  fis <- (he-ho)/he
  
  if(metric == "he") h <- he
  if(metric == "ho") h <- ho
  if(metric == "p") h <- p
  if(metric == "q") h <- q
  if(metric == "fis") h <- fis
  
  #cat("returning", metric)
  return(h)
}

# simulate population 
jarSim <- function(generations = 10, N = 100){
  
  pop1 <- matrix(nrow = 2, ncol = N, data = c(0,1))
  #print(pop1[,1:10])
  jar <- list(pop1)
  
  for(j in 1:generations){  
    probA <- sum(pop1)/(N*2)
    proba <- 1 - probA 
    
    for(i in 1:length(pop1)){
      pop1[i] <- sample(0:1, size = 1, replace = T, prob = c(proba, probA))  
      
    }
    jar[[j+1]] <- pop1
  }
  length(jar)
  names(jar) <- paste0("gen",0:generations)
  return(jar)
}




plot.fis.fst <- function(fixation = T){
  
  for(k in 1:50){
    jar1 <- jarSim(100, 50)
    jar2 <- jarSim(100, 50)
    jars <- list()
    for(i in 1:length(jar1)){
      jars[[i]] <- cbind(jar1[[i]], jar2[[i]])  
    }
    jars[[i]]
    
    heGen <- sapply(jars, cal.heterozygosity)
    hoGen <- sapply(jars, function(x) cal.heterozygosity(x, metric = "ho"))
    
    hsGen <- (sapply(jar1, cal.heterozygosity)+sapply(jar2, cal.heterozygosity))/2
    
    fis <- (heGen-hoGen)/heGen
    fst <- (heGen - hsGen)/heGen
    
    # plot(x = c(1, length(jars)), y = c(-1,1), col = "white")
    # points(fis, col = "red", pch = 16, cex = 2) 
    # points(fst, col = "black", pch = 16)
    
    dffis <- data.frame(gen = names(jar1), stat = "Fis", 
                        value = fis ,row.names = NULL, psize = 4,
                        i = as.numeric(sub("gen", "", names(jar1))))
    dffst <- data.frame(gen = names(jar1), stat = "Fst", value = fst ,
                        row.names = NULL,  psize = 2,
                        i = as.numeric(sub("gen", "", names(jar1))))
    df <- rbind(dffis, dffst)
    
    dftest <- df[complete.cases(df$value),]
    
    finalv <- dftest$value[nrow(dftest)]
    
    test <- ifelse(fixation, finalv > 0.5, finalv < 0.5)
    if(test){
      simPlot <- ggplot(df, aes(x = i, y = value, colour = stat))+
        geom_hline(yintercept = c(-1, 0, 1), 
                   linetype = c("solid", "dashed", "solid"),
                   colour = c("grey","black", "grey"))+
        geom_point(size = df$psize)+
        scale_color_manual(values = c("red", "black"))+
        theme_classic()+
        ylim(-1,1)+
        theme(legend.position = c(0.1, 0.8))+
        xlab("generation")
      print(simPlot)
      ggsave(file = "./figures/simulate_fis_fst.png", simPlot)
      
      stop("Not an Error ~ plot saved")
    }
  }
  
}

plot.fis.fst(T)
#plot.fis.fst(F)



