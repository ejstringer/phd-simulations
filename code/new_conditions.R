
# starting values ----------------

# phases and lengths
phaseNo <- paste0(rep(c('L','I', 'D'), 3), rep(1:3, each = 3))[-1]
fstch2 <- c(0.017, 0.034, 0.027, 0.005, 0.012, 0.02, 0.01, 0.03)
phaselength <- c(1, 2, 2, 4, 2, 4, 2, 2)
# migration 

mPhase <- c(L = em.mRate(mean(c(0.02,0.027)), 25), # low
            # D = em.mRate(mean(c(0.034,0.03, 0.012)), 25), # decrease
            I = em.mRate(mean(c(0.01,0.005)), 50))# increase

data.frame(fst = fstch2, n1 = c(100, 50, 25, 50,50,25,50,50),
           n = c(100, 50, 25, 100,50,25,100,50)) %>% 
  mutate(m_rate,
         m = em.mRate(fst, n),
         gens = c(1,2,2,4,2,4,2,2),
         fstdiff = fst - lag(fst),
         fst2 = fst + fstdiff,
         #fst2 = ifelse(fst < lag(fst), fst - fstdiff, fst + fstdiff),
         # fst2 = ifelse(fst2 <= 0, 0.001, fst2),
         m2 = em.mRate(fst2, n1),
         # Nm = m*n,
         # m_adjust = m*c(2,2),
         # m_adjust = ifelse(m_adjust > 1, 0.99, m_adjust),
         # npp = npp[-1,-1],
         tm = log(1/2)/log((1-m)^2*(1-(1/(2*n1)))),
         t =log(1/2)/log((1-m2)^2*(1-(1/(2*n1))))
         #t2 = lag(t)
         #tadjust = log(1/2)/log((1-m_adjust)^2*(1-(1/(2*n))))
  )

em.mRate(0.017, 100)*100

data.frame(fst = fstch2, npp = npp$npp[-1],
           phaseNo = paste0(rep(c('L', 'I', 'D'), 3),rep(1:3, each = 3))[-1],
           N = c(100, 50, 25, 50,50,25,50,50)) %>% 
  mutate(start = ifelse(fst == 0.017, 14.46, 0),
         phase = str_sub(phaseNo, 1,1),
         phase = ifelse(phaseNo == 'I1', 'int', phase),
         # phase = ifelse(phaseNo == 'D2', 'D2', phase),
         x = 14.46/60.71,
         Nm = (npp*x),
         m = Nm/N,
         m_rate,
         mfst = em.mRate(fst, N)
  ) %>% 
  group_by(phase, N) %>% 
  summarise(fst = mean(fst),
            npp = mean(npp)) %>% 
  mutate(m = em.mRate(fst, N)) %>% 
  filter(phase != 'int') %>% 
ggplot(aes(fst, m))+
  geom_bar(stat = 'identity')+
  theme_classic()


fst1 = 0.27
mfst = 0.005


N = 50
m2 = seq(0.1,0.5, 0.01)


data.frame(N = rep(c(25, 50, 100), each = length(m2)),
           m2) %>%
  mutate(t  = log(1/2)/log((1-m2)^2*(1-(1/(2*N))))) %>% 
  ggplot(aes(m2, t, colour = N, group = N))+
  theme_classic()+
  geom_hline(yintercept = c(2,4), lty = 2)+
  geom_line()



m.time <- function(N, t){
  N2 = 2*N
  tN <- (exp(log(1/2)/t))/ (1 - (1/(N2)))
  d = 1-tN
  b = 2
  a = 1
  
  ((-b)+sqrt(b^2-4*a*d))/2*a
  ((-b)-sqrt(b^2-4*a*d))/2*a
  
  
}


m = c()
log(1/2)/log((1-m)^2*(1-(1/(2*n))))

em.mRate(0.017, 50)
npp <- read.csv('../phd-analysis2/output/npp_means.csv')
mPhase
m_rate <- c(0.15, 0.09, 0.42, 0.66,0.3, 0.42, 0.66, 0.09)

ngrids <- c(37, 25, 17, 23, 16, 33, 27, 32)
## define simulation conditions ---------


conditions <- data.frame(gen = 1:sum(phaselength), 
                         period = rep(1:8, phaselength),
                         fstch2 = rep(fstch2, phaselength),
                         phaseNo = rep(phaseNo, phaselength), 
                         migration = rep(m_rate, phaselength),
                         rep = rep(c(1,1,1,2,2,2,3,3), phaselength),
                         year = 0.5) %>% 
  mutate(phase = str_sub(phaseNo, 1,1),
         N = case_when(
           phase == 'L' ~ 25,
           phase == 'I' ~ 100,
           phase == 'D' ~ 25,
         ),
         N = ifelse(phaseNo == 'I2', 250, N),
         offspring = 4) %>% 
  left_join(data.frame(grids = phase_nGrids, phaseNo = names(phase_nGrids)))

for(i in 2:nrow(conditions)){
  x <- ifelse(conditions$phase[i] == 'L', 1, 0.5)
  conditions$year[i] <- conditions$year[i-1] + x
}

conditions <- conditions %>% 
  rowwise() %>% 
  mutate(yearsSince = case_when(
    rep == '1' ~ year - 0,
    rep == '2' ~ year - 3.5,
    rep == '3' ~ year - 10.5,
  ))


conditions
