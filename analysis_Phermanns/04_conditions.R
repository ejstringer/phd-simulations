
# starting values ----------------

# phases and lengths
phaseNo <- paste0(rep(c('L','I', 'D'), 3), rep(1:3, each = 3))[-1]
fstch2 <- c(0.017, 0.034, 0.027, 0.005, 0.012, 0.02, 0.01, 0.03)
phaselength <- c(1, 2, 2, 4, 2, 4, 2, 2)
ngrids <- c(37, 25, 17, 23, 16, 33, 27, 32)

# migration 

N_sim <- c(100,      50,    25,   50,     50,   25,   50, 50)
mBased <- apply(data.frame(n = N_sim, f = fstch2),1, function(x) em.mRate(x[2], x[1]))

data.frame(n = N_sim, fst = fstch2, m = round(mBased,2)) %>% 
  mutate(Nm = round(n*m),
         s = phaseNo)

sapply(c(0.003, 0.007, 0.003, 0.001, NA, 0.005, 0, 0),
       function(x) em.mRate(x, 100))

## old migration estimates
mPhase <- c(L = em.mRate(mean(c(0.02,0.027)), 25), # low
            # D = em.mRate(mean(c(0.034,0.03, 0.012)), 25), # decrease
            I = em.mRate(mean(c(0.01,0.005)), 50))# increase
mPhase
m_rate <- c(0.15, 0.09, 0.42, 0.66,0.3, 0.42, 0.66, 0.09)


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
           phase == 'D' ~ 50,
         ),
         N = ifelse(phase == 'D' & duplicated(phaseNo), 25, N),
         N = ifelse(phaseNo == 'I2', 200, N),
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
