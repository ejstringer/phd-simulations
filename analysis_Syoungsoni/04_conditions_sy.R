
# starting values ----------------

# phases and lengths
phaseNo <- paste0(rep(c('L','I', 'D'), 3), rep(1:3, each = 3))[-1]
fstch2 <- c(0.003, 0.007, 0.003, 0.001, NA, 0.005, 0, 0)
phaselength <- c(1, 1, 1, 2, 1, 3, 1, 1)
# migration 
phaselength %>% sum

mPhase <- c(L = em.mRate(mean(fstch2[c(2,3,6,8)]), 100), # low
            D = em.mRate(mean(fstch2[c(2,3,6,8)]), 100), # decrease (same)
            I = em.mRate(mean(fstch2[c(4,7)]), 100))# increase
mPhase
em.mRate(0.003, 100)
m_rate <- c(0.83, 0.66, 0.66, 0.99,0.66, 0.66, 0.99, 0.66)
ngridsSy <- rowSums(table(ph@other$ind.metrics$phaseNo,
                          ph@other$ind.metrics$gridId)>0)
ngrids <- ngridsSy[-1]
## define simulation conditions ---------


conditions <- data.frame(gen = 1:sum(phaselength), 
                          period = rep(1:8, phaselength),
                          fstch2 = rep(fstch2, phaselength),
                          phaseNo = rep(phaseNo, phaselength), 
                          migration = rep(m_rate, phaselength),
                          rep = rep(c(1,1,1,2,2,2,3,3), phaselength),
                          year = 1) %>% 
  mutate(phase = str_sub(phaseNo, 1,1),
         N = case_when(
           phase == 'L' ~ 25,
           phase == 'I' ~ 100,
           phase == 'D' ~ 25,
         ),
         N = ifelse(phaseNo == 'I2', 250, N),
         N = 100,
         offspring = 6) %>% 
  left_join(data.frame(grids = phase_nGrids, phaseNo = names(phase_nGrids)))

for(i in 2:nrow(conditions)){
  x <- ifelse(conditions$phase[i] == 'L', 1, 1)
  conditions$year[i] <- conditions$year[i-1] + x
}
# 
# conditions <- conditions %>% 
#   rowwise() %>% 
#   mutate(yearsSince = case_when(
#     rep == '1' ~ year - 0,
#     rep == '2' ~ year - 3.5,
#     rep == '3' ~ year - 10.5,
#   ))


conditions
