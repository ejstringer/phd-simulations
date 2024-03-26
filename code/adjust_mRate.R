
#conditions adjust

mAdjust <- c(I1 = 0,D1 = -0.02,L2 = 0.02,I2 =  -0.02,D2 = -0.02,L3 =  0, I3 = 0.02,D3 = -0.02)

mAdjust <- c(I1 = 0,D1 = -0.04,L2 = 0.04,I2 =  -0.04,D2 = -0.04,L3 =  0, I3 = 0.04,D3 = -0.04)

mAdjust <- c(I1 = 0,D1 = -0.02,L2 = 0,I2 =  0,D2 = -0.02,L3 =  0, I3 = 0.0,D3 = -0.02)

m_rate <- c(0.15, 0.09, 0.42, 0.66,0.3, 0.42, 0.66, 0.09)
x <- m_rate


for (i in 1:10) {
  x<- x+mAdjust
  
  x <- ifelse(x <0 , 0, x)
  x <- ifelse(x >= 1, 1, x)
  print(unname(round(x,2)))
}

