first_sin<-sin(.2*seq(100))
second_sin<-sin(.2*seq(100)+.5)
plot(first_sin,type='l')
lines(second_sin,col='red')

res<-ccf(second_sin, first_sin)
plot(res)
res$lag[res$acf==max(res$acf)]
