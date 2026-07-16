#Heeger model

theta = c(0,20,45,70,90)
c_s = .7
c_c = .7
sigma = 1.5
a_c = 1
w = 20

overall_a  = 1

R = overall_a * (a_c * c_c / ((a_c * c_c) + c_s*exp((-theta) /w) + sigma))

plot(theta,R,'l')


# Brouwer & Heeger Cross-orientation suppression in human visual cortex
rmax = 1
baseline = 1
contrast = seq(0,2,.1)
n = 3
sigma = 1

R = rmax * (contrast^n /( contrast^n + sigma )) + baseline
plot(contrast,R)

sigma = 5
R = rmax * (contrast^n /( contrast^n + sigma )) + baseline
lines(contrast,R)
