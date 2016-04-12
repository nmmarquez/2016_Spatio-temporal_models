#Week 3 - temporal models
#Lecture R code
#Noble Hendrix
#March 18, 2016

# simulate Gaussian White Noise process
set.seed(123)
y = rnorm(250)
ts.plot(y,main="Gaussian White Noise Process",xlab="time",ylab="y(t)",
        col="blue", lwd=2)
abline(h=0)

# simulate random walk 
set.seed(321)
e = rnorm(250)
y.rw = cumsum(e)
ts.plot(y.rw, lwd=2, col="blue", main="Random Walk")
abline(h=0)

# simulate AR(1) process: alpha = 0.9
ar1.model = list(ar=0.9)
mu = 0
set.seed(123)
ar1.sim = arima.sim(model=ar1.model,n=250)
#compute the theoretical ACF for the AR(1) with alpha = 0.9
ar1.acf = ARMAacf(ar=0.9, ma=0, lag.max=10)

par(mfrow=c(1,2))
ts.plot(ar1.sim,main="AR(1) Process:  alpha=0.9",
xlab="time",ylab="y(t)", col="blue", lwd=2)
abline(h=0)
abline(h=mu, lty = 2)
# ACF for AR(1) model
plot(0:10, ar1.acf,type="h", col="blue", lwd=2,
main="Theoretical ACF for AR(1): alpha=0.9",xlab="lag",ylab="rho(j)", ylim = c(0,1))
abline(h=0)


# simulate MA(1) process with beta 0.9 and e(t) ~ N(0,1)
ma1.model = list(ma=0.9)
set.seed(123)
ma1.sim = arima.sim(model=ma1.model,n=250)
ma1.acf = ARMAacf(ar=0, ma=0.9, lag.max=10)

par(mfrow=c(1,2))
ts.plot(ma1.sim,main="MA(1) Process: beta = 0.9",
xlab="time",ylab="y(t)", col="blue", lwd=2)
abline(h=c(0))
plot(0:10, ma1.acf,type="h", col="blue", lwd=2,
main="Theoretical ACF for MA(1): beta=0.9",xlab="lag",ylab="rho(j)")

abline(h=0)
par(mfrow=c(1,1)


