################################################################################
#Section 4.3
################################################################################

library(haven)
library(xts)
library(zoo)
library(haven)
library(dplyr)
library(lubridate)
library(DFA)
library(CTRE)
library(pracma)
#install.packages("DFA")

#Clearly no pattern found in the empirical copula
ctre_d <- episode_details[c("t", "negDst")]
ctre_r <- ctre(ctre_d)
empcopula(ctre_r, OCTRE = FALSE)
ctre_d %>% ctre() %>% interarrival() %>% mlqqplot(tail = 0.9, log = 'xy')
###################################################################### Surrogate

#devtools::install_github("mlaib/MFDFA")
library(MFDFA)
library(nonlinearTseries)
second_series <- FFTsurrogate(sample(episode_details$negDst), 1)
second_series <- as.numeric(second_series)

#third_series <- high_intensity_events$negDst

########################################################################## MFDFA
#Code inspired by https://search.r-project.org/CRAN/refmans/MFDFA/html/MFXDFA.html
q<--10:10
m<-1
scale=10:100
mfdfa_first <- MFDFA(episode_details$negDst, scale, m, q)
mfdfa_second <- MFDFA(second_series, scale, m, q)

dev.new(width=12, height=8)

layout(matrix(c(1,2,3,4,5,6,7,8), 2, 4, byrow = TRUE), widths=c(4, 1, 4, 1), heights=c(4, 4))

poly_fit<-function(x,y,n){
  formule<-lm(as.formula(paste('y~',paste('I(x^',1:n,')', sep='',collapse='+'))))
  res1<-coef(formule)
  poly.res<-res1[length(res1):1]
  allres<-list(polyfit=poly.res, model1=formule)
  return(allres)}

#First plot: Fluctuation function Fq
par(mar=c(5, 5, 4, 2))
p1 <- c(1, which(q == 0), which(q == q[length(q)]))
plot(log2(scale), log2(mfdfa_first$Fqi[,1]), pch=16, col=1, axes = F, xlab = "s (days)",
     ylab=expression('log'[2]*'(F'[q]*')'), cex=1, cex.lab=1.6, cex.axis=1.6,
     main= "Fluctuation function Fq",
     ylim=c(min(log2(mfdfa_first$Fqi[,c(p1)]), log2(mfdfa_second$Fqi[,c(p1)])),
            max(log2(mfdfa_first$Fqi[,c(p1)]), log2(mfdfa_second$Fqi[,c(p1)]))))
lines(log2(scale), mfdfa_first$line[,1], type="l", col=1, lwd=2)
grid(col="midnightblue")
axis(2)
lbl <- scale[c(1, floor(length(scale)/8), floor(length(scale)/4),
               floor(length(scale)/2), length(scale))]
att <- log2(lbl)
axis(1, at=att, labels=lbl)
for (i in 2:3) {
  k <- p1[i]
  points(log2(scale), log2(mfdfa_first$Fqi[,k]), col=i, pch=16)
  lines(log2(scale), mfdfa_first$line[,k], type="l", col=i, lwd=2)
}

#Add the second time series to the first plot
for (i in 1:3) {
  k <- p1[i]
  points(log2(scale), log2(mfdfa_second$Fqi[,k]), col=i+3, pch=17) # Use different symbols
  lines(log2(scale), mfdfa_second$line[,k], type="l", col=i+3, lwd=2, lty=2) # Use different line types
}
par(mar=c(0, 0, 0, 0))
plot(0, 0, type = "n", xlab = "", ylab = "", axes = FALSE)
legend("left", legend=c(paste('q','=',q[p1], "First"), paste('q','=',q[p1], "Second")),
       col=c(1,2,3,4,5,6), pch=c(16,16,16,17,17,17), lty=c(1,1,1,2,2,2), cex=1.2, bty="n")

#Second plot: Hurst exponent
par(mar=c(5, 5, 4, 2))
plot(q, mfdfa_first$Hq, col=1, axes= F, ylab=expression('h'[q]), pch=16, cex.lab=1.8,
     cex.axis=1.8, main="Hurst exponent", ylim=c(min(c(mfdfa_first$Hq, mfdfa_second$Hq)), max(c(mfdfa_first$Hq, mfdfa_second$Hq))))
points(q, mfdfa_second$Hq, col=4, pch=17)
grid(col="midnightblue")
axis(1, cex=4)
axis(2, cex=4)
par(mar=c(0, 0, 0, 0))
plot(0, 0, type = "n", xlab = "", ylab = "", axes = FALSE)
legend("left", legend=c("Original data (Dst)", "Surrogate Original data"), col=c(1,4), pch=c(16,17), cex=1.2, bty="n")

#Third plot: Mass exponent
par(mar=c(5, 5, 4, 2))
plot(q, mfdfa_first$tau_q, col=1, axes=F, cex.lab=1.8, cex.axis=1.8,
     main="Mass exponent", pch=16, ylab=expression(tau[q]), ylim=c(min(c(mfdfa_first$tau_q, mfdfa_second$tau_q)), max(c(mfdfa_first$tau_q, mfdfa_second$tau_q))))
points(q, mfdfa_second$tau_q, col=4, pch=17)
grid(col="midnightblue")
axis(1, cex=4)
axis(2, cex=4)
par(mar=c(0, 0, 0, 0))
plot(0, 0, type = "n", xlab = "", ylab = "", axes = FALSE)
legend("left", legend=c("Original data (Dst)", "Surrogate Original data"), col=c(1,4), pch=c(16,17), cex=1.2, bty="n")

#Fourth plot: Multifractal spectrum
par(mar=c(5, 5, 4, 2))
plot(mfdfa_first$spec$hq, mfdfa_first$spec$Dq, col=1, axes=F, pch=16, ylab=bquote("f ("~alpha~")"), cex.lab=1.8, cex.axis=1.8, xlab=bquote(~alpha), xlim=c(0, 4), ylim=c(min(c(mfdfa_first$spec$Dq, mfdfa_second$spec$Dq)), max(c(mfdfa_first$spec$Dq, mfdfa_second$spec$Dq))))
points(mfdfa_second$spec$hq, mfdfa_second$spec$Dq, col=4, pch=17)
grid(col="midnightblue")
axis(1, cex=4)
axis(2, cex=4)

x1 <- mfdfa_first$spec$hq
y1 <- mfdfa_first$spec$Dq
x2 <- mfdfa_second$spec$hq
y2 <- mfdfa_second$spec$Dq

rr <- poly_fit(x1, y1, 4)
mm1 <- rr$model1
mm <- rr$polyfit
x_fit <- seq(min(c(x1, x2)), max(c(x1, x2)), 0.01)
curv <- mm[1]*x_fit^4 + mm[2]*x_fit^3 + mm[3]*x_fit^2 + mm[4]*x_fit + mm[5]
lines(x_fit, curv, col="red", lwd=2)

rr_second <- poly_fit(x2, y2, 4)
mm1_second <- rr_second$model1
mm_second <- rr_second$polyfit
curv_second <- mm_second[1]*x_fit^4 + mm_second[2]*x_fit^3 + mm_second[3]*x_fit^2 + mm_second[4]*x_fit + mm_second[5]
lines(x_fit, curv_second, col="blue", lwd=2, lty=2)
par(mar=c(0, 0, 0, 0))
plot(0, 0, type = "n", xlab = "", ylab = "", axes = FALSE)
legend("left", legend=c("Original data (Dst)", "Surrogate Original data"), col=c("red", "blue"), lwd=2, lty=c(1,2), cex=1.2, bty="n")

