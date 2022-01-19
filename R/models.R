SSgauss <- selfStart(~ h*exp(-(x-mu)^2/(2*sigma^2)), function(mCall, data, LHS, ...) {

    xy <- sortedXyData(mCall[["x"]], LHS, data)

    len <- dim(xy)[1]
    xyarea <- sum((xy[2:len,2]+xy[1:(len-1),2])*(xy[2:len,1]-xy[1:(len-1),1]))/2
    maxpos <- which.max(xy[,2])

    mu <- xy[maxpos,1]
    h <- xy[maxpos,2]
    sigma <- xyarea/(h*sqrt(2*pi))

    value <- c(mu, sigma, h)
    names(value) <- mCall[c("mu", "sigma", "h")]
    value

}, c("mu", "sigma", "h"))

etg <- function(x, H, t1, tt, k1, kt, lambda1, lambdat, alpha, beta)
    2*H*exp(0.5)/((1+lambda1*exp(k1*(t1-x)))^alpha + (1+lambdat*exp(kt*(x-tt)))^beta - 1)
