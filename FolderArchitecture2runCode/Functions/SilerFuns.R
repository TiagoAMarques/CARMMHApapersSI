#these are funtions related to the Siler model
#Based on funtions provided by Lori way back!

lx <- function(x,a1,a2,a3,b1,b3) {
  exp(-a2*x) * exp((-a1)*(1-exp(-b1*x))) * exp(a3*(1-exp(b3*x)))
}


# Age specific mortality rate
px <- function(x,a1,a2,a3,b1,b3,scaling=1) {
  lx(x+scaling,a1,a2,a3,b1,b3)/lx(x,a1,a2,a3,b1,b3)
}