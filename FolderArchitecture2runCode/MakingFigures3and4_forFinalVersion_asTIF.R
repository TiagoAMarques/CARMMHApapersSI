#Making Figure 3 high res

tiff(filename = "Figure3.tif",
     width = 1000, height = 800, units = "px", pointsize = 12,
     compression = "none",
     bg = "white")

cex.text <- 1.5

par(mfrow=c(1,1),mar=c(4.5,5,1.5,0.5))
#Value of A for other whales
# Pelagic dolphins
l1 <- 0.75
u1 <- 0.99
a1 <- 4.370439
b1 <- 3.580464
x <- seq(0, 1, length = 1000)
y1 <- dbeta(x, a1, b1)/(u1-l1)
plot(x* (u1-l1) + l1, y1, type = "l", xlim = c(0.6, 1),lwd=3, ylim = c(0, 16), 
     ylab = "Density", xlab = "Survival reduction factor",cex.axis=cex.text,cex.lab=cex.text,col="deeppink") 

#get realizations for reporting the mean below
n<-100000
x1s<-rbeta(10000,a1,b1)*(u1-l1)+l1
m1 <- mean(x1s)
sd1 <- sd(x1s)

# Mesopelagic divers
l3 <- 0.70
u3 <- 0.99
a3 <- 4.18
b3 <- 3.61
y3 <- dbeta(x, a3, b3)/(u3-l3)
lines(x* (u3-l3) + l3, y3, type = "l", col = "darkorchid1",lwd=3)

x3s<-rbeta(10000,a3,b3)*(u3-l3)+l3
m3 <- mean(x3s)
sd3 <- sd(x3s)

# sperm & beaked whales
l4 <- 0.63
u4 <- 0.99
a4 <- 4.16
b4 <- 2.96
y4 <- dbeta(x, a4, b4)/(u4-l4)
lines(x* (u4-l4) + l4, y4, type = "l", col = "dodgerblue",lwd=3)

x4s<-rbeta(10000,a4,b4)*(u4-l4)+l4
m4 <- mean(x4s)
sd4 <- sd(x4s)

legend("topleft",inset=0.03,
       legend = c("Pelagic dolphins", "Mesopelagic divers", "Bathypelagic divers"), 
       col = c("deeppink", "darkorchid1", "dodgerblue"), lty = c(1,1,1),lwd=3,cex=cex.text)

dev.off()


tiff(filename = "Figure4.tif",
     width = 1000, height = 800, units = "px", pointsize = 12,
     compression = "none",
     bg = "white")

cex.text <- 1.5


#offshore species
l2 <- 0
u2 <- 0.6
a2 <- 1.44
b2 <- 3.76
y2 <- dbeta(x, a2, b2)/(u2-l2)
plot(x* (u2-l2) + l2, y2, type = "l", col = "thistle4",lwd=3,ylab = "Density", xlab = "Recovery Proportion",cex.axis=cex.text,cex.lab=cex.text) 

dev.off()