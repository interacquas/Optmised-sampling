library(sp)
library(ggplot2)


accum <- read_excel("tempi-raggiungimento.xlsx", sheet = "accumulo")
accum <- as.data.frame(accum)


x <- accum$plotz
y1 <- accum$n1
y2 <- accum$n2
y3 <- accum$n3
y4 <- accum$n4
y5 <- accum$n5

plot(x, y1, type="s", col = "red",  lwd = 3, xlim= c(1, 12),
     main="Accumulation", ylab="Species Richness", xlab="Number of plots")

plot(x, y1, type="s", col = "black",  lwd = 1,
     main="Accumulation", ylab="Species Richness", xlab="Number of plots")
lines(x, y2, type="s", col = "red",  lwd = 3)

plot(x, y1, type="s", col = "black",  lwd = 1,
     main="Accumulation", ylab="Species Richness", xlab="Number of plots")
lines(x, y2, type="s", col = "black",  lwd = 1)
lines(x, y3, type="s", col = "red",  lwd = 3,)

plot(x, y1, type="s", col = "black",  lwd = 1,
     main="Accumulation", ylab="Species Richness", xlab="Number of plots")
lines(x, y2, type="s", col = "black",  lwd = 1)
lines(x, y3, type="s", col = "black",  lwd = 1)
lines(x, y4, type="s", col = "red",  lwd = 3)

plot(x, y1, type="s", col = "black",  lwd = 1,
     main="Accumulation", ylab="Species Richness", xlab="Number of plots")
lines(x, y2, type="s", col = "black",  lwd = 1)
lines(x, y3, type="s", col = "black",  lwd = 1)
lines(x, y4, type="s", col = "black",  lwd = 1)
lines(x, y5, type="s", col = "red",  lwd = 3)
