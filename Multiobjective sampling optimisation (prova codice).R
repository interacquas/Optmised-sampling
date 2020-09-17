library("plot3D")
library("scales")

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
} # funzione per normalizzare

alfa <- 0.2
beta <- 1
gamma <- 0



ndvi <- normalize(out1$`Aggregated matrix`$Variance)
ign <- normalize(out1$`Aggregated matrix`$`Mean Ignorance`)
dist <- normalize(out1$`Aggregated matrix`$`Mean Dist`)
final <- (alfa * ndvi) + (beta * ign) + (gamma * dist)

index <- match(max(final), final)
index


scatter3D(alfa*ndvi, beta*ign, gamma*dist, bty = "b2", colvar=final, zlab = "Miles gallon" ,clab = c("Multiobjective", "Sampling Optimisation"))
scatter3D(x = alfa*ndvi[index], y = beta*ign[index], z = gamma*dist[index], add = TRUE, colkey = FALSE, 
          pch = 18, cex = 3, col = "black")



plot(out1$`Aggregated matrix`$Try, final, type="o", col="black", pch="o", lty=1, ylim=c(0,3) )
lines(out1$`Aggregated matrix`$Try, final, col="black", pch="*")
#abline(h = max(final), col="red", lwd=3, lty=2)
#abline(V = index, col="red", lwd=3, lty=2)
points(index,max(final), col="red", cex=1)




