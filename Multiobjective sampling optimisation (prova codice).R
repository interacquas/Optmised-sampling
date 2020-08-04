out1$`Aggregated matrix`


normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}



ndvi <- normalize(out1$`Aggregated matrix`$Variance)
ign <- normalize(out1$`Aggregated matrix`$`Mean Ignorance`)
dist <- normalize(out1$`Aggregated matrix`$`Mean Dist`)
final <- ndvi * ign * dist

plot(out1$`Aggregated matrix`$Try, ndvi, type="o", col="blue", pch="o", lty=1, ylim=c(0,1) )
lines(out1$`Aggregated matrix`$Try, ndvi, col="blue", pch="*")


points(out1$`Aggregated matrix`$Try, ign, col="red", pch="*")
lines(out1$`Aggregated matrix`$Try, ign, col="red", pch="*")

points(out1$`Aggregated matrix`$Try, dist, col="green",pch="+")
lines(out1$`Aggregated matrix`$Try, dist, col="green", lty=3)


plot(out1$`Aggregated matrix`$Try, final, type="o", col="black", pch="o", lty=1, ylim=c(0,1) )
lines(out1$`Aggregated matrix`$Try, final, col="black", pch="*")
abline(h = max(final), col="red", lwd=3, lty=2)


