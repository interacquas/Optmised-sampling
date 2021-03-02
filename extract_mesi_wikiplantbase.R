library(sp)
library(ggplot2)

segntos <- read.csv(file = "C:/Users/Giuseppe Antonelli/Desktop/tesi/segntos.csv")

date <- segntos[,32]

date <- as.Date(date, format = '%Y%m%d')
date <- na.omit(date)
mesi <- format(date, '%m')
mesi <- as.numeric(mesi)

hist(mesi, main="Segnalazioni per mese")
