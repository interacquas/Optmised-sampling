library(sp)
library(ggplot2)

segntos <- read.csv(file = "C:/Users/Giuseppe Antonelli/Desktop/tesi/segntos.csv")

date <- segntos[,32]

date <- as.Date(date, format = '%Y%m%d')
date <- na.omit(date)
mesi <- format(date, '%m')
mesi <- as.numeric(mesi)

segnalaz <- table(mesi)

barplot(segnalaz, main="Segnalazioni per mese", ylab="segnalazioni", xlab="mesi",
        names.arg=c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"), cex.names=1.8)
     