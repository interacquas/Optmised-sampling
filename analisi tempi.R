library(readxl)

foglio2 <- read_excel("tempi-raggiungimento.xlsx", 
                     sheet = "target")

shapiro.test(foglio2$t)
shapiro.test(foglio2$n)

cor.test(foglio2$t, foglio2$n)

tsun <- lm(n~t, data=foglio2)
summary(tsun)
