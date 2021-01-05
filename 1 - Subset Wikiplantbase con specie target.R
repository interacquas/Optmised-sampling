wiki <- read.csv("Wikiplantbase/Wiki_agosto2020_ignorance.csv")

library(readxl)
target <- read_excel("Wikiplantbase/Target species.xlsx", col_names = TRUE)
head(target)

is.data.frame(wiki) #controllo che sia un dataframe

target <- target$lista_aggancio_nomiwiki #creo il vettore per fare il subset
is.vector(target)


wiki_new <- wiki[wiki$Taxon %in% target, ] # faccio il subset
wiki_final <- wiki_new[, 2:ncol(wiki_new)]
write.csv(wiki_final, "Wikiplantbase/wiki_final.csv", row.names = FALSE) # salvo il file wiki_final, cioè il dataframe da dare in pasto alla funzione di ignoranza 

wiki_list <- unique(wiki_final$Taxon) # creo la lista delle specie presenti nel dataframe

wiki_list <- sort(wiki_list)

setdiff(wiki_list, target) # le specie che sono nel dataframe ma non nel vettore: ovviamente il risultato deve essere un insieme vuoto

setdiff(target,wiki_list) # le specie che sono nella lista delle specie target ma non hanno fatto il match. Specie per specie c'è da capire perché

