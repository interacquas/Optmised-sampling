# Carico manualmente il file 'di Wiki' Wiki_ago2020_ignorance.csv", e lo chiamo wiki

# Carico manualmente il file "Target species.xlsx" e lo chiamo 'target'

target <- Target_species
wiki <- Wiki_agosto2020_ignorance

is.data.frame(wiki) #controllo che sia un dataframe

target <- target$...1 #creo il vettore per fare il subset
is.vector(target)


wiki_new <- wiki[wiki$Taxon %in% target, ] # faccio il subset

wiki_list <- unique(wiki_new$Taxon) # creo la lista delle specie presenti nel dataframe


setdiff(wiki_list, target) # le specie che sono nel dataframe ma non nel vettore: ovviamente il risultato deve essere un insieme vuoto

setdiff(target,wiki_list) # le specie che sono nella lista delle specie target ma non hanno fatto il match. Specie per specie c'è da capire perché




