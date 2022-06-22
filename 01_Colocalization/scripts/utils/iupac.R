iupac.map <- matrix(
    c(
        "A", "M", "R", "W",
        "M", "C", "S", "Y",
        "R", "S", "G", "K",
        "W", "Y", "K", "T"
    ), 
    nrow=4
)
rownames(iupac.map) <- diag(iupac.map)
colnames(iupac.map) <- diag(iupac.map)

alleles.iupac <- Vectorize(function(x, y) {
 
    return(iupac.map[x, y])
})