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
 
    if (x %in% rownames(iupac.map) && y %in% colnames(iupac.map)) {
        return(iupac.map[x, y])
    } else {
        return(NA)
    }
})