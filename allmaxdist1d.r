#Univariate version of allmaxdist(), a function of the convevol package (Stayton, 2015).
#v0.2

#Stayton CT. 2015. The definition, recognition, and interpretation of convergent evolution, and two new measures for quantifying and assessing the significance of convergence. Evolution 69:2140–2153. doi:10.1111/evo.12729
#convevol package: https://CRAN.R-project.org/package=convevol

allmaxdist1d<-function (phyl, phendata, mat = TRUE) 
{
    require(convevol)
    
    if (class(phyl) != "phylo") 
        stop("your tree must be class 'phylo.'")
    if (nrow(phendata) != length(phyl$tip)) 
        stop("your data matrix must have the same number of rows as tips in the tree.")
    dims <- dim(phyl$edge)
    ntips <- (dims[1]/2) + 1
    allmaxes <- matrix(0, ntips, ntips)
    maxeslist <- c()
    t1 <- 1
    while (t1 < ntips) {
        t2 <- t1 + 1
        while (t2 <= ntips) {
            m <- maxdist(phyl, phendata, t1, t2)
            allmaxes[t1, t2] <- m
            maxeslist <- c(maxeslist, m)
            t2 <- t2 + 1
        }
        t1 <- t1 + 1
    }
    if (mat == TRUE) {
        answer <- allmaxes
    }
    else {
        answer <- maxeslist
    }
}
