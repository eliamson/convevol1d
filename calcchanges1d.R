#Univariate version of calcchanges(), a function of the convevol package (Stayton, 2015).
#v0.2

#Stayton CT. 2015. The definition, recognition, and interpretation of convergent evolution, and two new measures for quantifying and assessing the significance of convergence. Evolution 69:2140–2153. doi:10.1111/evo.12729
#convevol package: https://CRAN.R-project.org/package=convevol

calcchanges1d<-function (phyl, phendata) 
{
    require(convevol)
    
    if (class(phyl) != "phylo") 
        stop("your tree must be class 'phylo.'")
    if (length(phendata) != length(phyl$tip)) 
        stop("your data matrix must have the same number of rows as tips in the tree.")
    if (is.null(names(phendata))) {
        warning("no row names for data.  Assuming that the rows are in the same order as tips.")
        #names(X) <- phyl$tip.label
    }
    allancs <- fastAnc(phyl, phendata)
    allvals<-c(phendata,allancs)

    pdims <- dim(phyl$edge)
    changes <- rep(0, pdims[1])
    j <- 1
    for (j in 1:pdims[1]) {
        ancnode <- phyl$edge[j, 1]
        desnode <- phyl$edge[j, 2]
        ancval <- allvals[ancnode]
        desval <- allvals[desnode]
        change <- abs(ancval - desval)
        changes[j] <- change
    }
    answer <- changes
}
