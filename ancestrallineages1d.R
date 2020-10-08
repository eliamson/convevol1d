#Univariate version of ancestrallineages(), a function of the convevol package (Stayton, 2015).
#v0.2

#Stayton CT. 2015. The definition, recognition, and interpretation of convergent evolution, and two new measures for quantifying and assessing the significance of convergence. Evolution 69:2140–2153. doi:10.1111/evo.12729
#convevol package: https://CRAN.R-project.org/package=convevol

ancestrallineages1d<-function (phyl, phendata, t1, t2) 
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
    if (is.finite(t1)) {
    }
    else {
        t1 <- labelstonumbers(phyl, t1)
    }
    if (is.finite(t2)) {
    }
    else {
        t2 <- labelstonumbers(phyl, t2)
    }
    if (t1 > length(phyl$tip)) 
        stop("your first tip isn't in the phylogeny.")
    if (t2 > length(phyl$tip)) 
        stop("your second tip isn't in the phylogeny.")
    allancs <- fastAnc(phyl, phendata)
    alldata<-c(phendata,allancs)
    anctimes <- node.depth.edgelength(phyl)
    combineddata <- cbind(anctimes, alldata)
    mrcas <- mrca(phyl)
    mrcat1t2 <- mrcas[t1, t2]
    t1path <- combineddata[t1, ]
    anc <- findanc(phyl, t1)
    t1path <- rbind(t1path, combineddata[anc[1], ])
    while (anc[1] != mrcat1t2) {
        anc <- findanc(phyl, anc[1])
        t1path <- rbind(t1path, combineddata[anc[1], ])
    }
    t2path <- combineddata[t2, ]
    anc <- findanc(phyl, t2)
    t2path <- rbind(t2path, combineddata[anc[1], ])
    while (anc[1] != mrcat1t2) {
        anc <- findanc(phyl, anc[1])
        t2path <- rbind(t2path, combineddata[anc[1], ])
    }
    t1pathtraits <- t1path[, 2]
    t2pathtraits <- t2path[, 2]
    answer <- list(t1pathtraits, t2pathtraits)
}
