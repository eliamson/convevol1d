#Univariate version of convratsig(), a function of the convevol package (Stayton, 2015).
#v0.2

#Bug fix 23.03.20: call convrat1dMod(), see bug fix there.

#Stayton CT. 2015. The definition, recognition, and interpretation of convergent evolution, and two new measures for quantifying and assessing the significance of convergence. Evolution 69:2140–2153. doi:10.1111/evo.12729
#convevol package: https://CRAN.R-project.org/package=convevol

convratsig1dMod<-function (phyl, phendata, convtips, nsim) 
{
    require(convevol)
    
    if (class(phyl) != "phylo") 
        stop("your tree must be class 'phylo.'")
    if (length(phendata) != length(phyl$tip)) 
        stop("your data matrix must have the same number of rows as tips in the tree.")
    obs <- convrat1dMod(phyl, phendata, convtips)
    pvar<-ratematrix(phyl,phendata)
    rootval<-fastAnc(phyl,phendata)[1]
    k <- 1
    C1s <- c()
    C2s <- c()
    C3s <- c()
    C4s <- c()
    while (k <= nsim) {
        simphendata <- sim.char(phyl, pvar, 1, model = c("BM"), root = rootval)
        rsimphendata <- simphendata[, , 1]
        simresults <- convrat1dMod(phyl, rsimphendata, convtips)
        C1s <- c(C1s, simresults[1])
        C2s <- c(C2s, simresults[2])
        C3s <- c(C3s, simresults[3])
        C4s <- c(C4s, simresults[4])
        k <- k + 1
    }
    C1greater <- 0
    C2greater <- 0
    C3greater <- 0
    C4greater <- 0
    for (i in 1:nsim) {
        if (C1s[i] >= obs[1]) {
            C1greater <- C1greater + 1
        }
        if (C2s[i] >= obs[2]) {
            C2greater <- C2greater + 1
        }
        if (C3s[i] >= obs[3]) {
            C3greater <- C3greater + 1
        }
        if (C4s[i] >= obs[4]) {
            C4greater <- C4greater + 1
        }
    }
    C1cutoff <- sort(C1s)[round(nsim * 0.95)]
    C2cutoff <- sort(C2s)[round(nsim * 0.95)]
    C3cutoff <- sort(C3s)[round(nsim * 0.95)]
    C4cutoff <- sort(C4s)[round(nsim * 0.95)]
    ObservedCs <- obs
    names(ObservedCs) <- c("C1", "C2", "C3", "C4")
    Pvals <- c(C1greater/(nsim + 1), C2greater/(nsim + 1), C3greater/(nsim + 
        1), C4greater/(nsim + 1))
    names(Pvals) <- c("C1", "C2", "C3", "C4")
    Cutoffs <- c(C1cutoff, C2cutoff, C3cutoff, C4cutoff)
    names(Cutoffs) <- c("C1", "C2", "C3", "C4")
    AllSimCs <- cbind(C1s, C2s, C3s, C4s)
    return(list(ObservedCs = ObservedCs, Pvals = Pvals, AllSimCs = AllSimCs, 
        Cutoffs = Cutoffs))
}
