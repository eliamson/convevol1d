#Univariate version of convrat(), a function of the convevol package (Stayton, 2015).
#v0.2
#Bug fix 23.03.20: C1 and C2 can't be negative, even when  ntax > 2

#Stayton CT. 2015. The definition, recognition, and interpretation of convergent evolution, and two new measures for quantifying and assessing the significance of convergence. Evolution 69:2140–2153. doi:10.1111/evo.12729
#convevol package: https://CRAN.R-project.org/package=convevol


convrat1dMod<-function (phyl, phendata, convtips) 
{
  require(convevol)
  
    if (class(phyl) != "phylo") 
        stop("your tree must be class 'phylo.'")
    if (length(phendata) != length(phyl$tip)) 
        stop("your data matrix must have the same number of rows as tips in the tree.")
    ntax <- length(convtips)
    if (ntax == 2) {
        tipsdist <- abs(phendata[convtips[1]] - phendata[convtips[2]])
        mxdist <- maxdist1d(phyl, phendata, convtips[1], convtips[2])
        if (tipsdist >= mxdist) {
            C1 <- 0
            C2 <- 0
        }
        if (mxdist > tipsdist) {
            C1 <- 1 - (tipsdist/mxdist)
            C2 <- mxdist - tipsdist
        }
        lineagepaths <- ancestrallineages1d(phyl, phendata, convtips[1], convtips[2])
        t1totalevolution <- 0
        t2totalevolution <- 0

	  for (i in 2:length(lineagepaths[[1]])) {
		t1totalevolution<-t1totalevolution + abs(lineagepaths[[1]][i]-lineagepaths[[1]][i-1])
        	}

 	  for (i in 2:length(lineagepaths[[2]])) {
		t2totalevolution<-t2totalevolution + abs(lineagepaths[[2]][i]-lineagepaths[[2]][i-1])
        	}


        totallineageevolution <- t1totalevolution + t2totalevolution
        C3 <- C2/totallineageevolution
        cMRCA <- findMRCA(phyl, tips = convtips, type = "node")
        subphyl <- extract.clade(phyl, cMRCA)
        wholephylchanges <- sum(calcchanges1d(subphyl, phendata[subphyl$tip.label 
            ]))
        C4 <- C2/wholephylchanges
    }
    if (ntax > 2) {
        C1s <- c()
        C2s <- c()
        C3s <- c()
        C4s <- c()
        for (i in 1:ntax) {
            j <- i + 1
            while (j <= ntax) {
                tipsdist <- abs(phendata[convtips[i]] - phendata[convtips[j]])
                mxdist <- maxdist1d(phyl, phendata, convtips[i], convtips[j])

                if (tipsdist >= mxdist) {
                  C1 <- 0
                  C2 <- 0
                }
                if (mxdist > tipsdist) {
                  C1 <- 1 - (tipsdist/mxdist)
                  C2 <- mxdist - tipsdist
                }
                
                lineagepaths <- ancestrallineages1d(phyl, phendata, convtips[i], convtips[j])
                t1totalevolution <- 0
                t2totalevolution <- 0
               
	          for (ii in 2:length(lineagepaths[[1]])) {
		        t1totalevolution<-t1totalevolution + abs(lineagepaths[[1]][ii]-lineagepaths[[1]][ii-1])
        	        }

 	          for (ii in 2:length(lineagepaths[[2]])) {
		        t2totalevolution<-t2totalevolution + abs(lineagepaths[[2]][ii]-lineagepaths[[2]][ii-1])
                    }

                totallineageevolution <- t1totalevolution + t2totalevolution
                C3 <- C2/totallineageevolution
                cMRCA <- findMRCA(phyl, tips = convtips, type = "node")
                subphyl <- extract.clade(phyl, cMRCA)
                wholephylchanges <- sum(calcchanges1d(subphyl, phendata[subphyl$tip.label]))
                wholephylchanges <- sum(calcchanges1d(phyl, phendata))
                C4 <- C2/wholephylchanges
                C1s <- c(C1s, C1)
                C2s <- c(C2s, C2)
                C3s <- c(C3s, C3)
                C4s <- c(C4s, C4)
                j <- j + 1
            }
        }
        C1 <- mean(C1s)
        C2 <- mean(C2s)
        C3 <- mean(C3s)
        C4 <- mean(C4s)
    }
    answer <- c(C1, C2, C3, C4)
    names(answer) <- c("C1", "C2", "C3", "C4")
    answer
}
