#convnum() function of the convevol package (Stayton, 2015) altered to fix Bug (06.05.20): use predict to plot ellipsoid. 
#v0.2

#Stayton CT. 2015. The definition, recognition, and interpretation of convergent evolution, and two new measures for quantifying and assessing the significance of convergence. Evolution 69:2140–2153. doi:10.1111/evo.12729
#convevol package: https://CRAN.R-project.org/package=convevol

convnumMod <- function (phyl, phendata, convtips, plot = TRUE, ellipse = NULL, 
          plotellipse = NULL) 
{
  require(convevol)
  
  if (class(phyl) != "phylo") 
    stop("your tree must be class 'phylo.'")
  if (nrow(phendata) != length(phyl$tip)) 
    stop("your data matrix must have the same number of rows as tips in the tree.")
  if (is.list(convtips) == TRUE) {
    convtips <- unlist(convtips)
  }
  if (length(convtips) <= ncol(phendata)) 
    stop("You must have fewer variables than putatively convergent taxa")
  phendata <- as.matrix(phendata)
  phendata <- phendata[phyl$tip.label, ]
  convtaxa <- phendata[convtips, ]
  if (is.null(ellipse)) {
    convell <- ellipsoidhull(convtaxa, 1e-04)
  }
  else {
    convell <- ellipse
  }
  if (is.null(plotellipse)) {
    plotell <- ellipsoidhull(convtaxa[, 1:2], 1e-04)
  }
  else {
    plotell <- plotellipse
  }
  allancstates <- apply(phendata, 2, fastAnc, tree = phyl)
  alldata <- rbind(phendata, allancstates)
  if (plot == TRUE) {
    phyl$node.label <- NULL
    phylomorphospace(phyl, phendata[, 1:2], label = "off")
    points(phendata[unlist(convtips), 1:2], col = "white")
    #plotellipse(plotell)
    lines(predict(convell), col="red")
  }
  cross <- 0
  nbran <- dim(phyl$edge)
  crossedges <- c()
  i <- 1
  while (i <= nbran[1]) {
    isinanc = FALSE
    isindes = FALSE
    anc <- phyl$edge[i, 1]
    des <- phyl$edge[i, 2]
    ancval <- alldata[anc, ]
    desval <- alldata[des, ]
    ancdev <- ancval - convell$loc
    desdev <- desval - convell$loc
    cutoff <- 1.05 * convell$d2
    if (t(ancdev) %*% ginv(convell$cov) %*% ancdev <= cutoff) {
      isinanc = TRUE
    }
    else {
      isinanc = FALSE
    }
    if (t(desdev) %*% ginv(convell$cov) %*% desdev <= cutoff) {
      isindes = TRUE
    }
    else {
      isindes = FALSE
    }
    if (isinanc != TRUE & isindes == TRUE) {
      cross <- cross + 1
      crossedges <- c(crossedges, i)
      if (plot == TRUE) {
        arrows(ancval[1], ancval[2], desval[1], desval[2], 
               length = 0.1, angle = 30, code = 2, col = "red")
      }
    }
    i <- i + 1
  }
  return(list(cross, convell, plotell))
  crossedges
}  