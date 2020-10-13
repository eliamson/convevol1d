#Convergence direction
#Based on the function convnum() (Stayton, 2015) of the convevol package, this function assesses the number of tips having increased or decreased their trait value (this function is for univariate data) when compared to the reconstructed value (as give by phytool's function fastAnc()) of the ancestral node or that of the root of the tree.
#In addition to the total number of increase(s) or decrease(s), it also tests whether the difference exceeds the 95% confidence interval (95CI) of the reconstructed node values (either ancestor or root).

#Abbreviations:
# i+ : strong increase (exceeding 95CI) 
# i: weak increase (not exceeding 95CI)
# d: weak decrease (not exceeding 95CI)
# d+: strong decrease (exceeding 95CI) 

#Reference: Stayton, C. T. 2015. The definition, recognition, and interpretation of convergent evolution, and two new measures for quantifying and assessing the significance of convergence. Evolution 69:2140–2153

#Tested with R version 3.6.2, phytools version 0.7-20, and convevol version 1.3

require(phytools)
require(convevol) 

#### Simulate a dataset ####

phyl_S<-pbtree(n=50)
phendata_S <- fastBM(phyl_S)
convtips_S <- names(phendata_S)[c(1,2)] #Two converging tips are chosen

#### Running the function with the simulated dataset ####

convDir(phyl_S, phendata_S, convtips_S) 

#### Sourcing the function ####

convDir<-function (phyl, phendata, convtips) 
{
  if (class(phyl) != "phylo") 
    stop("your tree must be class 'phylo.'")
  if (length(phendata) != length(phyl$tip)) 
    stop("your data matrix must have the same number of rows as tips in the tree.")
  
  
  nConv <- length(convtips) # Number on convergent tips
  allancs <- fastAnc(phyl, phendata,CI=T) # Estimation of all ancestral values
  rootNum <- nodepath(phyl)[[1]][1] # Locate the root
  rootVal <- c()
  rootVal["val"]<- allancs$ace[as.character(rootNum)];rootVal["CIsL"]<- allancs$CI95[as.character(rootNum),1]; rootVal["CIsU"]<- allancs$CI95[as.character(rootNum),2] # Concatenate root values in one array
  AncDir <- data.frame()
  RootDir <- data.frame()
    for (i in 1:nConv) { # Looking at each converging tip
      AncDir[i,"convtip"] <- RootDir[i,"convtip"] <- convtips[i]
      t1 <- labelstonumbers(phyl, convtips[i]) # Retrieve tip name
      anc <- findanc(phyl,t1) # Locate direct ancestral node
      tipVal<- phendata[convtips[i]] # Retrieve tip value
      ancVal <- c()
      ancVal["val"]<- allancs$ace[as.character(anc[1])];ancVal["CIsL"]<- allancs$CI95[as.character(anc[1]),1]; ancVal["CIsU"]<- allancs$CI95[as.character(anc[1]),2] # Concatenate ancestral node values in one array
      if (tipVal>ancVal["CIsU"]){ #Tip value strongly increase from ancestor
        AncDir[i,"Dir"] <- "i+"
      } 
      if (tipVal<ancVal["CIsU"] & tipVal>ancVal["val"]){ #Tip value weakly increase from ancestor
        AncDir[i,"Dir"] <- "i" 
      }
      if (tipVal<ancVal["val"]& tipVal>ancVal["CIsL"] ){ #Tip value weakly decrease from ancestor
        AncDir[i,"Dir"] <- "d"
      }
      if (tipVal<ancVal["CIsL"] ){ #Tip value strongly decrease from ancestor
        AncDir[i,"Dir"] <- "d+"
      }
      
      if (tipVal>rootVal["CIsU"]){ #Tip value strongly increase from root
        RootDir[i,"Dir"] <- "i+"  
      } 
      if (tipVal<rootVal["CIsU"] & tipVal>rootVal["val"]){ #Tip value weakly increase from root
        RootDir[i,"Dir"] <- "i"  
      }
      if (tipVal<rootVal["val"]& tipVal>rootVal["CIsL"] ){ #Tip value weakly decrease from root
        RootDir[i,"Dir"] <- "d"
      }
      if (tipVal<rootVal["CIsL"] ){ #Tip value strongly decrease from root
        RootDir[i,"Dir"] <- "d+"
      }
    }
    
  list(RootDir=RootDir,AncDir=AncDir) # Append results
}
