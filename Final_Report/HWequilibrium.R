
require(HardyWeinberg)

#Link:
#https://cran.r-project.org/web/packages/HardyWeinberg/HardyWeinberg.pdf



# HWChisqMat executes the Chisquare test for HWE for each row in a matrix.
HWChisqMat(X, ...)
X <- HWData(100,10)
colnames(X) <- c("MM","MN","NN")
Results <- HWChisqMat(X)
Output <- cbind(X,Results$chisqvec,Results$pvalvec)
print(Output)


# Test for an autosomal blood group marker
# HWChisq performs the chi-square test for Hardy-Weinberg equilibrium both for autosomal and Xchromosomal markers.
x <- c(MM=120919,MN=156574,NN=50775)
HW.test <- HWChisq(x,verbose=TRUE)

x <- c(MM=3375,MN=3958,NN=1231)
HW.test <- HWChisq(x,verbose=TRUE)

x <- c(MM=120919+3375,MN=156574+3958,
       NN=50775+1231)
HW.test <- HWChisq(x,verbose=TRUE)
