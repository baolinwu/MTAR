# MTAR
  - an R package implementing various methods for genome-wide association test of multiple traits
  - installation within R: install_github("baolinwu/MTAR")  (assuming the "devtools" R package installed)

-----
## Fast and accurate genome-wide association test of multiple quantitative traits
  - Reference:  
     - Wu,B. and Pankow,J.S. (2017) Fast and accurate genome-wide association test of multiple continuous traits. *tech report*
  - Jointly test the association of a SNP with multiple continuous traits.
     - accurate calculation of P-values.
     - very efficient and extremely scalable to genome-wide association test.
  - Sample codes
```r
  library(VGAM)
  library(MTAR)
  Z = rbinom(1000,1,0.5)
  G = rbinom(1000,2,0.25)
  ##
  X = rnorm(1000)
  e = rnorm(1000)
  Y1 = Z+X + 0.15*G + rnorm(1000)+e
  Y2 = Z+X + 0.1*G + rnorm(1000)+e
  Y = cbind(Y1,Y2)
  obj = MLM.null(Y,cbind(Z,X))
  MQTAc(obj, G)
  ##
  X1 = rnorm(1000)
  X2 = rnorm(1000)
  Y1 = Z+X1 + 0.15*G + rnorm(1000)+e
  Y2 = Z+X2 + 0.1*G + rnorm(1000)+e
  Y = cbind(Y1,Y2)
  YX = list(Y, cbind(X1,Z), cbind(X2,Z))
  objd = MDM.null(YX,pux=2)
  MQTAd(objd, G)
```
