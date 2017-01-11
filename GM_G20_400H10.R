
install.packages("glmnet")
library("glmnet")
install.packages("dplyr")
library("dplyr")
install.packages("ggplot2")
library("ggplot2")
install.packages("broom")
library("broom")
install.packages("ggfortify")
library(ggfortify)
install.packages("grid")
library(grid)
install.packages("gridExtra")
library(gridExtra)
install.packages("penalizedSVM")
library(penalizedSVM)
install.packages("BayesLogit")
library(BayesLogit)

ptmGM <- proc.time()
data <- read.table("G2M1000S0400H10.txt" ,header=TRUE)
num=1000
xdata = data[,1:num]
xdata[1,2*num]
x<- array(0,dim=c(400,num*2))
y2 = data$CLASS
for (i in 1:num) {
for (j in 1:400){
  k=(i*2)
  l=(i*2)-1
  if(xdata[j,i]==0) {
    x[j,l] <- 1
    x[j,k] <- -0.5
  }
  else if (xdata[j,i]==1) 
  { x[j,l] <- 0
    x[j,k] <- 0.5
  } 
  else if (xdata[j,i]==2) 
  { x[j,l] <- -1
    x[j,k] <- -0.5
  } 
 }
}

x= data.frame(x)
x2= as.matrix(x)

n <- nrow(data)
p <- ncol(data)-1
test.ratio <- .2 # ratio of test(20%)/train(80%) samples
n.test <- ceiling(n*test.ratio)
testindex <- sample(1:n,n.test)
trainindex <- setdiff(1:n,testindex)

data.train <- data[trainindex,]
data.test <- data[testindex,]

####### create polinomial interaction ######
data2 = data.frame(x2 , y = y2)
xInter = model.matrix(y~(.)^2,data2)
ccut<- array(0,dim=c(1,num+1))
 for(i in 0:num){
   if(i==0)
   {
     ccut[i+1]<-1
     ccut[i+2]<-ccut[i+1]+(2*num)+1 
   }
   else if(num>i)
   {
     ccut[i+2]<-ccut[i+1]+(4*(num-i))+1 
   }
   else if(num==i)
   {
     break
   }
 }

xInter2 = xInter[,-ccut]
TimeGM<-proc.time()-ptmGM


ptmLasso <- proc.time()
Glasso<-cv.glmnet(xInter2[trainindex,] , y2[trainindex] ,alpha=1, family = "binomial")

a=coef(Glasso,Glasso$lambda.1se)
Gnylasso1<-predict(Glasso, s=Glasso$lambda.1se, newx=xInter2[testindex,],type = "class")
table(Gnylasso1 , data.test$CLASS)
mean(Gnylasso1 == data.test$CLASS)
TimeLasso<-proc.time()-ptmLasso

ptmElas <- proc.time()
Gelnet<-cv.glmnet(xInter2[trainindex,] , y2[trainindex] ,alpha=0.5, family = "binomial")
b=coef(Gelnet,Gelnet$lambda.1se)
GnyElastic1<-predict(Gelnet, s=Gelnet$lambda.1se, newx=xInter2[testindex,],type = "class")
table(GnyElastic1 , data.test$CLASS)
mean(GnyElastic1 == data.test$CLASS)
TimeElas<-proc.time()-ptmElas

ptmRid <- proc.time()
Gridge<-cv.glmnet(xInter2[trainindex,] , y2[trainindex] ,alpha=0, family = "binomial")
c=coef(Gridge,Gridge$lambda.1se)
GnyRidge1<-predict(Gridge, s=Gridge$lambda.1se, newx=xInter2[testindex,],type = "class")
table(GnyRidge1 , data.test$CLASS)
mean(GnyRidge1 == data.test$CLASS)
Timerid<-proc.time()-ptmRid 




sink('output/GM1000_G2M1000S0400H10_s.txt')
cat("\nridge =\n")
c
mean(GnyRidge1 == data.test$CLASS)
cat("time = \n")
Timerid + TimeGM
sink()

sink('output/GM1000_G2M1000S0400H10_s.txt', append = TRUE)
cat("\nelastic =\n")
b
mean(GnyElastic1 == data.test$CLASS)
cat("time = \n")
TimeElas + TimeGM
sink()

sink('output/GM1000_G2M1000S0400H10_s.txt', append = TRUE)
cat("\nlasso =\n")
a
mean(Gnylasso1 == data.test$CLASS)
cat("time = \n")
TimeLasso + TimeGM
sink()

#####################scad l2########################
ptmScad <- proc.time()
yt=ifelse(y2[trainindex]==0,-1,1)
system.time( scad<- svm.fs(xInter2[trainindex,], y=yt, fs.method="scad+L2",cross.outer= 0, grid.search = "interval", maxIter = 10,inner.val.method = "cv", cross.inner= 5, maxevals=500, parms.coding = "log2", show="none", verbose=FALSE ) )
ytest=ifelse(y2[testindex]==0,-1,1)
test.error.scad<-predict(scad, newdata=xInter2[testindex,],newdata.labels=ytest)
testn=ifelse(test.error.scad$fitted<0,-1,1)
mean(testn == ytest)
TimeScad <- proc.time()-ptmScad



sink('output/GM1000_G2M1000S0400H10_s.txt', append = TRUE)
cat("\nScad L2 =\n")
scad$model$w
mean(testn == ytest)
cat("time = \n")
TimeScad + TimeGM
sink()

