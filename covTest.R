library("glmnet")
library("dplyr")
library("ggplot2")
library("broom")
library(ggfortify)
library(grid)
library(gridExtra)
library(penalizedSVM)
library(BayesLogit)
library("covTest")

ptmGM <- proc.time()
data <- read.table("Data#1_3-3_1-1_20_1600_1.txt" ,header=TRUE )
num=10
row=nrow(data)
xdata = data[,1:num]
xdata[1,2*num]
x<- array(0,dim=c(row,num*2))
y2 = data$Class
for (i in 1:num) {
  for (j in 1:row){
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
Covlasso <- lars.glm(xInter2[trainindex,] , y2[trainindex] , family = "binomial")
cov <- covTest(Covlasso,xInter2[trainindex,] , y2[trainindex])


sink('output/cov10_m1_s16.txt')
Covlasso 
cov
sink()

