####  正态分布参数检验 ######### 

#### 单总体 正态分布均值检验  ######### 

mu.test.known=function(data, mu0, Sigma0, alpha=0.05)   ## ith sample is put in the ith line in data
{
data=as.matrix(data)
n=nrow(data)
p=ncol(data)

X.bar=apply(data, 2, mean)
T1=n*t(X.bar-mu0)%*%solve(Sigma0)%*%(X.bar-mu0)

a2=qchisq(1-alpha, p)

reject=matrix(c(T1, a2), nrow=1)
rownames(reject)=c("Reject")
colnames(reject)=c("Obs", ">1-alpha")

pv=1-pchisq(T1, p)
return(list(Reject.area=reject, p.value=pv))
}


mu.test=function(data, mu0)   ## ith sample is put in the ith line in data
{
data=as.matrix(data)
n=nrow(data)
p=ncol(data)

X.bar=apply(data, 2, mean)
A=(n-1)*var(data)

T2=(n-1)*n*t(X.bar-mu0)%*%solve(A)%*%(X.bar-mu0)
F=(n-p)/((n-1)*p)*T2

p.two=1-pf(F, p, n-p)
return(list(p.value=p.two))
}


#### 多总体 正态分布均值检验  ######### 

#### 两总体等方差 均值检验  ######### 
two.mu.test=function(data1, data2)   ## ith sample is put in the ith line in data
{
data1=as.matrix(data1)
data2=as.matrix(data2)
n1=nrow(data1)
n2=nrow(data2)
p=ncol(data1)

X.bar=apply(data1, 2, mean) 
A1=(n1-1)*var(data1)
Y.bar=apply(data2, 2, mean)
A2=(n2-1)*var(data2) 
A=(A1+A2)/(n1+n2-2)

T2=(n1*n2/(n1+n2))*t(X.bar-Y.bar)%*%solve(A)%*%(X.bar-Y.bar)
F=(n1+n2-2-p+1)/((n1+n2-2)*p)*T2

p.two=1-pf(F, p, (n1+n2-p-1))
return(list(p.value=p.two))
}

#### 多总体等方差 均值检验 --多元方差分析 ######### 
multi.mu.test=function(data, k)            ## data 
{
ind=data$ind

n=nrow(data)
p=ncol(data)-1

data=data[ ,1:p]
T=(n-1)*var(data)
  
A=0
for (i in 1:k)                                
{
datai=data[ind==i, ]
ni=nrow(datai)                                 
A=A+(ni-1)*var(datai)
}

Lambda=det(A)/det(T)
n1=n-k
n2=k-1
r=n1-(p-n2+1)/2
Chi=(-1)*r*log(Lambda)

p.value=1-pchisq(Chi, p*n2)
return(list(p.value=p.value))
}


#### 单总体 正态分布方差检验  ######### 
var.test=function(data, Sigma0)
{
n=nrow(data)
p=ncol(data)
A=(n-1)*var(data)
S=A%*%solve(Sigma0)
W=p^p*det(S)/(sum(diag(S)))^p
T5=-(n-1-(2*p^2+p+2)/(6*p))*log(W)
p.value=1-pchisq(T5, p*(p+1)/2-1)
return(list(p.value=p.value))
}

#### 多总体 正态分布方差检验  ######### 
multi.var.test=function(data, k)
{
ind=data$ind

n=nrow(data)
p=ncol(data)-1
data=data[ ,1:p]

A=0
for (i in 1:k)                                
{
datai=data[ind==i, ]
ni=nrow(datai)                                 
A=A+(ni-1)*var(datai)
}
  
det.A=0
for (i in 1:k)                                
{
datai=data[ind==i, ]
ni=nrow(datai)                                 
det.A=det.A+(ni-1)*log(det(var(datai)))
}

M=(n-k)*log(det(A/(n-k)))-det.A
d=(2*p^2+3*p-1)*(k+1)/(6*(p+1)*(n-k))
f=p*(p+1)*(k-1)/2

T6=(1-d)*M

p.value=1-pchisq(T6, f)
return(list(p.value=p.value))
}


### 多总体 正态分布均值方差同时检验  ######### 
multi.mean.var.test=function(mydata, k)
{
  ind=mydata$ind
  
  n=nrow(mydata)
  p=ncol(mydata)-1
  mydata=mydata[ ,1:p]
  bi=0
  
  A=1
  for (i in 1:k){
    mydatai=mydata[ind==i, ]
    ni=nrow(mydatai)
    bi=bi+1/(ni-1)
    A=A*(det((ni-1)*var(mydatai))/(ni-1)^p)^((ni-1)/2)
  }
  TT=(det((n-1)*var(mydata))/((n-k)^p))^((n-k)/2)
  
  b=(bi-1/(n-k))*((2*p^2+3*p-1)/(6*(p+3)*(k-1)))-(p-k+2)/((n-k)*(p+3))
  f=p*(p+3)*(k-1)/2
  T5=-2*(1-b)*log(A/TT)
  
  p.value=pchisq(T5, f, lower.tail = FALSE)
  return(list(p.value=p.value))
}


#### 多元正态 独立性检验  ######### 
norm.independent.test=function(mydata, subvector, k)
{
  n=nrow(mydata)
  p=ncol(mydata)
  pa = rep(0, k)
  
  Aii=1
  A=(n-1)*var(mydata)
  for (i in 1:k){
    datai=subset(mydata, select = (subvector==i))
    pa[i]=ncol(datai)
    Aii=Aii*det(as.matrix(A[subvector==i, subvector==i]))
  }
  
  
  V=det(A)/Aii
  b=n-3/2-(p^3-sum(pa^3))/(3*(p^2-sum(pa^2)))
  f=0.5*(p*(p+1)-sum(pa*(pa+1)))
  
  T7=-b*log(V)
  
  p.value=pchisq(T7, f, lower.tail = FALSE)
  return(list(t=T7, p.value=p.value))
}