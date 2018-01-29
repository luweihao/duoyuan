### 多总体的判别分析 ###
## TrnX--训练样本
## TrnG--训练样本 分类因子 
## TstX---待判样本  

distinguish.distance=function(TrnX, TrnG, TstX=NULL, var.equal=FALSE)
{
   flag=0
   if (is.null(TstX)==TRUE){
     TstX=TrnX
     flag=1
   }
   if (is.vector(TstX)==TRUE)  TstX=t(as.matrix(TstX))
   else if (is.matrix(TstX)!=TRUE) TstX=as.matrix(TstX)   
   if (is.matrix(TrnX)!=TRUE) TrnX=as.matrix(TrnX)

   nx=nrow(TstX)
   blong=matrix(rep(0, nx), nrow=1, dimnames=list("blong", 1:nx))
   
   g=length(levels(TrnG))
   mu=matrix(0, nrow=g, ncol=ncol(TrnX))
   
   for (i in 1:g)  
   {
   mu[i,]=colMeans(TrnX[TrnG==i,]) 
   }
    
   D=matrix(0, nrow=g, ncol=nx)
   
   if (var.equal==TRUE)
   {
      for (i in 1:g)
         D[i,]= mahalanobis(TstX, mu[i,], var(TrnX))
   }else
   {
      for (i in 1:g)
         D[i,]= mahalanobis(TstX, mu[i,], var(TrnX[TrnG==i,]))
   }
   
   for (j in 1:nx)
   {
      dmin=Inf
      for (i in 1:g)
          if (D[i,j]<dmin)
          {
          dmin=D[i,j]; blong[j]=i
          }
   }
   if(flag){
     origin = TrnG
     table1 <- as.data.frame( t(rbind(origin, blong)) )
     tab1=xtabs(~ origin+blong, data = table1)
     print(tab1)
     tab2=prop.table(tab1)
     print(tab2, digits = 3)
     wrong <- which((as.numeric(origin)-as.numeric(blong))!=0)
     cat("The index of wrong data are:", wrong, "\n")
     #return(table1)
   }else{return(blong)}
}
