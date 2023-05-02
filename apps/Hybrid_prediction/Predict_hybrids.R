rm(list=ls())

library(BGLR)

#setwd("~/Documents/Documentos Paulino/CIMMyT/Advanta-2019/3. Summary GS and GBLUP/Part3/examples/")

load("cornHybrid.RData")

ls()

str(cornHybrid)
names(cornHybrid)

A=cornHybrid$K
dim(A)
colnames(A)

tmp1=cornHybrid$hybrid
str(tmp1)
#write.csv(tmp1,file="h1.csv")

colnames(tmp1)

y=tmp1$Yield
Loc=tmp1$Location
boxplot(y~Loc)

class(tmp1)
colnames(tmp1)
table(tmp1$Location)
boxplot(tmp1$Yield~tmp1$Location)

#Male
unique(tmp1$GCA1)
levels(tmp1$GCA1)
nlevels(tmp1$GCA1)

#Female
unique(tmp1$GCA2)
levels(tmp1$GCA2)
nlevels(tmp1$GCA2)

#Hibs
unique(tmp1$SCA)
levels(tmp1$SCA)
nlevels(tmp1$SCA)

ZE=model.matrix(~Location-1, data=tmp1)
dim(ZE)
write.csv(ZE,file="ZE.csv")
Z1=model.matrix(~GCA1-1, data=tmp1)
dim(Z1)
#write.csv(Z1,file="Z1.csv")
Z2=model.matrix(~GCA2-1, data=tmp1)

Z3=model.matrix(~SCA-1, data=tmp1)
dim(Z3)

dim(A)
rownames(A)

K1=A[levels(tmp1$GCA1), levels(tmp1$GCA1)]
dim(K1)
rownames(K1)
K2=A[levels(tmp1$GCA2), levels(tmp1$GCA2)]
dim(K2)
rownames(K2)
K3=kronecker(K1,K2)
dim(K3)

K1star=Z1%*%K1%*%t(Z1)
K2star=Z2%*%K2%*%t(Z2)
K3star=Z3%*%K3%*%t(Z3)

ETA=list(list(X=ZE,model="FIXED"),
	     list(K=K1star,model="RKHS"),
         list(K=K2star,model="RKHS"),
         list(K=K3star,model="RKHS"))

fm=BGLR(y=y,ETA=ETA,nIter=10000,burnIn = 5000,thin=10,
        verbose=FALSE)

plot(fm$y,fm$yHat)
plot(fm)
summary(fm)

unlink("*.dat")

#Results
results=data.frame(tmp1,yieldHat=fm$yHat)
write.csv(results,file="predictions_hybrids.csv")

cor(fm$y,fm$yHat,use="complete.obs")

fm$ETA[[2]]$varU #Male
fm$ETA[[3]]$varU #Female
fm$ETA[[4]]$varU #Hybrid

#Tamanho pob. prueba
sum(is.na(y))

#Tamanho pob. entrenamiento
sum(!is.na(y))

#Correlación entre observados y predichos
entrenamiento=(!is.na(y))
prueba=is.na(y)

y_entrenamiento=y[entrenamiento]
y_entrenamiento_hat=fm$yHat[entrenamiento]
Localidad_entrenamiento=tmp1$Location[entrenamiento]

plot(y_entrenamiento,y_entrenamiento_hat)
cor(y_entrenamiento,y_entrenamiento_hat)

result=data.frame(Loc=Localidad_entrenamiento,
                  y=y_entrenamiento,
                  yHat=y_entrenamiento_hat)
                  
#Componente de varianza para padre
fm$ETA[[2]]$varU

#Componente de varianza para mama
fm$ETA[[3]]$varU

#Componente de varianza para hijo
fm$ETA[[4]]$varU

#Error
fm$varE

predicciones=data.frame(Loc=Loc, yObs=y,yPred=fm$yHat,tmp1$SCA)
#write.csv(predicciones,file="~/Desktop/predicciones.csv")

#BLUPs
#Padre
plot(fm$ETA[[2]]$u)

#Madre
fm$ETA[[3]]$u

#Cruza
fm$ETA[[4]]$u
