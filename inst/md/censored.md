### Censored Regression

BGLR supports right, left and interval censoring. For censored outcomes the response is represented with a triplet $a_i, y_i, b_i$, where $a_i$ and $b_i$ are the lower and upper bounds for the phenotype $(y_i)$. The following table describes the configuration of the tripled for un-censored, right, left and interval censored data

Case               | $a_i$  |  $y_i$  | $b_i$
-------------------|-----|------|----
 Not-censored      | NA  |  $y_i$  | NA
 Right-Censored    | $a_i$  |  NA  | Inf
 Left-Censored     |-Inf |  NA  | $b_i$ 
 Interval-Censored | $a_i$  |  NA  | $b_i$


```R
 library(BGLR); data(wheat)
 y=wheat.Y[,1]
 
 # Introducing censoring (for simplicity done at fixed points -1 and 1)
 yNA=y
 a=rep(NA,599)
 b=rep(NA,599)
 type=sample(1:4,prob=c(2/3,1/9,1/9,1/9),size=599,replace=T)
# left-censored
 a[type==2]=-Inf
 b[type==2]=1
 yNA[type==2]=NA


# right-censored
 a[type==2]=-1
 b[type==2]=Inf
 yNA[type==2]=NA
 
# right-censored
 a[type==2]=-1
 b[type==2]=1
 yNA[type==2]=NA 
 
 fm=BGLR(y=yNA,a=a,b=b)
 plot(y,fm$yHat)

```

[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
