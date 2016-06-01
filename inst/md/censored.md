### Censored Regression

BGLR supports right, left and interval censoring. For censored outcomes the response is represented with a triplet ai,yi,bi, where ai and bi are the lower and upper bounds for the phenotype (yi). The following table describes the configuration of the tripled for un-censored, right, left and interval censored data

Case               | ai  |  yi  | bi
-------------------|-----|------|----
 Not-censored      | NA  |  yi  | NA
 Right-Censored    | ai  |  NA  | NA
 Left-Censored     | NA  |  NA  | bi 
 Interval-Censored | ai  |  NA  | bi

