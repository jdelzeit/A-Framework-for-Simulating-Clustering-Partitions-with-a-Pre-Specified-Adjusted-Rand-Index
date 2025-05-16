library(MASS)
library(mclust)
library(knitr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(kableExtra)
library(stringr)


mcreps = 500

n.vec = c(seq(from = 50, to = 500, by = 50), seq(from =600, to = 1000, by = 100))

ari.vec = seqari = seq(from= 0.1, to = 0.9, by = 0.1)

alpha = 0.05
effect.size = seq(from = 0.05, to = 0.2, by = 0.05)

RejectH0_Theo = data.frame(matrix(nrow = mcreps, ncol = length(effect.size))) 

for(i in 1:length(effect.size)){
  colnames(RejectH0_Theo)[i] = paste0("Theo",effect.size[i])
}                                                

Coverage = NULL
Result = NULL
Summary = NULL

for (setari in 1:length(seqari)){
  set.adj.ri = seqari[setari]
  
  adjri = NULL
  
  for (nsamp in 1:length(n.vec)){

    n = n.vec[nsamp]
    nc2 = choose(n,2)
    
    for (mc in 1:mcreps){
      
      
      neg = 0
      while (neg ==0 ) {  
        
        # Randomly simulate data for c1  (crosstab row sums)
        c1 = sample(1:c1numgroup, size=n, replace = TRUE)
        n1 = table(c1)

        prob = rep(0, c2numgroup)
        
        if(set.adj.ri >= 0 && set.adj.ri <= 0.1){
          prob = rep(1/c2numgroup, c2numgroup)
        }else if(set.adj.ri > 0.1 && set.adj.ri <= 0.7){
          probloop= rep(0, c2numgroup)
          for(c2num in 1:c2numgroup){
            probloop[c2num] = mean(rbinom(n, 1, p=ifelse(c1 != c2num, 0.05, .95)))
          }
          prob = probloop
        }else{
          probloop= rep(0, c2numgroup)
          for(c2num in 1:c2numgroup){
            probloop[c2num] = mean(rbinom(n, 1, p=ifelse(c1 != c2num, 0.001, .999)))
          }
          prob = probloop
        }
        
        prob = ifelse(prob <= 0 & set.adj.ri > 0.1 & set.adj.ri < 0.7, 0.05, 
                      ifelse(prob <= 0 & set.adj.ri >= 0.7, 0.01, prob))
        
        # # Randomly simulate data for c2
        c2 = sample(1:c2numgroup, size=n, replace = TRUE, prob = prob)
        n2 = table(c2)
        
        # # Randomly simulate data for c2
        while(length(n2) != c2numgroup){
          c2 = sample(1:c2numgroup, size=n, replace = TRUE, prob = prob)
          n2 = table(c2)
          
        }

        
        
        Exp.Value = sum(choose(n1,2))*sum(choose(n2,2))/nc2
        Max.Value = 0.5*(sum(choose(n1,2))+ sum(choose(n2,2)))
        
        
        # get n(1,1) value
        
        n11 = round(( ( 2* (set.adj.ri* (Max.Value - Exp.Value) + Exp.Value) - n1[1]^2 + n*n1[1] - n2[1]^2 + n*n2[1] - 3/4*n^2 + n )^(1/2) + n1[1] + n2[1] - n/2 )/2,0)
        
        n12 = n1[1] - n11
        n21 = n2[1] - n11
        n22 = n - n1[1] - n2[1] + n11
        
        neg = min(n11>0, n12 > 0, n21 > 0, n22 >0)
      }
      
      # c1 = c(rep(1,n1[1]), rep(2,n1[2]))
      
      c1.g1.index = which(c1 == 1)
      c1.g2.index = which(c1 == 2)
      
      c2 = c1-1
      

      #c2.switchtog2 = ifelse(length(c1.g1.index) > 0 & b != 0, sample(c1.g1.index, size = b, replace = FALSE), NULL) 
      
      c2.switchtog2 = sample(c1.g1.index, size = n12, replace = FALSE) 
      
      
      # Finally, c - These are c1 g2 that are c2 g1
      
      c2.switchtog1 = sample(c1.g2.index, size = n21, replace = FALSE)
      
      
      index.needs.switched = c(c2.switchtog2, c2.switchtog1)
      index.needs.switched = index.needs.switched[which(index.needs.switched !=0)]
      
      if(sum(index.needs.switched) > 0){
        for(i in 1:length(index.needs.switched)){
          
          c2[index.needs.switched[i]] = abs(c2[index.needs.switched[i]]-1)
          
        }
      }
      c2 = c2 + 1
      
      #c1; c2; rand.index(c1, c2)
      
      adjri[mc] = adjustedRandIndex(c1,c2)
      
      tabsq = n11^2 + n12^2 + n21^2 + n22^2
      a = (tabsq - n)/2
      b = (sum(n1^2) - tabsq)/2
      c = (sum(n2^2) - tabsq)/2
      d = (tabsq + n^2 - sum(n1^2) - sum(n2^2))/2
      
      denom = (a + b)*(a+c) + (b+d)*(c+d)
      
      e = 2*sum(n1^2) - (n+1)*n
      f = 2*sum(n2^2) - (n+1)*n
      g = 4*sum(n1^3) - 4*(n+1)* sum(n1^2) + (n+1)^2*n
      h = n*(n-1)
      i = 4*sum(n2^3) - 4*(n+1)* sum(n2^2) + (n+1)^2*n
      
      
      Var.aplusd = 1/16* (2*n*(n-1)- ( e*f/(n*(n-1)))^2 + 4*(g-h)*(i-h)/(n*(n-1)*(n-2))) + 
        1/16*( ((e^2 - 4*g + 2*h) * (f^2 - 4*i + 2*h))/(n*(n-1)*(n-2)*(n-3)) )
      Var.ARI = (choose(n,2))^2 * Var.aplusd / (((choose(n,2))^2 - denom)^2)
      Theo_LCL = adjri[mc] - 1.96*sqrt(Var.ARI)
      Theo_UCL = adjri[mc] + 1.96*sqrt(Var.ARI)
      
      
 
      
      CritValue=adjri[mc] - 1.645*sqrt(Var.ARI)
      
      
      
      for (es in 1:length(effect.size)){
        RejectH0_Theo[mc,es] = ifelse((set.adj.ri - effect.size[es]) < CritValue, 1, 0 )
        

      }
      
      Coverage[mc] = ifelse(Theo_LCL < set.adj.ri &  set.adj.ri > Theo_UCL, 0, 1)
      
      

      
      
      
    } # end mc loop
    
    
    
    Result.sampsize = as.data.frame(cbind(n, set.adj.ri,t(colMeans(RejectH0_Theo)), mean(adjri), sd(adjri), quantile(adjri,0.025), quantile(adjri, 0.5), quantile(adjri, 0.975)), nrow = 1, ncol = 7+1*ncol(RejectH0_Theo))
    
    
    
    Result = rbind(Result, Result.sampsize)  
    

  } #end sampsize loop
  
  
  
} #end ari loop

kable(Result)