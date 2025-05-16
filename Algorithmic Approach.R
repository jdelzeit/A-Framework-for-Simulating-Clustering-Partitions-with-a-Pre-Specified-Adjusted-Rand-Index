library(mclust)
library(knitr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(kableExtra)



n.vec = c(seq(from = 50, to = 500, by = 50), seq(from =600, to = 1000, by = 100))
ari.vec = seqari = seq(from= 0.1, to = 0.9, by = 0.1)

alpha = 0.05
effect.size = seq(from = 0.05, to = 0.2, by = 0.05)


c1numgroup = 2 #number of labels in clustering group 1
c2numgroup = 2 #number of labels in clustering group 2

num.mc.loops = 500
alpha = 0.05


RejectH0_Theo = data.frame(matrix(nrow = num.mc.loops, ncol = length(effect.size))) 

for(i in 1:length(effect.size)){
  colnames(RejectH0_Theo)[i] = paste0("Theo",effect.size[i])
}                                                

Result = NULL
Summary = NULL


for(sampari in 1:length(ari.vec)){
  
  set.adj.ri = ari.vec[sampari]
  adjri = NULL
  
  
  for (nsamp in 1:length(n.vec)){
    
    n = n.vec[nsamp]
    nc2 = choose(n,2)
    
    for (mc in 1:num.mc.loops){
      
      
      # 
      #   # Randomly simulate data for c1  (crosstab row sums)
      #   c1 = sample(1:c1numgroup, size=n, replace = TRUE)
      #   n1 = table(c1)
      #   
      #   # Randomly simulate data for c2
      #  c2 = sample(1:c2numgroup, size=n, replace = TRUE)
      # n2 = table(c2)
      # 
      # 
      
      #   
      # item1 = if(set.adj.ri >= 0 && set.adj.ri <= 0.1){
      #           rbinom(2, size = n, prob = 0.5)
      #       # }else if(set.adj.ri>0 && set.adj.ri <= 0.6) {
      #       #   rhyper(2,n*10*set.adj.ri, n*10*(1-set.adj.ri), n)
      #       }else if(set.adj.ri > 0.1 && set.adj.ri <= 0.9){
      #          i11 =  rbinom(n, 1, prob = set.adj.ri) 
      #          i12 = rbinom(n, 1, p=ifelse(i11==1, .95, .05))
      #          c(sum(i11),sum(i12))
      #       } else{
      #         i11 =  rbinom(n, 1, prob = set.adj.ri) 
      #          i12 = rbinom(n, 1, p=ifelse(i11==1, 0.999, 0.001))
      #          c(sum(i11),sum(i12))
      #       }
      
      
      # n1 = c(item1[1], n-item1[1])
      # n2 = c(item1[2], n-item1[2])
      # 
      # c1 = c(rep(1,n1[1]), rep(2,n1[2]))
      
      # Randomly simulate data for c1  (crosstab row sums)
      c1 = sample(1:2, size=n, replace = TRUE)
      n1 = table(c1)
      
      # Randomly simulate data for c2
      # c2 = sample(1:2, size=n, replace = TRUE)
      # n2 = table(c2)
      
      
      # c1 = NULL
      # for (i in 1:c1numgroup){
      #   c1 = c(c1, rep(i, n1[i]))
      # }
      
      
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
      
      # c2 = NULL
      # for (i in 1:c2numgroup){
      #   c2 = c(c2, rep(i, n2[i]))
      # }
      
      
      # n1;n2;
      
      Exp.Value = sum(choose(n1,2))*sum(choose(n2,2))/nc2
      Max.Value = 0.5*(sum(choose(n1,2))+ sum(choose(n2,2)))
      
      # Need sum(choose(nij,2)) value
      
      sum.choose.crosstab = set.adj.ri*(Max.Value - Exp.Value) + Exp.Value
      
      
      lower = sum.choose.crosstab *0.975
      upper = sum.choose.crosstab *1.025
      
      sumchoose_check = 0
      num.repeats = 1
      
      
      while(sumchoose_check == 0){
        
        
        
        crosstab = matrix(0, nrow =c1numgroup, ncol = c2numgroup)
        
        maxn1 = n1
        maxn2 = n2
        
        for (i in 1:(c1numgroup-1)) {
          for (j in 1:(c2numgroup-1)) {
            maxvalue = min(maxn1[i], maxn2[j])
            crosstab[i, j] = sample(0:maxvalue, 1)
            
            maxn1[i] = maxn1[i]-crosstab[i,j]
            maxn2[j] = maxn2[j]-crosstab[i,j]
            
          }
        }
        
        crosstab[1:(c1numgroup-1),c2numgroup] = n1[1:(c1numgroup-1)]-rowSums(crosstab)[1:(c1numgroup-1)]
        crosstab[c1numgroup,1:c2numgroup] = n2 - colSums(crosstab)
        
        
        # Do the squared sum of these = sumchoosecrosstab?
        
        sumchoose = sum(choose(crosstab,2))
        
        
        # lower; sumchoose; upper;
        
        sumchoose_check = ifelse(sumchoose >= lower & sumchoose <= upper & (sum(crosstab < 0) == 0), 1, 0)
        
        num.repeats = num.repeats + 1
        
        if(num.repeats > 100){
          
          # Randomly simulate data for c1  (crosstab row sums)
          #  c1 = sample(1:c1numgroup, size=n, replace = TRUE)
          #  n1 = table(c1)
          # 
          # # 
          #  c2 = sample(1:c2numgroup, size=n, replace = TRUE)
          #  n2 = table(c2)
          # 
          
          #     
          # item1 = if(set.adj.ri >= 0 && set.adj.ri <= 0.1){
          #         rbinom(2, size = n, prob = 0.5)
          #     # }else if(set.adj.ri>0 && set.adj.ri <= 0.6) {
          #     #   rhyper(2,n*10*set.adj.ri, n*10*(1-set.adj.ri), n)
          #     }else if(set.adj.ri > 0.1 && set.adj.ri <= 0.8){
          #        i11 =  rbinom(n, 1, prob = set.adj.ri) 
          #        i12 = rbinom(n, 1, p=ifelse(i11==1, .95, .05))
          #        c(sum(i11),sum(i12))
          #     } else{
          #       i11 =  rbinom(n, 1, prob = set.adj.ri) 
          #        i12 = rbinom(n, 1, p=ifelse(i11==1, 0.999, 0.001))
          #        c(sum(i11),sum(i12))
          #     }
          
          
          
          # n1 = c(item1[1], n-item1[1])
          # n2 = c(item1[2], n-item1[2])
          # 
          # c1 = c(rep(1,n1[1]), rep(2,n1[2]))
          # 
          
          # Randomly simulate data for c1  (crosstab row sums)
          c1 = sample(1:2, size=n, replace = TRUE)
          n1 = table(c1)
          
          # Randomly simulate data for c2
          # c2 = sample(1:2, size=n, replace = TRUE)
          # n2 = table(c2)
          
          #         c1 = NULL
          # for (i in 1:c1numgroup){
          #   c1 = c(c1, rep(i, n1[i]))
          # }
          # 
          
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
          
          # c2 = NULL
          # for (i in 1:c2numgroup){
          #   c2 = c(c2, rep(i, n2[i]))
          # }
          # 
          # 
          
          Exp.Value = sum(choose(n1,2))*sum(choose(n2,2))/nc2
          Max.Value = 0.5*(sum(choose(n1,2))+ sum(choose(n2,2)))
          
          # Need sum(choose(nij,2)) value
          
          sum.choose.crosstab = set.adj.ri*(Max.Value - Exp.Value) + Exp.Value
          
          
          lower = sum.choose.crosstab *0.975
          upper = sum.choose.crosstab *1.025
          
          sumchoose_check = 0
          num.repeats = 1
          
        }        
      }
      
      
      c2 = c1
      
      for (i in 1:c1numgroup){
        switchindex = which(c1 == i)
        for (j in 1:c2numgroup){
          if( i != j && crosstab[i,j] > 0){
            indicestoswitch = sample(switchindex, size = crosstab[i,j], replace = FALSE)
            c2[indicestoswitch] = j
            
            switchindex = switchindex[-which(switchindex %in% indicestoswitch)]
          }
          else{c2=c2}
        }
      }
      
      #c1; c2; rand.index(c1, c2)
      observed.ari = adjustedRandIndex(c1,c2)
      adjri[mc] = observed.ari
      # observed.ari
      
      
      tabsq = sum(crosstab^2)
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
      
      
      
      
      
      CritValue=adjri[mc] - 1.645*sqrt(Var.ARI)
      
      
      
      for (es in 1:length(effect.size)){
        RejectH0_Theo[mc,es] = ifelse((set.adj.ri - effect.size[es]) < CritValue, 1, 0 )
        
      }
      
      
    } # end mc loop
    
    
    
    Result.sampsize = as.data.frame(cbind(n, set.adj.ri,t(colMeans(RejectH0_Theo)), mean(adjri), sd(adjri), quantile(adjri,0.025), quantile(adjri, 0.5), quantile(adjri, 0.975)), nrow = 1, ncol = 7+1*ncol(RejectH0_Theo))
    
    
    
    Result = rbind(Result, Result.sampsize)  
    
  } #end sampsize loop
  
  
  
} #end ari loop

kable(Result)
