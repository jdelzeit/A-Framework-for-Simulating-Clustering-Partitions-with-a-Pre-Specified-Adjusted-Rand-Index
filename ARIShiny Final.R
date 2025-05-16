

# Load packages
library(shiny)
library(shinythemes)
library(mclust)
library(knitr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(kableExtra)
library(shinyjs)

library(DataCombine)
library(data.table)
library(tidyverse)
library(shinyWidgets)
library(formattable)
library(readr)
library(DT)


ui <- fluidPage(
  withMathJax(),
  
  theme = shinytheme("lumen"), 
  titlePanel("Adjusted Rand Index Power and Sample Size Input"),
  sidebarLayout(
    sidebarPanel(
      
      selectInput(inputId = "criteria", label = "Select Analysis",
                  choices = list( 
                    "Determine Necessary Sample Size" = "CritN", 
                    "Determine Power" = "CritPower"
                  ), selected = NULL, multiple = FALSE), 
      
      
      checkboxInput("upload", "Upload Cluster Labels for Cluster 1", value = FALSE),
      conditionalPanel(
        condition = "input.upload == true",
        fileInput("file1", "Choose CSV File", accept = ".csv"),
        uiOutput("column_ui")),
      
      conditionalPanel(
        condition = "input.upload == false",
        numericInput(inputId = "c1num", label = "Select Number of Labels for Cluster 1", value = 2, min = 1, max = 20, step = 1)),
      
      numericInput(inputId = "c2num", label = "Select Number of Labels for Cluster 2", value = 2, min = 1, max = 20, step = 1),
      

      
      # conditionalPanel(
      #   condition = "input.upload == false || input.upload == true && input.criteria == 'CritN'",
      # uiOutput("dynamic_input_crit")),

      # uiOutput("dynamic_input_crit"),
      
      conditionalPanel(
        condition = "input.criteria == 'CritN'",
      numericInput("minpower", "Enter Minimum Required Power:", value = 0.8, min = 0, max = 1),
      numericInput("nstart", "Enter Starting Sample Size:", value = 50, min = 40, max = 10000),
      numericInput("nstep", "Enter Value to Increase Sample Size:", value = 25, min = 1, max = 1000)),
      
      conditionalPanel(
        condition = "input.criteria == 'CritPower'",
        numericInput("givensampsize", "Enter Given Sample Size:", value = 50, min = 30, max = 10000)),
        
      
      numericInput(inputId = "ARIHave", label = "Input Actual ARI: ", value = 0.4, min = 0, max = 1, step = 0.05),
    
      
      uiOutput("dynamic_input_hypoth"),
      
      numericInput(inputId = "NumMCReps", label = "Input Number of Repetitions: ", value = 100, min = 0, max = 1000, step = 50),
      
      
      h4("Calculating:", style="color: red;"),
      # h4(uiOutput("hypothesis_text"), style="color: red;"),
      h4(textOutput("hypothesis_text"), style="color: red;"),
      
      
      # Action button to start the calculation
      actionButton("calc_button", "Start Calculation"),
      
      
      hr(),
      h5("Author:", style="font-style: italic;"),
      h5("Jennifer Delzeit", style="color: blue;"),
      h5(today(), style="color: blue;")
      
      #, width = 3.5
      
    ),
    
    mainPanel(
      
      # Action button to start the calculation
      # actionButton("calc_button", "Start Calculation"),
      
      # "Calculating..." message 
      # div(id = "calculating_msg", "Calculating...", style = "font-size: 20px; color: red; display: none;"),
      
      
      textOutput("complete"),
      
      # verbatimTextOutput("selected_column") 
      dataTableOutput(outputId = "analysis"),
      
      conditionalPanel(
        condition = "input.criteria == 'CritN'",
      plotOutput(outputId = "powerplot")
      )
      
    )#mainpanel
  ) #sidebar
  
)




###############################################################
########### Define server function  ###########################
###############################################################

server <- function(session,input, output) {
  
  dataset <- reactive({
    req(input$upload)  # Only proceed if upload is checked
    req(input$file1)  # Ensure file input is available when upload is selected
    read.csv(input$file1$datapath)  # Read CSV file
  })
  
  output$column_ui <- renderUI({
    req(dataset())  
    column_names <- colnames(dataset())  
    selectInput("column", "Choose Column for Cluster 1 Labels", choices = column_names)
  })
  

  output$dynamic_input_hypoth = renderUI({
    numericInput(inputId = "ARIWant", label = "Input Test ARI: ", value = 0.3, min = 0, max = input$ARIHave-0.05, step = 0.05)
    
  })
  

  output$hypothesis_text = renderText({
    # paste0("P(", "\\(\\hat{ARI}\\)", intToUtf8(8805), input$ARIWant, "| ARI = ", input$ARIHave, ")") 
    paste0("P(ARI", intToUtf8(8805), input$ARIWant, "| ARI = ", input$ARIHave, ")") 
    
    
  })
  
  observeEvent(input$calc_button, {
    showModal(modalDialog("Calculating...", footer=NULL))
    
    SetUpSimC1Num = reactive({
      if (input$upload == TRUE) {
        req(input$column)  
        req(dataset())  
        selected_col = dataset()[[input$column]] 
        c1num = length(unique(selected_col))  
      } else {
        c1num = input$c1num 
      }
      
      return(c1num)  
    })
    

  
    
    SampleN1 = reactive({
      c1numgroup = SetUpSimC1Num()

      n1 = if (input$upload == TRUE) {
        req(input$column)  
        req(dataset())  
        selected_col = dataset()[[input$column]]  
        n1 = table(selected_col)
        return(n1)
      } else {
        c1 = sample(1:c1numgroup, size=n, replace = TRUE)
        n1 = table(c1)
        return(n1)
      }

      return(n1)
    })


    ResultTable = reactive({
      
      
      
      if (input$criteria == "CritN"){
        num.mc.loops = input$NumMCReps
        
        ari.vec = input$ARIHave
        
        alpha = 0.05
        effect.size = input$ARIHave - input$ARIWant
        n.vec = input$nstart
        
        
        c1numgroup = SetUpSimC1Num() #number of labels in clustering group 1
        c2numgroup = input$c2num #number of labels in clustering group 2
        
        statistic.theo = NULL
        Power_All = NULL
        
        RejectH0 = data.frame(matrix(nrow = num.mc.loops, ncol = length(effect.size))) 
        
        
        RejectH0_Theo = data.frame(matrix(nrow = num.mc.loops, ncol = 1)) 
        
        colnames(RejectH0_Theo) = "Result"
          
                                            
        
        Result = NULL
        Summary = NULL
        
        Power = 0
        
        
        sampari = 1
        set.adj.ri = input$ARIHave
        adjri = NULL
        
        

          while(Power < input$minpower){
            
            n = n.vec
            nc2 = choose(n,2)
            
            for (mc in 1:num.mc.loops){
              
              # Randomly simulate data for c1  (crosstab row sums)
              
        
            if (input$upload == TRUE) {
                req(input$column)  
                req(dataset())  
                selected_col = dataset()[[input$column]]  
                n1 = table(selected_col)
              } else {
                c1 = sample(1:c1numgroup, size=n, replace = TRUE)
                n1 = table(c1)
              }
              

              c1 = NULL
              for (i in 1:c1numgroup){
                c1 = c(c1, rep(i, n1[i]))
              }
              
              
              prob = rep(0, c2numgroup)
              
              if(set.adj.ri >= 0 && set.adj.ri <= 0.1){
                prob = rep(1/c2numgroup, c2numgroup)
              }else if(set.adj.ri > 0.1 && set.adj.ri <= 0.7){
                probloop= rep(0, c2numgroup)
                for(c2num in 1:c2numgroup){
                  probloop[c2num] = mean(rbinom(n, 1, p=ifelse(c1 != c2num, 0.1, .9)))
                }
                prob = probloop
              }else{
                probloop= rep(0, c2numgroup)
                for(c2num in 1:c2numgroup){
                  probloop[c2num] = mean(rbinom(n, 1, p=ifelse(c1 != c2num, 0.005, .995)))
                }
                prob = probloop
              }
              
              prob = ifelse(prob <= 0 & set.adj.ri > 0.1 & set.adj.ri <= 0.7, 0.1,
                            ifelse(prob <= 0 & set.adj.ri > 0.7, 0.005, prob))              
              # # Randomly simulate data for c2
              c2 = sample(1:c2numgroup, size=n, replace = TRUE, prob = prob)
              n2 = table(c2)
              
              # # Randomly simulate data for c2
              while(length(n2) != c2numgroup){
                c2 = sample(1:c2numgroup, size=n, replace = TRUE, prob = prob)
                n2 = table(c2)
                
              }
              
              c2 = NULL
              for (i in 1:c2numgroup){
                c2 = c(c2, rep(i, n2[i]))
              }
              
              
              

              Exp.Value = sum(choose(n1,2))*sum(choose(n2,2))/nc2
              Max.Value = 0.5*(sum(choose(n1,2))+ sum(choose(n2,2)))
              
              # Need sum(choose(nij,2)) value
              
              sum.choose.crosstab = set.adj.ri*(Max.Value - Exp.Value) + Exp.Value
              
              
              lower = sum.choose.crosstab *0.99
              upper = sum.choose.crosstab *1.01
              
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
                
                
                sumchoose_check = ifelse(sumchoose >= lower & sumchoose <= upper & (sum(crosstab < 0) == 0), 1, 0)
                
                num.repeats = num.repeats + 1
                
                if(num.repeats > 100){
                  
                  n1 = if (input$upload == TRUE) {
                    req(input$column) 
                    req(dataset())  
                    selected_col = dataset()[[input$column]] 
                    n1 = table(selected_col)
                  } else {
                    # Randomly simulate data for c1  (crosstab row sums)
                    c1 = sample(1:c1numgroup, size=n, replace = TRUE)
                    n1 = table(c1)
                  }

                  c1 = NULL
                  for (i in 1:c1numgroup){
                    c1 = c(c1, rep(i, n1[i]))
                  }
                  
                  prob = rep(0, c2numgroup)
                  
                  if(set.adj.ri >= 0 && set.adj.ri <= 0.1){
                    prob = rep(1/c2numgroup, c2numgroup)
                  }else if(set.adj.ri > 0.1 && set.adj.ri <= 0.7){
                    probloop= rep(0, c2numgroup)
                    for(c2num in 1:c2numgroup){
                      probloop[c2num] = mean(rbinom(n, 1, p=ifelse(c1 != c2num, 0.1, .9)))
                    }
                    prob = probloop
                  }else{
                    probloop= rep(0, c2numgroup)
                    for(c2num in 1:c2numgroup){
                      probloop[c2num] = mean(rbinom(n, 1, p=ifelse(c1 != c2num, 0.005, .995)))
                    }
                    prob = probloop
                  }
                  
                  prob = ifelse(prob <= 0 & set.adj.ri > 0.1 & set.adj.ri <= 0.7, 0.1,
                                ifelse(prob <= 0 & set.adj.ri > 0.7, 0.005, prob))
                  
                  c2 = sample(1:c2numgroup, size=n, replace = TRUE, prob = prob)
                  n2 = table(c2)
                  
                  # # Randomly simulate data for c2
                  while(length(n2) != c2numgroup){
                    c2 = sample(1:c2numgroup, size=n, replace = TRUE, prob = prob)
                    n2 = table(c2)
                    
                  }
                  
                  c2 = NULL
                  for (i in 1:c2numgroup){
                    c2 = c(c2, rep(i, n2[i]))
                  }
                  
                  
                  # 
                  Exp.Value = sum(choose(n1,2))*sum(choose(n2,2))/nc2
                  Max.Value = 0.5*(sum(choose(n1,2))+ sum(choose(n2,2)))
                  
                  # Need sum(choose(nij,2)) value
                  
                  sum.choose.crosstab = set.adj.ri*(Max.Value - Exp.Value) + Exp.Value
                  
                  
                  lower = sum.choose.crosstab *0.99
                  upper = sum.choose.crosstab *1.01
                  
                  sumchoose_check = 0
                  num.repeats = 1
                  
                }        
              }
              
              
              
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
              
              
              tabsq = sum(crosstab^2)
              a = (tabsq - n)/2
              b = (sum(n1^2) - tabsq)/2
              c = (sum(n2^2) - tabsq)/2
              d = (tabsq = n^2 - sum(n1^2) - sum(n2^2))/2
              
              denom = (a + b)*(a+c) + (b+d)*(c+d)
              
              e = 2*sum(n1^2) - (n+1)*n
              f = 2*sum(n2^2) - (n+1)*n
              g = 4*sum(n1^3) - 4*(n+1)* sum(n1^2) + (n+1)^2*n
              h = n*(n-1)
              i = 4*sum(n2^3) - 4*(n+1)* sum(n2^2) + (n+1)^2*n
              
              
              Var.aplusd = 1/16* (2*n*(n-1)- ( e*f/(n*(n-1)))^2 + 4*(g-h)*(i-h)/(n*(n-1)*(n-2))) + 
                1/16*( ((e^2 - 4*g + 2*h) * (f^2 - 4*i + 2*h))/(n*(n-1)*(n-2)*(n-3)) )
              Var.ARI = (choose(n,2))^2 * Var.aplusd / (((choose(n,2))^2 - denom)^2)
              
              
              
              
              
              for (es in 1:length(effect.size)){
                RejectH0_Theo[mc,es] = ifelse(pnorm(abs(adjri[mc]-(set.adj.ri - effect.size[es]))/sqrt(Var.ARI), lower.tail = FALSE) < alpha, 1, 0)
                
              }
              
              
              
              
            } # end mc loop
            
            
            Result.sampsize = as.data.frame(cbind(n, set.adj.ri,t(colMeans(RejectH0_Theo)), mean(adjri), quantile(adjri,0.025), quantile(adjri, 0.5), quantile(adjri, 0.975)), nrow = 1, ncol = 6+ncol(RejectH0_Theo))
            
            
            Result = rbind(Result, Result.sampsize)  
            
            Power = as.numeric(t(colMeans(RejectH0_Theo)))
            # Power_DF = as.data.frame(cbind(n, Power), nrow = 1, ncol = 2)
            # Power_All = rbind(Power_All, Power_DF)
            
            n.vec = n.vec + input$nstep
  
            
          } #end While loop for Power
          
          
          
      } else if (input$criteria == "CritPower") {
        
        
        n.vec = input$givensampsize
        
        alpha = 0.05
        ari.vec = input$ARIHave
        
        alpha = 0.05
        # effect.size = seq(from = 0.05, to = 0.2, by = 0.05)
        effect.size = input$ARIHave - input$ARIWant
        
        c1numgroup = SetUpSimC1Num() #number of labels in clustering group 1
        c2numgroup = input$c2num #number of labels in clustering group 2
        
        num.mc.loops = input$NumMCReps
        
        RejectH0 = data.frame(matrix(nrow = num.mc.loops, ncol = length(effect.size))) 
        
        
        RejectH0_Theo = data.frame(matrix(nrow = num.mc.loops, ncol = length(effect.size))) 
        
        for(i in 1:length(effect.size)){
          colnames(RejectH0_Theo)[i] = paste0("Theo",effect.size[i])
          
        }                                                
        
        Result = NULL
        Summary = NULL
        
        
        sampari = 1
        set.adj.ri = input$ARIHave
        adjri = NULL
        
        
        nsamp = 1
        n = input$givensampsize
        nc2 = choose(n,2)
        
        for (mc in 1:num.mc.loops){
          
          if (input$upload == TRUE) {
            req(input$column)  
            req(dataset())  
            selected_col = dataset()[[input$column]]  
            n1 = table(selected_col)
          } else {
            # Randomly simulate data for c1  (crosstab row sums)
            c1 = sample(1:c1numgroup, size=n, replace = TRUE)
            n1 = table(c1)
          }          
          
          

          c1 = NULL
          for (i in 1:c1numgroup){
            c1 = c(c1, rep(i, n1[i]))
          }
          
          prob = rep(0, c2numgroup)
          
          if(set.adj.ri >= 0 && set.adj.ri <= 0.1){
            prob = rep(1/c2numgroup, c2numgroup)
          }else if(set.adj.ri > 0.1 && set.adj.ri <= 0.7){
            probloop= rep(0, c2numgroup)
            for(c2num in 1:c2numgroup){
              probloop[c2num] = mean(rbinom(n, 1, p=ifelse(c1 != c2num, 0.1, .9)))
            }
            prob = probloop
          }else{
            probloop= rep(0, c2numgroup)
            for(c2num in 1:c2numgroup){
              probloop[c2num] = mean(rbinom(n, 1, p=ifelse(c1 != c2num, 0.005, .995)))
            }
            prob = probloop
          }
          
          prob = ifelse(prob <= 0 & set.adj.ri > 0.1 & set.adj.ri <= 0.7, 0.1,
                        ifelse(prob <= 0 & set.adj.ri > 0.7, 0.005, prob))
          
          c2 = sample(1:c2numgroup, size=n, replace = TRUE, prob = prob)
          n2 = table(c2)
          
          # # Randomly simulate data for c2
          while(length(n2) != c2numgroup){
            c2 = sample(1:c2numgroup, size=n, replace = TRUE, prob = prob)
            n2 = table(c2)
            
          }
          
          c2 = NULL
          for (i in 1:c2numgroup){
            c2 = c(c2, rep(i, n2[i]))
          }
          
          
          Exp.Value = sum(choose(n1,2))*sum(choose(n2,2))/nc2
          Max.Value = 0.5*(sum(choose(n1,2))+ sum(choose(n2,2)))
          
          # Need sum(choose(nij,2)) value
          
          sum.choose.crosstab = set.adj.ri*(Max.Value - Exp.Value) + Exp.Value
          
          
          lower = sum.choose.crosstab *0.99
          upper = sum.choose.crosstab *1.01
          
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
            
            
            sumchoose_check = ifelse(sumchoose >= lower & sumchoose <= upper & (sum(crosstab < 0) == 0), 1, 0)
            
            num.repeats = num.repeats + 1
            
            if(num.repeats > 100){
              
              if (input$upload == TRUE) {
                req(input$column) 
                req(dataset())  
                selected_col = dataset()[[input$column]]  
                n1 = table(selected_col)
              } else {
                # Randomly simulate data for c1  (crosstab row sums)
                c1 = sample(1:c1numgroup, size=n, replace = TRUE)
                n1 = table(c1)
              }
              
              
              # c1 = sample(1:c1numgroup, size=n, replace = TRUE)
              # n1 = table(c1)
              
              c1 = NULL
              for (i in 1:c1numgroup){
                c1 = c(c1, rep(i, n1[i]))
              }
              
              prob = rep(0, c2numgroup)
              
              if(set.adj.ri >= 0 && set.adj.ri <= 0.1){
                prob = rep(1/c2numgroup, c2numgroup)
              }else if(set.adj.ri > 0.1 && set.adj.ri <= 0.7){
                probloop= rep(0, c2numgroup)
                for(c2num in 1:c2numgroup){
                  probloop[c2num] = mean(rbinom(n, 1, p=ifelse(c1 != c2num, 0.1, .9)))
                }
                prob = probloop
              }else{
                probloop= rep(0, c2numgroup)
                for(c2num in 1:c2numgroup){
                  probloop[c2num] = mean(rbinom(n, 1, p=ifelse(c1 != c2num, 0.005, .995)))
                }
                prob = probloop
              }
              
              prob = ifelse(prob <= 0 & set.adj.ri > 0.1 & set.adj.ri <= 0.7, 0.1,
                            ifelse(prob <= 0 & set.adj.ri > 0.7, 0.005, prob))
              
              c2 = sample(1:c2numgroup, size=n, replace = TRUE, prob = prob)
              n2 = table(c2)
              
              # # Randomly simulate data for c2
              while(length(n2) != c2numgroup){
                c2 = sample(1:c2numgroup, size=n, replace = TRUE, prob = prob)
                n2 = table(c2)
                
              }
              
              c2 = NULL
              for (i in 1:c2numgroup){
                c2 = c(c2, rep(i, n2[i]))
              }
              
              
              
              # 
              Exp.Value = sum(choose(n1,2))*sum(choose(n2,2))/nc2
              Max.Value = 0.5*(sum(choose(n1,2))+ sum(choose(n2,2)))
              
              # Need sum(choose(nij,2)) value
              
              sum.choose.crosstab = set.adj.ri*(Max.Value - Exp.Value) + Exp.Value
              
              
              lower = sum.choose.crosstab *0.99
              upper = sum.choose.crosstab *1.01
              
              sumchoose_check = 0
              num.repeats = 1
              
            }        
          }
          
          
          
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
          
          
          tabsq = sum(crosstab^2)
          a = (tabsq - n)/2
          b = (sum(n1^2) - tabsq)/2
          c = (sum(n2^2) - tabsq)/2
          d = (tabsq = n^2 - sum(n1^2) - sum(n2^2))/2
          
          denom = (a + b)*(a+c) + (b+d)*(c+d)
          
          e = 2*sum(n1^2) - (n+1)*n
          f = 2*sum(n2^2) - (n+1)*n
          g = 4*sum(n1^3) - 4*(n+1)* sum(n1^2) + (n+1)^2*n
          h = n*(n-1)
          i = 4*sum(n2^3) - 4*(n+1)* sum(n2^2) + (n+1)^2*n
          
          
          Var.aplusd = 1/16* (2*n*(n-1)- ( e*f/(n*(n-1)))^2 + 4*(g-h)*(i-h)/(n*(n-1)*(n-2))) + 
            1/16*( ((e^2 - 4*g + 2*h) * (f^2 - 4*i + 2*h))/(n*(n-1)*(n-2)*(n-3)) )
          Var.ARI = (choose(n,2))^2 * Var.aplusd / (((choose(n,2))^2 - denom)^2)
          
          
          
          
          
          for (es in 1:length(effect.size)){
            RejectH0_Theo[mc,es] = ifelse(pnorm(abs(adjri[mc]-(set.adj.ri - effect.size[es]))/sqrt(Var.ARI), lower.tail = FALSE) < alpha, 1, 0)
            
          }
          
          
          
          
        } # end mc loop
        
        
        Result.sampsize = as.data.frame(cbind(n, set.adj.ri,t(colMeans(RejectH0_Theo)), mean(adjri), quantile(adjri,0.025), quantile(adjri, 0.5), quantile(adjri, 0.975)), nrow = 1, ncol = 6+ncol(RejectH0_Theo))
        
        
        Result = rbind(Result, Result.sampsize)  
        
        
      
        
        } else {}
      
      
      Result
    }) # Close Reactive DT
    
    
    output$analysis = renderDT({
      DT = ResultTable()
      
      PowerAchieve = round(as.numeric(DT[nrow(DT),3]),4)
      SampleSizeReq = as.numeric(DT[nrow(DT),1])
      
      DT_Display = as.data.table(rbind(SampleSizeReq, PowerAchieve), nrow = 2, ncol = 1)
      
      colnames(DT_Display) = "Value"
      rownames(DT_Display) = c("The Required Sample Size", "The Power Acheived")
      
      datatable(DT_Display, options(dom="t"))
      
    })
    
    output$powerplot = renderPlot({
      
      DF = ResultTable()
      
      ggplot(data = DF, aes(x = DF[,1], y = DF[,3])) +
        geom_line()+
        geom_hline(yintercept = input$minpower, linetype="dashed",  color = "firebrick3", linewidth=1.25) +
        labs(x = "Sample Size", y = "Power") +
        ggtitle("Power Curve")
    })
    

    removeModal()
    
    output$complete = renderText({
      "Calculation complete!"
    })
    
  }) #observe event close
  
} # close server

# Create Shiny object
shinyApp(ui = ui, server = server)













