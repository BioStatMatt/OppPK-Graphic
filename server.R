library('shiny')
library('rhandsontable')

source('Bayes.R')

# Define server logic required to draw plot
server <- function(input, output) {
  # App title with line break
  output$titleText <- renderUI({
    HTML(paste0("Therapeutic drug monitoring:", '<br/>',
                "Assessing individual pharmacokinetic heterogeneity"))
  })
  
  
  #rhandsontable for dosing information
  output$dosing <- renderRHandsontable({
    if(1 - input$common){
      vec <- numeric(10)
      doseDF = data.frame("Start" = vec, "End" = vec, "Infusion Rate" = vec,
                          check.names = FALSE)
      rhandsontable(doseDF)
    }
  })
  
  #rhandsontable for sample information
  output$sample <- renderRHandsontable({
    vec <- numeric(10)
    sampDF = data.frame("Time" = vec, "Concentration" = vec)
    rhandsontable(sampDF)
  })
  
  
  ##########################################################
  
  # These functions allow me to grab the data from rhandsontable to use in plotting
  sampTable <- eventReactive(input$go, {
    live_data = hot_to_r(input$sample)
    return(live_data)
  })
  doseTable <- eventReactive(input$go, {
    live_data = hot_to_r(input$dosing)
    return(live_data)
  })
  
  ##########################################################
  # Make table if a common dosing pattern was used
  comTab <- eventReactive(input$goT, {
    if(input$common){
      nDose <- as.numeric(input$num)
      begin <- seq(0, (nDose - 1) * input$freq, by = input$freq)
      end <- begin + input$duration
      kR <- rep(input$infRate, nDose)
      comPat <- data.frame("Start" = begin, "End" = end, "Infusion Rate" = kR,
                           check.names = FALSE)      
    }else{
      comPat <- NULL
    }
    return(comPat)
  })
  
  output$commonPattern <- renderTable({
    comTab()
  })
  
  ##########################################################
  
  # Create the plot
  output$plot <- renderPlot({
    # Get data from rHandsonTable
    stab <- sampTable()
    if(input$common){
      dtab <- comTab()
    }else{
      dtab <- doseTable()
    }
    
    # Won't produce plot unless both sample information and dosing schedule has been provided
    if(sum(stab > 0) > 0 & sum(dtab > 0) > 0){
      
      # SAMPLE INFORMATION
      datHot <- stab[apply(stab, MARGIN = 1, function(x) any(x > 0)),]
      # Required for compatibility with functions from Bayes.R
      names(datHot) <- c("time_h", "conc_mg_dl")
      
      # DOSING INFORMATION
      # Get data from rhandsontable
      ivtHot <- dtab[apply(dtab, MARGIN = 1, function(x) any(x > 0)),]
      
      # Required for compatibility with functions from Bayes.R
      names(ivtHot) <- c("begin", "end", "k_R")
      ivtHot <- apply(ivtHot, MARGIN = 1, function(x) list(begin = x[1], end = x[2], k_R = x[3])) 
      ivtData <- ivtHot
      
      # Functions from Bayes.R
      # To try default example, use:
      #    dat = data.frame(time_h = c(1,4,40), conc_mg_dl = c(82.7,80.4,60))
      #    ivt = ivt_d
      est <- optim(lpr_mean_d, log_posterior, ivt=ivtData,
                   dat=datHot, control = list(fnscale=-1), hessian=TRUE)
      plot_post_conc(est, ivtData, datHot)
      
      
    }else if(sum(stab > 0) == 0){
      if(sum(dtab > 0) > 0){
        dat <- data.frame("empty" = numeric(0))
        
        # DOSING INFORMATION
        # Get data from rhandsontable
        ivtHot <- dtab[apply(dtab, MARGIN = 1, function(x) any(x > 0)),]
        
        # Required for compatibility with functions from Bayes.R
        names(ivtHot) <- c("begin", "end", "k_R")
        ivtHot <- apply(ivtHot, MARGIN = 1, function(x) list(begin = x[1], end = x[2], k_R = x[3])) 
        ivtData <- ivtHot
        
        est <- optim(lpr_mean_d, log_posterior, ivt=ivtData,
                     dat=dat, control = list(fnscale=-1), hessian=TRUE)
        plot_post_conc(est, ivtData, dat)
      }
    }
  })
  
}
