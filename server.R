library('shiny')
library('rhandsontable')

source('Bayes.R')

# Define server logic required to draw plot
server <- function(input, output) {
  # App title with line break
  output$titleText <- renderUI({
    HTML("Piperacillin Therapeutic Drug Monitoring")
  })
  
  
  #rhandsontable for dosing information
  output$dosing <- renderRHandsontable({
    if(1 - input$common){
      if(is.null(input$dosing)){
        doseDF  <- data.frame(
          "Start (h)" = c(0, 8, 16, 24, 32, rep(0,5)),
          "End (h)" = c(0.5, 8.5, 16.5, 24.5, 32.5, rep(0,5)),
          "Rate (g/h)" = c(6, 6, 6, 6, 6, rep(0,5)), check.names=FALSE)
      }else{
        doseDF <- hot_to_r(input$dosing)
      }
      rhandsontable(doseDF, colWidths = c(65,65,70))
    }else{
      comPat <- hot_to_r(input$dosing)
      nDose <- as.numeric(input$num)
      begin <- seq(0, (nDose - 1) * input$freq, by = input$freq)
      end <- begin + input$duration
      kR <- rep(input$infRate, nDose)
      comPat[1:input$num, "Start (h)"] <- begin
      comPat[1:input$num, "End (h)"] <- end
      comPat[1:input$num, "Rate (g/h)"] <- kR
      if(input$num > 10){
        comPat[(input$num + 1):10] <- matrix(0, ncol = 3, nrow = 10 - input$num)
      }
      rhandsontable(comPat, colWidths = c(65,65,70))
    }
  })
  
  #rhandsontable for sample information
  output$sample <- renderRHandsontable({
    vec <- numeric(10)
    sampDF = data.frame("Time (h)" = vec, 
                        "Conc. (μg/ml)" = vec,
                        check.names=FALSE)
    rhandsontable(sampDF, colWidths = c(65,100))
  })
  
  
  ##########################################################
  
  # These functions allow me to grab the data from rhandsontable to use in plotting
  sampTable <- eventReactive(input$goPlot, {
    live_data = hot_to_r(input$sample)
    return(live_data)
  })
  
  doseTable <- eventReactive(input$goPlot, {
    live_data = hot_to_r(input$dosing)
    return(live_data)
  })
  
  
  ##########################################################
  
  # Create the plot
  output$plot <- renderPlot({
    # Get data from rHandsonTable
    stab <- sampTable()
    dtab <- doseTable()
    
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
        plot_post_conc(est, ivtData, dat, thres = input$thres)
        
        #Display coordinates when hovering over a point
        output$info <- renderText({
          paste("Time=", round(input$plot_hover$x,2), "h",
                 "\nConcentration=", round(input$plot_hover$y,2), "μg/ml", sep=" ")
        })
      }
    }
  })
  
}
