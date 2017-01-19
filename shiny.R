#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(rhandsontable)

source('Bayes.R')

# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("OppSampPK"), 
  
  # Sidebar 
  sidebarLayout(
    sidebarPanel(
      p("This application is designed to inform the user about an individual patient's 
        response to treatment over time. Currently, the application is based on the
        pharmacodynamics of piperacillin-tazobactam, used to treat infection in patients
        with acute kidney injury."),
      p("Below, users can enter information on the patient's dosing schedule and 
        measurements of drug concentration from available blood samples."),
      ####
      h4("Dosing Schedule"),
      # Set up checkbox and table for common dosing pattern
      checkboxInput("common", "Common Dosing Pattern"),
      conditionalPanel(
        condition = "input.common == true",
        numericInput("freq", "Dose Frequency (h)", value = 8, min = 0),
        numericInput("duration", "Infusion Duration (h)", value = .5, min = 0),
        numericInput("infRate", "Infusion Rate (g/h)", value = 6, min = 0),
        numericInput("num", "Number of Doses", value = 5, min = 0),
        actionButton("goT", "Update Table"),
        tableOutput("commonPattern")
      ),
      rHandsontableOutput("dosing"),
      ####
      
      # Set up area for rhandsontable output
      h4("Sample Information"),
      rHandsontableOutput("sample"),
      # Button for updating plot
      actionButton("go", "Update Plot")
      ),
    # Show a plot
    # Hide a warning that results from the odd way I'm grabbing the time/dose data from input (in server)
    mainPanel(plotOutput("plot"), tags$style(type="text/css",
                                             ".shiny-output-error { visibility: hidden; }",
                                             ".shiny-output-error:before { visibility: hidden; }"))
    )
)


# Define server logic required to draw plot
server <- function(input, output) {
  
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
    dtab <- doseTable()
    
    # Won't produce plot unless both sample information and dosing schedule has been provided
    if(sum(stab > 0) > 0 & sum(dtab > 0) > 0){
      
      # SAMPLE INFORMATION
      datHot <- stab[apply(stab, MARGIN = 1, function(x) any(x > 0)),]
      # Required for compatibility with functions from Bayes.R
      names(datHot) <- c("time_h", "conc_mg_dl")
      
      # DOSING INFORMATION
      if(input$common){
        # Common dosing pattern
        ivtData <- comTab()
        
        # Required for compatibility with functions from Bayes.R
        names(ivtData) <- c("begin", "end", "k_R")
        ivtData <- apply(ivtData, MARGIN = 1, function(x) list(begin = x[1], end = x[2], k_R = x[3]))  
      }else{ # Common dosing pattern not specified
        # Get data from rhandsontable
        ivtHot <- dtab[apply(dtab, MARGIN = 1, function(x) any(x > 0)),]
        
        # Required for compatibility with functions from Bayes.R
        names(ivtHot) <- c("begin", "end", "k_R")
        ivtHot <- apply(ivtHot, MARGIN = 1, function(x) list(begin = x[1], end = x[2], k_R = x[3])) 
        ivtData <- ivtHot
      }
      
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
        
        # Get dosing information
        if(input$common){
          # Common dosing pattern
          ivtData <- comTab()
          
          # Required for compatibility with functions from Bayes.R
          names(ivtData) <- c("begin", "end", "k_R")
          ivtData <- apply(ivtData, MARGIN = 1, function(x) list(begin = x[1], end = x[2], k_R = x[3]))  
        }else{ # Common dosing pattern not specified
          # Get data from rhandsontable
          ivtHot <- dtab[apply(dtab, MARGIN = 1, function(x) any(x > 0)),]
          
          # Required for compatibility with functions from Bayes.R
          names(ivtHot) <- c("begin", "end", "k_R")
          ivtHot <- apply(ivtHot, MARGIN = 1, function(x) list(begin = x[1], end = x[2], k_R = x[3])) 
          ivtData <- ivtHot
        }
        
        est <- optim(lpr_mean_d, log_posterior, ivt=ivtData,
                     dat=dat, control = list(fnscale=-1), hessian=TRUE)
        plot_post_conc(est, ivtData, dat)
      }
    }
  })
  
}


# Run the application 
shinyApp(ui = ui, server = server)

