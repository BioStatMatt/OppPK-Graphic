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
      # Consider replacing these with numeric inputs
      # Drop down to select number of available sample measurements
      selectInput("nMeas", "Number of Sample Measurements",
                  choices = 0:10,
                  selected = 0),
      # Drop down to select number of doses administered
      selectInput("nDoseAdmin", "Number of Doses Administered", 
                  choices = 0:10, # I'm not sure what a reasonable upper bound would be here
                  selected = 0),
      # Set up area for rhandsontable output
      helpText("Sample Information"),
      rHandsontableOutput("sample"),
      helpText("Dosing Schedule"),
      rHandsontableOutput("dosing"),
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
  
  #rhandsontable for sample information
  sampDat <- reactive({
    if (is.null(input$hotSamp)) {
      vec <- numeric(input$nMeas)
      sampDF = data.frame("Time" = vec, "Concentration" = vec)
    } else {
      sampDF = hot_to_r(input$hotSamp)
    }
    rhandsontable(sampDF)
  })
  
  output$sample <- renderRHandsontable({
    sampDat()
  })
  
  
  #rhandsontable for dosing information
  doseDat <- reactive({
    if (is.null(input$hotDose)) {
      vec <- numeric(input$nDoseAdmin)
      doseDF = data.frame("Start" = vec, "End" = vec, "Infusion_Rate" = vec)
    } else {
      doseDF = hot_to_r(input$hotDose)
    }
    rhandsontable(doseDF)
  })
  
  output$dosing <- renderRHandsontable({
    doseDat()
  })
  
  # These functions allow me to grab the data from rhandsontable to use in plotting
  sampTable <- eventReactive(input$go, {
    live_data = hot_to_r(input$sample)
    return(live_data)
  })
  doseTable <- eventReactive(input$go, {
    live_data = hot_to_r(input$dosing)
    return(live_data)
  })
  
  
  
  # Create the plot
  output$plot <- renderPlot({
    # Won't produce plot unless both sample information and dosing schedule has been provided
    if(input$nMeas > 0 & input$nDoseAdmin > 0){
      
      # Get data from rhandsontable
      datHot <- sampTable()
      ivtHot <- doseTable()
      
      # Required for compatibility with functions from Bayes.R
      names(datHot) <- c("time_h", "conc_mg_dl")
      names(ivtHot) <- c("begin", "end", "k_R")
      ivtHot <- apply(ivtHot, MARGIN = 1, function(x) list(begin = x[1], end = x[2], k_R = x[3]))
      
      # Functions from Bayes.R
      # To try default example, use:
      #    dat = data.frame(time_h = c(1,4,40), conc_mg_dl = c(82.7,80.4,60))
      #    ivt = ivt_d
      est <- optim(lpr_mean_d, log_posterior, ivt=ivtHot,
                   dat=datHot, control = list(fnscale=-1), hessian=TRUE)
      plot_post_conc(est, ivtHot, datHot)
      
      
    }else if(input$nMeas == 0){
      if(input$nDoseAdmin > 0){
        dat <- data.frame("empty" = numeric(0))
        
        ivt_d <- list(list(begin=0.0, end=0.5, k_R=6),
                      list(begin=8.0, end=8.5, k_R=6),
                      list(begin=16.0, end=16.5, k_R=6),
                      list(begin=24.0, end=24.5, k_R=6),
                      list(begin=32.0, end=32.5, k_R=6))
        
        ivtHot <- doseTable()
        names(ivtHot) <- c("begin", "end", "k_R")
        ivtHot <- apply(ivtHot, MARGIN = 1, function(x) list(begin = x[1], end = x[2], k_R = x[3]))
        
        est <- optim(lpr_mean_d, log_posterior, ivt=ivtHot,
                     dat=dat, control = list(fnscale=-1), hessian=TRUE)
        plot_post_conc(est, ivtHot, dat)
      }
    }
  })
  
  
}


# Run the application 
shinyApp(ui = ui, server = server)

