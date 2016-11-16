#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

source('Bayes.R')

# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("Opportunistic Sampling"),
  
  # Sidebar 
  sidebarLayout(
    sidebarPanel(
      # Drop down to select number of available sample measurements
      selectInput("nMeas", "Number of Measurements", 
                  choices = 0:10,
                  selected = 0
      ),
      # Will dynamically create time/dose inputs based on number of measurements (in server)
      uiOutput("measures"),
      # Set up area for rhandsontable output
      helpText("Sample Information"),
      rHandsontableOutput("sample"),
      helpText("Dosing Schedule"),
      rHandsontableOutput("dosing"),
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
  
  # Create input boxes given the number of available samples
  # Can remove if I can get data frames from rhandsontable for plot
  output$measures <- renderUI({
    if(input$nMeas > 0){
      times <- lapply(1:input$nMeas, function(i){
        timeId <- paste0("time", i)
        timeName <- paste("Time", i, sep = " ")
        numericInput(timeId, timeName, min = 0, value = 0)
      })
      doses <- lapply(1:input$nMeas, function(i){
        concId <- paste0("conc", i)
        concName <- paste("Concentration", i, sep = " ")
        numericInput(concId, concName, min = 0, value = 0)
      })
      do.call(tagList, rbind(times, doses))
    }
  })
  
  ###################
  
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
      vec <- numeric(input$nMeas)
      doseDF = data.frame("Start" = vec, "End" = vec, "Infusion_Rate" = vec)
    } else {
      doseDF = hot_to_r(input$hotDose)
    }
    rhandsontable(doseDF)
  })
  
  output$dosing <- renderRHandsontable({
    doseDat()
  })
  
  #Currently not sure how to pull info from rhandsontable and use in plot
  #The functions below may be useful
  sampTable <- eventReactive(input$go, { 
    live_data = hot_to_r(input$sample)
    return(live_data)
  })
  doseTable <- eventReactive(input$go, { 
    live_data = hot_to_r(input$dosing)
    return(live_data)
  })  
  
  
  
  #Create the plot
  output$plot <- renderPlot({
    if(input$nMeas > 0){
      #There must be a better way to grab only the input$timeX and input$doseX values I need
      dat <- data.frame(time_h = c(input$time1, input$time2, input$time3, input$time4, input$time5,
                                   input$time6, input$time7, input$time8, input$time9, input$time10),
                        conc_mg_dl = c(input$conc1, input$conc2, input$conc3, input$conc4, input$conc5,
                                       input$conc6, input$conc7, input$conc8, input$conc9, input$conc10))
      
      #Restrict to actual measurements
      dat <- dat[1:input$nMeas,]
      
      #Try to use rhandsontable
      datHot <- sampTable()
      ivtHot <- doseTable()
      
      #Current roadblock: `optim` function gives `Error in [: incorrect number of dimensions` when 
      #  I plug in datHot for dat
      est <- optim(lpr_mean_d, log_posterior, ivt=ivt_d,
                   dat=dat, control = list(fnscale=-1), hessian=TRUE)
      plot_post_conc(est, ivt_d, dat)
      
      
    }else if(input$nMeas == 0){
      #Put in code to produce plot based only on prior
    }
  })
  
}


# Run the application 
shinyApp(ui = ui, server = server)

