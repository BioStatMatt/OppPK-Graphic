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
                  choices = 1:10,
                  selected = 1
      ),
      # Will dynamically create time/dose inputs based on number of measurements (in server)
      uiOutput("measures")
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
  output$measures <- renderUI({
    times <- lapply(1:input$nMeas, function(i){
      timeId <- paste0("time", i)
      timeName <- paste("Time", i, sep = " ")
      numericInput(timeId, timeName, min = 0, value = 0)
    })
    doses <- lapply(1:input$nMeas, function(i){
      doseId <- paste0("dose", i)
      doseName <- paste("Dose", i, sep = " ")
      numericInput(doseId, doseName, min = 0, value = 0)
    })
    do.call(tagList, rbind(times, doses))
  })
  
 
  
  output$plot <- renderPlot({
    #There must be a better way to grab only the input$timeX and input$doseX values I need
    dat <- data.frame(time_h = c(input$time1, input$time2, input$time3, input$time4, input$time5,
                                 input$time6, input$time7, input$time8, input$time9, input$time10),
                      conc_mg_dl = c(input$dose1, input$dose2, input$dose3, input$dose4, input$dose5,
                                     input$dose6, input$dose7, input$dose8, input$dose9, input$dose10))

    #Restrict to actual measurements
    dat <- dat[1:input$nMeas,]
    
    
    # No changes made here
    est <- optim(lpr_mean_d, log_posterior, ivt=ivt_d,
                 dat=dat, control = list(fnscale=-1), hessian=TRUE)
    plot_post_conc(est, ivt_d, dat)
  })
}


# Run the application 
shinyApp(ui = ui, server = server)

