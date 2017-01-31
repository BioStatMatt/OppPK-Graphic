library('shiny')
library('shinydashboard')
library('rhandsontable')

ui <- dashboardPage(
  dashboardHeader(
    title = "Piperacillin TDM"
  ),
  dashboardSidebar(
    tags$head(tags$style(".wrapper {overflow: visible !important;}")),
    hr(),
    p("This application is designed to provide individualized estimates of drug exposure for
      critically ill patients with sepsis that have received one or more doses of piperacillin."),
    p("Enter patient dosing schedule (historical and/or proposed) and 
      measurements of drug concentration in blood (if any)."),
    hr(),
    
    checkboxInput("common", "Common Dosing Pattern"),
    em("Manual entry of infusion data is disabled while this box is checked."),
    
    conditionalPanel(
      condition = "input.common == true",
      sliderInput("num", "Number of Doses", min = 1, max = 10, value = 5, step = 1),
      numericInput("freq", "Dose Frequency (h)", value = 8, min = 0),
      numericInput("duration", "Infusion Duration (h)", value = .5, min = 0),
      numericInput("infRate", "Infusion Rate (g/h)", value = 6, min = 0)
    )
    
    ),
  dashboardBody(
    h3("Piperacillin Therapeutic Drug Monitoring"),
    
    fluidRow(
      column(width = 6,
             h4("Infusion Schedule"),
             rHandsontableOutput("dosing")
      ),
      column(width = 6,
             h4("Concentration Data"),
             rHandsontableOutput("sample")
      )
    ),
    
    HTML('<br/>'),
    actionButton("goPlot", "Update Plot"),
    HTML('<br/><br/>'),
    
    numericInput('thres', "Threshold (μg/ml)", value = 64, min = 0, width = '60%'),
    em("Threshold value is typically some multiple of the minimum inhibitory concentration 
       (MIC) for the target microorganism."),
    HTML('<br/><br/>'),
    em("'fT > threshold' in the plot legend provides an estimate of the fraction of time
       spent above the specified threshold."),
      
    plotOutput("plot", hover = "plot_hover"), 
    verbatimTextOutput("info"),
    
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }")
  )
)

