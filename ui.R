library('shiny')
library('rhandsontable')

# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("", "Piperacillin TDM"), 
  
  # Sidebar 
  sidebarLayout(
    sidebarPanel(
      h3("Piperacillin TDM"),
      hr(),
      p("This application is designed to provide individualized estimates of drug exposure for
        critically ill patients with sepsis that have received one or more doses of piperacillin."),
      p("Enter patient dosing schedule (historical and/or proposed) and 
         measurements of drug concentration in blood (if any)."),
      hr(),
      ####
      h4("Infusion Schedule"),
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
      h4("Concentration Data"),
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

