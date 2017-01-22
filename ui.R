library('shiny')
library('rhandsontable')

# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel(htmlOutput("titleText")), 
  
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

