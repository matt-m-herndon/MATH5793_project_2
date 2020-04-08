library(shiny)
library(shinydashboard)
library(shinyMatrix)
library(matrixcalc)
library(project_2)

# default sigma
sigma = matrix(c(100,40,40,400), 2, 2)

ui <- dashboardPage(
  dashboardHeader(title = paste("Simulation of Xbar and ", expression(theta)), titleWidth = 500),
  dashboardSidebar(
    fluidRow(
      box(width=12,
        solidHeader = FALSE,
        collapsible = FALSE,
        background = "black",
        title = "Simulation Controls",
        numericInput("nnumber", "Sample size", 50),
        numericInput("iternumber", "Iterations", 1000)
      )
    ),
    fluidRow(
      box(width=12,
        solidHeader = FALSE,
        collapsible = FALSE,
        background = "black",
        numericInput("mu1", withMathJax("$$\\mu_1$$"), 0),
        numericInput("mu2", withMathJax("$$\\mu_2$$"), 0)
      ),
    ),
    fluidRow(
      box(width=12,
        solidHeader = FALSE,
        collapsible = FALSE,
        background = "black",
        title = withMathJax("$$\\Sigma$$"),
        matrixInput(
          "sigma",
          class = "numeric",
          value = sigma
        )
      ),
    )
  ),
  dashboardBody(
    # Boxes need to be put in a row (or column)
    fluidRow(
      box(width=12, plotOutput("plot1", height = 500)),
    )
  )
)

server <- function(input, output, session) {

  output$plot1 <- renderPlot({
    lsigma = input$sigma
    diag_update = NA
    if (sigma[1,2] != lsigma[1,2]){
      diag_update = lsigma[1,2]
    }else if (sigma[1,2] != lsigma[2,1]){
      diag_update = lsigma[2,1]
    }
    if (!is.na(diag_update)) {
      sigma[1,2] = diag_update
      sigma[2,1] = diag_update
      updateMatrixInput(session, 'sigma', value=sigma)
    }
    if (is.positive.definite(sigma)){
      mui = c(input$mu1, input$mu2)
      xbarthetadist(n=input$nnumber, iter=input$iternumber, mu=mui, sigma=sigma)
    }else{
      print('Sigma is not positive definite!')
    }
  })
}

shinyApp(ui, server)
