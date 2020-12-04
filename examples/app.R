library(shiny)

ui <- fluidPage(
  mainPanel(
    plotOutput('plot')
  )
)

server <- function(input, output, session) {
  library(MASS)
  library(Group9LinearModel)
  data(Boston)

  fit = myLm(Boston$crim, Boston[c("age", "medv")])

  output$plot <- renderPlot({
    plot(fit)
  })
}

shinyApp(ui=ui, server=server)

