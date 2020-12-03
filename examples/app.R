library(shiny)

ui <- fluidPage(
  mainPanel(
    plotOutput('plot')
  )
)

#   pageWithSidebar(
#
#   # headerPanel('Group9LinearModel'),
#
#   # sidebarPanel(
#     # numericInput('clusters', 'Cluster count', 3, min = 1, max = 9)
#   # ),
#
#   mainPanel(
#     plotOutput('plot')
#   )
# )

server <- function(input, output, session) {
  library(MASS)
  library(Group9LinearModel)
  data(Boston)

  fit = myLm(Boston$crim, Boston[c("age", "medv")])


  # df <- iris[, c("Sepal.Length", "Sepal.Width")]

  # clusters <- reactive({kmeans(df, input$clusters)})

  output$plot <- renderPlot({
    plot(fit)
  })
}

shinyApp(ui=ui, server=server)

