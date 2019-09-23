library(shiny)
library(Seurat)
library(ggplot2)

load('experiment_merged.RData')
genes = experiment.merged@assays$RNA
meta_nums <- colnames(dplyr::select_if(experiment.merged@meta.data, is.numeric))
meta_cats <- colnames(dplyr::select_if(experiment.merged@meta.data, is.factor))
pcs <- list('PC_1','PC_2','PC_3','PC_4','PC_5','PC_6','PC_7','PC_8','PC_9')


server = function(input, output, session){
  outVar = reactive({
    if (input$dataset == 'Genes'){mydata=row.names(genes)}
    else if (input$dataset == 'Numeric Metadata') {mydata=meta_nums}
    else if (input$dataset == 'PCs') {mydata=pcs}
    mydata
  })
  observe({
    updateSelectInput(session, "numeric",
                      choices = outVar()
  )})
  observe({
    updateSelectInput(session, "numeric2",
                      choices = outVar()
  )})
  formulaText <- reactive({
    paste("Marker Gene ~", input$numeric)
  })
  formulaText2 <- reactive({
    paste("Identity ~", input$categorical)
  })
  formulaText3 <- reactive({
    paste("Identity vs Marker Genes ~", input$categorical, input$numeric, input$numeric2)
  })

  # Return the formula text for printing as a caption ----
  output$caption <- renderText({
    formulaText()
  })

  output$caption2 <- renderText({
    formulaText2()
  })

  output$caption3 <- renderText({
    formulaText3()
  })
  
  # Get the graphs
  output$MarkerGenePlot <- renderPlot({
    FeaturePlot(
      experiment.merged,
      input$numeric,
      cols = c("lightgrey", "blue")
    )
  })

  output$CategoricalPlot <- renderPlot({
    DimPlot(object = experiment.merged, group.by=input$categorical, pt.size=0.5, do.label = TRUE, reduction = "tsne", label = T)
  })

  output$ViolinPlot <- renderPlot({
    Idents(experiment.merged) <- input$categorical
    VlnPlot(object =  experiment.merged, features = c(input$numeric, input$numeric2), pt.size = 0.05)
  })
}

ui <- fluidPage(

  # App title
  titlePanel("Identifying Marker Genes"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      selectInput("dataset", "Numeric Analysis Type:",
                  c('Genes', 'Numeric Metadata','PCs')),
      selectInput("categorical", "Identity:",
                  c(meta_cats)),
      selectInput("numeric", "Primary Numeric:", ""),

      selectInput('numeric2', 'Secondary Numeric', "")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      #h3(textOutput("caption")),
      plotOutput("MarkerGenePlot"),
      #h3(textOutput("caption2")),
      plotOutput("CategoricalPlot"),
      #h3(textOutput("caption3")),
      plotOutput("ViolinPlot")
    )
  )
)

shinyApp(ui, server)