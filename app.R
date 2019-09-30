library(shiny)
library(Seurat)
library(ggplot2)

load('experiment_merged.RData')
genes = experiment.merged@assays$RNA
meta_nums <- colnames(dplyr::select_if(experiment.merged@meta.data, is.numeric))
meta_cats <- colnames(dplyr::select_if(experiment.merged@meta.data, is.factor))
pcs <- list('PC_1','PC_2','PC_3','PC_4','PC_5','PC_6','PC_7','PC_8','PC_9')
agg_cats <- colnames(dplyr::select_if(experiment.merged@meta.data, is.factor))


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
  observe({
    updateSelectInput(session, "numeric_b",
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
      c("Xkr4", "Tcea1"), blend=TRUE
    )
  })

  output$CategoricalPlot <- renderPlot({
    DimPlot(object = experiment.merged, group.by=input$categorical, pt.size=0.5, do.label = TRUE, reduction = "tsne", label = T)
  })

  output$ViolinPlot <- renderPlot({
    Idents(experiment.merged) <- input$categorical
    VlnPlot(object =  experiment.merged, features = c(input$numeric, input$numeric2), pt.size = 0.05)
  })
  
  output$MarkerSet <- renderPlot({
    Idents(experiment.merged) <- input$categorical_b
    markers = input$numeric_b
    expr.cutoff = 3
    widedat <- FetchData(experiment.merged, markers)
    widedat$Cluster <- Idents(experiment.merged)
    longdat <- gather(widedat, key = "Gene", value = "Expression", -Cluster)
    longdat$Is.Expressed <- ifelse(longdat$Expression > expr.cutoff, 1, 0)
    longdat$Cluster <- factor(longdat$Cluster)
    longdat$Gene <- factor(longdat$Gene)
    
    # Need to summarize into average expression, pct expressed (which is also an average)
    plotdat <- group_by(longdat, Gene, Cluster) %>% summarize(`Percentage of Expressed Cells` = mean(Is.Expressed), `Mean Expression` = mean(Expression))
    
    
    ggplot(plotdat, aes(x = Gene, y = Cluster)) +
      geom_point(aes(size = `Percentage of Expressed Cells`, col = `Mean Expression`)) +
      labs(size = "Percentage\nof Expressed\nCells", col = "Mean\nExpression", x = NULL) +
      scale_color_gradient(low = "grey", high = "slateblue4") + theme_grey(base_size = 15) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }, height = 1000, width = 900 )

  output$summary <- renderPrint({
    summary(d())
  })
  
  # Generate an HTML table view of the data ----
  output$table <- renderTable({
    d()
  })
  
}

ui <- fluidPage(
  
  # App title ----
  titlePanel("scRNA Seurat Analysis"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
        conditionalPanel(condition = "input.tabselected == -999",
                selectInput("dataset", "Numeric Analysis Type:",
                            c('Genes', 'Numeric Metadata','PCs')),
                selectInput("categorical", "Identity:",
                            c(meta_cats)),
                selectInput("numeric", "Primary Numeric:", ""),
                
                selectInput('numeric2', 'Secondary Numeric', "")
        ),
        
        conditionalPanel(condition = "input.tabselected == 2",
                         selectInput("categorical_b", "Identity:",
                                     c(agg_cats)),
                         selectInput("numeric_b", "Primary Numeric:", "", multiple=TRUE)
         ),
        conditionalPanel(condition = "input.tabselected == 3")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Marker Genes (TSNE)", value=-999,
                               #h3(textOutput("caption")),
                               plotOutput("MarkerGenePlot"),
                               plotOutput("ViolinPlot"),
                               #h3(textOutput("caption2")),
                               plotOutput("CategoricalPlot")
                               #h3(textOutput("caption3")),
                           ),
                  
                  tabPanel("Marker Set (Grid)", value=2,
                               plotOutput("MarkerSet")
                          ),

                  tabPanel("Documentation", value=3,
                           includeMarkdown("docs/testing.md"),
                           includeMarkdown("docs/todo.md")
                         ),
                  
                  id = "tabselected"
      )
    )
  )
)
shinyApp(ui, server)