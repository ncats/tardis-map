shinyUI(pageWithSidebar(
    # Application title
    headerPanel(HTML("<img src='tardis.png' width=40>TarDis - Explore Target-Disease Associations (v0.1.5)")),
    
    # Sidebar with a slider input for number of bins
    sidebarPanel(
        HTML("This is an interactive heatmap representation of target-disease associations, obtained by a text mining protocol developed as part of the <a rhef='http://diseases.jensenlab.org/Search'>DISEASES</a> database. The value of each association is a z-score. By default, targets are aggregated at various levels defined by the <a href='http://targetcentral.ws/index'>IDG KMC</a> (such as Target Development Level and protein family)<p>"),
        HTML("<br>You can provide a comma separated list of HGNC gene symbols in the text area to regenerate the visualization for those specific targets<p>"),
        selectInput("targetAggLevel", "Target Level", c("TDL", "Family"), selected='TDL'),
        selectInput("diseaseAggLevel", "Disease Level", c("Level1","Level2", "Level3", "Raw"), selected='Level1'),
        textAreaInput("syms", "Target Symbols", rows = 5)
        
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
        plotlyOutput("plotlyHM", width = "850px", height = "1600px")        
    )
    
))
