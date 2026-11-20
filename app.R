# ======== PACKAGES ========
suppressMessages({
  library(shiny)
  library(tidyverse)
  library(uwot)
  library(NMF)
  library(mclust)
})

# ======== LOAD DATA ========
TBnanostring <- readRDS("TBnanostring.rds")
x <- dplyr::select(TBnanostring, -TB_Status) %>%
  as.matrix()
y <- TBnanostring$TB_Status

# ======== DIMENSION REDUCTION (Computed Once) ========
# PCA: First two PC scores
pca <- prcomp(x)
pca_2d <- pca$x[ , 1:2]

# UMAP (Use same settings as Task 2.1)
set.seed(2026)
umap_2d <- umap(x, n_neighbors = 15)

# NMF (Use same settings as Tasks 3.1-3.2)
x_nmf <- x # Copy gene expression data
for (n in seq(1,ncol(x))) {
  x_nmf[ ,n] <- x_nmf[ ,n] - min(x[ ,n])    # Subtract minimum from corresponding gene column
}
set.seed(1952) # Set seed for reproducibility
nmf_res <- nmf(x_nmf, rank = 2, nrun = 10) # Perform NMF
nmf_2d <- basis(nmf_res) # Extract basis matrix

# Combine the reduced dimensions 
dims <- list(PCA = pca_2d, UMAP = umap_2d, NMF = nmf_2d)

# Apply same col names to each so serve can handle them same way
dims <- lapply(dims, function(m) {
  m <- as.data.frame(m)
  colnames(m) <- c("Dim1","Dim2")
  m
})

# Store axis labels for each method
axes <- list(PCA = c("PC1","PC2"),
             UMAP = c("UMAP1","UMAP2"),
             NMF = c("NMF1","NMF2"))



# ======== UI ========
ui <- fluidPage(
  
  titlePanel("TBnanostring: Dimension Reduction and K-Means Clustering"),
  
  sidebarLayout(
    # Inputs 
    sidebarPanel(
      selectInput("method", "Dimension Reduction Method", 
                  choices = c("PCA" = "PCA",
                              "UMAP" = "UMAP",
                              "NMF (Rank 2)" = "NMF")),
      numericInput("clusters", "Cluster Count (k)", 2, 
                   min = 1, max = 9)
    ),
    
    # Outputs
    mainPanel(
      plotOutput("plot1"),
      h4("Clusters vs. TB Status"),
      tableOutput("table1"),
      textOutput("ari")
    )
  )
)


# ======== SERVER ========
server <- function(input, output, session) {
  
  # Pin 2D coordinates for selected method
  selectedData <- reactive({
    dims[[input$method]]
  })
  
  # Run k-means on the reduced dimensions
  clusters <- reactive({
    set.seed(1999)
    kmeans(selectedData(), 
           centers = input$clusters,
           nstart = 25)
  })
  
  # Create plot of two dimensions, colored by cluster and shaped by TB status
  output$plot1 <- renderPlot({
    plot_df <- data.frame(selectedData(),
                          Cluster = factor(clusters()$cluster),
                          TB_Status = y)
    centers_df <- as.data.frame(clusters()$centers)
    xy_labs <- axes[[input$method]]
    
    ggplot(plot_df, aes(Dim1, Dim2, color = Cluster, shape = TB_Status)) +
      geom_point(size = 3) +
      geom_point(data = centers_df, aes(Dim1, Dim2),
                 inherit.aes = FALSE,
                 shape = 4, size = 6, stroke = 2) +
      coord_fixed(ratio = 1) +
      labs(title = paste0(input$method, " with k-means (k = ", input$clusters, ")"),
           x = xy_labs[1],
           y = xy_labs[2],
           color = "Clusters",
           shape = "TB Status") +
      theme_bw()
  })
  
  # Table: cross-tabulation of clusters vs. TB status
  output$table1 <- renderTable({
    tab <- table(clusters()$cluster, y) 
    tab_df <- as.data.frame.matrix(tab)
    cbind(Cluster = rownames(tab_df),tab_df)
  })
  
  # Compute ARI between clusters and TB status
  output$ari <- renderText({
    ari <- adjustedRandIndex(clusters()$cluster, y)
    paste("Adjusted Rand Index (Clusters vs. TB Status):", round(ari,3))
  })
  
}

shinyApp(ui = ui, server = server)