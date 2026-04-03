library(shiny)
library(ggplot2)
library(jsonlite)
library(TDA)
library(TDAmapper)
library(igraph)

# 1. User Interface
ui <- fluidPage(
  
  titlePanel(
    tagList(
      "TopologgeR: 2D Topology Explorer",
      br(), # Line break
      tags$span(
        style = "font-size: 16px; font-weight: normal; color: #555;",
        "New to Topological Data Analysis (TDA)? ",
        tags$a(href = "https://github.com/leliavski/tda-app", 
               target = "_blank", # Opens in a new tab
               "Read the beginner's guide here.")
      )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      # Reduced default N heavily to prevent crashing during VR computation
      sliderInput("n_points", "Number of points (N):", 
                  min = 20, max = 150, value = 60, step = 10),
      sliderInput("noise", "Noise parameter:", 
                  min = 0, max = 20, value = 5, step = 1),
      
      hr(),
      h4("Topology Parameters"),
      sliderInput("max_scale", "Max Filtration Scale (Max Epsilon):", 
                  min = 10, max = 200, value = 100, step = 10),
      sliderInput("epsilon", "Visualize VR Complex at Epsilon:", 
                  min = 0, max = 200, value = 30, step = 2),
      
      hr(),
      h4("Mapper Parameters"),
      selectInput("mapper_filter", "Filter Function (Lens):", 
                  choices = c("X-Coordinate (Left to Right)" = "x",
                              "Y-Coordinate (Bottom to Top)" = "y",
                              "Radial Distance (Center Outward)" = "radial",
                              "Eccentricity (Branch Extremes)" = "eccentricity")),
      sliderInput("mapper_intervals", "Number of Intervals (Resolution):", 
                  min = 2, max = 30, value = 10, step = 1),
      sliderInput("mapper_overlap", "Percent Overlap:", 
                  min = 10, max = 80, value = 50, step = 5),
      
      actionButton("submit", "Submit Shape", class = "btn-primary"),
      actionButton("clear", "Clear Canvas"),
      
      hr(),
      helpText("Draw a shape and click Submit.")
    ),
    
    mainPanel(
      fluidRow(
        column(6, 
               h4("Drawing Canvas"),
               tags$div(
                 tags$canvas(id = "drawCanvas", width = "350", height = "350",
                             style = "border: 2px solid #ccc; cursor: crosshair; touch-action: none;")
               )
        ),
        column(6, 
               h4("Data & VR Complex (1-skeleton)"),
               plotOutput("pointPlot", width = "350px", height = "350px")
        )
      ),
      fluidRow(
        column(6, 
               h4("Persistence Barcodes"),
               plotOutput("barcodePlot", height = "300px")
        ),
        column(6, 
               h4("Persistence Diagram"),
               plotOutput("diagramPlot", height = "300px")
        )
      ),
      fluidRow(
        column(6, 
               h4("Mapper Pre-images & Coverings"),
               plotOutput("mapperDataPlot", height = "400px")
        ),
        column(6, 
               h4("Mapper Graph"),
               plotOutput("mapperPlot", height = "400px")
        )
      ),
      
      # JavaScript Bridge (Unchanged from previous version)
      tags$script(HTML("
        var canvas = document.getElementById('drawCanvas');
        var ctx = canvas.getContext('2d');
        ctx.lineWidth = 2;
        ctx.lineCap = 'round';
        
        var drawing = false;
        var coords = [];
        var stroke_id = 0; // Track individual pen strokes
      
        canvas.addEventListener('mousedown', function(e) {
          drawing = true;
          stroke_id++; // Increment ID for every new stroke
          var rect = canvas.getBoundingClientRect();
          var x = e.clientX - rect.left;
          var y = e.clientY - rect.top;
          ctx.beginPath();
          ctx.moveTo(x, y);
          coords.push({x: x, y: y, stroke: stroke_id});
        });
      
        canvas.addEventListener('mousemove', function(e) {
          if(drawing) {
            var rect = canvas.getBoundingClientRect();
            var x = e.clientX - rect.left;
            var y = e.clientY - rect.top;
            ctx.lineTo(x, y);
            ctx.stroke();
            coords.push({x: x, y: y, stroke: stroke_id});
          }
        });
      
        canvas.addEventListener('mouseup', function(e) {
          drawing = false;
          Shiny.setInputValue('draw_coords', JSON.stringify(coords));
        });
      
        document.getElementById('clear').addEventListener('click', function() {
          ctx.clearRect(0, 0, canvas.width, canvas.height);
          coords = [];
          stroke_id = 0;
          Shiny.setInputValue('draw_coords', '');
        });
      "))
    )
  )
)

# 2. Server Logic
server <- function(input, output, session) {
  
  # Reactive block 1: Generate the point cloud grouped by strokes
  processed_data <- eventReactive(input$submit, {
    req(input$draw_coords)
    if (input$draw_coords == "") return(NULL)
    
    df <- fromJSON(input$draw_coords)
    if (nrow(df) < 2) return(NULL)
    
    # Flip Y-axis
    df$y <- 350 - df$y
    
    # Split the dataframe into a list of dataframes, one for each stroke
    strokes <- split(df, df$stroke)
    
    # Calculate the physical length of each stroke
    stroke_lengths <- sapply(strokes, function(s) {
      if(nrow(s) < 2) return(0)
      sum(sqrt(diff(s$x)^2 + diff(s$y)^2))
    })
    
    total_length <- sum(stroke_lengths)
    req(total_length > 0)
    
    # Interpolate each stroke independently
    sampled_list <- lapply(names(strokes), function(st_name) {
      s <- strokes[[st_name]]
      slen <- stroke_lengths[[st_name]]
      
      # Skip if it's just a single dot or zero length
      if (slen == 0 || nrow(s) < 2) return(NULL)
      
      # Apportion the N points based on how long this stroke is compared to the total
      # (Guarantee at least 2 points per valid stroke to prevent errors)
      n_this_stroke <- max(2, round(input$n_points * (slen / total_length)))
      
      # Interpolate just this specific stroke
      dist <- sqrt(diff(s$x)^2 + diff(s$y)^2)
      cumdist <- c(0, cumsum(dist))
      
      target_dists <- seq(0, slen, length.out = n_this_stroke)
      sampled_x <- approx(cumdist, s$x, target_dists)$y
      sampled_y <- approx(cumdist, s$y, target_dists)$y
      
      data.frame(x = sampled_x, y = sampled_y)
    })
    
    # Bind the independent strokes back into a single dataframe
    final_df <- do.call(rbind, sampled_list)
    
    # Apply global Gaussian noise
    final_df$x <- final_df$x + rnorm(nrow(final_df), mean = 0, sd = input$noise)
    final_df$y <- final_df$y + rnorm(nrow(final_df), mean = 0, sd = input$noise)
    
    return(final_df)
  })
  
  # Reactive block 2: Calculate Persistent Homology
  topology_data <- reactive({
    df <- processed_data()
    req(df)
    
    # Calculate Vietoris-Rips filtration
    diag <- ripsDiag(X = df, maxdimension = 1, maxscale = input$max_scale, 
                     library = "GUDHI", printProgress = FALSE)
    
    diag_mat <- diag$diagram
    
    # Safety check 1: If no topological features are found at all
    if (is.null(diag_mat) || length(diag_mat) == 0) {
      return(data.frame(dimension = factor(), Birth = numeric(), Death = numeric(), id = integer()))
    }
    
    # Safety check 2: If the matrix collapses into a vector (only 1 feature found)
    # We force it back into a 3-column format
    if (!is.matrix(diag_mat)) {
      diag_mat <- matrix(diag_mat, ncol = 3, byrow = TRUE)
    }
    
    # Safely construct the data frame column by column
    diag_df <- data.frame(
      dimension = as.factor(diag_mat[, 1]),
      Birth = as.numeric(diag_mat[, 2]),
      Death = as.numeric(diag_mat[, 3])
    )
    
    # Sort for neater barcode plotting
    diag_df <- diag_df[order(diag_df$dimension, diag_df$Birth), ]
    diag_df$id <- seq_len(nrow(diag_df))
    
    return(diag_df)
  })
  
  # Reactive block 3: Centralized Mapper Calculation
  mapper_model <- reactive({
    df <- processed_data()
    req(df)
    
    d_mat <- dist(df)
    full_d_mat <- as.matrix(d_mat)
    
    # Apply the selected Filter Function
    if (input$mapper_filter == "x") {
      filter_val <- df$x
    } else if (input$mapper_filter == "y") {
      filter_val <- df$y
    } else if (input$mapper_filter == "radial") {
      center_x <- mean(df$x)
      center_y <- mean(df$y)
      filter_val <- sqrt((df$x - center_x)^2 + (df$y - center_y)^2)
    } else if (input$mapper_filter == "eccentricity") {
      filter_val <- apply(full_d_mat, 1, max)
    }
    
    m1 <- tryCatch({
      mapper1D(
        distance_matrix = d_mat,
        filter_values = filter_val,
        num_intervals = input$mapper_intervals,
        percent_overlap = input$mapper_overlap,
        num_bins_when_clustering = 10
      )
    }, error = function(e) return(NULL))
    
    # Return everything the plots will need
    list(m1 = m1, df = df, filter_val = filter_val)
  })
  
  # Output 1: Data Plot + VR Complex 1-Skeleton
  output$pointPlot <- renderPlot({
    df <- processed_data()
    req(df)
    
    # Calculate distance matrix to find edges for the VR complex
    d_mat <- as.matrix(dist(df))
    
    # Find all pairs of points where distance <= input$epsilon
    edges <- which(d_mat <= input$epsilon & upper.tri(d_mat), arr.ind = TRUE)
    
    p <- ggplot(df, aes(x = x, y = y))
    
    # If there are edges to draw, add them under the points
    if (nrow(edges) > 0) {
      edge_df <- data.frame(
        x = df$x[edges[, 1]], y = df$y[edges[, 1]],
        xend = df$x[edges[, 2]], yend = df$y[edges[, 2]]
      )
      p <- p + geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend),
                            color = "skyblue", alpha = 0.4, size = 0.5)
    }
    
    # Add the actual data points
    p + geom_point(color = "#2c3e50", size = 2) +
      coord_fixed(xlim = c(0, 350), ylim = c(0, 350)) +
      theme_minimal()
  })
  
  # Output 2: Persistence Barcode
  output$barcodePlot <- renderPlot({
    df <- topology_data()
    req(df)
    
    ggplot(df, aes(x = Birth, xend = Death, y = id, yend = id, color = dimension)) +
      geom_segment(size = 1.5) +
      geom_vline(xintercept = input$epsilon, linetype = "dashed", color = "red", alpha = 0.6) +
      scale_color_manual(values = c("0" = "black", "1" = "red"), 
                         labels = c("H0 (Components)", "H1 (Loops)")) +
      labs(x = "Filtration Scale (Epsilon)", y = "Feature ID") +
      theme_minimal() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  })
  
  # Output 3: Persistence Diagram
  output$diagramPlot <- renderPlot({
    df <- topology_data()
    req(df)
    
    ggplot(df, aes(x = Birth, y = Death, color = dimension, shape = dimension)) +
      geom_abline(slope = 1, intercept = 0, color = "grey") + # The diagonal
      geom_point(size = 3, alpha = 0.7) +
      scale_color_manual(values = c("0" = "black", "1" = "red"),
                         labels = c("H0 (Components)", "H1 (Loops)")) +
      labs(x = "Birth", y = "Death") +
      coord_fixed(xlim = c(0, input$max_scale), ylim = c(0, input$max_scale)) +
      theme_minimal()
  })
  
  # Output 4: Mapper Graph
  output$mapperPlot <- renderPlot({
    model <- mapper_model()
    req(model)
    m1 <- model$m1
    
    # Dynamic Title
    plot_title <- switch(input$mapper_filter,
                         "x" = "Mapper Graph (Lens: X-Coordinate)",
                         "y" = "Mapper Graph (Lens: Y-Coordinate)",
                         "radial" = "Mapper Graph (Lens: Radial Distance)",
                         "eccentricity" = "Mapper Graph (Lens: Eccentricity)")
    
    if (is.null(m1) || length(m1$adjacency) == 0) {
      return(ggplot() + 
               annotate("text", x = 0, y = 0, label = "No graph generated.") + 
               theme_void())
    }
    
    g <- graph.adjacency(m1$adjacency, mode = "undirected")
    V(g)$size <- sapply(m1$points_in_vertex, length)
    
    node_levels <- m1$level_of_vertex
    color_palette <- colorRampPalette(c("lightblue", "darkblue"))(input$mapper_intervals)
    V(g)$color <- color_palette[node_levels]
    
    plot(g, layout = layout_with_kk, vertex.label = NA,
         vertex.frame.color = "white", edge.color = "gray50",
         main = plot_title)
  })
  
  # Output 5: Mapper Pre-images on Data
  output$mapperDataPlot <- renderPlot({
    model <- mapper_model()
    req(model)
    m1 <- model$m1
    df <- model$df
    
    # Base plot: just the grey points
    p <- ggplot(df, aes(x = x, y = y)) + 
      geom_point(color = "grey80", size = 1.5) +
      coord_fixed(xlim = c(0, 350), ylim = c(0, 350)) +
      theme_minimal()
    
    if (is.null(m1) || length(m1$adjacency) == 0) return(p)
    
    # 1. Extract all clustered points
    clustered_pts_list <- lapply(seq_along(m1$points_in_vertex), function(i) {
      pts <- df[m1$points_in_vertex[[i]], , drop = FALSE]
      pts$vertex <- i
      pts$level <- m1$level_of_vertex[i]
      return(pts)
    })
    clustered_pts_df <- do.call(rbind, clustered_pts_list)
    
    # 2. Calculate Convex Hulls for clusters with 3 or more points
    hull_list <- lapply(clustered_pts_list, function(pts) {
      if (nrow(pts) >= 3) {
        hull_idx <- chull(pts$x, pts$y)
        return(pts[hull_idx, ])
      }
      return(NULL)
    })
    hull_df <- do.call(rbind, hull_list)
    
    # 3. Color Palette matching the Graph
    color_palette <- colorRampPalette(c("lightblue", "darkblue"))(input$mapper_intervals)
    
    # 4. Add colored polygons (hulls) and highlighted points
    if (!is.null(hull_df) && nrow(hull_df) > 0) {
      p <- p + geom_polygon(data = hull_df, aes(x = x, y = y, group = vertex, fill = level),
                            alpha = 0.3, color = "black", linewidth = 0.2)
    }
    
    p <- p + geom_point(data = clustered_pts_df, aes(x = x, y = y, color = level), size = 2) +
      scale_fill_gradientn(colors = color_palette, limits = c(1, input$mapper_intervals), guide = "none") +
      scale_color_gradientn(colors = color_palette, limits = c(1, input$mapper_intervals), guide = "none") +
      labs(x = "", y = "")
    
    return(p)
  })
  
}

# 3. Run the App
shinyApp(ui = ui, server = server)