library(shiny)
library(ggplot2)

# Define UI for application
ui <- fluidPage(
  titlePanel("Michaelis-Menten Kinetics Simulation"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("vmax", "Vmax:", 1, min = 0),
      numericInput("km", "Km:", 1, min = 0),
      sliderInput("substrate_range", "Substrate Concentration Range (S):", 
                  min = 0, max = 100, value = c(0, 10)),
      numericInput("comp_inh", "Competitive Inhibitor Concentration (I):", 0, min = 0),
      numericInput("noncomp_inh", "Non-competitive Inhibitor Concentration (I):", 0, min = 0),
      numericInput("ki_comp", "Inhibition constant (Ki) for competitive:", 1, min = 0),
      numericInput("ki_noncomp", "Inhibition constant (Ki) for noncompetitive:", 1, min = 0),
      
      checkboxGroupInput("inhibitions", "Select Inhibitions to Display:",
                         choices = list("Non-Inhibited" = "none",
                                        "Competitive" = "competitive",
                                        "Non-Competitive" = "noncompetitive",
                                        "Both Combined" = "combined"),
                         selected = c("none", "competitive", "noncompetitive")),
      
      downloadButton("downloadPlot", "Download EPS")
    ),
    
    mainPanel(
      plotOutput("mm_plot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Function to calculate Michaelis-Menten velocity
  mm_velocity <- function(S, vmax, km) {
    (vmax * S) / (km + S)
  }
  
  # Competitive inhibition equation
  competitive_inhibition <- function(S, vmax, km, I, Ki) {
    (vmax * S) / (km * (1 + I/Ki) + S)
  }
  
  # Noncompetitive inhibition equation
  noncompetitive_inhibition <- function(S, vmax, km, I, Ki) {
    (vmax * S) / (km + S) * (1 / (1 + I/Ki))
  }
  
  # Combined competitive and noncompetitive inhibition
  combined_inhibition <- function(S, vmax, km, comp_I, comp_Ki, noncomp_I, noncomp_Ki) {
    (vmax * S) / ((km * (1 + comp_I / comp_Ki)) + S) * (1 / (1 + noncomp_I / noncomp_Ki))
  }
  
  # Reactive expression to generate the plot based on inputs
  output$mm_plot <- renderPlot({
    S <- seq(input$substrate_range[1], input$substrate_range[2], length.out = 100)
    
    # Initialize dataframe for plotting
    data <- data.frame(S = S)
    
    # Add non-inhibited curve
    if ("none" %in% input$inhibitions) {
      data$None <- mm_velocity(S, input$vmax, input$km)
    }
    
    # Add competitive inhibition curve
    if ("competitive" %in% input$inhibitions) {
      data$Competitive <- competitive_inhibition(S, input$vmax, input$km, input$comp_inh, input$ki_comp)
    }
    
    # Add noncompetitive inhibition curve
    if ("noncompetitive" %in% input$inhibitions) {
      data$Noncompetitive <- noncompetitive_inhibition(S, input$vmax, input$km, input$noncomp_inh, input$ki_noncomp)
    }
    
    # Add combined inhibition curve
    if ("combined" %in% input$inhibitions) {
      data$Combined <- combined_inhibition(S, input$vmax, input$km, input$comp_inh, input$ki_comp, input$noncomp_inh, input$ki_noncomp)
    }
    
    # Plotting
    plot_data <- reshape2::melt(data, id = "S", variable.name = "Condition", value.name = "Velocity")
    ggplot(plot_data, aes(x = S, y = Velocity, color = Condition)) +
      geom_line(size = 1) +
      labs(title = "Michaelis-Menten Kinetics", x = "Substrate Concentration (S)", y = "Reaction Velocity (V)") +
      theme_minimal()
  })
  
  # Function to download the plot as an EPS file
  output$downloadPlot <- downloadHandler(
    filename = function() { "mm_plot.eps" },
    content = function(file) {
      S <- seq(input$substrate_range[1], input$substrate_range[2], length.out = 100)
      
      # Initialize dataframe for plotting
      data <- data.frame(S = S)
      
      # Add non-inhibited curve
      if ("none" %in% input$inhibitions) {
        data$None <- mm_velocity(S, input$vmax, input$km)
      }
      
      # Add competitive inhibition curve
      if ("competitive" %in% input$inhibitions) {
        data$Competitive <- competitive_inhibition(S, input$vmax, input$km, input$comp_inh, input$ki_comp)
      }
      
      # Add noncompetitive inhibition curve
      if ("noncompetitive" %in% input$inhibitions) {
        data$Noncompetitive <- noncompetitive_inhibition(S, input$vmax, input$km, input$noncomp_inh, input$ki_noncomp)
      }
      
      # Add combined inhibition curve
      if ("combined" %in% input$inhibitions) {
        data$Combined <- combined_inhibition(S, input$vmax, input$km, input$comp_inh, input$ki_comp, input$noncomp_inh, input$ki_noncomp)
      }
      
      # Plotting
      plot_data <- reshape2::melt(data, id = "S", variable.name = "Condition", value.name = "Velocity")
      p <- ggplot(plot_data, aes(x = S, y = Velocity, color = Condition)) +
        geom_line(size = 1) +
        labs(title = "Michaelis-Menten Kinetics", x = "Substrate Concentration (S)", y = "Reaction Velocity (V)") +
        theme_minimal()
      
      # Save as EPS
      postscript(file)
      print(p)
      dev.off()
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
