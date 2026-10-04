

# Install these packages if you don't already have them:
# install.packages("shiny")
# install.packages("ggplot2")

library(shiny)
library(ggplot2)


ui <- fluidPage(
  
  titlePanel("Lotka–Volterra Two-Species Competition Model"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      # -------------------------
      # Species 1 parameters
      # -------------------------
      
      h4("Species 1:", style = "color: #0000FF;"),
      
      sliderInput(
        inputId = "r1",
        label = "Growth rate (r₁):",
        min = 0,
        max = 3,
        value = 1,
        step = 0.01
      ),
      
      sliderInput(
        inputId = "K1",
        label = "Carrying capacity (K₁):",
        min = 100,
        max = 200,
        value = 150,
        step = 1
      ),
      
      sliderInput(
        inputId = "alpha",
        label = "Competition coefficient (α):",
        min = 0,
        max = 3,
        value = 1,
        step = 0.01
      ),
      
      sliderInput(
        inputId = "N1",
        label = "Starting population (N₁):",
        min = 0,
        max = 400,
        value = 50,
        step = 1
      ),
      
      hr(),
      
      # -------------------------
      # Species 2 parameters
      # -------------------------
      
      h4("Species 2", style = "color: #FF0000;"),
      
      sliderInput(
        inputId = "r2",
        label = "Growth rate (r₂):",
        min = 0,
        max = 3,
        value = 1,
        step = 0.01
      ),
      
      sliderInput(
        inputId = "K2",
        label = "Carrying capacity (K₂):",
        min = 100,
        max = 200,
        value = 150,
        step = 1
      ),
      
      sliderInput(
        inputId = "beta",
        label = "Competition coefficient (β):",
        min = 0,
        max = 3,
        value = 1,
        step = 0.01
      ),
      
      sliderInput(
        inputId = "N2",
        label = "Starting population (N₂):",
        min = 0,
        max = 400,
        value = 50,
        step = 1
      ),
      
      hr(),
      
      # -------------------------
      # Time
      # -------------------------
      
      sliderInput(
        inputId = "time",
        label = "Simulation time:",
        min = 10,
        max = 100,
        value = 50,
        step = 1
      )
    ),
    
    # -------------------------
    # Graph
    # -------------------------
    
    mainPanel(
      plotOutput(
        outputId = "competitionPlot",
        height = "600px"
      )
    )
  )
)


server <- function(input, output) {
  
  output$competitionPlot <- renderPlot({
    
    # -------------------------
    # Initial population sizes
    # -------------------------
    
    initial_state <- c(
      N1 = input$N1,
      N2 = input$N2
    )
    
    
    # -------------------------
    # Define the model
    # -------------------------
    
    competition_model <- function(state, parameters) {
      
      N1 <- state[1]
      N2 <- state[2]
      
      r1 <- parameters["r1"]
      K1 <- parameters["K1"]
      alpha <- parameters["alpha"]
      
      r2 <- parameters["r2"]
      K2 <- parameters["K2"]
      beta <- parameters["beta"]
      
      # Species 1
      dN1 <- r1 * N1 * (
        1 - (N1 + alpha * N2) / K1
      )
      
      # Species 2
      dN2 <- r2 * N2 * (
        1 - (N2 + beta * N1) / K2
      )
      
      return(c(dN1, dN2))
    }
    
    
    # -------------------------
    # Model parameters
    # -------------------------
    
    parameters <- c(
      r1 = input$r1,
      K1 = input$K1,
      alpha = input$alpha,
      r2 = input$r2,
      K2 = input$K2,
      beta = input$beta
    )
    
    
    # -------------------------
    # Time points
    # -------------------------
    
    times <- seq(
      from = 0,
      to = input$time,
      by = 0.1
    )
    
    
    # -------------------------
    # Runge-Kutta 4th-order solver
    # -------------------------
    
    rk4 <- function(initial_state, times, model, parameters) {
      
      n <- length(times)
      n_states <- length(initial_state)
      
      # Matrix to store results
      result <- matrix(
        0,
        nrow = n,
        ncol = n_states + 1
      )
      
      # Column 1 = time
      # Remaining columns = state variables
      result[, 1] <- times
      
      # Initial state
      result[1, 2:(n_states + 1)] <- initial_state
      
      
      # Step through time
      for (i in 1:(n - 1)) {
        
        # Current time step
        h <- times[i + 1] - times[i]
        
        # Current state
        current_state <- result[i, 2:(n_states + 1)]
        
        
        # RK4 calculations
        k1 <- model(
          current_state,
          parameters
        )
        
        k2 <- model(
          current_state + (h / 2) * k1,
          parameters
        )
        
        k3 <- model(
          current_state + (h / 2) * k2,
          parameters
        )
        
        k4 <- model(
          current_state + h * k3,
          parameters
        )
        
        
        # Calculate next state
        next_state <- current_state +
          (h / 6) * (
            k1 +
              2 * k2 +
              2 * k3 +
              k4
          )
        
        
        # Store result
        result[i + 1, 2:(n_states + 1)] <- next_state
      }
      
      
      # Convert to data frame
      result <- as.data.frame(result)
      
      names(result) <- c(
        "time",
        "N1",
        "N2"
      )
      
      return(result)
    }
    
    
    # -------------------------
    # Solve the differential equations
    # -------------------------
    
    output_data <- rk4(
      initial_state = initial_state,
      times = times,
      model = competition_model,
      parameters = parameters
    )
    
    
    # -------------------------
    # Convert to long format
    # -------------------------
    
    plot_data <- rbind(
      data.frame(
        Time = output_data$time,
        Population = output_data$N1,
        Species = "Species 1"
      ),
      data.frame(
        Time = output_data$time,
        Population = output_data$N2,
        Species = "Species 2"
      )
    )
    
    
    # -------------------------
    # Make the graph
    # -------------------------
    
    ggplot(
      plot_data,
      aes(
        x = Time,
        y = Population,
        color = Species
      )
    ) +
      geom_line(linewidth = 0.6) +
      
      scale_color_manual(
        values = c(
          "Species 1" = "blue",
          "Species 2" = "red"
        )
      ) +
      
      labs(
        x = "Time",
        y = "Population Size",
        color = ""
      ) +
      
      theme_minimal(base_size = 14) +
      
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom"
      ) +
      
      coord_cartesian(ylim = c(0, 400))
  })
}


shinyApp(ui = ui, server = server)
