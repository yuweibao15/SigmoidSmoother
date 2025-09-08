# plot_SmoothSkygrid_G_H.R
# Interactive Shiny app for G(t,s) and H(t) (step function version)
# Run:
#   library(shiny)
#   shiny::runApp("/Users/ybao2/Github_Repo/SmoothSkygrid/R/plot_SmoothSkygrid_G_H.R")

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(gridExtra)
})

# -------------------------
# Math helpers
# -------------------------

# G(t, s) = 1/theta1 + sum_k (1/theta_{k+1} - 1/theta_k) / (1 + exp(-s (t - x_k)))
G_fun <- function(t, thetas, x, s) {
  stopifnot(length(thetas) == length(x) + 1L, s > 0)
  inv_thetas <- 1 / thetas
  if (length(x) == 0L) return(rep(inv_thetas[1], length(t)))
  delta <- inv_thetas[-1] - inv_thetas[-length(inv_thetas)]
  # plogis(u) = 1/(1+exp(-u)) for numerical stability
  L <- plogis(outer(t, x, function(tt, xx) s * (tt - xx)))
  inv_thetas[1] + as.vector(L %*% delta)
}

# H(t) = 1/theta1 + sum_k (1/theta_{k+1} - 1/theta_k) * I{ t >= x_k }
H_fun <- function(t, thetas, x) {
  stopifnot(length(thetas) == length(x) + 1L)
  inv_thetas <- 1 / thetas
  if (length(x) == 0L) return(rep(inv_thetas[1], length(t)))
  delta <- inv_thetas[-1] - inv_thetas[-length(inv_thetas)]
  kstar <- findInterval(t, x)       # number of knots <= t (0..M)
  cumdelta <- c(0, cumsum(delta))   # length M+1, index by kstar+1
  inv_thetas[1] + cumdelta[kstar + 1]
}

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
  shiny::withMathJax(),
  titlePanel(HTML("Interactive \\(G(t,s)\\) and \\(H(t)\\)")),
  helpText(
    "\\(G(t, s) = \\frac{1}{\\theta_1} + \\sum_{k=1}^M 
      \\frac{\\tfrac{1}{\\theta_{k+1}} - \\tfrac{1}{\\theta_k}}{1 + e^{-s (t - x_k)}}\\)",
    br(),
    "\\(H(t) = \\frac{1}{\\theta_1} + \\sum_{k=1}^M
      (\\tfrac{1}{\\theta_{k+1}} - \\tfrac{1}{\\theta_k})\\,\\mathbf{1}\\{t\\ge x_k\\}\\)"
  ),
  
  sidebarLayout(
    sidebarPanel(
      # Core controls
      sliderInput("M", "Number of transitions (M):", min = 0, max = 10, value = 4, step = 1),
      sliderInput("s", "Steepness (s):", min = 0, max = 1000, value = 10, step = 0.1),
      sliderInput("tminmax", "Time range [t_min, t_max]:", min = 0, max = 10, value = c(0, 3), step = 0.1),
      
      checkboxInput("logH", "Plot H(t) on log10 scale", value = FALSE),
      tags$hr(),
      
      actionButton("autospace_x", "Auto-space knots x_k over time range"),
      actionButton("smooth_plateaus", "Auto-smooth θ plateaus"),
      tags$hr(),
      
      h4("θ (M+1 values) and x (M values)"),
      helpText("θ must be > 0. x must be strictly increasing."),
      uiOutput("theta_inputs"),
      uiOutput("x_inputs"),
      width = 4
    ),
    
    mainPanel(
      plotOutput("plots", height = "520px"),
      tags$hr(),
      verbatimTextOutput("diag"),
      width = 8
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output, session) {
  
  # Defaults that adapt to M
  default_thetas <- reactive({
    M <- input$M
    base <- 8000; span <- 3000
    val <- base + span * sin(seq(0, pi, length.out = M + 1))
    round(pmax(val, 1e-6), 3)
  })
  
  default_x <- reactive({
    M <- input$M
    if (M == 0) return(numeric(0))
    rng <- input$tminmax
    xs <- seq(rng[1], rng[2], length.out = M + 2)
    xs <- xs[-c(1, length(xs))]
    round(xs, 3)
  })
  
  # Buttons to auto-fill
  observeEvent(input$autospace_x, {
    M  <- input$M
    xs <- default_x()
    for (k in seq_len(M)) {
      updateNumericInput(session, paste0("x_", k), value = xs[k])
    }
  })
  
  observeEvent(input$smooth_plateaus, {
    th <- default_thetas()
    for (j in seq_along(th)) {
      updateNumericInput(session, paste0("theta_", j), value = th[j])
    }
  })
  
  # Dynamic UI for θ and x
  output$theta_inputs <- renderUI({
    M  <- input$M
    th <- default_thetas()
    tagList(lapply(seq_len(M + 1), function(j) {
      numericInput(
        inputId = paste0("theta_", j),
        label   = sprintf("theta[%d]", j),
        value   = th[j],
        min     = 1e-9, step = 1
      )
    }))
  })
  
  output$x_inputs <- renderUI({
    M <- input$M
    if (M == 0) return(NULL)
    xs <- default_x()
    tagList(lapply(seq_len(M), function(k) {
      numericInput(
        inputId = paste0("x_", k),
        label   = sprintf("x[%d]", k),
        value   = xs[k],
        step    = 0.001
      )
    }))
  })
  
  # Collect parameter vectors
  thetas_reac <- reactive({
    M <- input$M
    vals <- sapply(seq_len(M + 1), function(j) input[[paste0("theta_", j)]], USE.NAMES = FALSE)
    validate(need(all(is.finite(vals)) && all(vals > 0), "All θ must be positive and finite."))
    vals
  })
  
  x_reac <- reactive({
    M <- input$M
    if (M == 0) return(numeric(0))
    vals <- sapply(seq_len(M), function(k) input[[paste0("x_", k)]], USE.NAMES = FALSE)
    validate(
      need(all(is.finite(vals)), "All x must be finite."),
      need(all(diff(vals) > 0), "x must be strictly increasing (x1 < x2 < ... < xM).")
    )
    vals
  })
  
  # Time grid (fixed resolution; keep UI simple)
  t_grid <- reactive({
    rng <- input$tminmax
    seq(rng[1], rng[2], length.out = 600)
  })
  
  # Compute curves
  curves <- reactive({
    t  <- t_grid()
    th <- thetas_reac()
    xs <- x_reac()
    s  <- input$s
    
    Gvals <- G_fun(t, th, xs, s)
    Hvals <- H_fun(t, th, xs)
    list(t = t, G = Gvals, H = Hvals, x = xs)
  })
  
  # Plots
  output$plots <- renderPlot({
    crv <- curves()
    df1 <- data.frame(t = crv$t, val = crv$G)
    df2 <- data.frame(t = crv$t, val = crv$H)
    
    p1 <- ggplot(df1, aes(t, val)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = crv$x, linetype = "dashed", color = "grey70") +
      labs(title = expression(G(t, s)), x = "t", y = "G(t,s)") +
      theme_minimal(base_size = 13)
    
    p2 <- ggplot(df2, aes(t, val)) +
      geom_step(linewidth = 1) +
      geom_vline(xintercept = crv$x, linetype = "dashed", color = "grey70") +
      labs(
        title = expression(H(t)), x = "t",
        y = if (isTRUE(input$logH)) "log10 H(t)" else "H(t)"
      ) +
      theme_minimal(base_size = 13)
    
    if (isTRUE(input$logH)) {
      p2 <- p2 + scale_y_continuous(trans = "log10")
    }
    
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  }, res = 120)
  
  # Diagnostics
  output$diag <- renderPrint({
    th <- thetas_reac()
    xs <- x_reac()
    list(
      M = input$M,
      thetas = signif(th, 6),
      x = if (length(xs)) signif(xs, 6) else numeric(0),
      s = input$s,
      t_range = input$tminmax
    )
  })
}

shinyApp(ui, server)