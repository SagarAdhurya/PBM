####Model declaration####
# Load necessary library
if (require(deSolve) == F) {
  install.packages('deSolve',repos='http://cran.r-project.org', dependencies = T);
}

# Define the observed values for r as [time, r]
r_observed <- data.frame(
  time = c(0, 10, 17, 27, 36, 45, 61, 70, 84, 100),
  r = c(0.5, 1.2, 1.9, 2.1, 0.9, 0.6, 0.1, 0.2, 0.7, 0.6)
)

# Create a time-dependent function for r using interpolation
r_func <- approxfun(x = r_observed$time, y = r_observed$r, rule = 2)

# Define the Lotka-Volterra system with time-dependent r
prey_predator_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Interpolated r value at time t
    r_t <- r_func(t)  # Prey growth rate from observed values
    a_t <- a  # Constant predation rate
    f_t <- f  # Predator's efficiency to convert food into predator's offspring
    q_t <- q  # Constant predator death rate
    
    # Lotka-Volterra equations
    dPrey <- r_t * Prey - a_t * Prey * Predator
    dPredator <- f_t * a_t * Prey * Predator - q_t * Predator
    
    # Return rate of change
    list(c(dPrey, dPredator))
  })
}

# Initial state (Prey and Predator populations)
initial_state <- c(Prey = 40, Predator = 9)

# Time sequence for the simulation
time <- seq(0, 100, by = 1)  # Simulating for 100 time units

# Parameters for the constant rates (excluding r)
parameters <- list(
  a = 0.2,
  f = 0.1, 
  q = 0.3   
)

# Run the simulation with original parameters
results <- ode(
  y = initial_state,
  times = time,
  func = prey_predator_model,
  parms = parameters,
  method = 'rk4'
)

# Convert the results to a data frame for easier plotting
results_df <- as.data.frame(results)
write.csv(results_df, file = "sim_result.csv", row.names = F)

#Plot and export the graph as PNG
png("prey_predator_plots.png", width = 2000, height = 1500, units = "px", res = 300)
{
  # Plot prey and predator populations over time
  plot(results_df$time, results_df$Prey, type = "l", col = "blue", lwd = 2, ylab = "Population", xlab = "Time",
       main = "Prey-Predator Dynamics")
  lines(results_df$time, results_df$Predator, col = "red", lwd = 2)
  legend("top", legend = c("Prey", "Predator"), col = c("blue", "red"), lwd = 2)
  grid()
}

dev.off()  # Close the graphical device

####Sensitivity analysis####
# Sensitivity Analysis Function
calculate_sensitivity <- function(param_name, perturbation_factor) {
  sensitivity <- list()
  
  for (perturbation in c(-perturbation_factor, perturbation_factor)) {
    # Perturb the parameter
    perturbed_params <- parameters
    perturbed_params[[param_name]] <- parameters[[param_name]] * (1 + perturbation)
    
    # Simulate the model with perturbed parameter using rk4 method
    perturbed_results <- ode(
      y = initial_state,
      times = time,
      func = prey_predator_model,
      parms = perturbed_params,
      method = "rk4"
    )
    perturbed_results_df <- as.data.frame(perturbed_results)
    
    # Calculate sensitivity for each SV
    for (state_var in c("Prey", "Predator")) {
      original <- results_df[[state_var]]
      perturbed <- perturbed_results_df[[state_var]]
      
      # Sensitivity function: S_ij = (dX_j / X_j) / (dP_i / P_i)
      dX_j <- (perturbed - original) / original
      dP_i <- perturbation
      S_ij <- dX_j / dP_i
      
      # Initialisation and aggregation of sensitivities
      if (is.null(sensitivity[[state_var]])) {
        sensitivity[[state_var]] <- list(mean = 0, msqr = 0, mabs = 0)
      }
      sensitivity[[state_var]]$mean <- sensitivity[[state_var]]$mean + mean(S_ij) / 2
      sensitivity[[state_var]]$msqr <- sensitivity[[state_var]]$msqr + sqrt(mean(S_ij^2)) / 2
      sensitivity[[state_var]]$mabs <- sensitivity[[state_var]]$mabs + mean(abs(S_ij)) / 2
    }
  }
  return(sensitivity)
}

# Perform sensitivity analysis for each parameter
perturbation_factor <- 0.1  # 10% perturbation
sensitivity_results <- lapply(names(parameters), function(param_name) {
  calculate_sensitivity(param_name, perturbation_factor)
})

# Format and display results
names(sensitivity_results) <- names(parameters)
for (param_name in names(sensitivity_results)) {
  cat("\nSensitivity for parameter:", param_name, "\n")
  print(sensitivity_results[[param_name]])
}


# Calculate the overall sensitivity for each state variable (∆(X_j))
overall_sensitivity <- list()
for (state_var in c("Prey", "Predator")) {
  # Average of mean sensitivities across all parameters
  f_mean <- sapply(sensitivity_results, function(res) res[[state_var]]$mean)
  overall_sensitivity[[state_var]] <- mean(f_mean)  # ∆(X_j)
}

# Print overall sensitivity results
cat("\nOverall Sensitivity (∆X_j) for each state variable:\n")
print(overall_sensitivity)

####APC####
# Define observed data for Prey and Predator (at 10 specific time points)
observed_data <- data.frame(
  time = c(0, 10, 17, 27, 36, 45, 61, 70, 84, 100),
  Prey = c(40, 90, 38, 70, 35, 10, 25, 12, 18, 26),    # Observed data for Prey
  Predator = c(9, 13, 4, 10, 7, 4, 3, 3, 4, 6)    # Observed data for Predator
)

# Define the cost function to minimize
cost_function <- function(params) {
  # Update parameters
  parameters$a <- params[1]
  parameters$f <- params[2]
  parameters$q <- params[3]
  
  # Run the simulation with updated parameters
  cal <- ode(
    y = initial_state,
    times = time,
    func = prey_predator_model,
    parms = parameters,
    method = 'rk4'
  )
  
  # Convert results to a data frame
  cal_df <- as.data.frame(cal)
  
  # Interpolate observed data to match simulation time points
  observed_prey <- approx(x = observed_data$time, y = observed_data$Prey, xout = cal_df$time)$y
  observed_predator <- approx(x = observed_data$time, y = observed_data$Predator, xout = cal_df$time)$y
  
  # Calculate cost for Prey and Predator
  cost_prey <- sum((cal_df$Prey - observed_prey)^2 / (observed_prey^2), na.rm = TRUE)
  cost_predator <- sum((cal_df$Predator - observed_predator)^2 / (observed_predator^2), na.rm = TRUE)
  
  # Return the total cost (sum of both state variables' costs)
  return(cost_prey + cost_predator)
}

# Perform parameter optimisation
library(ecolMod)
optim <- pricefit(par = c(0.2, 0.1, 0.3), func = cost_function, minpar = c(0.2, 0.1, 0.001), maxpar = c(0.2, 0.1, 1), numiter = 10000)
optim$par

# Extract optimised parameters
optimized_parameters <- list(
  a = optim$par[1],
  f = optim$par[2],
  q = optim$par[3]
)

# Run the simulation with optimised parameters (calibrated)
optimized_results <- ode(
  y = initial_state,
  times = time,
  func = prey_predator_model,
  parms = optimized_parameters,
  method = 'rk4'
)

# Run the simulation with uncalibrated parameters (uncalibrated model)


# Convert results to data frames for plotting
optimized_results_df <- as.data.frame(optimized_results)
uncalibrated_results_df <- as.data.frame(results)

# Plot the results
png("calibrated_preypredator_plots_with_observed_interpolated.png", width = 3000, height = 3000, res = 300, units = 'px')
{
  par(mfrow = c(2, 1))
  
  # Determine xlim and ylim for Prey plot
  xlim_prey <- range(c(observed_data$time, uncalibrated_results_df$time, optimized_results_df$time))
  ylim_prey <- range(c(observed_data$Prey, uncalibrated_results_df$Prey, optimized_results_df$Prey))
  
  # Determine xlim and ylim for Predator plot
  xlim_predator <- range(c(observed_data$time, uncalibrated_results_df$time, optimized_results_df$time))
  ylim_predator <- range(c(observed_data$Predator, uncalibrated_results_df$Predator, optimized_results_df$Predator))
  
  # Plot Prey population dynamics
  plot(observed_data$time, observed_data$Prey, type = "p", col = "black", pch = 16, ylab = "Population", xlab = "Time",
       main = "Prey-Predator Dynamics (Optimized vs Unoptimized)", xlim = xlim_prey, ylim = ylim_prey)
  lines(uncalibrated_results_df$time, uncalibrated_results_df$Prey, col = "blue", lwd = 2)  # Uncalibrated simulation
  lines(optimized_results_df$time, optimized_results_df$Prey, col = "green", lwd = 2)  # Calibrated simulation
  legend("topright", legend = c("Observed Prey", "Uncalibrated Prey", "Calibrated Prey"), 
         col = c("black", "blue", "green"), pch = c(16, NA, NA), lwd = c(NA, 2, 2))
  
  # Plot Predator population dynamics
  plot(observed_data$time, observed_data$Predator, type = "p", col = "black", pch = 16, ylab = "Population", xlab = "Time",
       main = "Prey-Predator Dynamics (Optimized vs Unoptimized)", xlim = xlim_predator, ylim = ylim_predator)
  lines(uncalibrated_results_df$time, uncalibrated_results_df$Predator, col = "blue", lwd = 2)  # Uncalibrated simulation
  lines(optimized_results_df$time, optimized_results_df$Predator, col = "green", lwd = 2)  # Calibrated simulation
  legend("topright", legend = c("Observed Predator", "Uncalibrated Predator", "Calibrated Predator"), 
         col = c("black", "blue", "green"), pch = c(16, NA, NA), lwd = c(NA, 2, 2))
}
dev.off()
