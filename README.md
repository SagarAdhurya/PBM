---

# Process-based Modelling (PBM)

This repository contains R code for simulating prey-predator dynamics using a time-dependent Lotka-Volterra model. The model incorporates observed data for prey growth rate (`r`) through interpolation and performs sensitivity analysis and parameter calibration using Advanced Parameter Calibration (APC).

## Features

- **Time-dependent Lotka-Volterra Model:** Utilises observed data for prey growth rate (`r`) using interpolation.
- **Sensitivity Analysis:** Evaluates the impact of parameters (`a`, `f`, and `q`) on state variables (Prey and Predator) using mean, mean square root, and mean absolute sensitivity indices.
- **Advanced Parameter Calibration (APC):** Calibrates model parameters to observed data, minimising the cost function.
- **Visualisation:** Generates dynamic plots for simulated prey and predator populations compared with observed data.

## Requirements

- `R` (version 4.0 or higher)
- `deSolve`: For solving differential equations
- `ecolMod`: For parameter optimisation

Install necessary packages in R:

```r
install.packages(c('deSolve', 'ecolMod'))
```

## Usage

1. Clone this repository:

```sh
git clone https://github.com/SagarAdhurya/PBM
```

2. Open the R script and run it in RStudio or any compatible R environment.
3. The simulation results will be saved as `sim_result.csv`.
4. Plots will be exported as PNG files in the working directory.

## Structure

- `prey_predator_model.R`: Main R script containing model declaration, sensitivity analysis, and parameter calibration.
- `sim_result.csv`: Simulation results (automatically generated).
- `prey_predator_plots.png`: Plot showing Prey-Predator dynamics.
- `calibrated_preypredator_plots_with_observed_interpolated.png`: Plots comparing calibrated and uncalibrated models with observed data.

## Methodology

1. **Model Declaration:** Defines a time-dependent Lotka-Volterra model using observed prey growth rate (`r`) interpolated over time.
2. **Sensitivity Analysis:** Perturbs model parameters (`a`, `f`, and `q`) by ±10% to evaluate their effects on state variables.
3. **Automatic Parameter Calibration (APC):** Utilises a cost function to minimise discrepancies between simulated and observed data, adjusting parameters accordingly.
4. **Visualisation:** Generates dynamic plots to compare simulated results with observed values.

## Results

The model effectively simulates prey-predator dynamics with time-dependent growth rate (`r`). The sensitivity analysis reveals parameter influences, and APC enhances model accuracy by calibrating parameters.

## License

Project Title is released under the MIT License. See the **[LICENSE](https://github.com/SagarAdhurya/PBM/blob/main/LICENSE)** file for details.

## Authors

Project was created and maintained by **[Dr. Sagar Adhurya](https://github.com/SagarAdhurya)**.

## **Contact**

If you have any questions or comments about Project Title, please contact **[Dr. Sagar Adhurya](mailto:sagaradhurya@gmail.com)**.

---
