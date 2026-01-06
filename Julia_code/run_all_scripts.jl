# Hint: The following 10 substeps can be run independently


# 1. Convergence reults (Figures 2 and 3)

include("./Convergence/get_convergence_data.jl") # Computation

include("./Convergence/plot_convergence.jl") # Plotting

# 2. Asymptotic analysis (Figure 4)

include("./Asymptotic_convergence/paper_asymptotic_convergence.jl") # Computation/Plotting

# 3. Condition numbers (Table 1)

include("./Heat_simulation/condition_number.jl") # Computation

# 4. Simulation results for the implicit midpoint RK scheme
#    applied to the heat equation (Figure 5)

include("./Heat_simulation/paper_data_heat_simulation.jl") # Computation

include("./Heat_simulation/paper_plot_heat_simulation.jl") # Plotting