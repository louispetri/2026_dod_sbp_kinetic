using LinearAlgebra, PolynomialBases, DelimitedFiles
include("../functions.jl")
include("../functionsRK.jl")
include("../flux_operators.jl")
include("../analysis_tools.jl")


cut_cells = [0.49, 0.3, 0.1, 10^-2, 10^-7]
epsilons = [0.5, 0.1, 0.01]
a = 1
T = Float64
exprange = [4,8]
Tmax = 1.0
CFL_prefac = [0.5, 0.3, 0.15]
flux_tuple = [("altlr", "upwind_diss_symm", "upwind_diss_symm"),
                ("central", "central_symm", "upwind_diss_symm")]
fix_eta = true
c = [1.0, 0.55, 0.45]
include_b_h = true
TMMs = ["ARS3"]
degrees = 0:2
for eq_type in ["telegraph", "heat"]
    if eq_type in ["heat"]
        stepdim = exprange[2]-exprange[1]+1
    else
        stepdim = (exprange[2]-exprange[1]+1)*2
    end
    for TMM in TMMs
        for basis_type in [GaussLegendre]
            for (fluxtype, J1_type, J1_type_A) in flux_tuple
                if typeof(TMM) == String
                    TMM_str = TMM
                else
                    TMM_str = "$(TMM)"
                end
                mkpath(joinpath(@__DIR__,"./data/" * eq_type * "/"  * TMM_str * "/" * "/$(basis_type)/" * J1_type))
                for epsilon in epsilons
                    output_matrix = zeros(T, stepdim, length(degrees)+1) # store the error values for every degree
                    for (ideg, deg) in enumerate(degrees)
                        steps, errors, solution_output = convergence_test(;
                                                T=T,
                                                eq_type = eq_type, epsilon = epsilon, a = a,
                                                exprange = exprange, Tmax = Tmax, deg = deg, CFL_prefac = CFL_prefac[ideg], basis_type = basis_type,
                                                fluxtype = fluxtype, J1_type = J1_type, J1_type_A = J1_type_A,
                                                cut_cells = cut_cells, fix_eta = fix_eta, c = c[ideg], include_b_h = include_b_h,
                                                TMM = TMM)
                        output_matrix[:, ideg + 1] = errors
                        if ideg == 1
                            output_matrix[:, 1] = steps
                        elseif ideg != 1 && steps != output_matrix[:, 1]
                            throw(ArgumentError("ACHTUNG! ZELLENANZAHL UNTERSCHEIDET SICH BEI UNTERSCHEIDLICHEM GRAD"))
                        end
                        writedlm(joinpath(@__DIR__,"./data/" * eq_type * "/" * TMM_str * "/" * "/$(basis_type)/" * J1_type * "/simulation_data_epsilon=$(epsilon)_deg=$(deg).txt"), solution_output)
                    end
                    writedlm(joinpath(@__DIR__,"./data/" * eq_type * "/"  * TMM_str * "/" * "/$(basis_type)/" * J1_type * "/conv_error_epsilon=$(epsilon).txt"), output_matrix)
                end
            end
        end
    end
end