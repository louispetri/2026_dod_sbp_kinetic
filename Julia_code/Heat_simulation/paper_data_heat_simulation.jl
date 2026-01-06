using LinearAlgebra, CairoMakie, PolynomialBases, DelimitedFiles, LaTeXStrings, Colors
include("../functions.jl")
include("../functionsRK.jl")
include("../flux_operators.jl")
include("../analysis_tools.jl")


id_set = ("heat", "telsin" , -pi, pi)
# irrelevant, but has to be defined
epsilon = 0.5
a = 1.0
T = Float64


N = 2^5
Tmax = 5.0

c_upwind_diss_symm = [1.0, 0.511, 0.337, 0.247, 0.183, 0.148]
c_central_symm = [1.0, 0.563, 0.365, 0.264, 0.191, 0.150]


for do_stabilize in [false,true]
    for include_cut_cells in [false,true]
        if include_cut_cells == true
            if do_stabilize == false
                titlestring = "Unstabilized"
                ysc = identity
                settingstr = "unstabil"
            else
                titlestring = "DoD (stabilized)"
                ysc = identity
                settingstr = "stabil"
            end
        else
            titlestring = "Background"
            ysc = identity
            settingstr = "backgr"
        end
        for (fluxtype, J1_type, cs) in [("altlr", "upwind_diss_symm", c_upwind_diss_symm)]
            for (ideg, deg) in enumerate([0, 1, 2])
                CFL = 0.1/(2*deg+1)/N
                J1_type_A = "upwind_diss_symm"
                c = cs[deg+1]

                problem = setup_problem_eq(id_set[1], id_set[2], id_set[3], id_set[4], N, Tmax = Tmax, a = a, CFL = CFL, bcs = "periodic", epsilon = epsilon);
                if include_cut_cells == true
                    problem = include_cut_cell(problem, 0.001, 3);
                    problem = include_cut_cell(problem, 0.4, 7);
                    problem = include_cut_cell(problem, 0.0000001, 10);
                    problem = include_cut_cell(problem, 0.25, 14);
                    problem = include_cut_cell(problem, 0.49, 20);
                    problem = include_cut_cell(problem, 0.1, 26);
                end
                RHS_mat, problem, SBP_storage = DGsemidiscretization_DoD_telegraph(problem, deg, GaussLegendre,  do_stabilize = do_stabilize, fix_eta = true,
                                                                                c = cs[deg+1], fluxtype = fluxtype, include_b_h = true, ext_test_func = true,
                                                                                J1_type = J1_type, J1_type_A = J1_type_A);
                                                                                x_d = problem["x_d"]
                v = problem["v"]
                metadata = zeros(length(v) + 2)
                deg = problem["deg"]
                cellnumber = problem["cellnumber"]
                metadata[1] = deg
                metadata[2] = cellnumber
                metadata[3:end] = v
                if problem["nodes"] != GaussLegendre
                    throw(ArgumentError("Basistyp nicht implementiert"))
                end
                u0 = problem["u0"]
                Nx_plot = determine_Nx_plot(problem)
                Tableau = get_RK_tableau("ImplMidpoint")
                sol, u_exact = ImEx(problem, RHS_mat, Tableau, only_explicit = false)
                solution_output = zeros(T, length(sol[:, end]), 3)
                solution_output[:, 1] = x_d
                solution_output[:, 2] = sol[:, end]
                solution_output[:, 3] = u_exact[:, end]
                mkpath(joinpath(@__DIR__,"./data/N=$(N)/"))
                writedlm(joinpath(@__DIR__,"./data/N=$(N)/$(settingstr)_deg$(deg).txt"), solution_output)
                writedlm(joinpath(@__DIR__,"./data/N=$(N)/$(settingstr)_metadata_deg$(deg).txt"), metadata)
            end
        end
    end
end



