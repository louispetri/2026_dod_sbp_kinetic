using LinearAlgebra, PolynomialBases, DelimitedFiles
include("../functions.jl")
include("../functionsRK.jl")
include("../flux_operators.jl")
include("../analysis_tools.jl")

##################################
###### set up equation ###########
##################################
id_set = ("heat", "telsin" , -pi, pi)

# equation parameters
epsilon = 10^-2
a = 1.0

##################################
###### Numerical parameters ######
##################################
N = 2^7
Tmax = 5.0

c_upwind_diss_symm = [1.0, 0.511, 0.337, 0.247, 0.183, 0.148]
c_central_symm = [1.0, 0.563, 0.365, 0.264, 0.191, 0.150]

kappas_M = zeros(6,6)

for (fluxtype, J1_type, cs) in [("altlr", "upwind_diss_symm", c_upwind_diss_symm), ("central", "central_symm", c_central_symm)]
    for deg in [0, 1, 2, 3, 4, 5]
        for include_cut_cells in [true, false]
            for do_stabilize in [true, false]
                CFL = 0.05/(2*deg+1)/N
                J1_type_A = "upwind_diss_symm"
                c = cs[deg+1]

                ##################################
                ###### set up problem ############
                ##################################
                problem = setup_problem_eq(id_set[1], id_set[2], id_set[3], id_set[4], N, Tmax = Tmax, a = a, CFL = CFL, bcs = "periodic", epsilon = epsilon);
                if include_cut_cells == true
                    problem = include_cut_cell(problem, 0.001, 3);
                    problem = include_cut_cell(problem, 0.4, 7);
                    problem = include_cut_cell(problem, 0.0000001, 10);
                    problem = include_cut_cell(problem, 0.25, 14);
                    problem = include_cut_cell(problem, 0.49, 20);
                    problem = include_cut_cell(problem, 0.1, 26);
                end
                #
                ##################################
                ###### discretize in space #######
                ##################################
                RHS_mat, problem, SBP_storage = DGsemidiscretization_DoD_telegraph(problem, deg, GaussLegendre,  do_stabilize = do_stabilize, fix_eta = true,
                                                                                c = c, fluxtype = fluxtype, include_b_h = true, ext_test_func = true,
                                                                                J1_type = J1_type, J1_type_A = J1_type_A);
                dt=problem["t_d"][2]-problem["t_d"][1]
                invM_im=inv(I-dt*RHS_mat)
                if include_cut_cells == false
                    ind1 = 0
                elseif do_stabilize == false # cut cells ohne stabilisierung
                    ind1 = 2
                else # cut cells mit stabilisierung
                    ind1 = 4
                end
                if J1_type == "upwind_diss_symm"
                    ind2 = 1
                elseif J1_type == "central_symm"
                    ind2 = 2 
                end

                #condition in M-norm
                kappas_M[deg+1, ind1 + ind2] = calc_op_norm(problem, I-dt*RHS_mat)*calc_op_norm(problem, invM_im)
            end
        end
    end
end

writedlm(joinpath(@__DIR__,"./condition_numbers_Mnorm.txt"), round.(kappas_M, digits=4))
