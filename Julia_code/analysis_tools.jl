##########################################################################################
#############################      Hilfsfunktionen         ###############################
##########################################################################################

function construct_id_set_for_analysis(eq_type)
    if eq_type in ["telegraph", "telegraph_symm", "heat"]
        id_set = (eq_type, "telsin", -pi, pi)
    elseif eq_type in ["wave", "wave_symm", "transport"]
        id_set = (eq_type, "sin", 0, 1)
    end
    return id_set
end

function take_care_of_cut_cells(N, cut_cells)
    if 3*length(cut_cells) > N - 2 + length(cut_cells)
        i = 0
        while 3*length(cut_cells) > N - 2 + length(cut_cells)
            i = i+1
            cut_cells = cut_cells[1:end-1]
        end
        println("Warning: The last $(i) cut cells were not implemented, due to insufiicient cells!")
    end
    return cut_cells
end

function get_n_comp(eq_type)
    eq_types_1D = ["transport", "heat"]
    eq_types_2D = ["telegraph", "telegraph_symm", "wave", "wave_symm"]

    if eq_type in eq_types_1D
        n_comp = 1
    elseif eq_type in eq_types_2D
        n_comp = 2
    end
    return n_comp
end

#############################################################################################
###############################        Konvergenztest        ################################
#############################################################################################
function convergence_test(;
    T = Float64,
    eq_type = "telegraph", epsilon = 0.5, a = 1,
    exprange = [3,8], Tmax = 3.0, deg = 1, CFL_prefac = 0.5, basis_type = GaussLegendre,
    fluxtype = "altlr", J1_type = "upwind_diss_symm", J1_type_A = "upwind_diss_symm",
    cut_cells = [], do_stabilize = true, fix_eta = true, c = 0.4, include_b_h = true,
    TMM = SSPRK3
                        )
    if !(eq_type == "heat") deriv_fac = 0 else deriv_fac = 1 end
    id_set = construct_id_set_for_analysis(eq_type)

    n_comp = get_n_comp(eq_type)
    steps = 2 .^range(exprange[1],exprange[2])
    errors = zeros(T, length(steps)*n_comp)
    solution_output = Matrix{T}
    for (istep, step) in enumerate(steps)
        if !(eq_type == "heat")
            CFL = CFL_prefac*epsilon/(2*deg+1) # one factor of 1/N gets included in setup_problem_eq -> this is dt~epsilon/dx
        else
            # refined deg == 2 case to reduce the CFL that strong just in the required cases. Experimentally obtained.
            if deg == 2
                if basis_type == GaussLegendre
                    CFL = CFL_prefac/(2*deg+1)/(step) * 0.25
                elseif basis_type == LobattoLegendre
                    CFL = CFL_prefac/(2*deg+1)/(step) * 0.5
                end
            else
                CFL = CFL_prefac/(2*deg+1)/(step)  # one factor of 1/N gets included in setup_problem_eq -> this is dt~1/dx^2
            end
        end
        problem = setup_problem_eq(id_set[1], id_set[2], id_set[3], id_set[4], step, Tmax = Tmax, a = a, CFL = CFL, bcs = "periodic", epsilon = epsilon);
        cut_cells = take_care_of_cut_cells(step, cut_cells)
        for (i_cc, cut_cell) in enumerate(cut_cells)
            problem = include_cut_cell(problem, cut_cell, 3*i_cc)
        end

        RHS_mat, problem, SBP_storage = DGsemidiscretization_DoD_telegraph(problem, deg, basis_type, do_stabilize = do_stabilize, fix_eta = fix_eta,
                                                                    c = c, fluxtype = fluxtype, include_b_h = include_b_h, ext_test_func = true,
                                                                    J1_type = J1_type, J1_type_A = J1_type_A);

        #### Time stepping
        if typeof(TMM) == String
            Tableau = get_RK_tableau(TMM)
            sol, u_exact = ImEx(problem, RHS_mat, Tableau, only_explicit = false)
        else
            sol, u_exact = TMM(problem, RHS_mat)
        end
        #### Error calculation
        x_d = problem["x_d"]
        v = problem["v"]
        deg = problem["deg"]
        nodes = problem["nodes"]
        comprange = problem["comprange"]
        for comp in 1:n_comp
            errors[istep + (comp-1)*length(steps)] = 0
            for i in 1:problem["cellnumber"]
                basis1 = nodes(deg)
                for j in 1:deg+1
                    # Inserting these nodes into the structure basis1, to adapt the field basis1.interpolation_matrix to the used nodes.
                    basis1.nodes[j] = x_d[(deg+1)*(i-1) + j] 
                end
                value_range = (deg+1)*(i-1)+1 + (comp-1)*comprange:(deg+1)*i + (comp-1)*comprange
                errors[istep + (comp-1)*length(steps)] += integrate((sol[value_range, end] - u_exact[value_range, end]).^2, basis1.weights.*(v[i+1]-v[i])/2)
            end
            errors[istep + (comp-1)*length(steps)] = sqrt(errors[istep + (comp-1)*length(steps)])
        end
        if istep == length(steps)
            solution_output = zeros(T, length(sol[:, end]), 4)
            if eq_type in ["transport", "heat"]
                solution_output[:, 1] = x_d
            else 
                solution_output[:, 1] = vcat(x_d, x_d)
                steps = vcat(steps, steps)
            end
            solution_output[:, 2] = sol[:, end]
            solution_output[:, 3] = u_exact[:, end]
        end
    end
    return steps, errors, solution_output
end



#############################################################################################
################################        Operatornorm        #################################
#############################################################################################


function calc_op_norm(setup, RHS_mat)
    T = setup["T"]
    x_d = setup["x_d"]
    M_global = setup["M"]
    dim_M = length(M_global[1,:])
    sqrtM = sqrt.(M_global)
    invsqrtM = zeros(T, dim_M, dim_M)
    for ism in 1:dim_M
        for jsm in 1:dim_M
            if M_global[ism,jsm] != 0
                invsqrtM[ism,jsm] = 1 ./(sqrt.(M_global[ism,jsm]))
            end
        end
    end
    return opnorm(sqrtM * RHS_mat * invsqrtM)
end


function op_norm_vs_cut_fac(;
    T = Float64,
    eq_type = "telegraph", epsilon = 0.5, a = 1,
    N = 2^4, Tmax = 3.0, deg = 1, CFL_prefac = 0.5, basis_type = GaussLegendre,
    fluxtype = "altlr", J1_type = "upwind_diss_symm", J1_type_A = "upwind_diss_symm",
    cut_cell_refinement = 100, fix_eta = true, c = 0.4, include_b_h = true,
    include_cut_cells = true, do_stabilize = true,
                        )
    id_set = construct_id_set_for_analysis(eq_type)
    n_comp = get_n_comp(eq_type)
    cut_cell_factors = range(0.01, 0.49, length = cut_cell_refinement)
    opnorms = zeros(T, length(cut_cell_factors))
    CFL = 1 # cfl irrelevant here, but needs to be defined
    for (ialpha,alpha) in enumerate(cut_cell_factors)
        problem = setup_problem_eq(id_set[1], id_set[2], id_set[3], id_set[4], N, Tmax = Tmax, a = a, CFL = CFL, bcs = "periodic", epsilon = epsilon);
        if include_cut_cells == true
            for i in 1:6
                problem = include_cut_cell(problem, alpha, 3*i)
            end
        end
        RHS_mat, problem, SBP_storage = DGsemidiscretization_DoD_telegraph(problem, deg, basis_type, T = T, do_stabilize = do_stabilize, fix_eta = fix_eta,
                                                                    c = c, fluxtype = fluxtype, include_b_h = include_b_h, ext_test_func = true,
                                                                    J1_type = J1_type, J1_type_A = J1_type_A);
        opnorms[ialpha] = calc_op_norm(problem, RHS_mat)
    end
    return cut_cell_factors, opnorms
end

function minimize_operatornorm(;
                            T = Float64,
                            eq_type = "telegraph", epsilon = 0.5, a = 1,
                            N = 2^4, deg = 1, basis_type = GaussLegendre,
                            fluxtype = "altlr", J1_type = "upwind_diss_symm", J1_type_A = "upwind_diss_symm",
                            cut_cell_refinement = 51, fix_eta = true, include_b_h = true,
                            include_cut_cells = true, do_stabilize = true,
                            verbose = true
                                                )

    # CFL has no impact, but needs to be defined here
    CFL = 1.0
    
    # initializing
    best_norm = 0
    best_c = 0
    new_cs = zeros(T, 25)
    ######################
    #first step
    ######################
    cs = range(0.01, 1, length = 51) 
    norm_matrix_vs_c = zeros(T, cut_cell_refinement, length(cs))
    # calculate all operatornorms
    if verbose == true
        println("Time for first 100 evaluations:")
        @time for (ic, c) in enumerate(cs)
            norm_matrix_vs_c[:, ic] = op_norm_vs_cut_fac(;  T = T, eq_type = eq_type, epsilon = epsilon, a = a,
                                                            N = N, deg = deg, CFL_prefac = CFL_prefac, basis_type = basis_type,
                                                            fluxtype = fluxtype, J1_type = J1_type, J1_type_A = J1_type_A,
                                                            cut_cell_refinement = cut_cell_refinement, fix_eta = fix_eta, c = c, include_b_h = include_b_h,
                                                            include_cut_cells = true, do_stabilize = true,
                                                            )[2]
        end
    else
        for (ic, c) in enumerate(cs)
            norm_matrix_vs_c[:, ic] = op_norm_vs_cut_fac(;  T = T, eq_type = eq_type, epsilon = epsilon, a = a,
                                                            N = N, deg = deg, CFL_prefac = CFL_prefac, basis_type = basis_type,
                                                            fluxtype = fluxtype, J1_type = J1_type, J1_type_A = J1_type_A,
                                                            cut_cell_refinement = cut_cell_refinement, fix_eta = fix_eta, c = c, include_b_h = include_b_h,
                                                            include_cut_cells = true, do_stabilize = true,
                                                            )[2]
        end
    end
    # calculate c for all alphas
    (best_norm, ind_best_c) = findmin(maximum(norm_matrix_vs_c, dims = 1))
    best_c = cs[ind_best_c[2]]
    if verbose == true
        println("First approximation result: norm(L) = $(best_norm), c = $(best_c)")
    end
    ######################
    # iterate to more decimals
    ######################
    for exp in [-2, -3, -4]
        new_cs = range(best_c - 10.0^exp, best_c + 10.0^exp, length = 25)
        norm_matrix_vs_c = zeros(T, cut_cell_refinement, length(new_cs))
        if verbose == true
            println("Time for step $(-exp):")
            @time for (ic, c) in enumerate(new_cs)
                norm_matrix_vs_c[:, ic] = op_norm_vs_cut_fac(;  T = T, eq_type = eq_type, epsilon = epsilon, a = a,
                                                                N = N, deg = deg, CFL_prefac = CFL_prefac, basis_type = basis_type,
                                                                fluxtype = fluxtype, J1_type = J1_type, J1_type_A = J1_type_A,
                                                                cut_cell_refinement = cut_cell_refinement, fix_eta = fix_eta, c = c, include_b_h = include_b_h,
                                                                include_cut_cells = true, do_stabilize = true,
                                                                )[2]
            end
        else
            for (ic, c) in enumerate(new_cs)
            norm_matrix_vs_c[:, ic] = op_norm_vs_cut_fac(;  T = T, eq_type = eq_type, epsilon = epsilon, a = a,
                                                            N = N, deg = deg, CFL_prefac = CFL_prefac, basis_type = basis_type,
                                                            fluxtype = fluxtype, J1_type = J1_type, J1_type_A = J1_type_A,
                                                            cut_cell_refinement = cut_cell_refinement, fix_eta = fix_eta, c = c, include_b_h = include_b_h,
                                                            include_cut_cells = true, do_stabilize = true,
                                                            )[2]
            end
        end
       (best_norm, ind_best_c) = findmin(maximum(norm_matrix_vs_c, dims = 1))
        best_c = new_cs[ind_best_c[2]]
        if verbose == true
            println("Approximation result step $(-exp): norm(L) = $(best_norm), c = $(best_c)")
        end
    end
    return best_norm, best_c
end

#############################################################################################
############################        Asymptotic analysis        ##############################
#############################################################################################

function asymptotic_error(;
    T = Float64,
    epsilon_refinement = 30, a = 1,
    N = 2^5, Tmax = 3.0, deg = 1, CFL_prefac = 0.5, basis_type = GaussLegendre,
    fluxtype = "altlr", J1_type = "upwind_diss_symm", J1_type_A = "upwind_diss_symm",
    cut_cells = [], fix_eta = true, c = 0.4, include_b_h = true,
    TMM = SSPRK3,
    CFL_type = "eps/dx^2",
                        )
    epsilons = (1/2) .^(1:epsilon_refinement)
    errors = zeros(T, epsilon_refinement)
    sol_telegraph = Vector{Matrix{T}}()
    sol_heat = Vector{Matrix{T}}()
    for (ieps, eps) in enumerate(epsilons)
        if CFL_type == "eps/dx^2"
            CFL = CFL_prefac*eps/(2*deg+1)/N
        elseif CFL_type == "eps/dx"
            CFL = CFL_prefac*eps/(2*deg+1)
        elseif CFL_type == "1/dx^2"
            CFL = CFL_prefac/(2*deg+1)/N
        end
        # calulation of telegraph equation
        problem = setup_problem_eq("telegraph", "telegraph_well_prep" , -pi, pi, N, Tmax = Tmax, a = a, CFL = CFL, bcs = "periodic", epsilon = eps, T = T,);
        cut_cells = take_care_of_cut_cells(N, cut_cells)
        for (i_cc, cut_cell) in enumerate(cut_cells)
            problem = include_cut_cell(problem, cut_cell, 3*i_cc)
        end
""
        RHS_mat, problem, SBP_storage = DGsemidiscretization_DoD_telegraph(problem, deg, basis_type, do_stabilize = true, fix_eta = fix_eta,
                                                                    c = c, fluxtype = fluxtype, include_b_h = include_b_h, ext_test_func = true,
                                                                    J1_type = J1_type, J1_type_A = J1_type_A);

        if typeof(TMM) == String
            Tableau = get_RK_tableau(TMM)
            push!(sol_telegraph, ImEx(problem, RHS_mat, Tableau, only_explicit = false, optimized_ARS = true)[1])
        else
            push!(sol_telegraph, TMM(problem, RHS_mat)[1])
        end

        # calculation of heat equation
        problem = setup_problem_eq("heat", "telsin" , -pi, pi, N, Tmax = Tmax, a = a, CFL = CFL, bcs = "periodic", epsilon = eps, T = T);
        cut_cells = take_care_of_cut_cells(N, cut_cells)
        for (i_cc, cut_cell) in enumerate(cut_cells)
            problem = include_cut_cell(problem, cut_cell, 3*i_cc)
        end

        RHS_mat, problem, SBP_storage = DGsemidiscretization_DoD_telegraph(problem, deg, basis_type, do_stabilize = true, fix_eta = fix_eta,
                                                                    c = c, fluxtype = fluxtype, include_b_h = include_b_h, ext_test_func = true,
                                                                    J1_type = J1_type, J1_type_A = J1_type_A);

        #### Time stepping
        if typeof(TMM) == String
            Tableau = get_RK_tableau(TMM)
            push!(sol_heat, ImEx(problem, RHS_mat, Tableau, only_explicit = false)[1])
        else
            push!(sol_heat, TMM(problem, RHS_mat)[1])
        end

        #### Error calculation
        x_d = problem["x_d"]
        v = problem["v"]
        deg = problem["deg"]
        nodes = problem["nodes"]
        errors[ieps] = 0
        # calculate error just for the first component, as this is the one to compare with the heat equation
        for i in 1:problem["cellnumber"]
            basis1 = nodes(deg)
            for j in 1:deg+1
                basis1.nodes[j] = x_d[(deg+1)*(i-1) + j] 
            end
            value_range = (deg+1)*(i-1)+1:(deg+1)*i
            errors[ieps] += integrate((sol_telegraph[ieps][value_range, end] - sol_heat[ieps][value_range, end]).^2, basis1.weights.*(v[i+1]-v[i])/2)
        end
        errors[ieps] = sqrt(errors[ieps])
    end
    return epsilons, errors
end


#############################################################################################
###################################        Error        #####################################
#############################################################################################


function calc_error(problem, u_sol, u_exact; t_dest = length(u_sol[1, :]))
        x_d = problem["x_d"]
        v = problem["v"]
        deg = problem["deg"]
        nodes = problem["nodes"]
        error = 0
        for i in 1:problem["cellnumber"]
            basis1 = nodes(deg)
            for j in 1:deg+1
                basis1.nodes[j] = x_d[(deg+1)*(i-1) + j] 
            end
            value_range = (deg+1)*(i-1)+1:(deg+1)*i
            error+= integrate((u_sol[value_range, t_dest] - u_exact[value_range, t_dest]).^2, basis1.weights.*(v[i+1]-v[i])/2)
        end
        error = sqrt(error)
    return error
end
