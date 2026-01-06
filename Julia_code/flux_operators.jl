function stabilize_condition(setup, i)
    do_stabilize = setup["do_stabilize"]
    h = setup["h"]
    v = setup["v"]
    return do_stabilize && 2*(v[i+1]-v[i])/h<1
end

function determine_eta(setup, i)
    fix_eta = setup["fix_eta"]
    deg = setup["deg"]
    v = setup["v"]
    t_d = setup["t_d"]
    a = setup["a"]
    h = setup["h"]
    c = setup["c"]
    if fix_eta == false
        eta = 1-minimum([1/(2*deg+1)*(v[i+1]-v[i])/((t_d[2]-t_d[1])*a), 1.0])
    else
        eta = 1-minimum([(v[i+1]-v[i])/h*1/c, 1.0])
    end
    return eta
end

function get_Vol_interpol_matrices(setup, i)
    T = setup["T"]
    deg = setup["deg"]
    x_d = setup["x_d"]
    nodes = setup["nodes"]
    basisleft = nodes(deg, T)
    basisright = nodes(deg, T)
    for j in 1:deg+1
        basisleft.nodes[j] = x_d[(deg+1)*(i-2) + j] 
        basisright.nodes[j] = x_d[(deg+1)*(i) + j] 
    end
    Vol_interpol_left = interpolation_matrix(x_d[(deg+1)*(i-1) + 1:(deg+1)*(i)] , basisleft)
    Vol_interpol_right = interpolation_matrix(x_d[(deg+1)*(i-1) + 1:(deg+1)*(i)] , basisright)
    return Vol_interpol_left, Vol_interpol_right 
end

###############################################################################
###############################################################################

function get_Dminus(setup; for_bh = false)
    deg = setup["deg"]
    nodes = setup["nodes"]
    T = setup["T"]
    do_stabilize = setup["do_stabilize"]
    c = setup["c"]
    fix_eta = setup["fix_eta"]
    fluxtype = setup["fluxtype"]
    include_b_h = setup["include_b_h"]
    ext_test_func = setup["ext_test_func"]
    J1_type = setup["J1_type"]
    J1_type_A = setup["J1_type_A"]
    cellnumber = setup["cellnumber"]
    v = setup["v"]
    a = setup["a"]
    u0_eval = setup["u0_eval"]
    t_d = setup["t_d"]
    h = setup["h"]
    init_cond = setup["init_cond"]
    eq_type = setup["eq_type"]
    epsilon = setup["epsilon"]
    basis = setup["basis"]
    basis1 = setup["basis1"]
    x_d = setup["x_d"]
    comprange = setup["comprange"]
    n_comp = setup["n_comp"]
    xrange = setup["xrange"]
    V = setup["V"]
    TranspV = setup["TranspV"]
    DM = setup["DM"]

    if for_bh == false
        J1_choice = J1_type
    else
        J1_choice = J1_type_A
    end
    if J1_choice == "upwind_diss_symm"
        return get_D_upwind_diss_symm(setup, "minus")
    end

    B = zeros(T, cellnumber*2, cellnumber*2)
    DM_J1 = zeros(T, comprange, comprange)
    B_J0 = zeros(T, cellnumber*2, cellnumber*2)
    V_J0E1 = zeros(T, cellnumber*2, comprange)
    V_J0E2 = zeros(T, cellnumber*2, comprange)
    B_J0E1 = zeros(T, cellnumber*2, cellnumber*2)
    B_J0E2 = zeros(T, cellnumber*2, cellnumber*2)
    # Upwind Flux, periodic boundary condition
    B[1, 2*cellnumber] = -1/2 *2 # lila
    for i in 1:cellnumber
        # Upwind Flux
        if i != 1
            B[2*i-1, 2*(i-1)] = -1/2 *2 #lila
        end
        B[2*i, 2*i] = 1/2 *2 # pink
        #DoD-Terms
        if stabilize_condition(setup, i)
            eta = determine_eta(setup, i)
            #Extrapolation to correct boundaries for J0-Terms
            for j in 1:deg+1
                basis1.nodes[j] = x_d[(deg+1)*(i-2) + j] 
            end
            V_J0E1[2*(i-1)+2, (deg+1)*(i-2) + 1:(deg+1)*(i-1)] = interpolation_matrix(v[i+1], basis1) 
            V_J0E2[2*(i-1)+3, (deg+1)*(i-2) + 1:(deg+1)*(i-1)] = interpolation_matrix(v[i+1], basis1)
            # J0-Terms                 
            B_J0[2*i, 2*i] = -eta
            B_J0[2*i + 1, 2*i] = eta
            B_J0E1[2*i, 2*i] = eta
            B_J0E2[2*i + 1, 2*i + 1] = -eta
            # J1-Terms
            ac = (deg+1)*(i-1) + 1:(deg+1)*i
            aL = (deg+1)*(i-2) + 1:(deg+1)*(i-1)
            aR = (deg+1)*i + 1:(deg+1)*(i+1)
            Vol_interpol_left, Vol_interpol_right = get_Vol_interpol_matrices(setup, i) 
            loc_DM = basis.D' * diagm(basis.weights) * eta
            if J1_choice == "upwind_classic" || J1_choice == "central_classic"
                # flux to right
                DM_J1[ac, ac] += - loc_DM
                DM_J1[ac, aL] += loc_DM * Vol_interpol_left
                # extra stability terms
                if ext_test_func == true
                    # flux to right
                    DM_J1[aL, ac] += Vol_interpol_left' *  loc_DM
                    DM_J1[aL, aL] +=  - Vol_interpol_left' * loc_DM * Vol_interpol_left
                end
            elseif J1_choice == "upwind_symm"
                # Extra Factor (-1) because its defined by LHS 
                loc_DM = loc_DM * (-1)
                DM_J1[ac, aL] += - loc_DM * Vol_interpol_left
                DM_J1[ac, ac] += loc_DM
                if ext_test_func == true
                    DM_J1[aR, aL] += 1/2 * Vol_interpol_right' * loc_DM * Vol_interpol_left
                    DM_J1[aR, aR] += - 1/2 * Vol_interpol_right' * loc_DM * Vol_interpol_right
                    DM_J1[aL, aL] += 1/2 * Vol_interpol_left' * loc_DM * Vol_interpol_left
                    DM_J1[aL, ac] += - Vol_interpol_left' * loc_DM
                    DM_J1[aL, aR] += 1/2 * Vol_interpol_left' * loc_DM * Vol_interpol_right
                end
            end
        end

    end
    Dminus = V'* (B*V + B_J0*V + B_J0E1*V_J0E1 + B_J0E2*V_J0E2) - DM - DM_J1
    return Dminus
end
###############################################################################
###############################################################################

function get_Dplus(setup; for_bh = false)
    deg = setup["deg"]
    nodes = setup["nodes"]
    T = setup["T"]
    do_stabilize = setup["do_stabilize"]
    c = setup["c"]
    fix_eta = setup["fix_eta"]
    fluxtype = setup["fluxtype"]
    include_b_h = setup["include_b_h"]
    ext_test_func = setup["ext_test_func"]
    J1_type = setup["J1_type"]
    J1_type_A = setup["J1_type_A"]
    cellnumber = setup["cellnumber"]
    v = setup["v"]
    a = setup["a"]
    u0_eval = setup["u0_eval"]
    t_d = setup["t_d"]
    h = setup["h"]
    init_cond = setup["init_cond"]
    eq_type = setup["eq_type"]
    epsilon = setup["epsilon"]
    basis = setup["basis"]
    basis1 = setup["basis1"]
    x_d = setup["x_d"]
    comprange = setup["comprange"]
    n_comp = setup["n_comp"]
    xrange = setup["xrange"]
    V = setup["V"]
    TranspV = setup["TranspV"]
    DM = setup["DM"]

    if for_bh == false
        J1_choice = J1_type
    else
        J1_choice = J1_type_A
    end
    
    if J1_choice == "upwind_diss_symm"
        return get_D_upwind_diss_symm(setup, "plus")
    end
    
    B = zeros(T, cellnumber*2, cellnumber*2)
    DM_J1 = zeros(T, comprange, comprange)
    B_J0 = zeros(T, cellnumber*2, cellnumber*2)
    V_J0E1 = zeros(T, cellnumber*2, comprange)
    V_J0E2 = zeros(T, cellnumber*2, comprange)
    B_J0E1 = zeros(T, cellnumber*2, cellnumber*2)
    B_J0E2 = zeros(T, cellnumber*2, cellnumber*2)
    # Upwind Flux, periodic boundary condition
    B[2*cellnumber, 1] = 1 #blau
    for i in 1:cellnumber
        # Upwind Flux
        B[2*i-1, 2*i-1] = -1  # grün
        if i != cellnumber
            B[2*i, 2*i+1] = 1  # blau
        end
        #DoD-Terms
        if stabilize_condition(setup, i)
            eta = determine_eta(setup, i)
            #Extrapolation to correct boundaries for J0-Terms
            for j in 1:deg+1
                basis1.nodes[j] = x_d[(deg+1)*(i) + j]
            end
            V_J0E1[2*i-1, (deg+1)*i + 1:(deg+1)*(i+1)] = interpolation_matrix(v[i], basis1) 
            V_J0E2[2*i-2, (deg+1)*i + 1:(deg+1)*(i+1)] = interpolation_matrix(v[i], basis1)
            # J0-Terms                 
            B_J0[2*i - 1, 2*i - 1] = eta
            B_J0[2*i - 2, 2*i - 1] = -eta
            B_J0E1[2*i - 1, 2*i - 1] = -eta
            B_J0E2[2*i - 2, 2*i - 2] = eta
            if deg>0
                basisleft=nodes(deg)
                basisright=nodes(deg)
                for j in 1:deg+1
                    basisleft.nodes[j] = x_d[(deg+1)*(i-2) + j] 
                    basisright.nodes[j] = x_d[(deg+1)*(i) + j] 
                end
                # J1-Terms
                ac = (deg+1)*(i-1) + 1:(deg+1)*i
                aL = (deg+1)*(i-2) + 1:(deg+1)*(i-1)
                aR = (deg+1)*i + 1:(deg+1)*(i+1)
                Vol_interpol_left, Vol_interpol_right = get_Vol_interpol_matrices(setup, i) 
                loc_DM = basis.D' * diagm(basis.weights) *eta
                if J1_choice == "upwind_classic" || J1_choice == "central_classic"
                    # flux to left
                    DM_J1[ac, ac] += - loc_DM
                    DM_J1[ac, aR] += loc_DM * Vol_interpol_right
                    # extra stability terms
                    if ext_test_func == true
                    # flux to left
                    DM_J1[aR, ac] += Vol_interpol_right' *  loc_DM
                    DM_J1[aR, aR] += - Vol_interpol_right' * loc_DM * Vol_interpol_right
                    end
                elseif J1_choice == "upwind_symm"
                    # Extra Factor (-1) because its defined by LHS 
                    loc_DM = loc_DM * (-1)
                    DM_J1[ac, aR] += - loc_DM * Vol_interpol_right
                    DM_J1[ac, ac] += loc_DM
                    if ext_test_func == true
                        DM_J1[aL, aR] += 1/2 * Vol_interpol_left' * loc_DM * Vol_interpol_right
                        DM_J1[aL, aL] += - 1/2 * Vol_interpol_left' * loc_DM * Vol_interpol_left
                        DM_J1[aR, aR] += 1/2 * Vol_interpol_right' * loc_DM * Vol_interpol_right
                        DM_J1[aR, ac] += - Vol_interpol_right' * loc_DM
                        DM_J1[aR, aL] += 1/2 * Vol_interpol_right' * loc_DM * Vol_interpol_left
                    end
                end
            end
        end
    end
    Dplus = V'* (B*V + B_J0*V + B_J0E1*V_J0E1 + B_J0E2*V_J0E2) - DM - DM_J1
    return Dplus
end

function get_DM_J1_for_central_symm(setup)
    deg = setup["deg"]
    nodes = setup["nodes"]
    T = setup["T"]
    do_stabilize = setup["do_stabilize"]
    c = setup["c"]
    fix_eta = setup["fix_eta"]
    fluxtype = setup["fluxtype"]
    include_b_h = setup["include_b_h"]
    ext_test_func = setup["ext_test_func"]
    J1_type = setup["J1_type"]
    J1_type_A = setup["J1_type_A"]
    cellnumber = setup["cellnumber"]
    v = setup["v"]
    a = setup["a"]
    u0_eval = setup["u0_eval"]
    t_d = setup["t_d"]
    h = setup["h"]
    init_cond = setup["init_cond"]
    eq_type = setup["eq_type"]
    epsilon = setup["epsilon"]
    basis = setup["basis"]
    basis1 = setup["basis1"]
    x_d = setup["x_d"]
    comprange = setup["comprange"]
    n_comp = setup["n_comp"]
    xrange = setup["xrange"]
    V = setup["V"]
    TranspV = setup["TranspV"]
    DM = setup["DM"]

    DM_J1 = zeros(T, comprange, comprange)
    for i in 1:cellnumber
        if stabilize_condition(setup, i)
            eta = determine_eta(setup, i)
            if deg>0
                basisleft=nodes(deg)
                basisright=nodes(deg)
                for j in 1:deg+1
                    basisleft.nodes[j] = x_d[(deg+1)*(i-2) + j] 
                    basisright.nodes[j] = x_d[(deg+1)*(i) + j] 
                end
                # J1-Terms
                ac = (deg+1)*(i-1) + 1:(deg+1)*i
                aL = (deg+1)*(i-2) + 1:(deg+1)*(i-1)
                aR = (deg+1)*i + 1:(deg+1)*(i+1)
                Vol_interpol_left, Vol_interpol_right = get_Vol_interpol_matrices(setup, i)
                # Because of centralized version a factor 1/2, Factor (-1) because its defined by LHS 
                loc_DM = basis.D' * diagm(basis.weights) * 1/2 * eta * (-1)
                DM_J1[ac, ac] += 2 * loc_DM
                DM_J1[ac, aR] +=  - loc_DM * Vol_interpol_right
                DM_J1[ac, aL] +=  - loc_DM * Vol_interpol_left 
                if ext_test_func == true
                    DM_J1[aL, aR] += Vol_interpol_left' *  loc_DM * Vol_interpol_right 
                    DM_J1[aL, ac] += - Vol_interpol_left' *  loc_DM
                    DM_J1[aR, aL] += Vol_interpol_right' *  loc_DM * Vol_interpol_left 
                    DM_J1[aR, ac] += - Vol_interpol_right' *  loc_DM 
                end
            end
        end
    end
    return DM_J1
end


function get_D_central(setup)
    J1_type = setup["J1_type"]
    if J1_type == "central_classic"
        D_central = (get_Dminus(setup, for_bh = false) + get_Dplus(setup, for_bh = false))/2
    elseif J1_type == "central_symm"
        # For this case no J1-Terms get constructed inside get_Dplus/get_Dminus
        DM_J1 = get_DM_J1_for_central_symm(setup)
        D_central = (get_Dminus(setup, for_bh = false) + get_Dplus(setup, for_bh = false))/2 - DM_J1
    end
    return D_central
end


function get_sourceterm_telegraph(setup)
    deg = setup["deg"]
    T = setup["T"]
    cellnumber = setup["cellnumber"]
    basis = setup["basis"]
    comprange = setup["comprange"]
    v = setup["v"]
    scale(i) = (v[i+1]-v[i])/2

    M_source = zeros(T, comprange, comprange)
    for i in 1:cellnumber
        M_source[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(i-1) + 1:(deg+1)*i] = scale(i)*diagm(basis.weights)
    end
    return M_source
end

function get_D_upwind_diss_symm(setup, direction)
    deg = setup["deg"]
    nodes = setup["nodes"]
    T = setup["T"]
    do_stabilize = setup["do_stabilize"]
    c = setup["c"]
    fix_eta = setup["fix_eta"]
    fluxtype = setup["fluxtype"]
    include_b_h = setup["include_b_h"]
    ext_test_func = setup["ext_test_func"]
    J1_type = setup["J1_type"]
    J1_type_A = setup["J1_type_A"]
    cellnumber = setup["cellnumber"]
    v = setup["v"]
    a = setup["a"]
    u0_eval = setup["u0_eval"]
    t_d = setup["t_d"]
    h = setup["h"]
    init_cond = setup["init_cond"]
    eq_type = setup["eq_type"]
    epsilon = setup["epsilon"]
    basis = setup["basis"]
    basis1 = setup["basis1"]
    x_d = setup["x_d"]
    comprange = setup["comprange"]
    n_comp = setup["n_comp"]
    xrange = setup["xrange"]
    V = setup["V"]
    TranspV = setup["TranspV"]
    DM = setup["DM"]

    # construct D with central_symm flux
    DM_J1 = get_DM_J1_for_central_symm(setup)
    setup2 = copy(setup)
    setup2["J1_type"] = "central_symm"
    D_central = (get_Dminus(setup2, for_bh = false) + get_Dplus(setup2, for_bh = false))/2 - DM_J1
    left_interpolation = interpolation_matrix(-1, basis)
    right_interpolation = interpolation_matrix(1, basis)


    # adding the remaining dissipation background flux terms
    dissipation_flux_terms = zeros(T, comprange, comprange)
    dissipation_flux_terms[1:(deg+1), (deg+1)*(cellnumber-1) + 1:(deg+1)*cellnumber] += -1/2 * left_interpolation' * right_interpolation
    dissipation_flux_terms[(deg+1)*(cellnumber-1) + 1:(deg+1)*cellnumber, 1:(deg+1)] += -1/2 * right_interpolation' * left_interpolation
    for i in 1:cellnumber
        dissipation_flux_terms[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(i-1) + 1:(deg+1)*i] += 1/2 * right_interpolation' * right_interpolation
        dissipation_flux_terms[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(i-1) + 1:(deg+1)*i] += 1/2 * left_interpolation' * left_interpolation
        if i>1
            dissipation_flux_terms[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(i-2) + 1:(deg+1)*(i-1)] += -1/2 * left_interpolation' * right_interpolation
        end
        if i<cellnumber
            dissipation_flux_terms[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*i + 1:(deg+1)*(i+1)] += -1/2 * right_interpolation' * left_interpolation
        end
    end

    # adding remaining Stabilization terms fpr J0 and J1
    dissipation_J0_terms = zeros(T, comprange, comprange)
    dissipation_J0_terms1 = zeros(T, comprange, comprange)
    dissipation_J0_terms2 = zeros(T, comprange, comprange)
    dissipation_J1_terms = zeros(T, comprange, comprange)
    for i in 1:cellnumber
        if stabilize_condition(setup, i)
            eta = determine_eta(setup, i)
            aL = (deg+1)*(i-2) + 1:(deg+1)*(i-1)
            ac = (deg+1)*(i-1) + 1:(deg+1)*i
            aR = (deg+1)*i + 1:(deg+1)*(i+1)

            #Extrapolation to correct boundaries for J0-Terms
            for j in 1:deg+1
                basis1.nodes[j] = x_d[(deg+1)*(i-2) + j] 
            end
            B_J0cm1 = interpolation_matrix(v[i+1], basis1) 
            #Extrapolation to correct boundaries for J0-Terms
            for j in 1:deg+1
                basis1.nodes[j] = x_d[(deg+1)*(i) + j]
            end
            B_J0cp1 = interpolation_matrix(v[i], basis1)
            loc_fac = eta/2
            dissipation_J0_terms1[ac, aL] += loc_fac * right_interpolation' * B_J0cm1
            dissipation_J0_terms1[ac, ac] += loc_fac * -right_interpolation' * right_interpolation
            dissipation_J0_terms1[aR, aL] += loc_fac * -left_interpolation' * B_J0cm1
            dissipation_J0_terms1[aR, ac] += loc_fac * left_interpolation' * right_interpolation

            dissipation_J0_terms2[aL, ac] += loc_fac * right_interpolation' * left_interpolation
            dissipation_J0_terms2[aL, aR] += loc_fac * -right_interpolation' * B_J0cp1
            dissipation_J0_terms2[ac, ac] += loc_fac * -left_interpolation' * left_interpolation
            dissipation_J0_terms2[ac, aR] += loc_fac * left_interpolation' * B_J0cp1

            # J1-Terms
            Vol_interpol_left, Vol_interpol_right = get_Vol_interpol_matrices(setup, i) 
            loc_DM = basis.D' * diagm(basis.weights) * eta/2
            # constructed as direction minus -> for other direction, there is just a sign swap needed
            dissipation_J1_terms[ac, aR] += loc_DM * Vol_interpol_right
            dissipation_J1_terms[ac, aL] += - loc_DM * Vol_interpol_left 
            if ext_test_func == true
                dissipation_J1_terms[aR, ac] += Vol_interpol_right' * loc_DM 
                dissipation_J1_terms[aR, aR] += - Vol_interpol_right' * loc_DM * Vol_interpol_right
                dissipation_J1_terms[aL, aL] += Vol_interpol_left' * loc_DM * Vol_interpol_left
                dissipation_J1_terms[aL, ac] += - Vol_interpol_left' * loc_DM  
            end
        end
    end
    dissipation_J0_terms = dissipation_J0_terms1 + dissipation_J0_terms2
    
    # Symmetrizing Dissipation terms
    if direction == "minus"
        dissipation_terms_all =  dissipation_flux_terms + dissipation_J0_terms + dissipation_J1_terms
    elseif direction == "plus"
        dissipation_terms_all =  - (dissipation_flux_terms + dissipation_J0_terms + dissipation_J1_terms)
    else
        throw(ArgumentError("Wrong dircetions chosen in fluxtype upwind_diss_symm"))
    end
    dissipation_terms_symmetrized = 1/2*(dissipation_terms_all + dissipation_terms_all')

    if direction == "minus"
        global dissm = dissipation_terms_symmetrized
    elseif direction == "plus"
        global dissp = dissipation_terms_symmetrized
    end
    return D_central + dissipation_terms_symmetrized
end