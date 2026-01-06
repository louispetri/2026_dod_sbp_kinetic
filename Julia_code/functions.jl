function setup_problem_eq(eq_type, init_cond, xmin, xmax, cellnumber; T = Float64, Tmax = 1, a = 1, CFL = 0.9, bcs = "periodic", epsilon = 1.0)
    vertices = Array(range(xmin, xmax, length = cellnumber + 1))
    h = vertices[2]-vertices[1]
    t_d = range(0,Tmax, length = maximum([Integer(round(abs(Tmax/(CFL*h/a)))), 4]))
    # Setting initial condition as function
    random_ic = zeros(T, 3*cellnumber)
    if init_cond == "random"
        random_ic = randn(T, 3*cellnumber)
    end
    function u0_eval(x)
        if init_cond == "sin"        
            return sin(2*pi*x)
        elseif init_cond == "cos"
            return cos(4*pi*x)
        elseif init_cond == "squared"
            return (x-1/2)^2
        elseif init_cond == "const"
            return 1.0
        elseif init_cond == "linear"
            return x
        elseif init_cond == "hat"
            return hat(x)
        elseif init_cond == "discont"
            if x<1/2
                return -1.0
            else
                return 1.0
            end
        elseif init_cond == "plateau"
            if x >= 0.45 && x <= 0.55
                rho = 1.0
            else
                rho = 0.0
            end
        elseif init_cond in ["telsin", "telegraph_well_prep"]
                r = -2/(1+sqrt(1-4*epsilon^2))
            return 1/r*sin(x)
        elseif init_cond == "eq_data"
            return sin(x)
        elseif init_cond == "random"
            return random_ic[ceil(Int, abs(x/(xmax-xmin))*cellnumber)+1]
        end
    end
    
    problem = Dict("a" => a, "u0_eval" => u0_eval, "v" => vertices, "Tmax" => Tmax, "h" => h, "CFL" => CFL,
                    "cellnumber" => cellnumber, "t_d" => t_d, "init_cond" => init_cond, "xmin" => xmin,
                     "xmax" => xmax, "bcs" => bcs, "epsilon" => epsilon, "eq_type" => eq_type, "T" => T)
    return problem
end

function hat(x)
    if x > 0.4 && x < 0.6
        return 0.5-abs(x-0.5)
    else
        return 0.4
    end
end

function to_periodic(x)
    if x > 1 
        return x%1
    elseif x<0
        while x<0
            x = x+1
        end
        return x
    else 
        return x
    end
end

function get_exact_solution(setup)
    eq_type = setup["eq_type"]
    u0_eval = setup["u0_eval"]
    x_d = setup["x_d"]
    a = setup["a"]
    t_d = setup["t_d"]
    N = setup["cellnumber"]
    deg = setup["deg"]
    epsilon = setup["epsilon"]
    T = setup["T"]
    if setup["bcs"] == "periodic"
        if eq_type == "wave"
            u_wave_exact = zeros(T, (deg+1)*2*N, length(t_d));
            u_wave_exact[1:(deg+1)*N, :] = [1/2*(sin(2*pi*(x-1/epsilon*t))+epsilon*cos(2*pi*(x-1/epsilon*t))+sin(2*pi*(x+1/epsilon*t))-epsilon*cos(2*pi*(x+1/epsilon*t))) for x in x_d, t in t_d];
            u_wave_exact[(deg+1)*N+1:(deg+1)*2*N, :] = [1/2*(1/epsilon*sin(2*pi*(x-1/epsilon*t))+cos(2*pi*(x-1/epsilon*t))-1/epsilon*sin(2*pi*(x+1/epsilon*t))+cos(2*pi*(x+1/epsilon*t))) for x in x_d, t in t_d];
            return u_wave_exact
        elseif eq_type == "wave_symm"
            u_wave_exact = zeros(T, (deg+1)*2*N, length(t_d));
            u_wave_exact[1:(deg+1)*N, :] = [1/2*(sin(2*pi*(x-1/epsilon*t))+epsilon*cos(2*pi*(x-1/epsilon*t))+sin(2*pi*(x+1/epsilon*t))-epsilon*cos(2*pi*(x+1/epsilon*t))) for x in x_d, t in t_d];
            u_wave_exact[(deg+1)*N+1:(deg+1)*2*N, :] = [1/2*epsilon*(1/epsilon*sin(2*pi*(x-1/epsilon*t))+cos(2*pi*(x-1/epsilon*t))-1/epsilon*sin(2*pi*(x+1/epsilon*t))+cos(2*pi*(x+1/epsilon*t))) for x in x_d, t in t_d];
            return u_wave_exact
        elseif eq_type in ["telegraph", "telegraph_semistab"]
            if setup["init_cond"] != "telsin"
                println("Care: exact solution does not match inital conditions!")
            end
            r = -2/(1+sqrt(1-4*epsilon^2))
            u_tel_exact = zeros(T, (deg+1)*2*N, length(t_d));
            u_tel_exact[1:(deg+1)*N, :] = [1/r*exp(r*t)*sin(x) for x in x_d, t in t_d];
            u_tel_exact[(deg+1)*N+1:(deg+1)*2*N, :] = [exp(r*t)*cos(x) for x in x_d, t in t_d];
            return u_tel_exact
        elseif eq_type == "telegraph_symm"
            if setup["init_cond"] != "telsin"
                println("Care: exact solution does not match inital conditions!")
            end
            r = -2/(1+sqrt(1-4*epsilon^2))
            u_tel_exact = zeros(T, (deg+1)*2*N, length(t_d));
            u_tel_exact[1:(deg+1)*N, :] = [1/r*exp(r*t)*sin(x) for x in x_d, t in t_d];
            u_tel_exact[(deg+1)*N+1:(deg+1)*2*N, :] = [epsilon*exp(r*t)*cos(x) for x in x_d, t in t_d];
            return u_tel_exact
        elseif eq_type == "heat"
            r = -1
            u_tel_exact = zeros(T, (deg+1)*N, length(t_d));
            u_tel_exact[1:(deg+1)*N, :] = [1/r*exp(r*t)*sin(x) for x in x_d, t in t_d];
        elseif eq_type == "transport"
            u_exact = [u0_eval(to_periodic(x-(a*(t)))) for x in x_d, t in t_d]
            return u_exact
        else
            println("wrong eq type")
        end
    end
end

function f_RHS(u, t, RHS_mat, setup)
    return RHS_mat * u
end

function include_cut_cell(setup, cutcell_size, cutcell_pos)
    T = setup["T"]
    vert = setup["v"]
    create_cutcell = false
    if cutcell_size <=0.5
        if vert[cutcell_pos+1]-vert[cutcell_pos] > 1/2* setup["h"] && vert[cutcell_pos]-vert[cutcell_pos-1] > 1/2* setup["h"]
            create_cutcell = true
        end
    end
    if cutcell_size >0.5
        if vert[cutcell_pos+2]-vert[cutcell_pos + 1] > 1/2* setup["h"] && vert[cutcell_pos + 1]-vert[cutcell_pos] > 1/2* setup["h"]
            create_cutcell = true
        end
    end
    if create_cutcell
        vert_new = zeros(T, length(vert)+1)
        for i in range(1, length(vert)+1)
            if i < cutcell_pos
                vert_new[i] = vert[i] 
            elseif i == cutcell_pos
                vert_new[i] = vert[i] 
                vert_new[i+1] = ((1-cutcell_size)*vert[cutcell_pos] + (cutcell_size)*vert[cutcell_pos+1])
            elseif i > cutcell_pos + 1
                vert_new[i] = vert[i-1]
            end
        end
        setup["v"] = vert_new
        setup["cellnumber"] += 1
    else
        println("No new cut-cell created due to its neighbour being already a cut-cell")
    end
    return setup
end



####################################################################################################################################################
########################################################                          ##################################################################
########################################################    Telegraph equation    ##################################################################
########################################################                          ##################################################################
####################################################################################################################################################
function mymod(x, y)
    if x % y == 0
        return y
    else 
        return x % y
    end
end

function DGsemidiscretization_DoD_telegraph(setup, deg, nodes;
                                            do_stabilize = false, c = 1.0,
                                            fix_eta = true, fluxtype = "full", include_b_h = true,
                                            ext_test_func = true, J1_type = "upwind_classic", J1_type_A = "upwind_classic")
    # Parameter einlesen
    setup["deg"] = deg
    setup["nodes"] = nodes
    setup["do_stabilize"] = do_stabilize
    setup["c"] = c
    setup["fix_eta"] = fix_eta
    setup["fluxtype"] = fluxtype
    setup["include_b_h"] = include_b_h
    setup["ext_test_func"] = ext_test_func
    setup["J1_type"] = J1_type
    setup["J1_type_A"] = J1_type_A
    # Parameter auslesen
    T = setup["T"]
    fluxtype = setup["fluxtype"]
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

    setup["bcs"] = "periodic"

    splitRHS = Dict()
    SBP_storage = Dict()

    # basis defines the method, basis1 will be manipulated to serve interpolation matrices for other nodes
    basis = nodes(deg, T)
    basis1 = nodes(deg, T) # to be adapted for cut cell specification
    setup["basis"] = basis
    setup["basis1"] = basis1

    eq_types_1D = ["transport", "heat"]
    eq_types_2D = ["telegraph", "telegraph_symm", "wave", "wave_symm", "telegraph_semistab"]
    eq_types_tel = ["telegraph", "telegraph_symm", "telegraph_semistab"]

    # Some Parameters
    comprange = cellnumber*(deg+1)
    setup["comprange"] = comprange
    if eq_type in eq_types_1D
        n_comp = 1
    elseif eq_type in eq_types_2D
        n_comp = 2
    end
    setup["n_comp"] = n_comp
    xrange = n_comp*cellnumber*(deg+1)
    setup["xrange"] = xrange
    shift(i) = (v[i+1]+v[i])/2
    scale(i) = (v[i+1]-v[i])/2
    invscale(i) = 2/(v[i+1]-v[i])
    # Constructing the nodes
    x_d = zeros(T, comprange)
    for i in 1:cellnumber
        for j in 1:deg+1
            x_d[(deg+1)*(i-1) + j] = basis.nodes[j]*scale(i) + shift(i)
        end
    end
    setup["x_d"] = x_d

    Sc = zeros(T, xrange, xrange)
    M = zeros(T, xrange, xrange)
    M_global = zeros(T, xrange, xrange)
    # Inverse of the LHS mass-matrix 
    for i in 1:n_comp*cellnumber
        Sc[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(i-1) + 1:(deg+1)*i] = diagm(ones(T, deg+1)*invscale(mymod(i, cellnumber)))
        M[(deg+1)*(i-1)+1:(deg+1)*i, (deg+1)*(i-1)+1:(deg+1)*i] = diagm(basis.weights)
        M_global[(deg+1)*(i-1)+1:(deg+1)*i, (deg+1)*(i-1)+1:(deg+1)*i] = diagm(basis.weights) * scale(mymod(i, cellnumber))
    end
    invM_global = Sc*inv(M)
    setup["invM_global"] = invM_global
    # Cell boundary extrapolation operator V (for one component) and standard volume Term D'M
    left_interpolation = interpolation_matrix(-1 , basis)
    right_interpolation = interpolation_matrix(1, basis)
    V = zeros(T, cellnumber*2, comprange)
    TranspV = zeros(T, comprange, cellnumber*2)
    DM = zeros(T, comprange, comprange)
    for i in 1:cellnumber
        V[2*i-1, (deg+1)*(i-1) + 1:(deg+1)*i] = left_interpolation
        V[2*i, (deg+1)*(i-1) + 1:(deg+1)*i] = right_interpolation
        TranspV[(deg+1)*(i-1)+1:(deg+1)*i, 2*i-1:2*i] = V[2*i-1:2*i, (deg+1)*(i-1) + 1:(deg+1)*i]'
        DM[(deg+1)*(i-1) + 1:(deg+1)*i, (deg+1)*(i-1) + 1:(deg+1)*(i)] = basis.D' * diagm(basis.weights)
    end
    setup["V"] = V
    setup["TranspV"] = TranspV
    setup["DM"] = DM

    # RHS-Matrix anhand des Differentialgleichungssystems bestimmen
    RHS_mat = zeros(T, xrange, xrange)
    ex_RHS_mat = zeros(T, xrange, xrange)
    im_RHS_mat = zeros(T, xrange, xrange)
    if eq_type == "transport"
        if fluxtype == "upwind"
            Dminus = get_Dminus(setup, for_bh = false)
            Dplus = get_Dplus(setup, for_bh = false)
            SBP_storage["Dminus"] = Dminus
            SBP_storage["Dplus"] = Dplus
            if a>0
                RHS_mat = - invM_global * a*Dminus
            else
                RHS_mat = - invM_global * a*Dplus
            end
        end
        if fluxtype == "central"
            D_central = get_D_central(setup)
            SBP_storage["Dc"] = D_central
            RHS_mat = - invM_global * a*D_central
        end
    elseif eq_type == "wave"
        if fluxtype in ["altlr", "altrl"]
            Dminus = get_Dminus(setup, for_bh = false)
            Dplus = get_Dplus(setup, for_bh = false)
            SBP_storage["Dminus"] = Dminus
            SBP_storage["Dplus"] = Dplus
            if fluxtype == "altlr"
                ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -Dminus
                im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon^2*Dplus
            elseif fluxtype == "altrl"
                ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -Dplus
                im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon^2*Dminus
            end
        elseif fluxtype == "central"
            D_central = get_D_central(setup)
            SBP_storage["Dc"] = D_central
            ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -D_central
            im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon^2*D_central
        end
        ex_RHS_mat = invM_global*ex_RHS_mat
        im_RHS_mat = invM_global*im_RHS_mat
        RHS_mat = ex_RHS_mat + im_RHS_mat
    elseif eq_type == "wave_symm"
        if fluxtype in ["altlr", "altrl"]
            Dminus = get_Dminus(setup, for_bh = false)
            Dplus = get_Dplus(setup, for_bh = false)
            SBP_storage["Dminus"] = Dminus
            SBP_storage["Dplus"] = Dplus
            if fluxtype == "altlr"
                ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -1/epsilon*Dminus
                im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon*Dplus
            elseif fluxtype == "altrl"
                ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -1/epsilon*Dplus
                im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon*Dminus
            end
        elseif fluxtype == "central"
            D_central = get_D_central(setup)
            SBP_storage["Dc"] = D_central
            ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -1/epsilon*D_central
            im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon*D_central
        end
        ex_RHS_mat = invM_global*ex_RHS_mat
        im_RHS_mat = invM_global*im_RHS_mat
        RHS_mat = ex_RHS_mat + im_RHS_mat
    elseif eq_type in eq_types_tel
        Dminus_A = get_Dminus(setup, for_bh = true)
        Dplus_A = get_Dplus(setup, for_bh = true)
        A = zeros(T, comprange, comprange)
        M_source = get_sourceterm_telegraph(setup)
        if include_b_h== true
            A = 1/2*(Dplus_A-Dminus_A)
            setup["A"] = A
        end
        if eq_type == "telegraph"
            if fluxtype in ["altlr", "altrl"]
                Dminus = get_Dminus(setup, for_bh = false)
                Dplus = get_Dplus(setup, for_bh = false)
                SBP_storage["Dminus"] = Dminus
                SBP_storage["Dplus"] = Dplus
                setup["Dminus"] = Dminus
                setup["Dplus"] = Dplus
                setup["M_source"] = M_source
                if fluxtype == "altlr"
                    ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -Dminus
                    im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon^2*Dplus
                    ex_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = 1/epsilon*A
                    im_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = -1/epsilon^2*M_source
                elseif fluxtype == "altrl"
                    ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -Dplus
                    im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon^2*Dminus
                    ex_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = 1/epsilon*A
                    im_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = -1/epsilon^2*M_source
                end
            elseif fluxtype == "central"
                D_central = get_D_central(setup)
                SBP_storage["Dc"] = D_central
                setup["Dc"] = D_central
                ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -D_central
                im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon^2*D_central
                ex_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = 1/epsilon*A
                im_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = -1/epsilon^2*M_source
            end
            ex_RHS_mat = invM_global*ex_RHS_mat
            im_RHS_mat = invM_global*im_RHS_mat
            RHS_mat = ex_RHS_mat + im_RHS_mat
        elseif eq_type == "telegraph_semistab"
            setup_nonstab = copy(setup)
            setup_nonstab["do_stabilize"] = false
            if fluxtype in ["altlr", "altrl"]
                Dminus = get_Dminus(setup, for_bh = false)
                Dplus = get_Dplus(setup, for_bh = false)
                Dminus_nonstab = get_Dminus(setup_nonstab, for_bh = false)
                Dplus_nonstab = get_Dplus(setup_nonstab, for_bh = false)
                SBP_storage["Dminus"] = Dminus
                SBP_storage["Dplus"] = Dplus
                if fluxtype == "altlr"
                    ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -Dminus
                    im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon^2*Dplus_nonstab
                    ex_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = 1/epsilon*A
                    im_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = -1/epsilon^2*M_source
                elseif fluxtype == "altrl"
                    ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -Dplus
                    im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon^2*Dminus_nonstab
                    ex_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = 1/epsilon*A
                    im_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = -1/epsilon^2*M_source
                end
            elseif fluxtype == "central"
                D_central = get_D_central(setup)
                D_central_nonstab = get_D_central(setup_nonstab)
                SBP_storage["Dc"] = D_central
                ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -D_central
                im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon^2*D_central_nonstab
                ex_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = 1/epsilon*A
                im_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = -1/epsilon^2*M_source
            end
            ex_RHS_mat = invM_global*ex_RHS_mat
            im_RHS_mat = invM_global*im_RHS_mat
            RHS_mat = ex_RHS_mat + im_RHS_mat
        elseif eq_type == "telegraph_symm"
            if fluxtype in ["altlr", "altrl"]
                Dminus = get_Dminus(setup, for_bh = false)
                Dplus = get_Dplus(setup, for_bh = false)
                SBP_storage["Dminus"] = Dminus
                SBP_storage["Dplus"] = Dplus
                if fluxtype == "altlr"
                    ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -1/epsilon*Dminus
                    im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon*Dplus
                    ex_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = 1/epsilon*A
                    im_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = -1/epsilon^2*M_source
                elseif fluxtype == "altrl"
                    ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -1/epsilon*Dplus
                    im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon*Dminus
                    ex_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = 1/epsilon*A
                    im_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = -1/epsilon^2*M_source
                end
            elseif fluxtype == "central"
                D_central = get_D_central(setup)
                SBP_storage["Dc"] = D_central
                ex_RHS_mat[1:comprange, comprange + 1:2*comprange] = -1/epsilon*D_central
                im_RHS_mat[comprange + 1:2*comprange, 1:comprange] = -1/epsilon*D_central
                ex_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = 1/epsilon*A
                im_RHS_mat[comprange + 1:2*comprange, comprange + 1:2*comprange] = -1/epsilon^2*M_source
            end
            ex_RHS_mat = invM_global*ex_RHS_mat
            im_RHS_mat = invM_global*im_RHS_mat
            RHS_mat = ex_RHS_mat + im_RHS_mat
        end
    elseif eq_type == "heat"
        if fluxtype in ["altlr", "altrl"]
            Dminus = get_Dminus(setup, for_bh = false)
            Dplus = get_Dplus(setup, for_bh = false)
            SBP_storage["Dminus"] = Dminus
            SBP_storage["Dplus"] = Dplus
            if fluxtype == "altlr"
                RHS_mat = invM_global * Dminus * invM_global * Dplus
            elseif fluxtype == "altrl"
                RHS_mat = invM_global * Dplus * invM_global * Dminus
            end
        elseif fluxtype == "central"
            D_central = get_D_central(setup)
            SBP_storage["Dc"] = D_central
            RHS_mat = invM_global * D_central * invM_global * D_central
        end
    end

    setup["ex_RHS_mat"] = ex_RHS_mat
    setup["im_RHS_mat"] = im_RHS_mat

    # Constructing the initial conition (special case for j coordinate)
    u0_rho = [u0_eval(x) for x in x_d]
    if eq_type == "heat"
        if init_cond == "telsin"
            u0 = [-sin(x) for x in x_d]
        elseif init_cond == "eq_data"
            u0 = [sin(x) for x in x_d]
        end
    elseif eq_type == "transport"
         u0 = u0_rho
    else
        if fluxtype == "altlr"
            D_ic = get_Dplus(setup, for_bh = false)
        elseif fluxtype == "altrl"
            D_ic = get_Dminus(setup, for_bh = false)
        elseif fluxtype == "central"
            D_ic = get_D_central(setup)
        end
        if eq_type in ["telegraph", "telegraph_semistab"]
            if init_cond == "telsin"
                r = 1                                                                      # IC from Jang paper
                u0_j = 1/r*cos.(x_d)                                                       
            elseif init_cond == "telegraph_well_prep"
                u0_j = -  invM_global[1:length(u0_rho), 1:length(u0_rho)]*D_ic*u0_rho       # well prepared
            elseif init_cond == "eq_data"
                u0_j = cos.(x_d)
            end
        elseif eq_type == "telegraph_symm"
            if init_cond == "telsin"
                u0_j = epsilon*cos.(x_d)
            elseif init_cond == "telegraph_well_prep"
                u0_j = -  epsilon*invM_global[1:length(u0_rho), 1:length(u0_rho)]*D_ic*u0_rho   
            elseif init_cond == "eq_data"
                u0_j = epsilon*cos.(x_d)
            end
        elseif eq_type == "wave"
            u0_j = cos.(2*pi*x_d)
        elseif eq_type == "wave_symm"
            u0_j = cos.(2*pi*x_d)*epsilon
        end
        u0 = vcat(u0_rho, u0_j);
    end



    setup["u0"] = u0
    setup["M"] = M_global
    setup["x_d"] = x_d
    setup["basis"] = basis
    return RHS_mat, setup, SBP_storage
end



#### for Simulation
function determine_Nx_plot(problem)
    u0 = problem["u0"]
    if problem["eq_type"] in ["telegraph", "wave", "telegraph_symm", "wave_symm", "telegraph_semistab"]
        Nx_plot =  Int(length(u0)/2)
    elseif problem["eq_type"] == "transport" || problem["eq_type"] == "heat"
        Nx_plot = length(u0)
    else 
        Nx_plot = nothing
    end
    return Nx_plot
end