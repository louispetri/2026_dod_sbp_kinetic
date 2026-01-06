function ARS_optimized(setup, RHS_mat, Tableau)
    println("Caution! usage is limited")
    T = setup["T"]
    # Instead of Nx, use the following notation (where dim = Nx for system, compdim = Nx for scalar equations)
    ex_RHS_mat = setup["ex_RHS_mat"]
    if setup["eq_type"] in ["heat", "transport"]
        dim = size(ex_RHS_mat/2)[1]
        compdim = Int(dim)
    else
        dim = size(ex_RHS_mat/2)[1]
        compdim = Int(dim/2)
    end
    A_ex = Tableau["A_ex"]
    A_im = Tableau["A_im"]
    b_ex = Tableau["b_ex"]
    b_im = Tableau["b_im"]
    s = length(b_im)
    # General Initialisation
    x_d = setup["x_d"]
    a = setup["a"]
    u0 = setup["u0"]
    t_d = setup["t_d"]
    epsilon = setup["epsilon"]
    eq_type = setup["eq_type"]
    include_b_h = setup["include_b_h"]
    fluxtype = setup["fluxtype"]
    u_exact = get_exact_solution(setup)
    #
    invM = setup["invM_global"][1:compdim, 1:compdim]
    if include_b_h == true
        A = invM*setup["A"]
    else
        A = zeros(T,compdim, compdim)
    end
    if fluxtype == "altlr"
        Drho = invM*setup["Dminus"]
        Dg = invM*setup["Dplus"]
    elseif fluxtype == "altrl"
        Drho = invM*setup["Dplus"]
        Dg = invM*setup["Dg"]
    elseif fluxtype == "central"
        Drho = invM*setup["Dc"]
        Dg = invM*setup["Dc"]
    end
    #u_exact = 0
    Nt = length(t_d)
    u_sol = zeros(T, dim, Nt)
    if eq_type == "telegraph"
        # Pseudo-ImEx specialisation
        u_sol_rho = zeros(T, compdim, Nt)
        u_sol_j = zeros(T, compdim, Nt)
        for it in range(1,Nt, step = 1)
            if it == 1
                u_sol_rho[:,it] = u0[1:compdim, it]
                u_sol_j[:,it] = u0[compdim + 1:dim, it]
                u_sol[:, it] = vcat(u_sol_rho[:, it], u_sol_j[:, it])
            else
                dt = (t_d[it]-t_d[it-1])
                # define d_i
                d = zeros(T,s)
                tk_rho = zeros(T, compdim, s)
                tk_j = zeros(T, compdim, s)
                d[1] = 1
                for k = 2:s
                    d[k] = epsilon^2 + dt*A_im[k,k]
                end
                for i in 1:s
                    tk_rho[:, i] = prod(d[1:i-1]) * u_sol_rho[:,it - 1]
                    if i>1
                        tk_j[:, i] = epsilon^2*prod(d[1:i-1]) * u_sol_j[:,it - 1]
                    else
                        tk_j[:, i] = prod(d[1:i-1]) * u_sol_j[:,it - 1]
                    end
                    for j in 1:i-1
                        tk_rho[:, i] += -dt*Drho*A_ex[i,j]*prod(d[j+1: i-1])*tk_j[:, j]
                        tk_j[:, i] += dt*(epsilon*A*A_ex[i,j]*prod(d[j+1: i-1])*tk_j[:, j] - Dg*A_im[i,j]*prod(d[j: i-1])*tk_rho[:, j]-A_im[i,j]*prod(d[j+1: i-1])*tk_j[:, j])
                    end
                    # in diesem Teil tritt kein Produkt über d auf (wäre sowieso prod(d[i, i-1]))
                    tk_j[:, i] += -dt*Dg*A_im[i,i]*tk_rho[:, i]
                end
                # u^{n+1} = u^{(s)} mit u^{(s)} = (1/prod(d[1:i-1])*tilde{rho^{(s)}}, 1/prod(d[1:i])*tilde{g^{(s)}})
                u_sol_rho[:, it] = tk_rho[:, s]/prod(d[1:s-1])
                u_sol_j[:, it] = tk_j[:, s]/prod(d[1:s])
            end
            u_sol[:, it] = vcat(u_sol_rho[:, it], u_sol_j[:, it])
        end
    elseif eq_type == "heat"
        for it in range(1,Nt, step = 1)
            if it == 1
                u_sol[:, it] = u0[:, it]
            else
                dt = (t_d[it]-t_d[it-1])
                u_sol[:, it] = u_sol[:, it - 1] + dt*RHS_mat*u_sol[:, it - 1]
            end
        end
    end
    return u_sol, u_exact
end


function ImEx(setup, RHS_mat, Tableau; only_explicit = false, optimized_ARS = false)
    eq_type = setup["eq_type"]
    T = setup["T"]
    # matrices
    A_ex = Tableau["A_ex"]
    A_im = Tableau["A_im"]
    b_ex = Tableau["b_ex"]
    b_im = Tableau["b_im"]
    s = length(b_im)
    if optimized_ARS == true && eq_type in ["telegraph", "telegraph_symm", "wave", "wave_symm", "telegraph_semistab"] && A_im[:, 1] == zeros(T, length(A_im[:, 1]))
        return ARS_optimized(setup, RHS_mat, Tableau)
    end
    # RHS_mat not used, because of splitting
    # This implementation is pseudo ImEx, as the kinetic model can be solved explicit
    ex_RHS_mat = setup["ex_RHS_mat"]
    im_RHS_mat = setup["im_RHS_mat"]
    # Instead of Nx, use the following notation (where dim = Nx for system, compdim = Nx for scalar equations)
        if setup["eq_type"] in ["heat", "transport"]
            dim = size(ex_RHS_mat/2)[1]
            compdim = Int(dim)
        else
            dim = size(ex_RHS_mat/2)[1]
            compdim = Int(dim/2)
        end
    # Checking, if the splitting is possible
    im_RHS_mat[1:compdim, :] == zeros(T, compdim, dim) || throw(ArgumentError("pseudo-ImEx splitting not possible: implicit component has unexpected entrys "))
    # General Initialisation
    x_d = setup["x_d"]
    a = setup["a"]
    u0 = setup["u0"]
    t_d = setup["t_d"]
    u_exact = get_exact_solution(setup)
    #u_exact = 0
    Nt = length(t_d)
    u_sol = zeros(T, dim, Nt)
    ### If we have a fully implicit time-discretization (just works for impleuler now)
    if A_ex == zeros(T, size(A_ex)) && b_ex == zeros(T, size(b_ex))
        dim = size(RHS_mat)[1]
        u_sol = zeros(T, dim, Nt)
        if A_im == [1]
            for it in range(1,Nt, step = 1)
                if it == 1
                    u_sol[:, it] = u0[:, it]
                else
                    dt = (t_d[it]-t_d[it-1])
                    u_sol[:, it] = (I-dt*RHS_mat)\u_sol[:,it - 1]
                end
            end
        elseif A_im == [1/2]
            for it in range(1,Nt, step = 1)
                if it == 1
                    u_sol[:, it] = u0[:, it]
                else
                    dt = (t_d[it]-t_d[it-1])
                    k = (I-dt/2*RHS_mat)\u_sol[:,it - 1]
                    u_sol[:, it] = u_sol[:,it - 1] + dt*RHS_mat*k
                end
            end
        elseif A_im == [0 0; 1/2 1/2]
            for it in range(1,Nt, step = 1)
                if it == 1
                    u_sol[:, it] = u0[:, it]
                else
                    dt = (t_d[it]-t_d[it-1])
                    u_sol[:, it] = (I-dt/2*RHS_mat)\(u_sol[:,it - 1]+dt/2*RHS_mat*u_sol[:,it - 1])
                end
            end
        end
        return u_sol, u_exact
    end

    if (eq_type in ["telegraph", "telegraph_symm", "wave", "wave_symm", "telegraph_semistab"]) && only_explicit == false
        # Pseudo-ImEx specialisation
        u_sol_rho = zeros(T, compdim, Nt)
        u_sol_j = zeros(T, compdim, Nt)
        ex_RHS_j_to_rho = ex_RHS_mat[1:compdim, compdim + 1:dim]
        ex_RHS_j_to_j = ex_RHS_mat[compdim + 1:dim, compdim + 1:dim]
        im_RHS = im_RHS_mat[compdim + 1:dim, :]
    
        # as we restrict ourselves to periodic boundary conditions, calling f_RHS is not needed
        for it in range(1,Nt, step = 1)
            if it == 1
                u_sol_rho[:,it] = u0[1:compdim, it]
                u_sol_j[:,it] = u0[compdim + 1:dim, it]
                u_sol[:, it] = vcat(u_sol_rho[:, it], u_sol_j[:, it])
            else
                dt = (t_d[it]-t_d[it-1])
                k_rho = zeros(T, compdim, s)
                k_j = zeros(T, compdim, s)
                for i in 1:s
                    k_rho[:, i] = u_sol_rho[:,it - 1]
                    k_j[:, i] = u_sol_j[:,it - 1]
                    for l = 1:i-1
                        k_rho[:, i] += dt*A_ex[i, l]*ex_RHS_j_to_rho*k_j[:, l]
                        k_j[:, i] += dt*A_ex[i, l]*ex_RHS_j_to_j*k_j[:, l] + dt*A_im[i, l]*im_RHS*vcat(k_rho[:, l], k_j[:, l])
                    end
                    k_j[:, i] += dt*A_im[i,i]*im_RHS[:, 1:compdim]*k_rho[:, i]
                    k_j[:, i] = (I-dt*A_im[i,i]*im_RHS[:, compdim + 1:dim])\k_j[:, i]
                end
                u_sol_rho[:, it] = k_rho[:, s]
                u_sol_j[:, it] = k_j[:, s]
                
                u_sol_rho[:, it] = u_sol_rho[:, it - 1]
                u_sol_j[:, it] = u_sol_j[:, it - 1]
                for i in 1:s
                    u_sol_rho[:, it] += dt*b_ex[i]*ex_RHS_j_to_rho*k_j[:, i]
                    u_sol_j[:, it] += dt*b_ex[i]*ex_RHS_j_to_j*k_j[:, i] + dt*b_im[i]*im_RHS*vcat(k_rho[:, i], k_j[:, i])
                end
                
            end
            u_sol[:, it] = vcat(u_sol_rho[:, it], u_sol_j[:, it])
        end
    elseif (eq_type in ["heat", "transport"]) || only_explicit == true
        for it in range(1,Nt, step = 1)
            if it == 1
                u_sol[:, it] = u0[:, it]
            else
                dt = (t_d[it]-t_d[it-1])
                k = zeros(T, dim, s)
                for i in 1:s
                    k[:, i] = u_sol[:,it - 1]
                    for l = 1:i-1
                        k[:, i] += dt*A_ex[i, l]*RHS_mat*k[:, l]
                    end
                end
                u_sol[:, it] = u_sol[:,it - 1]
                for i in 1:s
                    u_sol[:,it] += dt*b_ex[i]*RHS_mat*k[:, i]
                end
            end
        end
    end
    return u_sol, u_exact
end

function reset_array!(arr)
    arr .= zero(eltype(arr))
end

function put_in!(c, a, b)
    half_len = size(a,1)
    full_length = size(c,1)
    for i in 1:half_len
        c[i] = a[i]
    end
    for i in half_len+1:full_length
        c[i] = b[i]
    end
end

function ImEx_reformulated(setup, RHS_mat, Tableau; only_explicit = false)
    T = setup["T"]
    # matrices
    A_ex = Tableau["A_ex"]
    A_im = Tableau["A_im"]
    b_ex = Tableau["b_ex"]
    b_im = Tableau["b_im"]
    s = length(b_im)
    # RHS_mat not used, because of splitting
    # This implementation is pseudo ImEx, as the kinetic model can be solved explicit
    ex_RHS_mat = setup["ex_RHS_mat"]
    im_RHS_mat = setup["im_RHS_mat"]
    # Instead of Nx, use the following notation (where dim = Nx for system, compdim = Nx for scalar equations)
        if setup["eq_type"] in ["heat", "transport"]
            dim = size(ex_RHS_mat/2)[1]
            compdim = Int(dim)
        else
            dim = size(ex_RHS_mat/2)[1]
            compdim = Int(dim/2)
        end
    # Checking, if the splitting is possible
    im_RHS_mat[1:compdim, :] == zeros(T, compdim, dim) || throw(ArgumentError("pseudo-ImEx splitting not possible: implicit component has unexpected entrys "))
    # General Initialisation
    x_d = setup["x_d"]
    a = setup["a"]
    u0 = setup["u0"]
    t_d = setup["t_d"]
    eq_type = setup["eq_type"]
    u_exact = get_exact_solution(setup)
    #u_exact = 0
    Nt = length(t_d)
    u_sol = zeros(T, dim, Nt)

    k_rho = zeros(T, compdim, s)
    k_j = zeros(T, compdim, s)
    rho_subst = zeros(T, compdim, s)
    j_subst = zeros(T, compdim, s)

    ### If we have a fully implicit time-discretization (just works for impleuler now)
    if A_ex == zeros(T, size(A_ex)) && b_ex == zeros(T, size(b_ex))
        dim = size(RHS_mat)[1]
        u_sol = zeros(T, dim, Nt)
        if A_im == [1]
            for it in range(1,Nt, step = 1)
                if it == 1
                    u_sol[:, it] = u0[:, it]
                else
                    dt = (t_d[it]-t_d[it-1])
                    u_sol[:, it] = (I-dt*RHS_mat)\u_sol[:,it - 1]
                end
            end
        elseif A_im == [1/2]
            for it in range(1,Nt, step = 1)
                if it == 1
                    u_sol[:, it] = u0[:, it]
                else
                    dt = (t_d[it]-t_d[it-1])
                    k = (I-dt/2*RHS_mat)\u_sol[:,it - 1]
                    u_sol[:, it] = u_sol[:,it - 1] + dt*RHS_mat*k
                end
            end
        elseif A_im == [0 0; 1/2 1/2]
            for it in range(1,Nt, step = 1)
                if it == 1
                    u_sol[:, it] = u0[:, it]
                else
                    dt = (t_d[it]-t_d[it-1])
                    u_sol[:, it] = (I-dt/2*RHS_mat)\(u_sol[:,it - 1]+dt/2*RHS_mat*u_sol[:,it - 1])
                end
            end
        end
        return u_sol, u_exact
    end

    if (eq_type in ["telegraph", "telegraph_symm", "wave", "wave_symm", "telegraph_semistab"]) && only_explicit == false
        # Pseudo-ImEx specialisation
        u_sol_rho = zeros(T, compdim, Nt)
        u_sol_j = zeros(T, compdim, Nt)
        ex_RHS_j_to_rho = ex_RHS_mat[1:compdim, compdim + 1:dim]
        ex_RHS_j_to_j = ex_RHS_mat[compdim + 1:dim, compdim + 1:dim]
        im_RHS = im_RHS_mat[compdim + 1:dim, :]
        # the following if-else cant be used in general, but just ensures that besides of RK-schemes
        # with non-invertible A_im's, also the ARS-schemes work
        if A_im[:, 1] == zeros(T, length(A_im[:, 1]))
            d_im = hcat(0, b_im[2:s]'*inv(A_im[2:s, 2:s]))
        else
            d_im = b_im*inv(A_im)
        end
        # as we restrict ourselves to periodic boundary conditions, calling f_RHS is not needed
        for it in range(1,Nt, step = 1)
            if it == 1
                u_sol_rho[:,it] = u0[1:compdim, it]
                u_sol_j[:,it] = u0[compdim + 1:dim, it]
                u_sol[:, it] = vcat(u_sol_rho[:, it], u_sol_j[:, it])
            else
                dt = (t_d[it]-t_d[it-1])
                # k_rho = zeros(T, compdim, s)
                # k_j = zeros(T, compdim, s)
                # rho_subst = zeros(T, compdim, s)
                # j_subst = zeros(T, compdim, s)
                reset_array!.((k_rho, k_j, rho_subst, j_subst))
                @views for i in 1:s
                    for l = 1:i-1
                        @. k_rho[:, i] += A_ex[i, l]*ex_RHS_j_to_rho*(dt*k_j[:, l] + u_sol_j[:,it - 1])
                        @. k_j[:, i] += A_ex[i, l]*ex_RHS_j_to_j*(dt*k_j[:, l] + u_sol_j[:,it - 1]) + A_im[i, l]*im_RHS*vcat((dt*k_rho[:, l] + u_sol_rho[:,it - 1]), (dt*k_j[:, l] + u_sol_j[:,it - 1]))
                    end
                    @. k_j[:, i] += A_im[i,i]*im_RHS[:, 1:compdim]*(dt*k_rho[:, i] + u_sol_rho[:,it - 1]) + A_im[i,i]*im_RHS[:, compdim + 1:dim] * u_sol_j[:,it - 1]
                    @. k_j[:, i] = (I-dt*A_im[i,i]*im_RHS[:, compdim + 1:dim])\k_j[:, i]
                end
                @views for i = 1:s
                    for l = 1:i-1
                        @. rho_subst[:, i] += A_ex[i, l]*ex_RHS_j_to_rho*(dt*k_j[:, l] + u_sol_j[:,it - 1])
                        @. j_subst[:, i] += A_ex[i, l]*ex_RHS_j_to_j*(dt*k_j[:, l] + u_sol_j[:,it - 1])
                    end
                end

                @views for i in 1:s
                    #classic final update
                    #u_sol_rho[:, it] += dt*b_ex[i]*ex_RHS_j_to_rho*(dt*k_j[:, i] + u_sol_j[:, it - 1])
                    #u_sol_j[:, it] += dt*b_ex[i]*ex_RHS_j_to_j*(dt*k_j[:, i] + u_sol_j[:, it - 1]) + dt*b_im[i]*im_RHS*vcat((dt*k_rho[:, i] + u_sol_rho[:,it - 1]), (dt*k_j[:, i] + u_sol_j[:,it - 1]))
                    # first try changing final update
                    #u_sol_rho[:, it] += dt*b_ex[i]*ex_RHS_j_to_rho*(dt*k_j[:, i] + u_sol_j[:, it - 1]) - d_im[i]*rho_subst[:, i]
                    #u_sol_j[:, it] += dt*b_ex[i]*ex_RHS_j_to_j*(dt*k_j[:, i] + u_sol_j[:, it - 1]) + d_im[i]*im_RHS*vcat(dt*k_rho[:, i], dt*k_j[:, i]) - d_im[i]*j_subst[:, i]
                    # correcting the problematic term
                    @. u_sol_rho[:, it] += dt*b_ex[i]*ex_RHS_j_to_rho*(dt*k_j[:, i] + u_sol_j[:, it - 1]) + dt*d_im[i]*k_rho[:, i] - dt*d_im[i]*rho_subst[:, i]
                    @. u_sol_j[:, it] += dt*b_ex[i]*ex_RHS_j_to_j*(dt*k_j[:, i] + u_sol_j[:, it - 1]) + dt*d_im[i]*k_j[:, i] - dt*d_im[i]*j_subst[:, i]
                end
                @. u_sol_rho[:, it] += u_sol_rho[:, it - 1]
                @. u_sol_j[:, it] += u_sol_j[:, it - 1]
            end
            # @views u_sol[:, it] .= vcat(u_sol_rho[:, it], u_sol_j[:, it])
            @views put_in!(u_sol[:,it], u_sol_rho[:, it], u_sol_j[:, it])
        end
    elseif (eq_type in ["heat", "transport"]) || only_explicit == true
        for it in range(1,Nt, step = 1)
            if it == 1
                @views u_sol[:, it] = u0[:, it]
            else
                dt = (t_d[it]-t_d[it-1])
                k = zeros(T, dim, s) # ARPIT: PRE ALLOCATE
                @views for i in 1:s
                    @. k[:, i] = u_sol[:,it - 1]
                    for l = 1:i-1
                        @. k[:, i] += dt*A_ex[i, l]*RHS_mat*k[:, l]
                    end
                end
                @. u_sol[:, it] = u_sol[:,it - 1]
                for i in 1:s
                    @. u_sol[:,it] += dt*b_ex[i]*RHS_mat*k[:, i]
                end
            end
        end
    end
    return u_sol, u_exact
end


function get_RK_tableau(tableau_string)
    Tableau = Dict()
    if tableau_string == "ImExEuler" #FSAL, SA, GSA
        Tableau["A_ex"] = [0 0; 1 0]
        Tableau["A_im"] = [0 0; 0 1]
        Tableau["b_ex"] = [1 0]
        Tableau["b_im"] = [0 1]
    elseif tableau_string == "ARS2" #FSAL, SA, GSA
        gamma = 1-1/sqrt(2)
        delta = 1-1/(2*gamma)
        Tableau["A_ex"] = [0 0 0; gamma 0 0; delta 1-delta 0]
        Tableau["A_im"] = [0 0 0; 0 gamma 0; 0 1-gamma gamma]
        Tableau["b_ex"] = [delta 1-delta 0]
        Tableau["b_im"] = [0 1-gamma gamma]
    elseif tableau_string == "ARS3" #FSAL, SA, GSA
        Tableau["A_ex"] = [0 0 0 0 0;
        1/2 0 0 0 0;
        11/18 1/18 0 0 0;
        5/6 -5/6 1/2 0 0;
        1/4 7/4 3/4 -7/4 0]
        Tableau["A_im"] = [0 0 0 0 0;
                    0 1/2 0 0 0;
                    0 1/6 1/2 0 0;
                    0 -1/2 1/2 1/2 0;
                    0 3/2 -3/2 1/2 1/2]
        Tableau["b_ex"] = [1/4 7/4 3/4 -7/4 0]
        Tableau["b_im"] = [0, 3/2, -3/2, 1/2, 1/2]
    elseif tableau_string == "AGSA342" #FSAL, SA, GSA (order 2, Type I ImEx)
        Tableau["A_ex"] = [ 0 0 0 0;
                    -139833537/38613965 0 0 0;
                    85870407/49798258 -121251843/1756367063 0 0;
                    1/6 1/6 2*1/3 0]
        Tableau["A_im"] = [ 168999711/74248304 0 0 0;
                    44004295/24775207 202439144/118586105 0 0;
                    -6418119/169001713 -748951821/1043823139 12015439/183058594 0;
                    -370145222/355758315 1/3 0 202439144/118586105]
        Tableau["b_ex"] = [1/6 1/6 2/3 0]
        Tableau["b_im"] = [-370145222/355758315 1/3 0 202439144*1/118586105]
    elseif tableau_string == "IGSA2" #FSAL, SA, GSA (order 2, Type I ImEx)
        Tableau["A_ex"] = [ 0 0 0 0;
                    1/3 0 0 0;
                    7/24 3/8 0 0;
                    1/2 -1/2 1 0]
        Tableau["A_im"] = [ 1/4 0 0 0;
                    0 1/4 0 0;
                    1/16 3/16 1/4 0;
                    1/4 1/4 1/4 1/4]
        Tableau["b_ex"] = [1/2 -1/2 1 0]
        Tableau["b_im"] = [1/4 1/4 1/4 1/4]
    elseif tableau_string == "SSP2ImEx332" # not FSAL, SA, not GSA (order 2, Type I ImEx) 
        Tableau["A_im"] = [1/4 0 0;
                    0 1/4 0;
                    1/3 1/3 1/3]
        Tableau["b_im"] = [1/3 1/3 1/3]
        Tableau["A_ex"] = [0 0 0;
                    1/2 0 0;
                    1/2 1/2 0]
        Tableau["b_ex"] = [1/3 1/3 1/3]
    elseif tableau_string == "GSA1" # FSAL, SA, GSA (order 1, Type I ImEx) 
        gamma = 0.5 # free choice for gamma > 0
        a_ex = 0.5 # free choice for gamma > 0
        w1_ex = 0.5 # free choice for gamma > 0
        Tableau["A_im"] = [gamma 0 0;
                    0 gamma 0;
                    1-gamma 0 gamma]
        Tableau["b_im"] = [1-gamma 0 gamma]
        Tableau["A_ex"] = [0 0 0;
                    a_ex 0 0;
                    w1_ex 1-w1_ex 0]
        Tableau["b_ex"] = [w1_ex 1-w1_ex 0]
    elseif tableau_string == "SSP3ImEx343" # not FSAL, not SA, not GSA (order 3, Type I ImEx)
        α = 0.24169426078821
        β = 0.06042356519705
        η = 0.12915286960590
        Tableau["A_im"] = [ α 0 0 0;
                   -α α 0 0;
                    0 1-α α 0;
                   β η 1/2−β−η−α α]
        Tableau["b_im"] = [0 1/6 1/6 2/3]
    
        Tableau["A_ex"] = [0 0 0 0;
                      0 0 0 0;
                      0 1 0 0;
                      0 1/4 1/4 0]
        Tableau["b_ex"] = [0 1/6 1/6 2/3]
    elseif tableau_string == "BPR343" # No Type 1, but FSAL, SA
        Tableau["A_im"] = [0 0 0 0 0;
                    1/2 1/2 0 0 0;
                    5/18 -1/9 1/2 0 0;
                    1/2 0 0 1/2 0;
                    1/4 0 3*1/4 -1/2 1/2]
        Tableau["b_im"] = [1/4, 0, 3/4, -1/2, 1/2]
        Tableau["A_ex"] = [0 0 0 0 0;
                    1 0 0 0 0;
                    4/9 2/9 0 0 0;
                    1/4 0 3/4 0 0;
                    1/4 0 3/4 0 0]
        Tableau["b_ex"] = [1/4, 0, 3/4, 0, 0]
    elseif tableau_string == "ImExEulernotGSA"
        Tableau["A_ex"] = [0]
        Tableau["A_im"] = [1]
        Tableau["b_ex"] = [1]
        Tableau["b_im"] = [1]
    elseif tableau_string == "ARK324L2SA_ERK" # not FSAL, SA, not GSA (3rd order type II ImEx)
        Tableau["A_ex"] = [0 0 0 0;
                    (1767732205903)/2027836641118 0 0 0;
                    (5535828885825)/10492691773637 (788022342437)/10882634858940 0 0;
                    (6485989280629)/16251701735622 (-4246266847089)/9704473918619 (10755448449292)/10357097424841 0]
        Tableau["A_im"] =  [0 0 0 0;
                    (1767732205903)/4055673282236 (1767732205903)/4055673282236 0 0;
                    (2746238789719)/10658868560708 (-640167445237)/6845629431997 (1767732205903)/4055673282236 0;
                    (1471266399579)/7840856788654 (-4482444167858)/7529755066697 (11266239266428)/11593286722821 (1767732205903)/4055673282236]
        Tableau["b_ex"] = [(1471266399579)/7840856788654 (-4482444167858)/7529755066697 (11266239266428)/11593286722821 (1767732205903)/4055673282236]
        Tableau["b_im"] = [(1471266399579)/7840856788654 (-4482444167858)/7529755066697 (11266239266428)/11593286722821 (1767732205903)/4055673282236]
    elseif tableau_string == "ImplEuler"
        Tableau["A_ex"] = [0]
        Tableau["A_im"] = [1]
        Tableau["b_ex"] = [0]
        Tableau["b_im"] = [1]
    elseif tableau_string == "ImplMidpoint"
        Tableau["A_ex"] = [0]
        Tableau["A_im"] = [1/2]
        Tableau["b_ex"] = [0]
        Tableau["b_im"] = [1]
    elseif tableau_string == "CrankNicholson"
        Tableau["A_ex"] = [0 0; 0 0]
        Tableau["A_im"] = [0 0; 1/2 1/2]
        Tableau["b_ex"] = [0, 0]
        Tableau["b_im"] = [1/2, 1/2]
    end
    return Tableau
end