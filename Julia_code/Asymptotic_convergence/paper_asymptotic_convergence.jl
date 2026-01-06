using LinearAlgebra, CairoMakie, PolynomialBases, DelimitedFiles, LaTeXStrings
include("../functions.jl")
include("../functionsRK.jl")
include("../flux_operators.jl")
include("../analysis_tools.jl")


TMMs = ["ARS3", "SSP2ImEx332"]
linestyles = [:solid, :dash, :dashdotdot]

include_b_h = true
epsilon_refinement = 38
T = Float64

degs = [0, 1, 2]
CFL_prefacs = [0.5, 0.15, 0.0375]
cs = [1.0, 0.55, 0.45]

flux_pairs = [("altlr", "upwind_diss_symm"), ("central", "central_symm")]


for (fluxtype, J1_type) in flux_pairs
    if J1_type == "upwind_diss_symm"
        flux_pair_str = L"D^-, D^+"
    elseif J1_type == "central_symm"
        flux_pair_str = L"D^z, D^z"
    end
    titlestring = L"num. flux pair: %$(flux_pair_str)"
    f = Figure(fontsize = 32)
    ax = Axis(f[1,1], xlabel = L"\varepsilon", ylabel = L"$L^2$-norm of $u_\varepsilon-u_0$", xlabelsize = 36, ylabelsize = 36,
            title = titlestring, titlesize = 21, xscale = log10, yscale = log10)
    icol = 1
    for (iTMM, TMM) in enumerate(TMMs)
        data = zeros(T, epsilon_refinement, length(degs)+1)
        for deg in degs
            for N in [2^4]
                steps, errors = asymptotic_error(;
                                T = T, epsilon_refinement = epsilon_refinement, a = 1,
                                N = N, Tmax = 0.5, deg = deg, CFL_prefac = CFL_prefacs[deg+1], basis_type = GaussLegendre,
                                fluxtype = fluxtype, J1_type = J1_type, J1_type_A = "upwind_diss_symm",
                                cut_cells = [10^-7, 10^-3, 0.1, 0.3, 0.49], fix_eta = true, c = cs[deg+1], include_b_h = include_b_h,
                                TMM = TMM,
                                CFL_type = "1/dx^2",
                                )
                if TMM == "SSP2ImEx332"
                    TMM_label_str = "IMEX SSP2"
                else
                    TMM_label_str = "$(TMM)"
                end
                lines!(ax, steps, errors, label = "$(TMM_label_str), p=$(deg)", linestyle = linestyles[iTMM], linewidth = 3.5, color = Makie.wong_colors()[icol])
                icol += 1
                data[:, 1] = steps 
                data[:, deg + 2] = errors
            end
        end
        mkpath(joinpath(@__DIR__,"./paper/data/$(T)_optim"))
        writedlm(joinpath(@__DIR__,"./paper/data/$(T)_optim/asymp_conv_$(fluxtype)_$(TMM).txt"), data)
    end
    mkpath(joinpath(@__DIR__,"./paper/$(T)_optim"))
    save(joinpath(@__DIR__,"./paper/$(T)_optim/asymp_conv_$(fluxtype).pdf"), f)
    g = Figure(fontsize = 20)
    gleg = Legend(g,ax, framevisible = true, labelsize = 17)
    gleg.orientation = :horizontal
    gleg.nbanks = 2
    g[1,1] = gleg
    save(joinpath(@__DIR__,"./paper/$(T)_optim/asymp_conv_legend.pdf"), g)
end
