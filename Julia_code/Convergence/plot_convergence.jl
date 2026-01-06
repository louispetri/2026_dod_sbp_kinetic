using CairoMakie, DelimitedFiles, LaTeXStrings, PolynomialBases

epsilons = [0.5, 0.1, 0.01]
degrees = 0:2
TMM = "ARS3"
TMM_str = TMM
eq_types = ["telegraph", "heat"]
basis_type = GaussLegendre
J1_types = ["upwind_diss_symm", "central_symm"]

for eq_type in eq_types
    for J1_type in J1_types
        if J1_type == "upwind_diss_symm"
            flux_pair_str = L"D^-, D^+"
        elseif J1_type == "central_symm"
            flux_pair_str = L"D^z, D^z"
        end
        for epsilon in epsilons
            B = readdlm(joinpath(@__DIR__,"./data/" * eq_type * "/" * TMM_str * "/" * "/$(basis_type)/" * J1_type * "/conv_error_epsilon=$(epsilon).txt"),'\t', Float64,'\n')
            if eq_type == "telegraph"
                fds = Int(length(B[:,1])/2)
                titlestring = L"$\epsilon = %$(epsilon)$, num. flux pair: %$(flux_pair_str)"
            elseif eq_type == "heat"
                fds = Int(length(B[:,1]))
                titlestring = L"num. flux pair: %$(flux_pair_str)"
            end
            # reference values
            r0 = B[1:fds,1].^-1
            r1 = B[1:fds,1].^-2
            r2 = B[1:fds,1].^-3

            f = Figure(fontsize = 32)
            ax = Axis(f[1,1], xlabel = L"N", ylabel = L"error$$", xlabelsize = 36, ylabelsize = 36,
                title = titlestring, titlesize = 36, xscale = log10, yscale = log10)
            scatterlines!(ax, B[1:fds,1], B[1:fds,2], label = "p=0", linewidth = 3.5, markersize = 15, color = Makie.wong_colors()[1])
                        lines!(ax, B[1:fds,1], r0, label = "ref. ord. 1", linestyle = :dash, linewidth = 3.5, color = Makie.wong_colors()[5])
            scatterlines!(ax, B[1:fds,1], B[1:fds,3], label = "p=1", linewidth = 3.5, markersize = 15, color = Makie.wong_colors()[2])
                        lines!(ax, B[1:fds,1], r1, label = "ref. ord. 2", linestyle = :dash, linewidth = 3.5, color = Makie.wong_colors()[6])
            scatterlines!(ax, B[1:fds,1], B[1:fds,4], label = "p=2", linewidth = 3.5, markersize = 15,  color = Makie.wong_colors()[3])
                        lines!(ax, B[1:fds,1], r2, label = "ref. ord. 3", linestyle = :dash, linewidth = 3.5, color = :green)

            #axislegend(ax, position = :rb, fontsize = 5)
            mkpath(joinpath(@__DIR__,"./plots"))
            save(joinpath(@__DIR__,"./plots/conv_$(eq_type)_$(J1_type)_eps=$(epsilon)_plot.pdf"), f)
            g = Figure(fontsize = 20)
            gleg = Legend(g,ax, framevisible = true, labelsize = 17)
            gleg.orientation = :horizontal
            gleg.nbanks = 2
            g[1,1] = gleg
            save(joinpath(@__DIR__,"./plots/legend_convergence.pdf"), g)
        end
    end
end