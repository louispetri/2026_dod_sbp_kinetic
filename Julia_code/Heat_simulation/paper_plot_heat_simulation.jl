using LinearAlgebra, CairoMakie, PolynomialBases, DelimitedFiles, LaTeXStrings, Colors


id_set = ("heat", "telsin" , -pi, pi)
# irrelevant, but has to be defined
epsilon = 0.5
a = 1.0


N = 2^5
Tmax = 5.0

c_upwind_diss_symm = [1.0, 0.511, 0.337, 0.247, 0.183, 0.148]
c_central_symm = [1.0, 0.563, 0.365, 0.264, 0.191, 0.150]
linewidths = [3.25, 2.5, 2.5]
cols = [:red2, :cyan2, :gold]

function lines_discont!(problem, axis, x_d, sol; label = NaN, linestyle = NaN, linewidth = NaN, color = NaN, titlestring = NaN,
                        plot_type = "discont_solid")
    deg = Int(problem["deg"])
    v = problem["v"]
    nodes = GaussLegendre
    cellnumber = Int(problem["cellnumber"])
    basis = nodes(deg) # same basis as in simulation
    degf = deg + 1 # degrees of freedom per cell
    x_vec = x_d
    sol_vec = sol
    if nodes == GaussLegendre
        int_shift = 0 # dazu da, um index korrekt nach hinten zu verschieben
        for i in 1:cellnumber
            insert_left = degf*(i-1) +1 + int_shift:degf*(i-1) + int_shift
            splice!(x_vec, insert_left, v[i])
            splice!(sol_vec, insert_left, interpolation_matrix(-1, basis)*sol_vec[degf*(i-1)+1 + int_shift:degf*i + int_shift])
            int_shift += 1
            insert_right = degf*(i) +1 + int_shift:degf*(i) + int_shift
            splice!(x_vec, insert_right, v[i+1])
            splice!(sol_vec, insert_right, interpolation_matrix(1, basis)*sol_vec[degf*(i-1)+1 + int_shift:degf*i + int_shift])
            int_shift += 1
            if !(plot_type in ["discont_solid_connected", "discont_dashed_connected"])
                insert_right = degf*(i) +1 + int_shift:degf*(i) + int_shift
                splice!(x_vec, insert_right, v[i+1])
                splice!(sol_vec, insert_right, NaN)
                int_shift += 1
            end
        end
    else
        throw(ArgumentError("Basistyp nicht implementiert"))
    end
    lines!(axis, x_vec, sol_vec, label = label, linestyle = linestyle, linewidth = linewidth, color = color)
end

for plot_type in ["discont_solid_connected"]
    mkpath(joinpath(@__DIR__,"./heat_plots/N=$(N)/$(plot_type)"))
    for do_stabilize in [false,true]
        for include_cut_cells in [false,true]
            if include_cut_cells == true
                if do_stabilize == false
                    titlestring = "Unstabilized"
                    ysc = identity
                    settingstr = "unstabil"
                else
                    titlestring = "    DoD (stabilized)"
                    ysc = identity
                    settingstr = "stabil"
                end
            else
                titlestring = "    Background"
                ysc = identity
                settingstr = "backgr"
            end
            marker_legend_control = true
            f = Figure(fontsize = 32)
            ax = Axis(f[1,1], xlabel = L"x", ylabel = L"\rho", xlabelsize = 36, ylabelsize = 36,
                title = titlestring, titlesize = 36, xscale = identity, yscale = ysc)
            xlims!(ax, [-pi, pi])
            for (fluxtype, J1_type, cs) in [("altlr", "upwind_diss_symm", c_upwind_diss_symm)]
                for (ideg, deg) in enumerate([0, 1, 2])
                    if J1_type == "upwind_diss_symm"
                        flux_pair_str = L"(D^-, D^+),\; p=%$(deg)"
                        if plot_type in ["discont_solid", "discont_solid_connected"]
                            linestyle = [:solid, :solid, :dash]
                        else
                            linestyle = [:solid, :dashdot, :dash]
                        end
                    else
                        throw(ArgumentError("J1-Fluss nicht implementiert!"))
                    end
                    CFL = 0.1/(2*deg+1)
                    J1_type_A = "upwind_diss_symm"
                    c = cs[deg+1]

                    metadata = readdlm(joinpath(@__DIR__,"./data/N=$(N)/$(settingstr)_metadata_deg$(deg).txt"),'\t', Float64,'\n')
                    data = readdlm(joinpath(@__DIR__,"./data/N=$(N)/$(settingstr)_deg$(deg).txt"),'\t', Float64,'\n')
                    x_d = data[:, 1]
                    sol = data[:, 2]
                    u_exact = data[:, 3]
                    problem = Dict("deg" => metadata[1], "cellnumber" => metadata[2], "v" => metadata[3:end])
                    v = problem["v"]
                    if include_cut_cells == true
                        if do_stabilize == true
                            ypos = 0.00
                        else
                            ypos = - 0.1
                        end
                        alpha = 1
                        for cut_pos in [3, 7, 10, 14, 20, 26]
                            markersize = 13
                            if marker_legend_control == true
                                marker_legend_control = false
                                scatter!(ax, (v[cut_pos]+v[cut_pos+1])/2, ypos; marker = :circle, markersize = markersize,
                                    color = RGBA(1.0, 0.75, 0.80, alpha), strokecolor = (:pink, 0.8), strokewidth = 2, label = "cut cell position")
                            else
                                scatter!(ax, (v[cut_pos]+v[cut_pos+1])/2, ypos; marker = :circle, markersize = markersize,
                                    color = RGBA(1.0, 0.75, 0.80, alpha), strokecolor = (:pink, 0.8), strokewidth = 2)
                            end
                        end
                    end
                    if plot_type == "cont_dashed"
                        lines!(ax, x_d, sol, label = flux_pair_str, linestyle = linestyle[ideg], linewidth = linewidths[ideg], color = cols[ideg])
                    elseif plot_type in ["discont_dashed", "discont_solid", "discont_solid_connected", "discont_dashed_connected"]
                        if include_cut_cells == true && do_stabilize == false
                            lines_discont!(problem, ax, x_d, sol, label = flux_pair_str, linestyle = linestyle[ideg], linewidth = linewidths[ideg], color = cols[ideg], titlestring = titlestring, plot_type = plot_type)
                        else
                            lines_discont!(problem, ax, x_d, sol * 10^3, label = flux_pair_str, linestyle = linestyle[ideg], linewidth = linewidths[ideg], color = cols[ideg], titlestring = titlestring, plot_type = plot_type)
                            Label(f[1, 1, Top()], halign = :left, L"\times 10^{-3}", fontsize = 33, height = 45)
                        end
                    end
                end
            end
            #axislegend(ax, position = :rt, fontsize = 5)
            save(joinpath(@__DIR__,"./heat_plots/N=$(N)/$(plot_type)/heat_cutcells=$(include_cut_cells)_stabilize=$(do_stabilize).pdf"), f)
            
            g = Figure(fontsize = 20)
            gleg = Legend(g,ax, framevisible = true, labelsize = 17)
            gleg.orientation = :horizontal
            gleg.nbanks = 2
            g[1,1] = gleg
            save(joinpath(@__DIR__,"./heat_plots/N=$(N)/$(plot_type)/heat_implmidpoint_legend.pdf"), g)
        end
    end
end


#############################################################################################
##############        Discontinuous plotting (requires CairoMakie.jl)        ################
#############################################################################################

function lines_discont!(problem, axis, x_d, sol; label = NaN, linestyle = NaN, linewidth = NaN, color = NaN, titlestring = NaN,
                        plot_type = "discont_solid")
    deg = Int(problem["deg"])
    v = problem["v"]
    nodes = GaussLegendre
    cellnumber = Int(problem["cellnumber"])
    basis = nodes(deg) # same basis as in simulation
    degf = deg + 1 # degrees of freedom per cell
    x_vec = x_d
    sol_vec = sol
    if nodes == GaussLegendre
        int_shift = 0 # dazu da, um index korrekt nach hinten zu verschieben
        for i in 1:cellnumber
            insert_left = degf*(i-1) +1 + int_shift:degf*(i-1) + int_shift
            splice!(x_vec, insert_left, v[i])
            splice!(sol_vec, insert_left, interpolation_matrix(-1, basis)*sol_vec[degf*(i-1)+1 + int_shift:degf*i + int_shift])
            int_shift += 1
            insert_right = degf*(i) +1 + int_shift:degf*(i) + int_shift
            splice!(x_vec, insert_right, v[i+1])
            splice!(sol_vec, insert_right, interpolation_matrix(1, basis)*sol_vec[degf*(i-1)+1 + int_shift:degf*i + int_shift])
            int_shift += 1
            if !(plot_type in ["discont_solid_connected", "discont_dashed_connected"])
                insert_right = degf*(i) +1 + int_shift:degf*(i) + int_shift
                splice!(x_vec, insert_right, v[i+1])
                splice!(sol_vec, insert_right, NaN)
                int_shift += 1
            end
        end
    else
        throw(ArgumentError("Basistyp nicht implementiert"))
    end
    lines!(axis, x_vec, sol_vec, label = label, linestyle = linestyle, linewidth = linewidth, color = color)
end