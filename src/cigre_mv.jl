using Graphs
using PowerDynamics # Download package from branch: https://github.com/JuliaEnergy/PowerDynamics.jl/tree/FluctuationNode
using OrdinaryDiffEq
using HDF5
using Interpolations
using CairoMakie, GraphMakie, ColorSchemes, NetworkLayout

datapath = joinpath(@__DIR__,"..","data")
figpath = joinpath(@__DIR__,"..","figures")

## Set Per Unit MV

const base_power = 1E6 # 1MW
const base_voltage = 20E3 # 20kV
const base_current = base_power / base_voltage # 50A
const base_admittance = base_power / base_voltage^2 # 0.0025Ω^-1
const ω = 2 * π * 50.0 # 314.1593rad/s

## Load Network Structure and Parameters from Data

include(joinpath(datapath,"cigre_mv/cigre_data.jl"))
buses, lines, ldata = CIGRE_static(close_switches=true)

## power flow solution

buses[1] = SlackAlgebraic(U=1)
pg = PowerGrid(buses, lines)

data, result = power_flow(pg);
v = [result["solution"]["bus"][string(k)]["vm"] for k in 1:length(pg.nodes)];
va = [result["solution"]["bus"][string(k)]["va"] for k in 1:length(pg.nodes)];

## initialize model

τ_Q = 8.0
K_P = 10.0
K_Q = 0.01
V_r = 1.0
τ_P =0.5

buses = Array{PowerDynamics.AbstractNode,1}([])

for bus in 1:11
    P_Gen = result["solution"]["gen"]["$bus"]["pg"]
    Q_Gen = result["solution"]["gen"]["$bus"]["qg"]
    push!(
        buses,
        VSIMinimal(τ_P=τ_P,τ_Q=τ_Q,K_P=K_P,K_Q=K_Q,V_r=V_r,P=P_Gen,Q=Q_Gen),
        )
end

pg = PowerGrid(buses, lines);
ic_guess = PowerDynamics.initial_guess(pg, v .* exp.(1im .* va))
x0 = State(pg, ic_guess);

## fluctuation time series

function time_interpolation(series,Δt)
  itp = interpolate(series, BSpline(Linear()))
  etp = extrapolate(itp, Interpolations.Flat()) #extrapolate first and last value to all times
  t -> etp(t/Δt .+1)
end

series = h5open(joinpath(datapath,"fluctuations_time_series.hdf"),"r") do file
    read(file, "sum_of_fluctuations")
end

series .-= sum(series)/length(series)
Fluct = time_interpolation(series,0.01)

## simulate flucutation

t_idx = Float64[]
tspan = (0.0, 100.0)
Ng = 10

for bus in 1:11

    acp = buses[bus].P # active power
    rep = buses[bus].Q # reactive power

    # add fluctuation to the load (exchange PQAlgebraic with FluctuationNode)
    buses[bus] = FluctuationNode(t -> acp + Fluct(t), t -> rep)
    pg = PowerGrid(buses, lines)

    # remove frequency variable from initial condition
    var_idx = deleteat!(collect(1:33),3*bus)
    ic = x0.vec[var_idx]

    # integrate with OrdinaryDiffEq
    ode = ODEProblem(rhs(pg), ic, tspan)
    dqsol = solve(ode, Rodas4())
    sol = PowerGridSolution(dqsol, pg);
    println("Simulate fluctuation at ",bus, ": ", sol.dqsol.retcode)

    # compute L2-norm of frequency deviation
    Δt = 0.01; T = 50
    bus_idx = deleteat!(collect(1:11),bus)
    sol = sol(0.0:Δt:T,bus_idx,:ω)
    norm = sum(abs2,sol)*Δt / T / Ng # 1/N ∑(1/T ∫ ωᵢ² dt)
    push!(t_idx,norm)

    buses[bus] = VSIMinimal(τ_P=τ_P,τ_Q=τ_Q,K_P=K_P,K_Q=K_Q,V_r=V_r,P=acp,Q=rep)

end

## plot fluctuation (using GraphMakie)

function flow_plot(op,t_idx,min,max)

    g = op.grid.graph
    nodes = op.grid.nodes
    dg = SimpleDiGraph(g)
    arrow_size = []
    rem_list = []

    linenumber = [(e,i) for (i,e) in enumerate(edges(g))] |> Dict
    lineset = edges(g) |> collect |> Set
    dg_edges = collect(dg |> edges)
    lines = op.grid.lines 

    for i in 1:ne(dg)

        e = dg_edges[i]
        src_bus = e.src
        dst_bus = e.dst
        u_src = op[src_bus,:u]
        u_dst = op[dst_bus,:u]

        if e in lineset
            line = lines[linenumber[e]]
        else
            e_inv = Edge(e.dst,e.src)
            line = lines[linenumber[e_inv]]
        end

        if typeof(line) == StaticLine
            y = line.Y
        else
            y = line.y
        end

        i = y*(u_src - u_dst)
        b = y |> imag |> abs
        p = b*sin(angle(u_src)-angle(u_dst))

        if p < 0
            push!(rem_list,e)
        else
            push!(arrow_size,10*p)
        end
    end
    
    map(e->rem_edge!(dg,e),rem_list)

    node_color = []
    node_shape = []
    cm = ColorSchemes.inferno
    ti = (t_idx .- min); ti /= (max-min);
    colorval = get(cm,ti)
    node_idx = collect(1:11)
    colordict = Dict(node_idx .=> colorval)

    for (idx,bus) in enumerate(nodes)
        if bus.P <= 0
            push!(node_color,colordict[idx])
            push!(node_shape,:rect)
        else
            push!(node_color,colordict[idx])
            push!(node_shape,'●')
        end
    end

    nlabels = repr.(1:nv(g))

    return dg, nlabels, arrow_size, node_color, node_shape

end

l2 = sqrt.(t_idx)

l2_min = 0.166 #minimum(l2) #0.083
l2_max = 0.177 #maximum(l2) #0.089
dg, nlabels, arrow_size, node_color, node_shape = flow_plot(x0,l2,l2_min,l2_max);
fig = Figure(resolution = (600, 700), font="CMU Serif");
ax = Axis(fig[1,1]);
#add_edge!(dg,4,11); add_edge!(dg,6,7)
graphplot!(dg, node_color=node_color,node_size=15,node_marker=node_shape,
            layout=Stress(seed=2),arrow_size=13,nlabels=["$i" for i in 1:11],nlabels_offset=Point(-0.05,0.1))
hidedecorations!(ax)
load_marker = MarkerElement(color = :black, marker = '□', markersize = 15)
gen_marker = MarkerElement(color = :black, marker = '○', markersize = 15)
flow_marker = MarkerElement(color = :black, marker = '➛', markersize = 20)
axislegend(ax,[load_marker,gen_marker,flow_marker],["Net Consumer","Net Producer","Power Flow"],position = :rb)
fig[2,1] = Colorbar(fig; vertical=false, flipaxis = false,
                   limits = (l2_min,l2_max),
                   colormap=:inferno, label=L"$L_2$ norm of the frequency response [Hz]",labelsize=18)
current_figure()
save(joinpath(figpath,"cigre_mv_closed.pdf"), fig)