using CSV
using DataFrames
using Graphs
using PowerDynamics # Download package from branch: https://github.com/JuliaEnergy/PowerDynamics.jl/tree/FluctuationNode
using OrdinaryDiffEq
using HDF5
using Interpolations
using CairoMakie, GraphMakie, ColorSchemes, NetworkLayout

datapath = joinpath(@__DIR__,"..","data")
figpath = joinpath(@__DIR__,"..","figures")

## Load Network Structure and Parameters from Data

BusData = CSV.read(joinpath(datapath,"rts96/Bus.csv"), DataFrame; types = [Int, Int, String, Int, Int])
LineData = CSV.read(joinpath(datapath,"rts96/Line.csv"),DataFrame)
    LineData[!, :source] = LineData.source .|> Int
    LineData[!, :dest] = LineData.dest .|> Int
GeneratorData = CSV.read(joinpath(datapath,"rts96/Generator.csv"),DataFrame)
LoadData = CSV.read(joinpath(datapath,"rts96/Load.csv"),DataFrame)
FlowData = CSV.read(joinpath(datapath,"rts96/Flow.csv"),DataFrame)

const N = nrow(BusData)
const L = nrow(LineData)
const G = nrow(GeneratorData)


g = SimpleGraph(N)
    for i = 1:L
        add_edge!(g, Int64(LineData[i, :source]), Int64(LineData[i, :dest]))
    end

node_df = outerjoin(BusData, GeneratorData, LoadData, FlowData, on=:ID, makeunique=true)

slack_idx = argmax(skipmissing(node_df.P_Gen))

## construct power grid

nodes = []
    for n in eachrow(node_df)
        if n.Number == slack_idx
            # in the data set, n.Vm is not exactly 1,
            # so our steady state will be slightly different
            push!(nodes, SlackAlgebraic(; U=complex(1.)) ) # n.Vm
        else
            if n.P_Gen |> ismissing
                push!(nodes, PQAlgebraic(; P=-n.P_Load, Q=0.0) )
            else
                push!(nodes, SwingEqLVS(; H=n.Inertia, P=n.P_Gen-n.P_Load, D=0.01, Ω=100π, Γ=10., V=1.) )
            end
        end
    end

lines = []
    for l in eachrow(LineData)
        if node_df[l.source, :Base_V] == node_df[l.dest, :Base_V] #normal line
            push!(lines, StaticLine(; from=l.source, to=l.dest, Y=inv(complex(l.r, l.x))))
        else
            t_ratio = 1.0 # the tap is already included via the p.u. system
            push!(lines, StaticLine(; from=l.source, to=l.dest, Y=inv(complex(l.r, l.x))))
        end
    end

pg = PowerGrid(g, nodes, lines);

## power flow solution

data, result = power_flow(pg);
v = [result["solution"]["bus"][string(k)]["vm"] for k in 1:length(pg.nodes)];
va = [result["solution"]["bus"][string(k)]["va"] for k in 1:length(pg.nodes)];

## exchange slack with swing equation

#=
We exchange the slack bus by a dynamic generator model with the same power dispatch.
This is gives a more realistic behavior of the dynamics, since a slack node would just
swallow any power imbalances instantaniously.
=#

Hs = node_df[slack_idx,:].Inertia
gen_idx = [gen for (gen,val) in data["gen"] if val["gen_bus"] == slack_idx][1]
P_Gen = result["solution"]["gen"][string(gen_idx)]["pg"]

SwEq = SwingEqLVS(; H=Hs, P=P_Gen, D=0.01, Ω=100π, Γ=10., V=1.)
nodes[slack_idx] = SwEq
pg = PowerGrid(g, nodes, lines);

## construct initial condition

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

loads = [idx for (idx,val) in enumerate(pg.nodes) if typeof(val) == PQAlgebraic];
gens = [idx for (idx,val) in enumerate(pg.nodes) if typeof(val) == SwingEqLVS];

t_idx = Float64[]

tspan = (0.0, 100.0)

#= 
We simulate single node fluctuations for every load node in the grid.
=#

for bus in loads

    acp = nodes[bus].P # active power
    rep = nodes[bus].Q # reactive power

    # add fluctuation to the load (exchange PQAlgebraic with FluctuationNode)
    nodes[bus] = FluctuationNode(t -> acp + Fluct(t), t -> rep)
    pg = PowerGrid(nodes, lines)

    timespan = (-5.,50.)
    ode = ODEProblem(rhs(pg), x0.vec, tspan)
    dqsol = solve(ode, Rodas4())
    sol = PowerGridSolution(dqsol, pg);
    println("Simulate fluctuation at ",bus, ": ", sol.dqsol.retcode)

    Δt = 0.01; T = 50

    Ng= length(gens)
    sol = sol(0.0:Δt:T,gens,:ω)

    norm = sum(abs2,sol)*Δt / T / Ng # 1/N ∑(1/T ∫ ωᵢ² dt)
    push!(t_idx,norm)

    nodes[bus] = PQAlgebraic(P=acp,Q=rep)

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

    loads = [idx for (idx,val) in enumerate(pg.nodes) if typeof(val) == PQAlgebraic];

    node_color = []
    node_shape = []
    cm = ColorSchemes.inferno
    ti = (t_idx .- min); ti /= (max-min);
    colorval = get(cm,ti)
    colordict = Dict(loads .=> colorval)

    for (idx,bus) in enumerate(nodes)

        if typeof(bus) == PQAlgebraic
            push!(node_color,colordict[idx])
            push!(node_shape,:rect)
        elseif typeof(bus) == FourthOrderEq
            push!(node_color,:black)
            push!(node_shape,'○')
        else
            push!(node_color,:black)
            push!(node_shape,'○')
        end
    end

    nlabels = repr.(1:nv(g))

    return dg, nlabels, arrow_size, node_color, node_shape

end

l2 = sqrt.(t_idx)

l2_min = minimum(l2)
l2_max = maximum(l2)
dg, nlabels, arrow_size, node_color, node_shape = flow_plot(x0,l2,l2_min,l2_max);
fig = Figure(resolution = (600, 700), font="CMU Serif");
ax = Axis(fig[1,1]);
graphplot!(dg, node_color=node_color,node_size=15,node_marker=node_shape,
            layout=Stress(),arrow_size=13)
hidedecorations!(ax)
load_marker = MarkerElement(color = :black, marker = '□', markersize = 15)
gen_marker = MarkerElement(color = :black, marker = '○', markersize = 15)
flow_marker = MarkerElement(color = :black, marker = '➛', markersize = 20)
axislegend(ax,[load_marker,gen_marker,flow_marker],["Load Bus","Generator Bus","Power Flow"],position = :rb)
fig[2,1] = Colorbar(fig; vertical=false, flipaxis = false,
                   limits = (l2_min,l2_max),
                   colormap=:inferno, label=L"$L_2$ norm of the frequency response [Hz]",labelsize=18)
current_figure()
save(joinpath(figpath,"rts96.pdf"), fig)