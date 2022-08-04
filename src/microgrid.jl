using NPZ
using Graphs
using CairoMakie, GraphMakie, ColorSchemes, NetworkLayout

function flow_plot(g,phases,power,t_idx,min,max)

    dg = SimpleDiGraph(g)
    rem_list = []
    dg_edges = collect(dg |> edges)

    for i in 1:ne(dg)
        e = dg_edges[i]
        if phases[e.src] < phases[e.dst]
            push!(rem_list,e)
        end
    end
    
    map(e->rem_edge!(dg,e),rem_list)

    node_color = []
    node_shape = []
    cm = ColorSchemes.inferno
    ti = (t_idx .- min); ti /= (max-min);
    colorval = get(cm,ti)
    nodes = dg |> vertices |> collect
    colordict = Dict(nodes .=> colorval)

    for i in 1:nv(dg)

        if power[i] > 0
            push!(node_color,colordict[i])
            push!(node_shape,:circle)
        else
            push!(node_color,colordict[i])
            push!(node_shape,'■')
        end
    end

    nlabels = repr.(1:nv(g))

    return dg, nlabels, node_color, node_shape

end

datapath = joinpath(@__DIR__,"..","data/microgrid")
figpath = joinpath(@__DIR__,"..","figures")

phases = npzread(joinpath(datapath,"phases.npy"))
power = npzread(joinpath(datapath,"power.npy"))
variance = npzread(joinpath(datapath,"variance_nl.npy"))
A = npzread(joinpath(datapath,"adjacency.npy"))
lat = npzread(joinpath(datapath,"lat.npy"))
lon = npzread(joinpath(datapath,"lon.npy"))

g = Graph(A)
pos = [Point2(lat[n],lon[n]) for n in 1:nv(g)]
layout = _ -> pos


t_idx = sum(variance,dims=2) .|> sqrt |> vec
idx_min = minimum(t_idx)
idx_max = maximum(t_idx)

dg, nlabels, node_color, node_shape = flow_plot(g,phases,power,t_idx,idx_min,idx_max)

fig = Figure(resolution = (600, 700), font="CMU Serif");
ax = Axis(fig[1,1:2]);

graphplot!(dg, node_color=node_color,node_size=15,node_marker=node_shape,arrow_size=13,
            layout=layout)
hidedecorations!(ax)
load_marker = MarkerElement(color = :black, marker = '□', markersize = 15)
gen_marker = MarkerElement(color = :black, marker = '○', markersize = 15)
flow_marker = MarkerElement(color = :black, marker = '➛', markersize = 20)
fig[2,1] = Colorbar(fig; vertical=false, flipaxis = false,
                   limits = (idx_min,idx_max),
                   colormap=:inferno, label=L"$L_2$ norm of the frequency response [Hz]",
                   labelsize=18,tellheight=true)
Legend(fig[2,2],[load_marker,gen_marker,flow_marker],["Net Consumer","Net Producer","Power Flow"],tellheight=false,tellwidth=true)#,orientation=:horizontal)

current_figure()
save(joinpath(figpath,"microgrid.pdf"), fig)