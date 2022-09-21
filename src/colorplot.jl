using NPZ
using CairoMakie, ColorSchemes

datapath = joinpath(@__DIR__,"..","data/microgrid")
figpath = joinpath(@__DIR__,"..","figures")

variance_nl = npzread(joinpath(datapath,"variance_nl.npy"))
variance_lin = npzread(joinpath(datapath,"variance_lin.npy"))
variance_blk = npzread(joinpath(datapath,"variance_blk.npy"))
variance_peak = npzread(joinpath(datapath,"variance_peak.npy"))

norm_nl = variance_nl .|> sqrt |> transpose
norm_lin = variance_lin .|> sqrt |> transpose
norm_blk = variance_blk .|> sqrt |> transpose
norm_peak = variance_peak .|> sqrt |> transpose

v_max = norm_lin |> maximum
v_min = norm_blk |> minimum

fig = Figure(resolution = (620, 725), font="CMU Serif")
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
ax3 = Axis(fig[2, 1])
ax4 = Axis(fig[2, 2])

heatmap!(ax1,norm_nl,colormap=:inferno,colorrange=(v_min,v_max))
ax1.xlabel = "Output Node"
ax1.ylabel = "Input Node"
ax1.title = "Nonlinear System"
ax1.titlesize = 18

heatmap!(ax2,norm_lin,colormap=:inferno,colorrange=(v_min,v_max))
ax2.xlabel = "Output Node"
ax2.ylabel = "Input Node"
ax2.title = "Linear System"
ax2.titlesize = 18

heatmap!(ax3,norm_peak,colormap=:inferno,colorrange=(v_min,v_max))
ax3.xlabel = "Output Node"
ax3.ylabel = "Input Node"
ax3.title = "Peak Approximation"
ax3.titlesize = 18

heatmap!(ax4,norm_blk,colormap=:inferno,colorrange=(v_min,v_max))
ax4.xlabel = "Output Node"
ax4.ylabel = "Input Node"
ax4.title = "Bulk Mode Contribution"
ax4.titlesize = 18

fig[3,1:2] = Colorbar(fig; limits = (v_min,v_max), colormap=:inferno,
                    vertical=false, flipaxis = false, width = 400, size = 20,
                    label=L"$L_2$ norm of the frequency response [Hz]", labelsize=16)
fig
save(joinpath(figpath,"colorplot.pdf"), fig)