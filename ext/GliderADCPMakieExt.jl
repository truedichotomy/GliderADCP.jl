module GliderADCPMakieExt

using GliderADCP
using Makie
using Dates

function GliderADCP.plot_sections(panels::AbstractVector;
                                  colorrange=(-0.5, 0.5), colormap=:balance,
                                  figsize=(1500, 60 + 300 * length(panels)))
    fig = Makie.Figure(size=figsize)
    local hm
    n = length(panels)
    for (row, (sec, field, title)) in enumerate(panels)
        A = getproperty(sec, field)
        nyo = length(sec.t)
        xt = round.(Int, range(1, nyo; length=min(8, nyo)))
        ax = Makie.Axis(fig[row, 1]; ylabel="depth (m)", yreversed=true, title=String(title),
            xticks=(xt, Dates.format.(sec.t[xt], "dd u")),
            xlabel=row == n ? "segment midpoints" : "")
        hm = Makie.heatmap!(ax, 1:nyo, sec.z, permutedims(A); colormap, colorrange)
    end
    Makie.Colorbar(fig[1:n, 2], hm; label="velocity (m/s)")
    return fig
end

end
