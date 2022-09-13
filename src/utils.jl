using VortexDistributions, GLMakie, Colors
function plot_iso(psi, X, visible=true, heal_2=false)
    density = abs2.(psi)
    pmax = maximum(density)
    density = density/pmax
    if heal_2
        scene = volume(X[1], X[2], X[3], density, algorithm = :iso, visible=visible, isovalue=0.65,isorange=0.075)
    else
        scene = volume(X[1], X[2], X[3], density, algorithm = :iso, visible=visible)
    end
    screen = display(scene)
    # resize!(screen, 2998, 1920)
end

function scatterVortsOnIso(vorts, markersize=200)
    vorts = vorts3DMatrix(vorts);
    
    scatter!(vorts[:, 1], vorts[:, 2], vorts[:, 3], color="red", markersize=markersize)
end

# function vorts3DMatrix(vorts)
#     vorts = vcat(vorts'...)
# end

function scatterClassifiedVortices(vortSets, vorts_3d, X, size, edges=false)
    colors = distinguishable_colors(length(vortSets),[RGB(0,0.35,0.25)],dropseed=true)
    v_matrix = vcat(vorts_3d'...)[:,1:3]'

    for i in 1:length(vortSets)
        vi = v_matrix[:, collect(vortSets[i])]
        if !edges
            vi = vi[:, [vortInBounds(vi[:, i], X) for i = 1:length(vi[1, :])]] # Filters vortices that aren't on the grid
        end
        scatter!(vi[1,:],vi[2,:],vi[3,:],markersize=size,color=colors[i])
    end
end

function vortInBounds(v, X)
    x = X[1]; y = X[2]; z = X[3];
    dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];
    if ((v[1] >= x[1]) && (v[1] <= x[end]) && 
        (v[2] >= y[1]) && (v[2] <= y[end]) && 
        (v[3] >= z[1]) && (v[3] <= z[end]))
        return true
    else 
        return false
    end
end

function plot_line(vort_sort, color, linewidth)
    vx = [vort_sort[i][1] for i in 1:length(vort_sort)]
    vy = [vort_sort[i][2] for i in 1:length(vort_sort)]
    vz = [vort_sort[i][3] for i in 1:length(vort_sort)]

    lines!(vx, vy, vz, linewidth = linewidth, color = color)
end

function periodicPlotting(vorts_sorted, X)
    x = X[1]; dx = x[2]-x[1];

    vorts_plotting = []

    for i in 1:length(vorts_sorted)
        vort_current = []
        vi = vorts_sorted[i]
        prev_break = 1
        for j in 1:length(vi)-1
            if euclid(vi[j], vi[j+1]) > dx
                push!(vort_current, vi[prev_break:j])
                prev_break = j+1
            end
        end
        push!(vort_current, vi[prev_break:end])
        if euclid(vort_current[1][1], vort_current[end][end]) < dx
            push!(vort_current[end], vort_current[1][1])
        end

        push!(vorts_plotting, vort_current)
    end
    
    colors = distinguishable_colors(length(vorts_plotting),[RGB(0,0.35,0.25)],dropseed=true);
    for i in 1:length(vorts_plotting)
        vi = vorts_plotting[i]
        try
            for j in 1:length(vi)
                plot_line(vi[j], colors[i], 5)
            end
        catch
        end
    end
end