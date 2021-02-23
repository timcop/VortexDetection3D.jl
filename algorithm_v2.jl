using Interpolations, LinearAlgebra

psi = xspace(psi_ring4, sim);
grad = gradient_3D_cent(psi, X);
grad_i = grad_itp(grad, X);
dz = z[2]-z[1];



vorts_3d = find_vortices3D_v2(psi, X);

# Find first zslice containing vortices
# Vort = [x y z charge label]

label = 1;
for zidx in 1:length(z)
    vorts_slice = vorts_3d[zidx];
    num_vorts = length(vorts_slice) # Number of vorts on z-slice
    if num_vorts != 0
        num_vorts = length(vorts_slice[:, 1])
        for v in 1:num_vorts
            if (vorts_slice[v, 5] == 0) # Check hasn't been labeled
                vorts_3d[zidx][v, 5] = label;
                label += 1;
            end
            
            vorts_sliceup = vorts_3d[zidx+1];
            
            if length(vorts_sliceup) != 0
                wps = wps_int(grad_i, [vorts_slice[v, 1], vorts_slice[v, 2], z[zidx]]);
                wps = wps .* (dz/abs(wps[3]));
                
                if vorts_slice[v, 4] > 0
                    approx_vort = vorts_slice[v, 1:2] + wps[1:2]
                else 
                    approx_vort = vorts_slice[v, 1:2] - wps[1:2]
                end
                vort_index = find_closest_tuple(vorts_sliceup[:, 1:2], [approx_vort[1], approx_vort[2]]);
                if vorts_sliceup[vort_index, 4] == vorts_slice[v, 4]
                    vorts_sliceup[vort_index, 5] = vorts_slice[v, 5] # Assign same label
                end

            end

        end
    end
end
vorts_label = Array{Float64, 2}[];
for i in 1:label-1
    current_vort = [1000 0 0 0 0];

    for zidx in 1:length(z)

        if length(vorts_3d[zidx]) != 0
            for v in 1:length(vorts_3d[zidx][:, 1])

                if vorts_3d[zidx][v, 5] == i
                    if current_vort[1] == 1000
                        current_vort = vorts_3d[zidx][v, :]';
                    else
                        current_vort = vcat(current_vort, vorts_3d[zidx][v, :]');
                    end
                end
            end 
        end

    end
    push!(vorts_label, current_vort)
end
vorts_label[10]
NUM_VORTS = length(vorts_label)


# For each v in vorts_label[i]
# If v[1, 3] != z[1] then vort not attached to bottom slice
# If v[end, 3] != z[end] then vort not attached to top slice
# If
x[1]
x[end]
dx = x[2] - x[1];
dy = y[2] - y[1]; 

function vortex_boundary_bottom(v, dx, dy, dz)
    at_boundary = true;
    if ((v[1, 1] - x[1] > dx) && (x[end] - v[1, 1] > dx)) # Then not at x Boundary
        if ((v[1, 2] - y[1] > dy) && (y[end] - v[1, 2] > dy)) # " y boundary 
            if ((v[1, 3] - z[1] > dz) && (z[end] - v[1, 3] > dz))
                at_boundary = false;
            end
        end
    end
    return at_boundary;
end

function vortex_boundary_top(v, dx, dy, dz)
    at_boundary = true;
    if ((v[end, 1] - x[1] > dx) && (x[end] - v[end, 1] > dx)) # Then not at x Boundary
        if ((v[end, 2] - y[1] > dy) && (y[end] - v[end, 2] > dy)) # " y boundary 
            if ((v[end, 3] - z[1] > dz) && (z[end] - v[end, 3] > dz))
                at_boundary = false;
            end
        end
    end
    return at_boundary;
end

function vortex_link_bottom(vidx, vorts_label)
    v = vorts_label[vidx];
    for i in 1:length(vorts_label)
        vi = vorts_label[i];
        if vi[1,5] != v[1,5] # different label
            # Find matching z level
            index = 0;
            for j in 1:length(vi[:, 1])
                if vi[j, 3] == v[1, 3]
                    index = j;
                end
            end
            # If matching z level exists
            if index != 0
                # If z level is at bottom of vi 
                if vi[index, 3] == vi[1, 3]
                    vorts_label[vidx][:, 5] .= vi[1, 5];
                    if !vortex_boundary_top(vi, dx, dy, dz)
                        vortex_link_top(i, vorts_label);
                    end
                elseif vi[index, 3] == vi[end, 3]
                    vorts_label[vidx][:, 5] .= vi[1, 5];
                    if !vortex_boundary_bottom(vi, dx, dy, dz)
                        vortex_link_bottom(i, vorts_label);
                    end
                end
            end
        end
    end
end

function vortex_link_top(vidx, vorts_label)
    v = vorts_label[vidx];
    for i in 1:length(vorts_label)
        vi = vorts_label[i];
        if vi[1,5] != v[1,5] # different label
            # Find matching z level
            index = 0;
            for j in 1:length(vi[:, 1])
                if vi[j, 3] == v[end, 3]
                    index = j;
                end
            end
            # If matching z level exists
            if index != 0
                # If z level is at bottom of vi 
                if vi[index, 3] == vi[1, 3]
                    vorts_label[vidx][:, 5] .= vi[1, 5];
                    if !vortex_boundary_top(vi, dx, dy, dz)
                        vortex_link_top(i, vorts_label);
                    end
                elseif vi[index, 3] == vi[end, 3]
                    vorts_label[vidx][:, 5] .= vi[1, 5];
                    if !vortex_boundary_bottom(vi, dx, dy, dz)
                        vortex_link_bottom(i, vorts_label);
                    end
                end
            end
        end
    end
end

function sort_vorts_label(vorts_label)
    
    temp = Array{Float64, 2}[]
    i = 1;
    for i in 1:length(vorts_label)
        for j in 1:NUM_VORTS
        
            if (vorts_label[i][1, 5] == j)
               #println(j)
                push!(temp, vorts_label[i]);
            end
        end
    end
    return temp;
end



vort_linked = Array{Float64, 2}[];
while (length(vorts_label) > 0)
    
    v = vorts_label[1]
    if !vortex_boundary_bottom(v, dx, dy, dz)
        vortex_link_bottom(1, vorts_label)
    end
    if !vortex_boundary_top(v, dx, dy, dz)
        vortex_link_top(1, vorts_label)
    end


    vorts_label = sort_vorts_label(vorts_label);
    reverse!(vorts_label);
    label = vorts_label[end][1, 5];
    while ((length(vorts_label) > 0) && (vorts_label[end][1, 5] == label))
        pop = pop!(vorts_label);
        push!(vort_linked, pop);
    end
    reverse!(vorts_label)
    println(length(vorts_label))
    # Pop all v in vorts_label with same label as v
end

vorts_label
vort_linked[6]








using Makie, AbstractPlotting
volume(dense(psi_ring4), algorithm = :iso, show_axis = true)


