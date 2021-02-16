using Interpolations, LinearAlgebra

psi = xspace(psi2, sim);
zidx = 32;
vorts_3d = find_vortices3D(psi, X)
grad = gradient_3D_cent(psi, X);
grad_i = grad_itp(grad, X);
dz = z[2]-z[1];

label = 1;

vorts_3d = find_vortices3D_v2(psi, X)

for i in 1:length(vorts_3d[1][:, 1])
    if vorts_3d[1][i, 3] > 0
        vorts_3d[1][i, 4] = label;
        label += 1;
    end
end

for i in 1:length(vorts_3d[32][:, 1])
    if vorts_3d[32][i, 3] < 0
        vorts_3d[32][i, 4] = label;
        label += 1;
    end
end


for zidx in 1:length(z)-1
    vorts = vorts_3d[zidx] 
    vorts_up = vorts_3d[zidx + 1]
    for v in 1:length(vorts[:, 1])
        if vorts[v, 3] > 0
            if vorts[v, 4] == 0
                vorts[v, 4] = label;
                label+=1;
            end
            wps = wps_int(grad_i, [vorts[v, 1], vorts[v, 2], z[zidx]]);
            wps = wps .* (dz/wps[3]);
            approx_v = vorts[v, 1:2] + wps[1:2];
            vort_index = find_closest_tuple(vorts_up[:, 1:2], [approx_v[1], approx_v[2]])
            vorts_up[vort_index, 4] = vorts[v, 4];
            
        end
    end    
end

for zidx in length(z):-1:2
    vorts = vorts_3d[zidx] 
    vorts_down = vorts_3d[zidx-1]
    for v in 1:length(vorts[:, 1])
        if vorts[v, 3] < 0
            if vorts[v, 4] == 0
                vorts[v, 4] = label;
                label+=1;
            end
            wps = wps_int(grad_i, [vorts[v, 1], vorts[v, 2], z[zidx]]);
            wps = wps .* (dz/wps[3]);
            approx_v = vorts[v, 1:2] + wps[1:2];
            vort_index = find_closest_tuple(vorts_down[:, 1:2], [approx_v[1], approx_v[2]])
            vorts_down[vort_index, 4] = vorts[v, 4];
        end
    end
end



