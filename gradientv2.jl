function gradient_3D_cent(arr, X)
    x = X[1]; y = X[2]; z = X[3];
    arr_r = real(arr);
    arr_i = imag(arr);

    dx = x[2] - x[1]; dy = y[2] - y[1]; dz = z[2] - z[1];
    dxarr_r = similar(arr_r);
    dyarr_r = similar(arr_r); 
    dzarr_r = similar(arr_r);
    dxarr_i = similar(arr_i);
    dyarr_i = similar(arr_i);
    dzarr_i = similar(arr_i); 

    # Periodic Boundary conditions
    dxarr_r[1, :, :] .= (arr_r[2, :, :] .- arr_r[end, :, :]) ./ (2*dx);
    dxarr_r[end, :, :] .= (arr_r[1, :, :] .- arr_r[end-1, :, :]) ./ (2*dx);
    dxarr_i[1, :, :] .= (arr_i[2, :, :] .- arr_i[end, :, :]) ./ (2*dx);
    dxarr_i[end, :, :] .= (arr_i[1, :, :] .- arr_i[end-1, :, :]) ./ (2*dx);

    for i in 2:length(x)-1
        dxarr_r[i,:,:] .= (arr_r[i+1, :, :] .- arr_r[i-1, :, :]) ./ (2*dx);
        dxarr_i[i,:,:] .= (arr_i[i+1, :, :] .- arr_i[i-1, :, :]) ./ (2*dx);
    end

    dyarr_r[:, 1, :] .= (arr_r[:, 2, :] .- arr_r[:, end, :]) ./ (2*dy);
    dyarr_r[:, end, :] .= (arr_r[:, 1, :] .- arr_r[:, end-1, :]) ./ (2*dy);
    dyarr_i[:, 1, :] .= (arr_i[:, 2, :] .- arr_i[:, end, :]) ./ (2*dy);
    dyarr_i[:, end, :] .= (arr_i[:, 1, :] .- arr_i[:, end-1, :]) ./ (2*dy);

    for i in 2:length(y)-1
        dyarr_r[:,i,:] .= (arr_r[:, i+1, :] .- arr_r[:, i-1, :]) ./ (2*dy);
        dyarr_i[:,i,:] .= (arr_i[:, i+1, :] .- arr_i[:, i-1, :]) ./ (2*dy);
    end

    dzarr_r[:, :, 1] .= (arr_r[:, :, 2] .- arr_r[:, :, end])./ (2*dz);
    dzarr_r[:, :, end] .= (arr_r[:, :, 1] .- arr_r[:, :, end-1])./ (2*dz);
    dzarr_i[:, :, 1] .= (arr_i[:, :, 2] .- arr_i[:, :, end])./ (2*dz);
    dzarr_i[:, :, end] .= (arr_i[:, :, 1] .- arr_i[:, :, end-1])./ (2*dz);

    for i in 2:length(z)-1
        dzarr_r .= (arr_r[:, :, i+1] - arr_r[:, :, i-1]) ./ (dz*2);
        dzarr_i .= (arr_i[:, :, i+1] - arr_i[:, :, i-1]) ./ (dz*2);
    end
    
    grad_r = [dxarr_r, dyarr_r, dzarr_r];
    grad_i = [dxarr_i, dyarr_i, dzarr_i];
    grad = [grad_r, grad_i];

    return grad;
end


# Testing wps at different z levels
zslice = 1;
psi = xspace(sol(50), sim);
grad = gradient_3D_cent(psi, X);

psi_torus = Torus(psi[:, :, zslice], x, y);
vorts = vortex_array(findvortices(psi_torus));
wps_arr = similar(real(psi[:, :, 1]))

for i in 1:length(x)
    for j in 1:length(y)
        wps_arr[i, j] = norm(wps_v(grad, [i,j, zslice]))
    end
end
wps_arr = abs.(wps_arr)
max = findmax(wps_arr)[1]
wps_arr = wps_arr ./max

using Plots
heatmap(x, y, wps_arr)
scatter!(vorts[:, 2], vorts[:, 1])


    


    
