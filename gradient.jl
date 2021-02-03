#########################################################

# Compute the gradient of 2D array, returns grad = [real, imag]

function gradient_2D(arr, X)
    arr_real = real(arr)
    arr_imag = imag(arr)
    x = X[1]; y = X[2];
    dx = x[2] - x[1];
    dy = y[2] - y[1];
    dxarr_real = similar(arr_real);
    dyarr_real = similar(arr_real); 
    dxarr_imag = similar(arr_imag);
    dyarr_imag = similar(arr_imag);

    for xin in 2:length(x)
        dxarr_real[xin, :] .= (arr_real[xin, :] .- arr_real[xin-1, :]) ./ dx;
        dxarr_imag[xin, :] .= (arr_imag[xin, :] .- arr_imag[xin-1, :]) ./ dx;
    end
    dxarr_real[1, :] .= 0;
    dxarr_imag[1, :] .= 0;
    
    for yin in 2:length(y)
        dyarr_real[:, yin] .= (arr_real[:, yin] .- arr_real[:, yin-1]) ./ dy;
        dyarr_imag[:, yin] .= (arr_imag[:, yin] .- arr_imag[:, yin-1]) ./ dy;
    end
    dyarr_real[:, 1] .= 0;
    dyarr_imag[:, 1] .= 0;
    GRADarr_real = [dxarr_real, dyarr_real];
    GRADarr_imag = [dxarr_imag, dyarr_imag];
    grad = [GRADarr_real, GRADarr_imag];
    return grad;
end

# Compute the gradient of 3D array, returns grad = [real, imag]

function gradient_3D(arr, X)

    arr_real = real(arr)
    arr_imag = imag(arr)
    x = X[1]; y = X[2]; z = X[3];
    dx = x[2] - x[1];
    dy = y[2] - y[1];
    dz = z[2] - z[1];
  
    dxarr_real = similar(arr_real);
    dyarr_real = similar(arr_real); 
    dzarr_real = similar(arr_real);
    dxarr_imag = similar(arr_imag);
    dyarr_imag = similar(arr_imag);
    dzarr_imag = similar(arr_imag); 

    for xin in 2:length(x)
        dxarr_real[xin, :, :] .= (arr_real[xin, :, :] .- arr_real[xin-1, :, :]) ./ dx;
        dxarr_imag[xin, :, :] .= (arr_imag[xin, :, :] .- arr_imag[xin-1, :, :]) ./ dx;
    end
    dxarr_real[1, :, :] .= 0;
    dxarr_imag[1, :, :] .= 0;

    for yin in 2:length(y)
        dyarr_real[:, yin, :] .= (arr_real[:, yin, :] .- arr_real[:, yin-1, :]) ./ dy;
        dyarr_imag[:, yin, :] .= (arr_imag[:, yin, :] .- arr_imag[:, yin-1, :]) ./ dy;
    end
    dyarr_real[:, 1, :] .= 0;
    dyarr_imag[:, 1, :] .= 0;

    for zin in 2:length(z)
        dzarr_real[:, :, zin] .= (arr_real[:, :, zin] .- arr_real[:, :, zin - 1]) ./ dz;
        dzarr_imag[:, :, zin] .= (arr_imag[:, :, zin] .- arr_imag[:, :, zin - 1]) ./ dz;
    end
    dzarr_real[:, :, 1] .= 0;
    dzarr_imag[:, :, 1] .= 0;

    GRADarr_real = [dxarr_real, dyarr_real, dzarr_real];
    GRADarr_imag = [dxarr_imag, dyarr_imag, dzarr_imag];
    GRADpsi = [GRADarr_real, GRADarr_imag];
    return GRADpsi;
end


# Pseudo-vorticity 3D, needs X = [x, y, z] for the grid and array gradient as parameters.
# Returns pseudo-vorticity vector 
function wps_v(gradpsi, X)
    x = X[1]; y = X[2]; z = X[3];
    grad_real = gradpsi[1]; grad_imag = gradpsi[2];
    dx_r = grad_real[1]; dy_r = grad_real[2]; dz_r = grad_real[3];
    dx_i = grad_imag[1]; dy_i = grad_imag[2]; dz_i = grad_imag[3];
    grad_v_r = [dx_r[x,y,z], dy_r[x,y,z], dz_r[x,y,z]];
    grad_v_i = [dx_i[x,y,z], dy_i[x,y,z], dz_i[x,y,z]];


    cross_v = cross(grad_v_r,  grad_v_i)
    wps = cross_v
    return wps;
end

# Pseudo-vorticity 2D, needs X = [x, y] for the grid and array gradient as parameters.
# Returns vector orthog to x, y. Good to sanity check as wps = 0 if no nearby vorticies.
function wps_v_2D(gradpsi, X)
    x = X[1]; y = X[2];
    grad_real = gradpsi[1]; grad_imag = gradpsi[2];
    dx_r = grad_real[1]; dy_r = grad_real[2];
    dx_i = grad_imag[1]; dy_i = grad_imag[2];
    grad_v_r = [dx_r[x,y], dy_r[x,y], 0];
    grad_v_i = [dx_i[x,y], dy_i[x,y], 0];

    cross_v = cross(grad_v_r, grad_v_i)
    wps = cross_v;
    return wps;
end

# Find closest element to x in array
function find_nearest(x, arr)
    nearest = arr[1];
    nearest_idx = 1;
    diff = abs(x-nearest);

    for i in 2:length(arr)
        if abs(x-arr[i]) < diff
            nearest = arr[i];
            nearest_idx = i;
            diff = abs(x-arr[i]);
        end
    end
    return [nearest,  Int(nearest_idx)];
end


#=
myarr = rand(128,128,32) + im*rand(128,128,32)
x = LinRange(-15, 15, 128);
y = x;
z = LinRange(-4, 4, 32)
X = [x, y, z]
grad = gradient_3D(myarr, X)
wps = wps_norm(grad, 10, 10, 10)
norm(wps)
grad_real = grad[1]; grad_imag = grad[2];
dx_r = grad_real[1]; dy_r = grad_real[2]; dz_r = grad_real[3];
dx_i = grad_imag[1]; dy_i = grad_imag[2]; dz_i = grad_imag[3];

=#
# Older slower way
#=
@time for zin in 1:length(z)
    for yin in 1:length(y)
        for xin in 2:length(x)
            ∇xarr_real[xin, yin, zin] = (arr_real[xin, yin, zin] - arr_real[xin-1, yin, zin]) / dx;
        end 
    end
end
#fix boundaries later
∇xarr_real[1, :, :] .= 0;

for zin in 1:length(z)
    for xin in 1:length(x)
        for yin in 2:length(y)
            ∇yarr_real[xin, yin, zin] = (arr_real[xin, yin, zin] - arr_real[xin, yin-1, zin]) / dy;
        end 
    end
end
∇yarr_real[:, 1, :] .= 0;

for xin in 1:length(x)
    for yin in 1:length(y)
        for zin in 2:length(z)
            ∇zarr_real[xin, yin, zin] = (arr_real[xin, yin, zin] - arr_real[xin, yin, zin-1]) / dz;
        end 
    end
end
∇zarr_real[:, :, 1] .= 0;

=#

# Returns [nearest value, index] 


