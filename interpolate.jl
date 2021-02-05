
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

function grad_itp(grad, X)
    x = X[1]; y = X[2]; z = X[3];
    x = LinRange(x[1], x[end], length(x));
    y = LinRange(y[1], y[end], length(y));
    z = LinRange(z[1], z[end], length(z));
    grad_r = grad[1]; grad_i = grad[2];
    grx = grad_r[1];
    gry = grad_r[2];
    grz = grad_r[3];
    gix = grad_i[1];
    giy = grad_i[2];
    giz = grad_i[3];
    grx_itp = interpolate(grx, BSpline(Cubic(Line(OnGrid()))));
    gry_itp = interpolate(gry, BSpline(Cubic(Line(OnGrid()))));
    grz_itp = interpolate(grz, BSpline(Cubic(Line(OnGrid()))));
    gix_itp = interpolate(gix, BSpline(Cubic(Line(OnGrid()))));
    giy_itp = interpolate(giy, BSpline(Cubic(Line(OnGrid()))));
    giz_itp = interpolate(giz, BSpline(Cubic(Line(OnGrid()))));
    sgrx_itp = scale(grx_itp, x, y, z);
    sgry_itp = scale(gry_itp, x, y, z);
    sgrz_itp = scale(grz_itp, x, y, z);
    sgix_itp = scale(gix_itp, x, y, z);
    sgiy_itp = scale(giy_itp, x, y, z);
    sgiz_itp = scale(giz_itp, x, y, z);
    gr_v = [sgrx_itp, sgry_itp, sgrz_itp];
    gi_v = [sgix_itp, sgiy_itp, sgiz_itp];
    grad = [gr_v, gi_v];
    return grad;
end

function wps_int(grad, Xin)
   
    grad_r = grad[1]; grad_i = grad[2];

    gr_v = [grad_r[1](Xin[1], Xin[2], Xin[3]), grad_r[2](Xin[1], Xin[2], Xin[3]), grad_r[3](Xin[1], Xin[2], Xin[3])];
    gi_v = [grad_i[1](Xin[1], Xin[2], Xin[3]), grad_i[2](Xin[1], Xin[2], Xin[3]),  grad_i[3](Xin[1], Xin[2], Xin[3])];

    wps = cross(gr_v, gi_v);
    return wps;
end

