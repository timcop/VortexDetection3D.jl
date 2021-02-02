#########################################################
function gradient_3D(sol, tslice)
    psi = xspace(sol(tslice), sim)
    psi_real = real(psi)
    psi_imag = imag(psi)

    dx = x[2] - x[1];
    dy = y[2] - y[1];
    dz = z[2] - z[1];
  
    dxpsi_real = similar(psi_real);
    dypsi_real = dxpsi_real; dzpsi_real = dxpsi_real;
    dxpsi_imag = similar(psi_imag);
    dypsi_imag = dxpsi_imag; dzpsi_imag = dxpsi_imag; 

    for xin in 2:length(x)
        dxpsi_real[xin, :, :] .= (psi_real[xin, :, :] .- psi_real[xin-1, :, :]) ./ dx;
        dxpsi_imag[xin, :, :] .= (psi_imag[xin, :, :] .- psi_imag[xin-1, :, :]) ./ dx;
    end
    dxpsi_real[1, :, :] .= 0;
    dxpsi_imag[1, :, :] .= 0;

    for yin in 2:length(y)
        dypsi_real[:, yin, :] .= (psi_real[:, yin, :] .- psi_real[:, yin-1, :]) ./ dy;
        dypsi_imag[:, yin, :] .= (psi_imag[:, yin, :] .- psi_imag[:, yin-1, :]) ./ dy;
    end
    dypsi_real[:, 1, :] .= 0;
    dypsi_imag[:, 1, :] .= 0;

    for zin in 2:length(z)
        dzpsi_real[:, :, zin] .= (psi_real[:, :, zin] .- psi_real[:, :, zin - 1]) ./ dz;
        dzpsi_imag[:, :, zin] .= (psi_imag[:, :, zin] .- psi_imag[:, :, zin - 1]) ./ dz;
    end
    dzpsi_real[:, :, 1] .= 0;
    dzpsi_imag[:, :, 1] .= 0;

    GRADpsi_real = [dxpsi_real, dypsi_real, dzpsi_real];
    GRADpsi_imag = [dxpsi_imag, dypsi_imag, dzpsi_imag];
    GRADpsi = [GRADpsi_real, GRADpsi_imag];
    return GRADpsi;
end

    ###### Lets try take the cross product of these
    #First lets say there's a vortex at (10, 10, 10)

    #grad at vortex
    using LinearAlgebra
function wps_norm(gradpsi, x, y, z)
    grad_real = gradpsi[1]; grad_imag = gradpsi[2];
    dx_r = grad_real[1]; dy_r = grad_real[2]; dz_r = grad_real[3];
    dx_i = grad_imag[1]; dy_i = grad_imag[2]; dz_i = grad_imag[3];
    grad_v_r = [dx_r[x,y,z], dy_r[x,y,z], dz_r[x,y,z]];
    grad_v_i = [dx_i[x,y,z], dy_i[x,y,z], dz_i[x,y,z]];


    cross_v = cross(grad_v_r,  grad_v_i)
    wps_v = cross_v

    normalize!(wps_v)
    return wps_v;
end




# Older slower way
#=
@time for zin in 1:length(z)
    for yin in 1:length(y)
        for xin in 2:length(x)
            ∇xpsi_real[xin, yin, zin] = (psi_real[xin, yin, zin] - psi_real[xin-1, yin, zin]) / dx;
        end 
    end
end
#fix boundaries later
∇xpsi_real[1, :, :] .= 0;

for zin in 1:length(z)
    for xin in 1:length(x)
        for yin in 2:length(y)
            ∇ypsi_real[xin, yin, zin] = (psi_real[xin, yin, zin] - psi_real[xin, yin-1, zin]) / dy;
        end 
    end
end
∇ypsi_real[:, 1, :] .= 0;

for xin in 1:length(x)
    for yin in 1:length(y)
        for zin in 2:length(z)
            ∇zpsi_real[xin, yin, zin] = (psi_real[xin, yin, zin] - psi_real[xin, yin, zin-1]) / dz;
        end 
    end
end
∇zpsi_real[:, :, 1] .= 0;

=#