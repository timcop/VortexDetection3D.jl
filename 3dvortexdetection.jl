using VortexDistributions

function find_vortices3D(sol, tslice)
    Nx = size(sol(1))[1];
    Ny = size(sol(1))[2];
    Nz = size(sol(1))[3];
    Lx = 2*Nx; Ly = 2*Ny; Lz = 2*Nz;
    x = LinRange(-Lx/2, Lx/2, Nx);
    y = LinRange(-Ly/2, Ly/2, Ny);
    z = LinRange(-Lz/2, Lz/2, Nz);

    psi = sol[tslice];
    psi = xspace(psi, sim);

    i = 1;
    psi_z = Torus(psi[:, :, i], x, y);
    vort_array = vortex_array(findvortices(psi_z));
    v = z[i]*ones(length(vort_array[:, 1]));
    

    for i = 1:Nz
        psi_z = Torus(psi[:,:,i], x, y);

        if i == 1
            vort_array = vortex_array(findvortices(psi_z));
            v = z[i]*ones(length(vort_array[:, 1]));
            vort_array = hcat(vort_array, v);
            vort_array = view(vort_array, :, [1,2,4,3]);
        
        else 
            temp = vortex_array(findvortices(psi_z));
            v = z[i]*ones(length(temp[:, 1]));
            temp = hcat(temp, v);
            temp = view(temp, :, [1, 2, 4, 3]);
            vort_array = vcat(vort_array, temp);
        end
    end

    return vort_array;
end


find_vortices3D(sol, 10)