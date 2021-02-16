using VortexDistributions

function find_vortices3D(psi, X)
    x = X[1]; y = X[2];
    Nx = length(x); Ny = length(y);
    if length(X) == 3
        z = X[3];
        Nz = length(z);
 
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
    else 
        
        vort_array = vortex_array(findvortices(psi));
    end
    temp = zeros(length(vort_array[:, 1]))
    vort_array = hcat(vort_array, temp);


    return vort_array;
end

function find_vortices3D_v2(psi, X)
    z = X[3];
    vorts_3d = Array{Float64, 2}[]

    for zidx in 1:length(z)
        vorts = vortex_array(findvortices(Torus(psi[:, :, zidx], x, y)));
        temp = zeros(length(vorts[:, 1]))
        vorts = hcat(vorts, temp);
        push!(vorts_3d, vorts);
    end
    return vorts_3d;
end
        

