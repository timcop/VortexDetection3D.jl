



function findidx_uniform(x, arr)
# For arrays where xi is between [-L, L] and uniform spacing
# So xi = -L/2 + (L/N)*i
# => i = (N/L)*(x + L/2)

    N = length(arr);
    L = arr[end]*2;
    idx = (N/L)*(x + L/2);
    if idx < 1
        println("x is not within the bounds of this array");
        return;
    else
        return idx;
    end
end








# Testing wps at different z levels
#=
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
=#
    

