################################
#wps array/heatmap
psi = xspace(sol(tslice), sim);
grad = gradient_3D_cent(psi, X);
grad_int = grad_itp(grad, X)

Nx = 10; Ny = 10;
A = LinRange(x[1], x[end], Nx); B = A;
zslice = 16;
tslice = 10;


wps_array = zeros(Nx, Ny)

for i in 1:Nx
    for j in 1:Ny
        wps_array[i, j] = norm(wps_int(grad, X, [A[i], B[j], z[zslice]]))
    end
end
wps_array
using Plots
heatmap(A, B, wps_array)

########################
using LinearAlgebra, Interpolations
x = LinRange(-15, 15, 128); y = x;
z = LinRange(-4, 4, 32);
X = [x, y, z]
myarr = rand(128, 128, 32) + im*rand(128, 128, 32);
grad = gradient_3D_cent(myarr, X)
grad_int = grad_itp(grad, X)
grad_int_r = grad_int[1];
####################


## Setup (params)
tslice = 20;
zslice = 16;
dz = z[2]-z[1]
psi = xspace(sol(tslice), sim);
grad = gradient_3D_cent(psi, X);
grad_int = grad_itp(grad, X)

# Vortices at current slice and slice above
vorts_zslice = vortex_array(findvortices(Torus(psi[:, :, zslice], x, y)))
vorts_zslice_up = vortex_array(findvortices(Torus(psi[:, :, zslice+1], x, y)))

# Pick a vort and check difference in coordinates, has to be positive charge
num_vort = 2
vort = vorts_zslice[num_vort, :]
vort_up = vorts_zslice_up[num_vort, :]; 
vort_diff = vort_up[1:2] - vort[1:2];

# Find pseudo vorticity and normalise it to wps_z = dz
wps = wps_int(grad_int, [vort[1], vort[2], z[zslice]])
wps_norm = wps .* (dz/abs(wps[3]))
vort_diff






