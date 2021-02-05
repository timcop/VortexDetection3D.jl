# Testing wps at different z levels
tslice = 3;
zslice = 16;
which_vort = 1;
psi = xspace(sol(tslice), sim);
grad = gradient_3D_cent(psi, X);
vorts_zslice = vortex_array(findvortices(Torus(psi[:, :, zslice], x, y)))
v1 = vorts_zslice[which_vort, :]
wps_v = wps_int(grad, X, [v1[2], v1[1], z[zslice]])

wps_v = wps_v ./ abs(wps_v[3]) #normalise so z component is 1
vorts_zslice_under = vortex_array(findvortices(Torus(psi[:, :, zslice-1], x, y)))
vorts_zslice_above = vortex_array(findvortices(Torus(psi[:, :, zslice+1], x, y)))

v1[1:2]
v_above = vorts_zslice_above[which_vort, 1:2]
v_under = vorts_zslice_under[which_vort, 1:2]
pseudo_vort = v1[1:2] + wps_v[1:2]

################################
#wps array/heatmap

Nx = 10; Ny = 10;
A = LinRange(x[1], x[end], Nx); B = A;
zslice = 16;
tslice = 10;
psi = xspace(sol(tslice), sim);
grad = gradient_3D_cent(psi, X);
grad_int = grad_itp(grad, X)

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
grad_int_r = grad_int[1






