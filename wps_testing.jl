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
using Interpolations, VortexDistributions, LinearAlgebra, Plots
zslice = 16;
dz = z[2]-z[1]
psi = xspace(psi2, sim);
grad = gradient_3D_cent(psi, X);
grad_int = grad_itp(grad, X);

# Vortices at current slice and slice above
vorts_zslice = vortex_array(findvortices(Torus(psi[:, :, zslice], x, y)))
vorts_zslice_up = vortex_array(findvortices(Torus(psi[:, :, zslice+1], x, y)))


# Pick a vort and check difference in coordinates, has to be positive charge
num_vort = 8
vort = vorts_zslice[num_vort, :]
vort_up = vorts_zslice_up[num_vort, :]
diff = vort_up[1:2] - vort[1:2]
# Find pseudo vorticity and normalise it to wps_z = dz
wps= wps_int(grad_int, [vort[1], vort[2], z[zslice]])
wps= wps .* (dz/abs(wps[3]))



#####################################
# Plots wps-heatmap and found vortices for certain zslice
using Interpolations, VortexDistributions, LinearAlgebra, Plots
psi = xspace(psi2, sim)
grad = gradient_3D_cent(psi, X);
grad_int = grad_itp(grad, X);
zslice = 16;
v1 = vortex_array(findvortices(Torus(psi[:, :, zslice], x, y)))
wps_int(grad_int, [v1[1, 1], v1[1, 2], z[zslice]])
wps_array_plot(psi, zslice)

gr()
anim = @animate for zidx = 1:length(z)
    wps_array_plot(psi, zidx)
end
gif(anim, "wps_heatmap.gif", fps = 8)

######################################
using Interpolations, LinearAlgebra, VortexDistributions
dz = z[2] - z[1];
ϵ = dz/10;
zslice = 16;
psi = xspace(psi1, sim);
psi_int = psi_itp(psi, X)
psi_2 = psi_itp_2D(psi_int, z[zslice] + 4*ϵ, X);
grad = gradient_3D_cent(psi, X);
grad_int = grad_itp(grad, X)

v1 = vortex_array(findvortices(Torus(psi[:, :, zslice], x, y)))
v2 = vortex_array(findvortices(Torus(psi_2, x, y)))

v1_1 = v1[3, :]
v2_1 = v2[3, :]
diff = v2_1 - v1_1
wps = wps_int(grad_int, [v1_1[1], v1_1[2], z[zslice]])
wps = wps .* ϵ/abs(wps[3])

diff = v2_1 - v1_1






