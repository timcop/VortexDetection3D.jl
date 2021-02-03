gr(xlabel="x",ylabel="y",legend=false)

# make a simple 2D test field
Nx = 400; Ny = Nx
Lx = 200; Ly = Lx
x = LinRange(-Lx / 2, Ly / 2, Nx); y = x
X = [x, y]
psi0 = one.(x*y') |> complex

# doubly periodic boundary conditions
psi = Torus(psi0,x,y)

# make a point vortex
#pv = PointVortex(x[15],x[15],1)

# make a scalar GPE vortex with exact core
#spv = ScalarVortex(pv);
#vortex!(psi,spv);

#=
for i in 1:10
    pv = PointVortex(x[i*10], x[i*10], 1)
    spv = ScalarVortex(pv)
    vortex!(psi, spv)
end
=#

# make some more random vortices
vort = rand_vortex(50,psi);
vortex!(psi,vort);

#vortex_array(pv)

vfound = vortex_array(findvortices(psi))
#find_vortices3D(psi, X)

@time vort = findvortices(psi)

@unpack ψ,x,y = psi

grad = gradient_2D(ψ, [x, y])
wps_v_2D(grad, [3, 3])

wps_arr = similar(real(ψ))

@time for i in 1:Nx
    for j in 1:Ny
        wps_arr[i, j] = wps_v_2D(grad, [i,j])[3]
    end
end
wps_arr = abs.(wps_arr)
max = findmax(wps_arr)[1]
wps_arr = wps_arr ./max
using Plots;

#heatmap(x, y, wps_arr)
heatmap(x, y, grad[2])
scatter!(vfound[:, 2], vfound[:, 1])
heatmap(x, y, angle.(ψ))
d = abs2.(ψ);
d = d ./ maximum(d)
heatmap(x, y, d)
scatter!(vfound[:, 2], vfound[:, 1]) 
