
using FourierGPE
|
#--- parameters
L = (30.,30.,8.)
N = (128,128,32)
sim = Sim(L,N)
@unpack_Sim sim;


#--- Initialize simulation
γ = 0.05
μ = 25.0
tf = 2.5/γ
Nt = 200
t = LinRange(0.,tf,Nt)

# potential
import FourierGPE.V
V(x,y,z,t) = 4*z^2


#--- random initial state
x,y,z = X
ψi = randn(N)+im*randn(N)
ϕi = kspace(ψi,sim)

@pack_Sim! sim;

#--- Evolve in k space
@time sol = runsim(sim); # will take a few minutes to run.


function dense(phi)
    ψm = xspace(phi,sim)
    density = abs2.(ψm)
    pmax = maximum(density)
    return density/pmax
end

psi = xspace(sol(100), sim);

grad = gradient_3D(psi, X)

zslice = 26;
v_array = vortex_array(findvortices(Torus(psi[:, :, zslice], x, y)))

wps_arr = similar(real(psi[:, :, zslice]));
for i in 1:128
    for j in 1:128
        wps_arr[i, j] = norm(wps_v(grad, [i, j, zslice]));
    end
end
wps_arr = abs.(wps_arr);
wps_arr = wps_arr ./ (findmax(wps_arr)[1])
heatmap(x, y, wps_arr)
scatter!(v_array[:, 2], v_array[:, 1])