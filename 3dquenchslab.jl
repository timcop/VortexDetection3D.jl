
using FourierGPE

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


#@time find_vortices3D(sol, 100)
gradpsi = gradient_3D(sol, 100)

wps = wps_norm(gradpsi, 10,10,10)
