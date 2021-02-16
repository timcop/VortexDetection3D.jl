using FourierGPE
using JLD2, FileIO

## parameters
    L = (30.,30.,8.)
    N = (128,128,32)
    sim = Sim(L,N)
    @unpack_Sim sim;


## Initialize simulation
    γ = 0.05
    μ = 25.0
    tf = 2.5/γ
    Nt = 100
    t = LinRange(0.,tf,Nt)

## potential
    import FourierGPE.V
    V(x,y,z,t) = 4*z^2


## random initial state
    x,y,z = X
    ψi = randn(N)+im*randn(N)
    ϕi = kspace(ψi,sim)

    @pack_Sim! sim;

## Evolve in k space
    @time sol = runsim(sim); # will take a few minutes to run.


function dense(phi)
    ψm = xspace(phi,sim)
    density = abs2.(ψm)
    pmax = maximum(density)
    return density/pmax
end

using JLD2 

psi1 = sol[50]
psi2 = sol[55]
psi3 = sol[60]

@save "3dquenchslab_data.jld2" psi1 psi2 psi3

## load test data
@load "3dquenchslab_data.jld2" psi1 psi2 psi3