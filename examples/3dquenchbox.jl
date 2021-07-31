using Plots, LaTeXStrings
gr(titlefontsize=12,size=(500,300),colorbar=false);

using FourierGPE

##Set simulation parameters
    L=(16.,16.,16.);
    N=(64,64,64);
    sim = Sim(L,N);
    @unpack_Sim sim;

## Initialse sim
    # parameters
    μ = 25.0;
    γ = 0.05;
    tf = 4/γ;
    Nt = 200;
    t = LinRange(0.,tf,Nt);

## Run sim
    x,y,z = X;
    ψi = randn(N)+im*randn(N);
    ϕi = kspace(ψi,sim);

    @pack_Sim! sim;

## Evolve in k-space
    #sol = runsim(sim); # will take a few minutes to run.


## Density Isosurface
    using Makie

    function dense(psi, sim)
        ψm = xspace(psi,sim)
        density = abs2.(ψm)
        pmax = maximum(density)
        return density/pmax
    end

    function densityfilm(Nt,saveto="media/3dquenchiso.gif")
        scene = Scene()
        tindex = Node(1)

        scene = volume(lift(i -> dense(i), tindex), algorithm = :iso,show_axis=false)

        record(scene, saveto, 1:Nt-10) do i
            tindex[] = i
        end
    end


## Plot density film
    #densityfilm(Nt)

## Plot time slice of density film
    Makie.volume(dense(psi_ring1, sim), algorithm = :iso, show_axis=true)

using JLD2 

## Some solutions with rings 
   #=
    psi_ring1 = sol[91]
    psi_ring2 = sol[105]
    psi_ring3 = sol[110]
    psi_ring4 = sol[65]
    psi_box_25 = sol(25)
    psi_tubes = sol(35)

    @save "3dquenchbox_data.jld2" psi_ring1 psi_ring2 psi_ring3
    @save "3dquenchbox_data2.jld2" psi_box_25 psi_tubes
    @save "3dquenchbox_data3.jld2" psi_ring4
    =#

## Load test data
    @load "src/3dquenchbox_data.jld2" psi_ring1 psi_ring2 psi_ring3
    @load "src/3dquenchbox_data2.jld2" psi_box_25 psi_tubes
    @load "src/3dquenchbox_data3.jld2" psi_ring4
    


## Testing area
    #=
    using Interpolations, LinearAlgebra
    psi = xspace(psi_ring3, sim);
    grad = gradient_3D_cent(psi, X);
    grad_i = grad_itp(grad, X);
    dz = z[2] - z[1]
    zslice = 40
    vorts1 = vortex_array(findvortices(Torus(psi[:, :, zslice], x, y)))
    vorts2 = vortex_array(findvortices(Torus(psi[:, :, zslice+1], x, y)))

    numvort = 2
    diff = vorts2[numvort, 1:2] - vorts1[numvort, 1:2]

    wps = wps_int(grad_i, [vorts1[numvort, 1], vorts1[numvort, 2], z[zslice]])
    wps = wps .* (dz/abs(wps[3]))

    diff
    =#

    



