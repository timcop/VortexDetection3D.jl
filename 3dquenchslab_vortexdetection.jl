# sol[t][x; y; z] where t ranges from 0 -> 200 and [x; y; z] is a 128x128x32 array in k space. First pick a time slice to examine using sol[t],
# then convert to position space using xspace(psi, sim)
using VortexDistributions, FourierGPE, FGPEexamples, Plots
#Changing something
#Params
Nx = 128; Ny = Nx;
Nz = 32;
Lx = 200; Ly = Lx;
Lz = 50;
x = LinRange(-Lx / 2, Ly / 2, Nx); y = x;
z = LinRange(-Lz/2, Lz/2, Nz);

psi_3d = sol[199];
psi_x = xspace(psi_3d, sim);
psi = Torus(psi_x[:,:,16], x,y);
vort_array = vortex_array(findvortices(psi))

vort_pos = vort_array[:, 1:2]
using Plotly
plotly()
Plots.heatmap(x, y, angle.(psi_x[:, :, 16]))

Plots.plot!(Plots.scatter!(vort_pos[:, 2], vort_pos[:, 1], color =:red))




function vortices_plot(psi_k, x, y, zidx)
    psi_xspc = xspace(psi_kspc, sim);
    psi = Torus(psi_xspc[:, :, zidx], x, y);
    vort_array = vortex_array(findvortices(psi));
    vort_pos = vort_array[:, 1:2];
 
    @views Plots.heatmap(x, y, angle.(psi_xspc[:, :, zidx]), title = "z level = $zidx", label = "Phase")
    Plots.scatter!(vort_pos[:, 2], vort_pos[:, 1], color =:white, label = "Vortices")
    Plots.xlabel!("x")
    Plots.ylabel!("y")

end

tslice = 80;
psi_kspc = sol[tslice];

using Plotly;
gr()
anim = @animate for zidx = 1:Nz
    #println("stage1")
    vortices_plot(psi_kspc, x, y, zidx)
    #println("stage2")
end
gif(anim, "anim_fps15.gif", fps = 8)

vortices_plot(psi_kspc, x, y, 32)


vorts = zeros(0);
for j in 1:20

    
    tslice = j*10;
    psi_3d = sol[tslice];
    psi_x = xspace(psi_3d, sim);
    psi = Torus(psi_x[:,:,1], x,y);

# Initialise the vortex_array at the first z level (z[i=1]), use hcat and vcat to concatenate the z position onto the predefined vortex_array for 2D.
# This is inefficent and would be better to ammend the PointVortex data type to allow x,y,z values
    vfound = findvortices(psi);
    vort_array = vortex_array(vfound);
    i = 1;
    v = z[i]*ones(length(vort_array[:, 1])); #vector of same length as vortices found at current z level, each element of vector is the current z level
    vort_array = hcat(vort_array, v);
    vort_array = view(vort_array, :, [1, 2, 4, 3]);
    total_vortices = length(vort_array[:, 1]); #initailise to the number of vortices at z[i=1]
    for i in 2:32

        println(i)
        psi = Torus(psi_x[:,:,i], x, y);
        temp = vortex_array(findvortices(psi));
        total_vortices += length(temp[:, 1]);

        v = z[i]*ones(length(temp[:, 1]));
        temp = hcat(temp, v)
        temp = view(temp, :, [1, 2, 4, 3])
        vort_array = vcat(vort_array, temp)

    end

# Compare length(vort_array[:, 1]) and total_vortices for sanity check
    append!(vorts, total_vortices);

end

using Plots;
x = 1:20;
gr()
plot(vorts)