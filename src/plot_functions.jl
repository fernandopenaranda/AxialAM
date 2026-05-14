function plot_am_formfactor(pam)
    N = 400
    # high_sym_line = (:Γ, :K, :M, :K, :Γ)
    n4as = as_fourth_neighbours(pam.a0)
    n4bs = bs_fourth_neighbours(pam.a0)
    kx = 1.2/pam.a0 .* collect(range(-π, π, length=N))
    ky = 1.2/pam.a0 .* collect(range(-π, π, length=N))
    Z = [am_order_parameter(n4as, n4bs, [kxi, kyi]) for kxi in kx, kyi in ky]
    # Plot
    fig = Figure(size = (400,530))
    ax = Axis(fig[1,1], xlabel="kx", ylabel="ky", title="ϕ(k)")
    heatmap!(ax, kx, ky, Z, colormap=:inferno)
    fig
end

function iwave_bands(sgnum, N, pam::AM_presets; high_sym_line = nothing)
    fig = Figure(size = (400,530))
    ax = Axis(fig[1,1], ylabel = "E (eV)")
    iwave_bands(ax, sgnum, N, pam)
    return fig
end

function iwave_bands(ax, sgnum, N, pam::AM_presets; high_sym_line = nothing)
    a0, R1, R2 = iwave_lattice_vectors(pam.a0)
    Rs = (R1,R2)
    dim = size(k_hamiltonian([0.,0.],pam),1)
    for i in 1:dim
        obs(k) = eigenvals_model(pam, k, i)
        if high_sym_line === nothing
            BZpaths.plot_observable_in_kpath!(ax, obs, Rs, sgnum, N, color = ifelse(i%2 == 1, :blue, :red))
        else
            BZpaths.plot_observable_in_kpath!(ax, obs, Rs, sgnum, N, high_sym_line)
        end
    end
end

function iwave_bands_3d(sgnum, N, pam::AM_presets_3d; high_sym_line = nothing)
    fig = Figure(size = (500,530))
    ax = Axis(fig[1,1], ylabel = "E (eV)", title = "a0 = $(pam.a0), μ = $(pam.μ), t1 = $(pam.t1), t2 = $(pam.t2), t4 = $(pam.t4), tperp = $(pam.tperp),
     Δsz = $(pam.Δsz),  Δtz = $(pam.Δtz), neel = $(pam.neel)")
    iwave_bands_3d(ax, sgnum, N, pam)
    return fig
end

function iwave_bands_3d(ax, sgnum, N, pam::AM_presets_3d; high_sym_line = nothing)
    a0, R1, R2 = iwave_lattice_vectors(pam.a0)
    Rs = (promote_3d(R1), promote_3d(R2), [0,0,a0])
    dim = size(k_hamiltonian_3d([0.,0., 0.],pam),1)
    obs(k, i) = eigenvals_model_3d(pam, k, i)
    obs_paths!(ax, dim, obs, Rs, sgnum, N, high_sym_line)
end

function iwave_Fijcuts(sgnum, N, pam::AM_presets_3d, p,  dirF1, dirF2; high_sym_line = nothing)
    fig = Figure(size = (400,530))
    
    ax = Axis(fig[1,1], ylabel = "F$(dirF1)$(dirF2) (Å³)")
    iwave_Fijcuts(ax, sgnum, N, pam, p,  dirF1, dirF2)
    return fig
end

function iwave_Fijcuts(ax, sgnum, N, pam::AM_presets_3d, p, dirF1, dirF2; high_sym_line = nothing)
    a0, R1, R2 = iwave_lattice_vectors(pam.a0)
    Rs = (promote_3d(R1), promote_3d(R2), [0,0,a0])
    dim = size(k_hamiltonian_3d([0.,0., 0.],pam),1)
    func = Optics_in_the_length_gauge.F
    obs(k, band) = func(p, k, dirF1, dirF2)[band]
    colorlist = [:red, :orange, :blue, :green]
    obs_paths!(ax, dim, obs, Rs, sgnum, N, high_sym_line, colorlist = colorlist)
end

function custom_fig_iwave_Fijcuts(sgnum, N, pam::AM_presets_3d, p; high_sym_line = nothing)
    fig = Figure(size = (1000,830))
    Label(
    fig[0, 1:3],
    "PS_QM = $(p.QM_switch), PS_orbital = $(p.PS_orbital_switch)",
    fontsize = 15
)
    tuple_list = [(1,1,:x,:x),(1,2,:x,:y),(1,3,:x,:z),(2,1,:y,:x),(2,2,:y,:y),(2,3,:y,:z),(3,1,:z,:x),(3,2,:z,:y), (3,3,:z,:z)]
    for (i,j, dirF1, dirF2) in tuple_list
        ax = Axis(fig[i,j], ylabel = "F$(dirF1)$(dirF2) (Å³)")
        iwave_Fijcuts(ax, sgnum, N, pam, p,  dirF1, dirF2)
    end
    return fig
end

function obs_paths!(ax, dim, obs, Rs, sgnum, N, high_sym_line; colorlist = nothing)
    for i in 1:dim
        band_obs(k) = obs(k, i)
        if high_sym_line === nothing
            if colorlist === nothing 
            BZpaths.plot_observable_in_kpath!(ax, band_obs, Rs, sgnum, N, color = ifelse(i%2 == 1, :blue, :red))
            else
                BZpaths.plot_observable_in_kpath!(ax, band_obs, Rs, sgnum, N, color = colorlist[i])
            end
        else
            BZpaths.plot_observable_in_kpath!(ax, band_obs, Rs, sgnum, N, high_sym_line)
        end
    end
end

function k_mesh_evals(func, pam, p; u0 = 0,  botbounds = [-0.5,-0.5], topbounds = [0.5,0.5], kpoints = 100)
    Gs = p.gs
    N = floor(Int, kpoints^(1/3))
    uas = range(botbounds[1], topbounds[1], length=N)
    ubs = range(botbounds[2], topbounds[2], length=N)
    kx = [Optics_in_the_length_gauge.transform_k([ua, ub], Gs)[1] for ua in uas, ub in ubs]
    ky = [Optics_in_the_length_gauge.transform_k([ua, ub], Gs)[2] for ua in uas, ub in ubs]
    f(ua,ub) = func(p, Optics_in_the_length_gauge.transform_k([ua,ub],Gs))
    Z = [f(ua,ub) for ua in uas, ub in ubs]
    return kx, ky, Z
end

function k_mesh_evals_3d(func, p; uz = 0,  botbounds = [-0.5,-0.5, -0.5], topbounds = [0.5,0.5,0.5], kpoints = 100)
    Gs = p.gs
    N = floor(Int, kpoints^(1/3))
    uas = range(botbounds[1], topbounds[1], length=N)
    ubs = range(botbounds[2], topbounds[2], length=N)
    kx = [Optics_in_the_length_gauge.transform_k([ua, ub, uz], Gs)[1] for ua in uas, ub in ubs]
    ky = [Optics_in_the_length_gauge.transform_k([ua, ub, uz], Gs)[2] for ua in uas, ub in ubs]
    f(ua,ub,uz) = func(p, Optics_in_the_length_gauge.transform_k([ua,ub, uz],Gs))
    Z = [f(ua,ub,uz) for ua in uas, ub in ubs]
    return kx, ky, Z
end

function plot_kresolved(kx, ky, Z; label = "f(a, b)", titlelab = "", colormap = missing, color = :black)
    fig = Figure(size = (600, 500))
    ax = Axis(fig[1, 1], xlabel="ka [1/Å]", ylabel="kb [1/Å]", title = titlelab)
    if colormap === missing
        cmap = cgrad([:white, color])#reverse(cgrad(:grays))  # CairoMakie built-in colormap, light-to-dark blues
    else cmap = colormap end
    lim = maximum(abs, Z)
    hm = surface!(ax, kx, ky, zeros(size(Z)),
                color=Z,
                shading=NoShading,
                colormap=cmap,
                colorrange = (-lim, lim))
    Colorbar(fig[1, 2], hm, label=label)
    # ylims!(ax, -2π-π/2, 2π+π/2)
    fig
end