#=
k-space model of an axial i-wave orbital AM
τ = orbital
σ = spin
=#

@with_kw struct AM_presets
    a0
    t1
    t4
    neel
end

am_presets() = AM_presets(1,1,1,[1,0,0])

k_hamiltonian(k, p::AM_presets) = k_hamiltonian(k, t1 = p.t1, t4 = p.t4, a0 = p.a0, neel = p.neel)

function k_hamiltonian(k; t1 = 1, t4 = 0.1, a0 = 1, neel = [1,1,1])
    a1, a2 = lattice_vectors(a0)
    ns = first_neighbours(a1, a2, [-0.5, -√3/2])
    n4as = as_fourth_neighbours(a0)
    n4bs = bs_fourth_neighbours(a0)
    h = zeros(ComplexF64, 4, 4)
    cos_dif = sum([cos(k'* na) for na in n4as] .- [ cos(k'* nb) for nb in n4bs])
    h[1:2,1:2]   .+= 2t4 * (cos_dif) .* τx
    h[3:4,3:4]   .+=  h[1:2,1:2]
    h[diagind(h)] .+= 2t1 * sum([cos(k' * n) for n in ns])
    h[1:2,3:4] .+=  (neel[1] -1im * neel[2])  .* τx
    h[3:4,1:2] .+= conj(h[1:2,3:4])
    h[1:2,1:2] .+= neel[3] .* τx
    h[3:4,3:4] .+= -neel[3] .* τx
    return h
end

function plot_am_formfactor(;  kws...)
    N = 200
    a1, a2 = lattice_vectors(a0)
    ns = first_neighbours(a1, a2, [-0.5, -√3/2])
    kx = collect(range(-π, π, length=N))
    ky = collect(range(-π, π, length=N))
    # Z =  [     sum([cos([kxi, kyi]' * n) for n in ns])                     for kxi in kx, kyi in ky]
    Z = [am_character([kxi, kyi]) for kxi in kx, kyi in ky]

    # Plot
    fig = Figure(resolution = (600,530))
    ax = Axis(fig[1,1], xlabel="kx", ylabel="ky", title="f(kx, ky) on [-π,π]^2")
    heatmap!(ax, kx, ky, Z, colormap=:viridis)
    # Colorbar(fig[1,2], ax; label="f(kx,ky)")

    fig
end


function am_character(k; a0 = 1)
    n4as = as_fourth_neighbours(a0)
    n4bs = bs_fourth_neighbours(a0)
    return sum([cos(k'* rot(-π, na)) for na in n4as] .+ [ cos(k'* rot(-π, nb)) for nb in n4bs])
end

first_neighbours(a1, a2, a3) = a1, -a1, a2, -a2, a3, -a3
as_fourth_neighbours(a0) =  fourth_neighbours(a0 .* [2.0, +1.7320508075688772])
bs_fourth_neighbours(a0) =  fourth_neighbours(a0 .* [2.5, -0.8660254037844386])

fourth_neighbours(r) = [rot(θ, r) for θ in 0:π/3:2π-π/3]



function lattice_vectors(a0)
    a1 = a0 .* [1.0,0]
    a2 = a0 .* [-1/2, √3/2]
    return a1, a2
end

rotmat(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]
rot120(v) = [cos(2π/3) -sin(2π/3); sin(2π/3) cos(2π/3)] * v

