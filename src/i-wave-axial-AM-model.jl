#=
k-space model of an axial i-wave orbital AM, τ = orbital, σ = spin
=#
@with_kw struct AM_presets
    a0::Float64
    μ::Float64
    t1::Float64
    t4::Float64
    Δsz::Float64
    Δtz::Float64
    neel::Vector{Float64}
end

am_presets() = AM_presets(1,0, 1,1,0,0,[1,0,0])

k_hamiltonian(k, p) = k_hamiltonian(k, p.a0, p.μ, p.t1, p.t4, p.Δsz, p.Δtz, p.neel)

function k_hamiltonian(k, a0, μ, t1, t4, Δsz, Δtz, neel)
    h = zeros(ComplexF64, 4, 4)
    k_hamiltonian!(h, promote_3d(k), a0, μ, t1, t4, Δsz, Δtz, neel)
    return h
end

function k_hamiltonian!(h, k, a0, μ, t1, t4, Δsz, Δtz, neel)
    a0, a1, a2 = iwave_lattice_vectors(a0)
    ns = promote_3d.(first_neighbours(a1, a2, a0 .* [-0.5, -√3/2]))
    n4as = as_fourth_neighbours(a0)
    n4bs = bs_fourth_neighbours(a0)
    cos_dif = am_order_parameter(n4as, n4bs, k)
    h[1:2,1:2]   .+= 2t4 * (cos_dif) .* τx
    h[3:4,3:4]   .+= 2t4 * (cos_dif) .* τx
    h[diagind(h)] .+= t1 * sum([cos(k' * n) for n in ns])
    h[1:2,3:4] .+=  (neel[1] -1im * neel[2])  .* τx
    h[3:4,1:2] .+= conj(h[1:2,3:4])
    h[1:2,1:2] .+= neel[3] .* τx
    h[3:4,3:4] .+= -neel[3] .* τx
    h .+=  Δsz .* kron([1 0;0 -1], I(2)) .+  Δtz .* kron(I(2), [1 0;0 -1]) .- μ * I(4)
    return h
end


function iwave_lattice_vectors(a0)
    a1 = a0 .* [1.0,0]
    a2 = a0 .* [-1/2, √3/2]
    return a0, a1, a2
end

am_order_parameter(n4as, n4bs, k) = 
    sum([cos(k' * n) for n in n4as[[1,2,3]]]) - sum([cos(k' * n) for n in n4bs[[1,2,3]]])

first_neighbours(a1, a2, a3) = a1, -a1, a2, -a2, a3, -a3
as_fourth_neighbours(a0) = promote_3d.(fourth_neighbours(a0 .* [-1/2, 3√3/2]))
bs_fourth_neighbours(a0) = promote_3d.(fourth_neighbours(a0 .* [ 1/2, 3√3/2]))
fourth_neighbours(r) = [rot(θ, r) for θ in 0:π/3:2π-π/3]

rot(θ, v) = rotmat(θ) * v
rotmat(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

#___________________________________________________________________
dh_k(k, pam, i) = dh_k(k, i, pam.a0, pam.t1, pam.t4)

function dh_k(k, i, a0, t1, t4)
    h = zeros(ComplexF64, 4, 4)
    dh_k!(h, promote_3d(k), i, a0, t1, t4)
    return h 
end

function dh_k!(h, k, i, a0, t1, t4)
    if i ≠ :z
        a0, a1, a2 = iwave_lattice_vectors(a0)
        ns = promote_3d.(first_neighbours(a1, a2, a0 .* [-0.5, -√3/2]))
        n4as = as_fourth_neighbours(a0)
        n4bs = bs_fourth_neighbours(a0)
        sin_dif = d_am_order_parameter(n4as, n4bs, k, i)
        h[1:2,1:2] .+= 2t4 * sin_dif .* τx
        h[3:4,3:4] .+= 2t4 * sin_dif .* τx
        h[diagind(h)] .+= t1 * sum([-sin(k' * n)*n[Optics_in_the_length_gauge.symb_to_ind(i)] 
            for n in ns])
    else nothing end
end

d_am_order_parameter(n4as, n4bs, k, i) = 
    sum([-sin(k' * n)*n[Optics_in_the_length_gauge.symb_to_ind(i)] for n in n4as[[1,2,3]]]) + 
    sum([sin(k' * n)*n[Optics_in_the_length_gauge.symb_to_ind(i)] for n in n4bs[[1,2,3]]])


ddh_k(k, pam, i,j) = ddh_k(k, i,j, pam.a0, pam.t1, pam.t4)
function ddh_k(k, i, j, a0, t1, t4)
    h = zeros(ComplexF64, 4, 4)
    ddh_k!(h, promote_3d(k), i, j, a0, t1, t4)
    return h
end


function ddh_k!(h, k, i, j, a0, t1, t4)
    if i ≠ :z && j ≠ :z
        a0, a1, a2 = iwave_lattice_vectors(a0)
        ns = promote_3d.(first_neighbours(a1, a2, a0 .* [-0.5, -√3/2]))
        n4as = as_fourth_neighbours(a0)
        n4bs = bs_fourth_neighbours(a0)
        cos_ddif = dd_am_order_parameter(n4as, n4bs, k, i, j)
        h[1:2,1:2] .+= 2t4 * (cos_ddif) .* τx
        h[3:4,3:4] .+= 2t4 * (cos_ddif) .* τx
        h[diagind(h)] .+= t1 * sum([-cos(k' * n)*n[Optics_in_the_length_gauge.symb_to_ind(i)]*n[Optics_in_the_length_gauge.symb_to_ind(j)] for n in ns])
    else nothing end
end

dd_am_order_parameter(n4as, n4bs, k, i, j) = 
    sum([-cos(k' * n)*n[Optics_in_the_length_gauge.symb_to_ind(i)]*n[Optics_in_the_length_gauge.symb_to_ind(j)] for n in n4as[[1,2,3]]]) + 
    sum([cos(k' * n)*n[Optics_in_the_length_gauge.symb_to_ind(i)]*n[Optics_in_the_length_gauge.symb_to_ind(j)] for n in n4bs[[1,2,3]]])