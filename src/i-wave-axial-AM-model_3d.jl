#= 3D k-space model of an axial i-wave orbital AM, τ = orbital, σ = spin =#
@with_kw struct AM_presets_3d
    a0::Float64
    μ::Float64
    t1::Float64
    t2::Float64
    t4::Float64
    tperp::Float64
    Δsz::Float64
    Δtz::Float64
    neel::Vector{Float64}
end

am_presets_3d() = AM_presets_3d(1,0, 1,1,1,1,0,0,[1,0,0])

k_hamiltonian_3d(k, p::AM_presets_3d) =
    k_hamiltonian_3d(k, p.a0, p.μ, p.t1,p.t2, p.t4, p.tperp, p.Δsz, p.Δtz, p.neel)

function k_hamiltonian_3d(k, a0, μ, t1,t2, t4, tperp, Δsz, Δtz, neel)
    h = zeros(ComplexF64, 4, 4)
    k_hamiltonian!(h, k, a0, μ, t1, t4, Δsz, Δtz, neel)
    k_hamiltonian_outofplane!(h, t2, tperp, k, a0)
    return h
end

dh_k_3d(k, pam, i) = dh_k_3d(k, i, pam.a0, pam.t1, pam.t2, pam.t4, pam.tperp)

function dh_k_3d(k, i, a0, t1, t2, t4, tperp)
    h = zeros(ComplexF64, 4, 4)
    dh_k!(h, k, i, a0, t1, t4)
    dk_hamiltonian_outofplane!(h, i, t2, tperp, k, a0)
    return h
end

ddh_k_3d(k, pam, i,j) = ddh_k_3d(k, i,j, pam.a0, pam.t1, pam.t2, pam.t4, pam.tperp)
function ddh_k_3d(k, i, j, a0, t1, t2, t4, tperp)
    h = zeros(ComplexF64, 4, 4)
    ddh_k!(h,k, i, j, a0, t1, t4)
    ddk_hamiltonian_outofplane!(h, i, j, t2, tperp, k, a0)
    return h
end

function k_hamiltonian_outofplane!(h, t2, tperp, k, a0)
    a0, a1, a2 = iwave_lattice_vectors(a0)
    ns = promote_3d.(first_neighbours(a1, a2, a0 .* [-0.5, -√3/2]))
    h .+= 2tperp * cos(k[3] * a0) .* I(4) # vertical hoppings
    h_inter = 2t2 * sum([cos(k' * (n+[0,0,a0])) for n in ns])  #inter-orbital first neighbor hopping out of plane
    h .+=  kron(I(2), [0 h_inter; conj(h_inter) 0])
end

function dk_hamiltonian_outofplane!(h, i, t2, tperp, k, a0)
    a0, a1, a2 = iwave_lattice_vectors(a0)
    ns = promote_3d.(first_neighbours(a1, a2, a0 .* [-0.5, -√3/2]))
    uz = [0,0,a0]
    dh = -2t2 * sum([sin(k' * (n+uz))*(n+uz)[Optics_in_the_length_gauge.symb_to_ind(i)] for n in ns]) 
    h[1:2,1:2] .+= [0 dh; conj(dh) 0]
    h[3:4,3:4] .+= [0 dh; conj(dh) 0]
    if i == :z
        h .+= -2a0* tperp * sin(k[3] * a0) .* I(4) 
    else nothing end
end

function ddk_hamiltonian_outofplane!(h, i, j, t2, tperp, k, a0)
    a0, a1, a2 = iwave_lattice_vectors(a0)
    ns = promote_3d.(first_neighbours(a1, a2, a0 .* [-0.5, -√3/2]))
    uz = [0,0,a0]
    dh = -2t2 * sum([cos(k' * (n+uz))*(n+uz)[Optics_in_the_length_gauge.symb_to_ind(i)]*(n+uz)[Optics_in_the_length_gauge.symb_to_ind(j)] for n in ns]) 
    h[1:2,1:2] .+= [0 dh; conj(dh) 0]
    h[3:4,3:4] .+= [0 dh; conj(dh) 0]
    if i == :z && j == :z
    h .+= -2a0^2* tperp * cos(k[3] * a0) .* I(4)
    else nothing end
end