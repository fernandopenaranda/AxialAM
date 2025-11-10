#=
k-space model of an axial i-wave orbital AM
τ = orbital
σ = spin
=#
function hamiltonian(k; t1 = 1, t4 = 0.1, N = 2, a0 = 1, neel = [1,1,1])
    a1, a2, a3 = lattice_vectors(a0)
    ns = first_neighbours(a1, a2, a3)
    n4as = as_fourth_neighbours(a1, a2, a3)
    n4bs = bs_fourth_neighbours(a1, a2, a3)

    h = zeros(ComplexF64, 4, 4)
    h[1:2,1:2]   .+= 2t4 * (sum([cos(k'* n4as[i]) - cos(k'* n4bs[i]) for i in 1:length(n4as)])) .* τx
    h[3:4,3:4]   .=  h[1:2,1:2]
    h[diagind(h)] .+= 2t1 * sum([cos(k' * n) for n in ns])
    h[1:2,3:4] .+=  neel' * [σx, σy, σz]
    h[3:4,1:2] .= conj(h[1:2,3:4])
    return h
end

first_neighbours(a1, a2, a3) = a1, -a1, a2, -a2, a3, -a3
as_fourth_neighbours(a1, a2, a3) = 2a1-a2, 2a2-a3, 2a3-a1, -2a1+a2, -2a2+a3, -2a3+a1
bs_fourth_neighbours(a1, a2, a3) = 2a1-a3, 2a2-a1, 2a3-a2, -2a1+a3, -2a2+a1, -2a3+a2

function lattice_vectors(a0)
    a1 = a0 .* [1,0]
    a2 = a0 .* [-1/2, √3/2]
    a3 = a0 .* [-1/2, -√3/2]
    return a1, a2, a3
end