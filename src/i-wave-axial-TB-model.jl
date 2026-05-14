#= 
    TB model for the i-wave orbital axial altermagnet on a triangular lattice
=#
"""
4 orbitals, inner is f orbital index 1, 2, outer is spin up down.
"""
function tb_hamiltonian(p::AM_presets)
    return hamiltonian(am_lat(p), am_onsite_model(p) + am_hopping_model(p), orbitals = 4)
end

"""
real space triangular lattice
"""
function am_lat(p::AM_presets)
    Rs = lattice_vectors(p.a0)# lattice vectors
    return lattice(sublat((0,0)), bravais = Tuple(Tuple.(Rs)))
end

# single valley model
function am_onsite_model(p::AM_presets)
    on_mat = SMatrix{4,4}([p.neel[3] .* τx 0I; 0I -p.neel[3] .* τx] + 
                [0I p.neel[1] .* τx; p.neel[1] .* τx  0I]+
                [0I -1im * p.neel[2].* τx; 1im * p.neel[2] .* τx  0I])
    return onsite(on_mat)
end

am_hopping_model(p::AM_presets) =
    hopping_1stneighbors(p) + hopping_4thneighbors(p)

hopping_1stneighbors(p) = hopping(p.t1*1I, range = neighbors(1)) 

"""
hopping 4 neighbours, t4 changes sign in the two sets of C6 related neighbors
"""
function hopping_4thneighbors(p) 
    Rs = lattice_vectors(p.a0)
    v= 2*Rs[1] - Rs[2] # position of a given 4th neighbor required to set the sign of t4
    return hopping((r, dr) -> t4_sign(r, v) * p.t4 .* SMatrix{4,4}([τx 0I; 0I τx]), 
        range = (neighbors(4), neighbors(4)))
end

# """ hoppings up to am. """
# function hopping_firstneighbours(p::Paramsθ)
#     t_HH = p.t_HH
#     t_TH = p.t_TH
#     model = hopping(t_HH, sublats = :XM => :MX, range = neighbors(2)) + 
#     hopping(t_HH, sublats = :MX => :XM, range = neighbors(2)) +
#     hopping((r, dr) -> -t_TH * conj(gk(-dr)), sublats = :MM => :MX, range = neighbors(2))+
#     hopping((r, dr) -> -t_TH * gk(dr), sublats = :MX => :MM, range = neighbors(2))  +
#     hopping((r, dr) -> t_TH* conj(gk(dr)), sublats = :XM => :MM, range = neighbors(2))  +
#     hopping((r, dr) -> t_TH* gk(-dr), sublats = :MM => :XM, range = neighbors(2))  
# end

"""
t4 picks a negative sign in one the two sets of 4th neighbors that are related by C6
"""
function t4_sign(r, v) 
    return ifelse(any(norm(v - rot(θ, 2 .* r)) < 1e-8 for θ in 0:π/3:2π-π/3) == true, 1, -1)
end



