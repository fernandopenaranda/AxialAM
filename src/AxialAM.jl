module AxialAM

    # using BZpaths
    using LinearAlgebra, Quantica, CairoMakie, StaticArrays
    using Revise

    const s0 = SA[1 0; 0 1]
    const σ0 = SA[1 0; 0 1]
    const σx = SA[0 1; 1 0]
    const σy = SA[0 -1im; 1im 0]
    const σz = SA[1 0; 0 -1]

    const τ0 = SA[1 0; 0 1]
    const τx = SA[0 1; 1 0]
    const τy = SA[0 -1im; 1im 0]
    const τz = SA[1 0; 0 -1]
    include("i-wave-axial-AM-model.jl")
    export  lattice_vectors, hamiltonian

end
