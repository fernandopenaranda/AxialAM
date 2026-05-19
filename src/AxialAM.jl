module AxialAM
    # using BZpaths
    using LinearAlgebra, StaticArrays
    using CairoMakie
    using SymScan
    using Revise, Parameters
    using BZpaths
    # using Brillouin
    using Crystalline
    using PhysicalConstants
    using PhysicalConstants.CODATA2018
    using Optics_in_the_length_gauge
    using Serialization

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
    include("i-wave-axial-AM-model_3d.jl")
    include("observables.jl")
    include("i-wave-axial-TB-model.jl")
    include("plot_functions.jl")
    include("wrappers.jl")
    include("cluster/cluster_tools.jl")
    include("cluster/create_bashfile.jl")

    export iwave_lattice_vectors, k_hamiltonian, iwave_bands, AM_presets, am_presets, sigma_abc_wrapper, ddh_k, dh_k
    export k_mesh_evals, plot_kresolved
    export AM_presets_3d, am_presets_3d, k_hamiltonian_3d, dh_k_3d, ddh_k_3d, sigma_abc_wrapper_3d
    export iwave_bands_3d, k_mesh_evals_3d, quantum_sigma_abc_wrapper_3d
    export iwave_Fijcuts, custom_fig_iwave_Fijcuts
    export slurm_conductivities, processing, create_bashfile
end
