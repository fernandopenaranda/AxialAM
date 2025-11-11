"""
bands computed from the k-space hamiltonian
usage: `bands(sgnum, N, pam; high_sym_line)`
    with `sgnum` the symmetry point group 17 for triangular lattice
        `N` number of k-points in the k-path
        `pam::AM_presets` the model presets
        the opt kwarg `high_sym_line` is a vector containing a path between high-symmetry-momenta
"""
function bands(sgnum, N, pam::AM_presets; high_sym_line = nothing)
    fig = Figure()
    Rs = lattice_vectors(pam.a0)
    ax = Axis(fig[1,1])
    dim = size(k_hamiltonian([0.,0.]),1)
    function eigenvals_model(k, i)
        es, _ = eigen(k_hamiltonian(k, pam))
        return real(es[i])
    end
    for i in 1:dim
        obs(k) = eigenvals_model(k, i)
        if high_sym_line === nothing
            BZpaths.plot_observable_in_kpath!(ax, obs, Rs, sgnum, N)
        else
            BZpaths.plot_observable_in_kpath!(ax, obs, Rs, sgnum, N, high_sym_line)
        end
    end
    return fig
end