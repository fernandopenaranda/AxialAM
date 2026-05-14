"""
bands computed from the k-space hamiltonian
usage: `bands(sgnum, N, pam; high_sym_line)`
    with `sgnum` the symmetry point group 17 for triangular lattice
        `N` number of k-points in the k-path
        `pam::AM_presets` the model presets
        the opt kwarg `high_sym_line` is a vector containing a path between high-symmetry-momenta
"""
function eigenvals_model(pam, k, i)
    es, _ = eigen(k_hamiltonian(k, pam))
    return real(es[i])
end

function eigenvals_model_3d(pam, k, i)
    es, _ = eigen(k_hamiltonian_3d(k, pam))
    return real(es[i])
end