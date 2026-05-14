"""
structure wrapper botbounds and topbounds in [-0.5,0.5]
"""
function sigma_abc_wrapper(pam; dirj=:x, dirE=:y, dirB=:z, T = 1, τ = 200, evals = 100, fermi_surface = false,
    integration_method = :montecarlo, botbounds = [-0.5,-0.5], topbounds = [0.5,0.5])
    a, R1, R2 = iwave_lattice_vectors(pam.a0)
    unit_convention_two_packages_t = 1e-15
    τ *= unit_convention_two_packages_t
    #hamiltonians and derivatives
    h(q) = k_hamiltonian(q, pam)
    dhx(q) = dh_k(q, pam, :x)
    dhy(q) = dh_k(q, pam, :y)
    dh(q) = [dhx(q), dhy(q)]
    didjh(q,i,j) = ddh_k(q, pam, i, j) 
    ddh(q) = [[didjh(q,:x,:x), didjh(q,:x,:y)],
              [didjh(q,:y,:x), didjh(q,:y,:y)]]

    # integral presets
    
    Rs = [R1,R2]
    Gs = dualbasis(Rs)
    computation = Transport_computation_3d_presets(botbounds,topbounds, evals, integration_method)
    return Optics_in_the_length_gauge.Classical_σijk_antisym(dirj, dirE, dirB, h, dh, ddh, Gs, τ, T, computation, fermi_surface)
end

function sigma_abc_wrapper_3d(pam; dirj=:x, dirE=:y, dirB=:z, T = 1, τ = 200, evals = 100, fermi_surface = false,
    integration_method = :montecarlo, botbounds = [-0.5,-0.5, -0.5], topbounds = [0.5,0.5, 0.5])
    a, R1, R2 = iwave_lattice_vectors(pam.a0)
    R3 = a*[0,0,1]
    unit_convention_two_packages_t = 1e-15
    τ *= unit_convention_two_packages_t
    #hamiltonians and derivatives
    h(q) = k_hamiltonian_3d(q, pam)
    dhx(q) = dh_k_3d(q, pam, :x)
    dhy(q) = dh_k_3d(q, pam, :y)
    dhz(q) = dh_k_3d(q, pam, :z)
    dh(q) = [dhx(q), dhy(q), dhz(q)]
    didjh(q,i,j) = ddh_k_3d(q, pam, i, j) 
    ddh(q) = [[didjh(q,:x,:x), didjh(q,:x,:y), didjh(q,:x,:z)],
              [didjh(q,:y,:x), didjh(q,:y,:y), didjh(q,:y,:z)],
              [didjh(q,:z,:x), didjh(q,:z,:y), didjh(q,:z,:z)]]
    # integral presets
    Rs = [promote_3d(R1), promote_3d(R2), R3]
    Gs = dualbasis(Rs)
    computation = Transport_computation_3d_presets(botbounds,topbounds, evals, integration_method)
    return Optics_in_the_length_gauge.Classical_σijk_antisym(dirj, dirE, dirB, h, dh, ddh, Gs, τ, T, computation, fermi_surface)
end

function quantum_sigma_abc_wrapper_3d(pam; dirj=:x, dirE=:y, dirB=:z, T = 1, τ = 200, evals = 100,
    integration_method = :montecarlo, botbounds = [-0.5,-0.5, -0.5], 
    topbounds = [0.5,0.5, 0.5], 
    omega_MM_switch = true, PS_switch = true, PS_orbital_switch =true, QM_switch = true, fermi_surface = false,
    which_mm = :orbital, epsilon = 1e-5)
    a, R1, R2 = iwave_lattice_vectors(pam.a0)
    R3 = a*[0,0,1]
    unit_convention_two_packages_t = 1e-15
    τ *= unit_convention_two_packages_t
    #hamiltonians and derivatives
    h(q) = k_hamiltonian_3d(q, pam)
    dhx(q) = dh_k_3d(q, pam, :x)
    dhy(q) = dh_k_3d(q, pam, :y)
    dhz(q) = dh_k_3d(q, pam, :z)
    dh(q) = [dhx(q), dhy(q), dhz(q)]
    didjh(q,i,j) = ddh_k_3d(q, pam, i, j) 
    ddh(q) = [[didjh(q,:x,:x), didjh(q,:x,:y), didjh(q,:x,:z)],
              [didjh(q,:y,:x), didjh(q,:y,:y), didjh(q,:y,:z)],
              [didjh(q,:z,:x), didjh(q,:z,:y), didjh(q,:z,:z)]]
    # integral presets
    Rs = [promote_3d(R1), promote_3d(R2), R3]
    Gs = dualbasis(Rs)
    computation = Transport_computation_3d_presets(botbounds,topbounds, evals, integration_method)
    return Optics_in_the_length_gauge.Quantum_correction_σijk_antisym(a, dirj, dirE, dirB, h, dh, ddh, Gs, τ, T, computation, which_mm, omega_MM_switch, PS_switch, PS_orbital_switch, QM_switch, fermi_surface, epsilon)
end

function promote_3d(v)
    if length(v) == 3
        return v
    else
        return [v[1], v[2], 0.0]
    end
end