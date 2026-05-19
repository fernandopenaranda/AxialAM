function compute_conductivities(partitions, evaluations; kws...)
    create_bashfile(partitions)
    slurm_conductivities(evals = evaluations; kws...)
end

evals = 20

kws = (a0 = 1.0, t1 = 0.3, t2 = 0, t4 = 1.0, tperp = 1.0, Deltasz = 1.0, Deltatauz = 3.0,
     neel= [0,0,1.0], integration_method = :montecarlo,
     T = 300,  mumin = -10, mumax = 10,  mupoints = 201)

function sweep_cond(dirj, dirE, dirB; evals = 20, subcubes = 150,
      omega_switch = false, ps_switch = true, ps_orbital_switch = true, qm_switch = false, kws...)
    compute_conductivities(subcubes, evals; 
         dirj= dirj, dirE=dirE, dirB = dirB, omega_switch =  omega_switch, ps_switch = ps_switch, 
         ps_orbital_switch = ps_orbital_switch, qm_switch = qm_switch, fermi_surface = false, kws...)
end



function axial_sweep_cond(dirj, dirE, dirB; evals = 20, subcubes = 20,
    omega_switch = false, ps_switch = true, ps_orbital_switch = true, qm_switch = false, kws...)
    
    computi(dirj, dirE, dirB) = sweep_cond(dirj, dirE, dirB; 
        evals = evals, subcubes = subcubes, omega_switch =  omega_switch, ps_switch = ps_switch, 
        ps_orbital_switch = ps_orbital_switch, qm_switch = qm_switch, fermi_surface = false, kws...)

    computi(:x,:y,:y)
    computi(:x,:z,:z)
    computi(:y,:x,:x)
    computi(:y,:z,:z)
    computi(:z,:x,:x)
    computi(:z,:y,:y)     
end