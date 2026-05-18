DEPOT_PATH[:] .= ["/scratch/ferpe/julia_depot"]
using Optics_in_the_length_gauge, CSV, JLD2, Serialization, AxialAM
job_id = parse(Int, ARGS[1]) # Get the task ID from SLURM
jobs_num = parse(Int, ARGS[2])
PID = ARGS[3]
dirj =  Symbol(ARGS[4])
dirE =  Symbol(ARGS[5])
dirB =  Symbol(ARGS[6])
T = parse(Float64, ARGS[7]) 
evals = parse(Int, ARGS[8])
omega_switch = parse(Bool, ARGS[9])
ps_switch = parse(Bool, ARGS[10])
ps_orbital_switch = parse(Bool, ARGS[11])
qm_switch = parse(Bool, ARGS[12])
fermi_surface = parse(Bool, ARGS[13])
epsilon = parse(Float64, ARGS[14])
which_mm = Symbol(ARGS[15])
integration_method = Symbol(ARGS[16])
a0 = parse(Float64, ARGS[17])
t1 = parse(Float64, ARGS[18])
t2 = parse(Float64, ARGS[19])
t4 = parse(Float64, ARGS[20])
tperp = parse(Float64, ARGS[21])
Deltasz = parse(Float64, ARGS[22])
Deltatauz = parse(Float64, ARGS[23])
neel1 = parse(Float64, ARGS[24])
neel2 = parse(Float64, ARGS[25])
neel3 = parse(Float64, ARGS[26])
mumin = parse(Float64, ARGS[27])
mumax = parse(Float64, ARGS[28])
mupoints= parse(Int, ARGS[29])

print("Starting...")
if mupoints == 0 || mupoints == 1
    muvec = [mumin]
else
    muvec = collect(range(mumin, mumax, length=mupoints))
end

subcubes_file = pwd() * "/subcubes.jls"
subcubes = deserialize(subcubes_file)
my_subcube = subcubes[job_id]
botbounds = [my_subcube[1][1], my_subcube[2][1], my_subcube[3][1]]
topbounds = [my_subcube[1][2], my_subcube[2][2], my_subcube[3][2]]

p = AM_presets_3d(am_presets_3d() , a0 =a0,  t1 = t1, t4 = t4, t2 = t2, 
    tperp = tperp, Δsz = Deltasz, Δtz = Deltatz, neel = [neel1, neel2, neel3])

keyws = (
    dirj = dirj,
    dirE = dirE,
    dirB = dirB,
    T = T,
    evals = evals,
    omega_MM_switch = omega_switch,
    PS_switch = ps_switch,
    PS_orbital_switch = ps_orbital_switch,
    QM_switch = qm_switch,
    fermi_surface = fermi_surface,
    epsilon = epsilon,
    which_mm = which_mm,
    integration_method = integration_method,
    botbounds = botbounds,
    topbounds = topbounds,
)

sijks = Float64[]
for mu in muvec  
    sijk_pres = quantum_sigma_abc_wrapper_3d(AM_presets_3d(p, μ = mu); keyws...) #computing struct for each mu
    push!(sijks,quantum_contribution(sijk_pres)) # compute
end
sijk_pres = quantum_sigma_abc_wrapper_3d(p; keyws...)

#_________________________________________________________________________________________
#store
#_________________________________________________________________________________________
data_folder = pwd() * "/Data/" * string(PID) * "/" * string(job_id)
mkpath(data_folder)
@save data_folder * "/presets.jld" sijk_pres
@save data_folder * "/calculation.jld" muvec sijks

str = pwd() * "/slurm-" * string(PID) * "." * string(job_id)
isfile(str * ".out") && mv(str * ".out", data_folder * "/output.out")
isfile(str * ".err") && mv(str * ".err", data_folder * "/error.err")