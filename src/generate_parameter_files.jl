using JSON3

model_parameters = Dict{String, Any}(
                        "L" => 5,                                   # size of LxL square array of atoms, atom number calculated in lattice mapper
                        "interaction-radius" => 7e-6,               # the interaction radius in meters
                        "C6" => 5.42e-24)                           # C6 interaction constant

sweep_parameters = Dict{String, Any}(
                        "delta-min" => -3.,                         # minimum delta/omega ratio
                        "delta-max" => 3.,                         # maximum delta/omega ratio
                        "switch-time" => 0.1,                       # NOT IMPLEMENTED (deprecated - protocol specific): time at which delta and omega scans switch (in microseconds)
                        "sweep-type" => 'r',                        # NOT IMPLEMENTED (deprecated - protocol specific): delta ramp type: options are 'l' for linear, or 'e' for exponential/logarithmic
                        "rydberg-ratio" => 1.2)                     # rydberg Rb/a ratio, implicitly determines maximum Omega value             
                    
simulation_parameters = Dict{String, Any}(
                        "sim-type" => "adiabatic_baseline",         # name of the sweep protocol (also generates a file by this name)
                        "cutoff" => 1e-10,                           # cutoff for SVD values in MPS evolution
                        "max-bond-dim" => 50,                       # maximum bond dimension for MPS
                        "shots" => 2000,                            # number of shots for sampling
                        "n-tau-steps" => 400)                        # time evolution step count (does this take precedence over above?)

measurement_parameters = Dict{String, Any}(
                        "compute-truncation-error" => false,         # compute truncation error during runtime
                        "compute-observables" => true,              # compute observables during runtime
                        "compute-fidelity-susceptibility" => true,  # compute fidelity susceptibility using raw value of tau as the perturbation step ALTER TO EXTRACT CHANGE IN DELTA
                        "sample-shots-runtime" => false)            # NOT YET IMPLEMENTED: samples shots at every timestep (potentially very expensive)

meta_parameters = Dict{String, Any}("generate-plots" => true)      # REDUNDANT: deprecated, generate plots after experiment is finished
#                        "path" => )                                 # store general working directory path

# PATHING PARAMETERS
# save path + files are determined by number of atoms (square number), max bond dimension, and sweep type ('r' or 'l')

L, chi, simtype, ratio = model_parameters["L"], simulation_parameters["max-bond-dim"], simulation_parameters["sim-type"], sweep_parameters["rydberg-ratio"];
ratio_str = replace(string(ratio), "." => "_")

function pathgen(path)
    if !isdir(path)
        mkdir(path)
    end
end

file_path = joinpath(dirname(@__DIR__), "data")
@info "Generating file structure in $file_path"
pathgen(file_path)

# NOTE: SWAPPED FROM ratio/L/type to type/ratio/L, SHOULD WORK - IF NOT, CHECK

file_path = joinpath(file_path, simtype)
pathgen(file_path)

file_path = joinpath(file_path, "ratio_"*ratio_str)
pathgen(file_path)

file_path = joinpath(file_path, "L_$L")
pathgen(file_path)

file_path = joinpath(file_path, "chi_$chi")
pathgen(file_path)

# generate core results and figures directories
pathgen(joinpath(file_path, "results"))
pathgen(joinpath(file_path, "figures"))

params = Dict("path" => file_path, # add path name to meta-params
            "model-parameters" => model_parameters,
            "sweep-parameters" => sweep_parameters,
            "simulation-parameters" => simulation_parameters,
            "measurement-parameters" => measurement_parameters,
            "meta-parameters" => meta_parameters)

open(joinpath(file_path, "input_params.json"), "w") do file
    JSON3.write(file, params)
end

@info "Saved input parameters to $file_path"