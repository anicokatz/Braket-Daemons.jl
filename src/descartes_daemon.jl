# Descartes unpacks relevant system/simulation parameters, reads the AHS object from memory (not from file) and executes the MPS evolution program
# Descartes does not save, post-process, or plot

using ITensors
using Dates
using Missings
using Random
using Logging
using JSON3

include("mps_utils.jl")

# merge and generate additional info
args_for_mps = merge(input_params["sweep-parameters"], 
                    input_params["measurement-parameters"],
                    input_params["model-parameters"], 
                    input_params["simulation-parameters"])
args_for_mps["number-of-atoms"] = input_params["model-parameters"]["L"] ^ 2
args_for_mps["experiment-path"] = results_path # NOTE: mps_utils.jl uses experiment-path so I've left in unchanged, but this is inconsistent with the daemons
args_for_mps["program-path"] = program_path

ahs_json = JSON3.read(read(program_path, String), Dict{String, Any})
# println(ahs_json)

# run ahs program, computing errors/observables/susceptibilities (if appropriate) at each timestep
results = run(ahs_json, args_for_mps)

#@info "Generating plots"
#if args["generate-plots"]
#    @info "Plotting results from $results_path"
#    plot_all(results_path)
#    @info "Plotting complete."
#end