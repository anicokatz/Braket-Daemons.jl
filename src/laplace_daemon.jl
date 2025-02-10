# program manager for the time-evolution aspects of the Quera project (integrate basic ITensor DMRG as well)
# takes in input_params.json file and:
#   - executes maxwell_daemon.jl to generate AnalogHamiltonianSimulation object with:
#       - params from model and sweep parameters from params.json
#       - time series defined by either adiabatic_program.jl or local_detuning.jl
#       - square lattice geometry
#   - executes descartes_daemon.jl with:
#       - output of maxwell_daemon.jl:
#       - merged system/measurement/simulation (truncation) parameters defined by params.json 
#   - NOT IMPLEMENTED: generate plots, saved to a subfolder, after data processing and backup

using JSON3

# TIDY UP IMPORTS AND PACKAGES - ESPECIALLY IN THE PROTOCOL FILES

@info "Parsing $(ARGS[1])..."
open(ARGS[1], "r") do file
    global input_json_string = JSON3.read(file)
end

input_params = Dict{String, Any}()

# manually parse JSON object into dict (I hate this, we should've used HDF5 or a better JSON package lol)
for (k, v) in input_json_string
    if String(k) == "path" # only the base path is not of type Dict{String, Any}
        input_params[String(k)] = v
    else # everything else we unpack, recasting symbols to Strings manually as we go (for some reason this can't happen automatically)
        tmp_dict = Dict{String, Any}()
        for (ks, vs) in v
            tmp_dict[String(ks)] = vs
        end
        input_params[String(k)] = tmp_dict
    end
end

# GO THROUGH AND SAVE OUTPUT PARAMETERS (Omega values, number of atoms, runtimes, experiment path etc.)
output_params = Dict()

# run maxwell_daemon.jl
@info "Executing maxwell_daemon.jl ..."
include("maxwell_daemon.jl")

# load into output parameters
output_params["omega-max"] = omega_max

# generate relevant time series (drive and shift) for the appropriate sweep type
@info "Generating ahs_object ..."
include("protocols/$(input_params["simulation-parameters"]["sim-type"]).jl") # set up a try/catch here to catch error re: improper sweep

# define compatible integration time step as a function of the step count and sweep protocol
input_params["simulation-parameters"]["tau"] = time_max/input_params["simulation-parameters"]["n-tau-steps"]

# generate filepaths for the json object and for the experiment output
# MUST OCCUR BEFORE EXECUTING DESCARTES
program_path = joinpath(input_params["path"], "ahs_object.json")
results_path = joinpath(input_params["path"], "results")

# process resulting ahs_object to a JSON object
ahs_json = JSON3.read(JSON3.write(ir(ahs_program)))

# Write the JSON object to a file
# REDUNDANCY, FOR SOME REASON mps_utils.jl SAVES THE JSON OBJECT AGAIN IN THE RESULTS FOLDER?
open(program_path, "w") do file
    JSON3.write(file, ahs_json)
end

# run Descartes: takes ahs program, merges parameter directories, and executes run() from mps_utils to run the simulation
@info "Executing descartes_daemon.jl ..."
include("descartes_daemon.jl")

# save results
# NOTE: results saving is dreadful, mixing filetypes and uses operator adjoints to data containing potentially complex errors?
# save_results function needs a total overhaul
@info "Saving results"
save_results(results, results_path)

# Add this line before using plot_all
include("plotter.jl")

# plot results
if input_params["meta-parameters"]["generate-plots"]
    figures_path = joinpath(input_params["path"], "figures")
    @info "Plotting results from $results_path to $figures_path"
    results_file = joinpath(results_path, "processed_results.h5")   
    plot_all(results_file, figures_path)
    @info "Plotting complete."
end

@info "Finished execution."