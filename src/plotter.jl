# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using CairoMakie
using CairoMakie: heatmap, scatter, lines, Colorbar
using Base
using ITensors.HDF5 

function plot_density(data, fig_path)

    vals = data["density"]
        
    num_atoms, num_times = size(vals)
    
    ts = [i for i in 0:num_times-1]
    xs = [i for i in 0:num_atoms-1]
    
    # all the data was hermitian conjugated here too? Ask Ale if this has been done anywhere else. / NB use permuteaxis for better extensibility
    fig, ax, hm = CairoMakie.heatmap(xs, ts, vals, colorrange = (0, 1), axis=(;title = "Density evolution: n_i(t)", xlabel = "atom index, i", ylabel="time step, t_j"))
    Colorbar(fig[:, end+1])  # equivalent
    plot_path = joinpath(fig_path, "density.png")
    save(plot_path, fig)

    fig, ax, plt = lines(xs, real(vals[:, end]),  axis=(;title = "Final density: n_i(T)", xlabel = "atom index, i", ylabel="density, n_i"))
    plot_path = joinpath(fig_path, "final_density.png")
    save(plot_path, fig)
    
end

function plot_staggered_magnetization(data, fig_path)

    vals = data["density"]
        
    mags, ts, Nx, Ny = calculate_staggered_magnetization(real(vals))

    fig, ax, plt = lines(ts, mags, xlabel = "time, t", ylabel="staggered magnetization, \$ M_{S} \$")
    plot_path = joinpath(fig_path, "staggered_magnetization.png")
    save(plot_path, fig)
    
end

function plot_fidelity_susceptibility(data, fig_path)

    vals = data["fidelity_susceptibility"]

    num_times = length(vals)
    
    ts = [i for i in 0:num_times-1]

    fig, ax, plt = lines(ts, real(vals), xlabel = "time, t", ylabel="fidelity susceptibility \$ \\chi_F \$")
    plot_path = joinpath(fig_path, "fidelity_susceptibility.png")
    save(plot_path, fig)

end

function plot_truncation_error(data, fig_path)

    vals = data["truncation_error"]
        
    num_times = length(vals)
    
    ts = [i for i in 0:num_times-1]

    fig, ax, plt = lines(ts, real(vals), xlabel = "time, t", ylabel="truncation_error, \$ \\mathcal{F}(\\psi, \\psi_{TRUE}) \$")
    plot_path = joinpath(fig_path, "truncation_error.png")
    save(plot_path, fig)
    
end

# function plot_atoms(experiment_path, fig_path)
#     csv_data = CSV.File(joinpath(experiment_path, "atom_coordinates.csv"))
#     df = DataFrame(csv_data)
#     coords = Matrix{Float64}(df)

#     fig, ax, hm = scatter(coords[1, :], coords[2, :], axis=(;title = "Atom coordinates: (x, y)", xlabel = "atom coordinate, x", ylabel="atom coordinate, y"))

#     plot_path = joinpath(fig_path, "atom_coordinates.png")
#     save(plot_path, fig)
# end

# function plot_bitstrings(experiment_path, fig_path, max_vals=20)
#     csv_data = CSV.File(joinpath(experiment_path, "bitstrings.csv"))
#     df = DataFrame(csv_data)
#     bit_values = Matrix{Int}(df)
#     bitstrings = [join(string.(row)) for row in eachrow(bit_values)]

#     # Create a counter (dictionary) from the list
#     counter = Dict{String, Int}()
#     for str in bitstrings
#         counter[str] = get(counter, str, 0) + 1
#     end    

#     # Create a list of counts corresponding to each string    
#     ks = collect(keys(counter))
#     vs = collect(values(counter))
    
#     # convert bitstrings to integer representation
#     # Sort by counts
#     sorted_indxs = sortperm(vs, rev = true)
#     num_vals = length(vs)

#     fig, ax, plt = barplot(1:num_vals, vs[sorted_indxs], axis = (xticks = (1:num_vals, ks[sorted_indxs]),
#                            title = "Sampled bitstrings: t=T", xlabelrotation=1), )
#     plot_path = joinpath(fig_path, "mps_samples.png")
#     save(plot_path, fig)
# end

function calculate_staggered_magnetization(mags)
    Nx = Ny = Int(sqrt(size(mags, 1)))
    num_times = size(mags, 2)
    ts = [i for i in 0:num_times-1]

    staggered_magnetization = zeros(Float64, num_times)
    for j in 1:Nx
        for i in 1:Ny
            idx_phys = Int(Ny*(j-1) + i)
            m = ( (-1)^(j+i - 2) ) * mags[idx_phys, :]
            staggered_magnetization += m
        end
    end
    return staggered_magnetization, ts, Nx, Ny
end

function plot_variance_staggered_magnetization(experiment_path, fig_path)
    # Get z_profile data for staggered magnetization
    csv_data = CSV.File(joinpath(experiment_path, "z_profile.csv"))
    df = DataFrame(csv_data)
    mags = Matrix{Float64}(df)
    staggered_magnetization, ts, Nx, Ny = calculate_staggered_magnetization(mags)
    variance_staggered_magnetization = -(staggered_magnetization .^ 2)

    # Get correlator data
    csv_data = CSV.File(joinpath(experiment_path, "correlator_zz.csv"))
    df = DataFrame(csv_data)
    corr_zz = Matrix{Float64}(df)
    
    # Calculate variance using the single correlator matrix
    corr_t = corr_zz  # Already in the right shape
    sum_corr = 0.0
    for j in 1:Nx
        for i in 1:Ny
            for jp in 1:Nx
                for ip in 1:Ny
                    # Subtract 1 from indices since we're 0-based indexing in the physical system
                    idx_phys_1 = Int(Ny*(j-1) + i) - 1
                    idx_phys_2 = Int(Ny*(jp-1) + ip) - 1
                    # Add 1 back for Julia's 1-based array indexing
                    m = ( (-1)^(j+i+jp+ip - 4) ) * corr_t[idx_phys_1 + 1, idx_phys_2 + 1]
                    sum_corr += m
                end
            end
        end
    end
    variance_staggered_magnetization[end] += sum_corr

    fig, ax, plt = lines(ts, variance_staggered_magnetization, axis=(;title = "Variance of staggered magnetization: <(Sz_i(T) - <Sz_i(T)>)^2>", xlabel = "time step, t_j", ylabel="variance of staggered magnetization, <(Sz_i(T) - <Sz_i(T)>)^2>"))
    plot_path = joinpath(fig_path, "variance_staggered_magnetization.png")
    save(plot_path, fig)
end

function plot_all(data, fig_path)
    # plot_correlator(experiment_path)
    # plot_atoms(experiment_path, fig_path)
    # plot_bitstrings(experiment_path, fig_path)
    # plot_z_profile(experiment_path)
    # plot_variance_staggered_magnetization(experiment_path)

    plot_density(data, fig_path)
    plot_staggered_magnetization(data, fig_path)

end

# Add main entry point to handle command line arguments
if abspath(PROGRAM_FILE) == @__FILE__

    if length(ARGS) != 1
        println("Usage: julia --project=. src/plotter.jl <experiment_path>")
        exit(1)
    end

    data_path = ARGS[1]

    h5open(data_path, "r") do fr
        global data = read(fr)
    end

    # @show data
    # @show data["density"]

    fig_path = joinpath(dirname(dirname(data_path)), "figures")

    plot_density(data, fig_path)
    plot_fidelity_susceptibility(data, fig_path)
    # plot_truncation_error(data, fig_path)
    plot_staggered_magnetization(data, fig_path)

end

# NOTE: ADD SUM OF DISCARDED EIGENVALUES TO THE mps_utils.jl CODE IF POSSIBLE