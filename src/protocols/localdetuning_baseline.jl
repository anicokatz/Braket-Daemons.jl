# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using ITensors
using CSV, DataFrames
using Dates
using Base.Filesystem
using Missings
using Random
using JSON3

using Markdown
using InteractiveUtils

using Braket
using Braket: AtomArrangement, AtomArrangementItem, TimeSeries, DrivingField, Pattern
using DataStructures, Statistics, Plots

#=
begin
	a = 6.8e-6
	L = 3  # number of atoms along x/y dimension
	register = AtomArrangement()
	
	for i in 0:(L-1)
		for j in 0:(L-1)
			push!(register, AtomArrangementItem((Float64(i), Float64(j)) .* a))
		end
	end
end
=#

global time_max = 2.5e-6


begin
	t1               		= 5e-7  # seconds
	t2			   			= 1e-6  # seconds
	t3						= 1.5e-6  # seconds
	t4						= 2e-6  # seconds
	t_max			        = time_max  # seconds
	Ω_max                   = 15800000.  # rad / sec
	Δ_start                 = 0 * Ω_max
	Δ_end                   = 0 * Ω_max
	
	Ω                       = TimeSeries()
	Ω[0.0]                  = 0.0
	Ω[t1]			        = 0.0
	Ω[t2]					= Ω_max
	Ω[t3]				 	= Ω_max
	Ω[t4		]	 		= 0.0
	Ω[t_max]				= 0.0
	
	Δ                       = TimeSeries()
	Δ[0.0]                  = Δ_start
	Δ[t_max]				= Δ_end

	Δ_loc                   = TimeSeries()
	Δ_loc[0.0]              = 0.
	Δ_loc[t1]              	= 125000000
	Δ_loc[t2]		        = 125000000
	Δ_loc[t3]				= 125000000
	Δ_loc[t4]				= 125000000
	Δ_loc[t_max]			= 0.

	ϕ           			= TimeSeries()
	ϕ[0.0]     				= 0.0
	ϕ[t_max]				= 0.0
end

begin
	drive = DrivingField(Ω, ϕ, Δ)
	# pt = Pattern([0. for i in 1:length(register)])
	# Below we change the pattern so as to make the local detuning prepare the checkerboard state
	# pt is an alternating pattern of 0s and 1s
	pt = Pattern([((-1.)^i+1.)/2 for i in 1:length(register)])
	shift = ShiftingField(Field(Δ_loc, pt))
	global ahs_program = AnalogHamiltonianSimulation(register, [drive, shift])
end

#=
json_str = JSON3.write(ir(ahs_program))
json_obj = JSON3.read(json_str)

# Define the file path
file_path = joinpath(dirname(@__DIR__), "local_det_prep", "N_$(length(register)).json")
# file_path = joinpath(dirname(@__DIR__), "examples", "ahs_program_default.json")

# Write the JSON object to a file
open(file_path, "w") do file
    JSON3.write(file, json_obj)
end
=#

