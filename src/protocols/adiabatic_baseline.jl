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

begin

	# hard-coded fundamental limits
	global time_max         = 4e-6  # seconds, IMPORTANT: global scope required here as its used to construct timestep tau
	
    # main ahs_program preparation - varies between sweep protocols

    time_ramp               = 1e-7  # seconds

	Ω_max                   = output_params["omega-max"]  # rad / sec, CHECK THIS WORKS IN SAVED PARAMS
	Δ_start                 = input_params["sweep-parameters"]["delta-min"] * Ω_max
	Δ_end                   = input_params["sweep-parameters"]["delta-max"] * Ω_max
	
	Ω                       = TimeSeries()
	Ω[0.0]                  = 0.0
	Ω[time_ramp]            = Ω_max
	Ω[time_max - time_ramp] = Ω_max
	Ω[time_max]             = 0.0
	
	Δ                       = TimeSeries()
	Δ[0.0]                  = Δ_start
	Δ[time_ramp]            = Δ_start
	Δ[time_max - time_ramp] = Δ_end
	Δ[time_max]             = Δ_end

	Δ_loc                   = TimeSeries()
	Δ_loc[0.0]              = 0.
	Δ_loc[time_max]         = 0.

	ϕ           = TimeSeries()
	ϕ[0.0]      = 0.0
	ϕ[time_max] = 0.0

end

begin
	
	drive = DrivingField(Ω, ϕ, Δ)

	pt = Pattern([0. for i in 1:length(register)])
	shift = ShiftingField(Field(Δ_loc, pt))

	global ahs_program = AnalogHamiltonianSimulation(register, [drive, shift])

end