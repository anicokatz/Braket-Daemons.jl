# Maxwell takes the model parameters and calculates simple important processed objects (Omega_max, etc.)
# Maxwell also generates the atomic arrangement
# These are then fed to the relevant backend time series generator

using Braket: AtomArrangement, AtomArrangementItem

begin

	L = input_params["model-parameters"]["L"]
	a = input_params["model-parameters"]["interaction-radius"]

    global register = AtomArrangement()

	for i in 0:(L-1)
		for j in 0:(L-1)
			push!(register, AtomArrangementItem((Float64(i), Float64(j)) .* a))
		end
	end

    ratio = input_params["sweep-parameters"]["rydberg-ratio"]
    C6 = input_params["model-parameters"]["C6"]
    a = input_params["model-parameters"]["interaction-radius"]
	Rb = a * ratio

	# calculate implicit target value of Omega
	global omega_max = C6/(Rb^6)

end
