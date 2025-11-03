using Printf
using HDF5
using Dumper
using DataInterpolations
using LinearAlgebra
using Statistics

function compute_ess(beta_final::AbstractVector{Float64}, path_file::AbstractString)
    #Load relevant data 
    f = h5open(path_file, "r")

    # ingridients to do it:
    betas = read(f["betas_traj"])[:] #get betas list 
    weights = read(f["norm_metts"])[:, :]  #get weights

    # Interpolate:
    weights_interp = interpolate_observable(beta_final, betas, weights)

    # Final array:
    data = zeros(length(beta_final))
    data2 = zeros(length(beta_final))

    #main loop:
    for (i, b) in enumerate(beta_final)
        wmax = maximum(weights_interp[:, i])
        w = exp.(weights_interp[:, i] .- wmax) # all weights are smaller than one!
        data[i] .= mean(w) 
        w = exp.(2.0*(weights_interp[:, i] .- wmax)) # all weights are smaller than one!
        data2[i] .= mean(w) 
    end

    return (data .* data) ./ data2 * 1.0
end