using Printf
using HDF5
using Dumper
using DataInterpolations
using LinearAlgebra
using Statistics



function reweight_observable(beta_final::AbstractVector{Float64},path_file::AbstractString,observable_tag::AbstractString)
    #Load relevant data 
    f = h5open(path_file, "r")
    
    # ingridients to do it:
    betas = read(f["betas_traj"])[:] #get betas list 
    weights = read(f["norm_metts"])[:,:]  #get weights
    observable = read(f[observable_tag])[:,:]
     
    # Interpolate:
    observable_interp = interpolate_observable(beta_final, betas, observable)
    weights_interp = interpolate_observable(beta_final, betas, weights)
    
    println("Before jackknife")

    # Final array:
    data = zeros(length(beta_final), 2)
    
    #main loop:
    for (i, b) in enumerate(beta_final)
        w =  exp.(weights_interp[:, i])
        ob = observable_interp[:, i]
        data[i, :] .= jackknife_weighted_ratio(ob, w) #get the error now!
    end
    return data
end


@doc raw"""
    jackknife_weighted_ratio(data::Vector{Float64}, weights::Vector{Float64}) -> (Float64, Float64)

Compute the jackknife estimate and standard error of a weighted average of the form:

    ⟨A⟩ = ∑ᵢ Aᵢ wᵢ / ∑ᵢ wᵢ

using the leave-one-out jackknife method.

# Arguments
- `data::Vector{Float64}`: The data values Aᵢ to be averaged.
- `weights::Vector{Float64}`: Corresponding weights wᵢ for each data point.

# Returns
- `μ::Float64`: The jackknife estimate of ⟨A⟩.
- `σ::Float64`: The standard error estimated via jackknife:
  
      σ² = (N - 1) * mean((Aᵢ^{(jack)} - μ)²)

# Details
This function is optimised for computing observables of the form:

    ⟨A⟩ = ∑ Aᵢ wᵢ / ∑ wᵢ
"""
function jackknife_weighted_ratio(data::AbstractVector{Float64}, weights::AbstractVector{Float64})
    
    N = length(data)
    
    @assert length(weights) == N "data and weights must be the same length"
    
    # Precompute full sums
    total_weighted = dot(data, weights)
    total_weights = sum(weights)
    
    # Allocate jackknife estimates
    estimates = similar(data)
    
    # Leave-one-out ratio estimates
    for i in 1:N
        numerator = total_weighted - data[i] * weights[i]
        denominator = total_weights - weights[i]
        estimates[i] = numerator / denominator
    end

    μ = mean(estimates)
    σ = sqrt((N - 1) * mean((estimates .- μ).^2))

    return μ, σ
end


function interpolate_observable(beta_final::AbstractVector{Float64}, betas::AbstractVector{Float64}, data::AbstractArray{Float64})
        
    dataf = zeros(size(data, 1), length(beta_final))

    for i in 1:size(data, 1)
        r = CubicSpline(data[i,:], betas; extrapolation = ExtrapolationType.Extension) # we will defntly compute out of the domain but it is not a problem!
        dataf[i, :] .= r.(beta_final)
    end
    
    return dataf
end