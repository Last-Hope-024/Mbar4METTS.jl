using Printf
using HDF5
using Dumper
using DataInterpolations
using LinearAlgebra
using Statistics


function reweight_first_dev_temp_observable(beta_final::AbstractVector{Float64}, path_file::AbstractString, observable_tag::AbstractString)
    #Load relevant data 
    f = h5open(path_file, "r")

    # ingridients to do it:
    betas = read(f["betas_traj"])[:] #get betas list 
    weights = read(f["norm_metts"])[:, :]  #get weights
    energy = read(f["energy"])[:, :]  #get energie
    observable = read(f[observable_tag])[:, :]

    # Interpolate:
    observable_interp = interpolate_observable(beta_final, betas, observable)
    dev_observable_interp = interpolate_observable_dev(beta_final, betas, observable)
    weights_interp = interpolate_observable(beta_final, betas, weights)
    energy_interp = interpolate_observable(beta_final, betas, energy)


    # Final array:
    data = zeros(length(beta_final), 2)

    #main loop:
    for (i, b) in enumerate(beta_final)
        w = exp.(weights_interp[:, i])
        ob = observable_interp[:, i]
        e = energy_interp[:,i]
        db_do = dev_observable_interp[:,i]
        data[i, :] .= jackknife_weighted_ratio_derivative(ob, db_do,e,w) #get the error now!
    end
    return data
end


function reweight_observable(beta_final::AbstractVector{Float64}, path_file::AbstractString, observable_tag::AbstractString)
    #Load relevant data 
    f = h5open(path_file, "r")

    # ingridients to do it:
    betas = read(f["betas_traj"])[:] #get betas list 
    weights = read(f["norm_metts"])[:, :]  #get weights
    observable = read(f[observable_tag])[:, :]

    # Interpolate:
    observable_interp = interpolate_observable(beta_final, betas, observable)
    weights_interp = interpolate_observable(beta_final, betas, weights)


    # Final array:
    data = zeros(length(beta_final), 2)

    #main loop:
    for (i, b) in enumerate(beta_final)
        w = exp.(weights_interp[:, i])
        ob = observable_interp[:, i]
        data[i, :] .= jackknife_weighted_ratio(ob, w) #get the error now!
    end
    return data
end


@doc raw"""
    jackknife_weighted_ratio(data::AbstractVector{Float64}, weights::AbstractVector{Float64})
    
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
    σ = sqrt((N - 1) * mean((estimates .- μ) .^ 2))

    return μ, σ
end

function jackknife_weighted_ratio_derivative(data::AbstractVector{Float64}, dev_data::AbstractVector{Float64}, energy::AbstractVector{Float64}, weights::AbstractVector{Float64})

    N = length(data)

    @assert length(weights) == N "data and weights must be the same length"

    # Precompute full sums
    data_weighted = dot(data, weights)  # Ob * weight

    data_product_energy = dot(data .* weights, energy) # Ob * E * weight

    total_dev_weighted = dot(dev_data, weights) # d_beta Ob * weight

    total_energy = dot(energy, weights) # E * weight

    total_weights = sum(weights) # weights

    # Allocate jackknife estimates
    estimates = similar(data)

    # Leave-one-out ratio estimates
    for i in 1:N
        dev_numerator = total_dev_weighted - dev_data[i] * weights[i]

        data_num = data_weighted - data[i] * weights[i]

        energy_num = total_energy - energy[i] * weights[i]

        data_product_energy_num = data_product_energy - data[i] * energy[i] * weights[i]

        denominator = total_weights - weights[i]

        energy_mean = energy_num / denominator


        estimates[i] = (data_num * energy_mean - data_product_energy_num + dev_numerator) / denominator
    end

    μ = mean(estimates)
    σ = sqrt((N - 1) * mean((estimates .- μ) .^ 2))

    return μ, σ
end


function interpolate_observable(beta_final::AbstractVector{Float64}, betas::AbstractVector{Float64}, data::AbstractArray{Float64})

    dataf = zeros(size(data, 1), length(beta_final))

    for i in 1:size(data, 1)
        r = CubicSpline(data[i, :], betas; extrapolation=ExtrapolationType.Extension) # we will defntly compute out of the domain but it is not a problem!
        dataf[i, :] .= r.(beta_final)
    end

    return dataf
end


function interpolate_observable_dev(beta_final::AbstractVector{Float64}, betas::AbstractVector{Float64}, data::AbstractArray{Float64})

    dataf = zeros(size(data, 1), length(beta_final))

    db = min(beta_final[end] - beta_final[end-1],0.01)

    for i in 1:size(data, 1)
        r = CubicSpline(data[i, :], betas; extrapolation=ExtrapolationType.Extension) # we will defntly compute out of the domain but it is not a problem!
        
        dataf[i, :] .= (r.(beta_final .- 2*db) + 8.0 .* r.(beta_final .+ db) - 8.0 .* r.(beta_final .- db) - r.(beta_final .+ 2*db)) ./(12*db)
    
    end

    return dataf
end