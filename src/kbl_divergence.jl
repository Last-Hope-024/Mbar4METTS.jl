
function kbl_divergence(beta_final::AbstractVector{Float64}, path_file::AbstractString; cut1::Int=1)
    #Load relevant data 
    f = h5open(path_file, "r")

    # ingridients to do it:
    betas = read(f["betas_traj"])[:] #get betas list 
    weights = read(f["norm_metts"])[:, cut1:end]  #get weights

    ########## Indicies that bound the values of the part actually computed ##############################
    idx1 = argmin(abs.(beta_final .-betas[1]))
    idx2 = argmin(abs.(beta_final .-betas[end]))

    # Interpolate:
    weights_interp = interpolate_observable(beta_final, betas, weights,idx1,idx2)

    # Final array:
    data = zeros(length(beta_final), 2)

    #main loop:
    for i in idx1:idx2
        data[i, :] .= jackknife_kbl(weights_interp[:, i]) #get the error now!
    end
    return data
end

function jackknife_kbl(weights::AbstractVector{Float64})

    N = length(weights)

    @assert length(weights) == N "data and weights must be the same length"

    # Precompute full sums using logsum
    total_weights = sum(weights)
    @show total_weights
    w_max = maximum(weights)
    
    pre_sum = sum(exp.(weights .- w_max))

    # Allocate jackknife estimates
    estimates = similar(weights)

    # Leave-one-out ratio estimates
    for i in 1:N
        partition =  log(pre_sum - exp(weights[i] - w_max))
        estimates[i] = -total_weights + weights[i] + w_max + partition
    end

    μ = mean(estimates)
    σ = sqrt((N - 1) * mean((estimates .- μ) .^ 2))

    return μ, σ
end